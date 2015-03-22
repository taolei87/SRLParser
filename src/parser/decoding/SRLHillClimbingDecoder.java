package parser.decoding;

import gnu.trove.list.array.TIntArrayList;

import java.util.Random;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import parser.*;
import utils.FeatureVector;
import utils.Utils;

public class SRLHillClimbingDecoder extends SRLDecoder {

	DependencyInstance pred, inst;
	int[][] goldlbids;
    
	int numframes;
	int len;
	
	SRLFeatureData sfd;
	boolean addLoss;
	
	double bestScore;	
	int unchangedRuns, totRuns;
	volatile boolean stopped;
	
	MaxMatchingDecoder decoder2;
	
    ExecutorService executorService;
	ExecutorCompletionService<Object> decodingService;
	HillClimbingTask[] tasks;
	
	public SRLHillClimbingDecoder(Options options) {
		this.options = options;
		this.nullWeight = options.nullWeight;
        executorService = Executors.newFixedThreadPool(options.numHcThreads);
		decodingService = new ExecutorCompletionService<Object>(executorService);
		tasks = new HillClimbingTask[options.numHcThreads];
		for (int i = 0; i < tasks.length; ++i) {
			tasks[i] = new HillClimbingTask();
			tasks[i].id = i;
			tasks[i].r = new Random(i);
		}
		
		decoder2 = new MaxMatchingDecoder(options);
	}

	@Override
	public void shutdown()
	{
		//System.out.println("shutdown");
		executorService.shutdownNow();
	}

	@Override
	public DependencyInstance decode(DependencyInstance inst,
			SRLFeatureData sfd, boolean addLoss) {
		this.inst = inst;
		this.sfd = sfd;
		this.addLoss = addLoss;
		
		pred = new DependencyInstance(inst);
		bestScore = Double.NEGATIVE_INFINITY;
		totRuns = 0;
		unchangedRuns = 0;
		stopped = false;
		
		numframes = inst.numframes;
		len = inst.length;
		
		// compute loss
		goldlbids = new int[numframes][len];
		if (addLoss) {
			for (int f = 0; f < numframes; ++f) {
				for (int i = 0, L = len; i < L; ++i) {
					// this is gold link
					goldlbids[f][i] = inst.frames[f].arglbids[i];
					
					// "project" gold links to the prediction
					// remove this link if it is not valid on the predicted tree
					if (goldlbids[f][i] >= 0 && sfd.isPruned(inst.frames[f].predid, i))
						goldlbids[f][i] = -1;
				}
			}
		}
		
		if (true) {
			DependencyInstance now = new DependencyInstance(inst);
			now = decoder2.decode(now, sfd, this.addLoss);
			TIntArrayList[] args = getArgs();
			hillClimbing(now, args);
			
        	double score = calcScore(now);
        	bestScore = score;
        	pred.frames = now.frames;
		}
		
		for (int i = 0; i < tasks.length; ++i) {
			decodingService.submit(tasks[i], null);			
		}
		
		for (int i = 0; i < tasks.length; ++i) {
			try {
				decodingService.take();
			} catch (InterruptedException e) {
				System.out.println("Semantic Hill climbing thread interupted!!!!");
			}
		}
		
        double goldScore = sfd.getScore(inst);
        double predScore = sfd.getScore(pred);
        double loss = predScore - goldScore;
        double dist = 0.0;
        
        if (addLoss) {
        	dist = sfd.getSRLCost(inst.frames, pred.frames);
        	loss += dist;
        }
		if (!addLoss && loss < -1e-6) {
			System.out.println("fail to find better result");
		}		
		return pred;		
	}
	
	private TIntArrayList[] getArgs() {
		TIntArrayList[] args = new TIntArrayList[numframes];
		for (int i = 0; i < numframes; ++i) {
			int pid = inst.frames[i].predid;
			args[i] = new TIntArrayList(5);
			for (int j = 0; j < len; ++j)
				if (!sfd.isPruned(pid, j))
					args[i].add(j);
		}
		return args;
	}
	
	private void hillClimbing(DependencyInstance now, TIntArrayList[] args) {
		// hill climbing
		boolean change = true;
		int loop = 0;
		while (change && loop < 100) {
			change = false;
			for (int i = 0; i < numframes; ++i) {
				boolean isChanged = findOptChange(goldlbids[i], now, sfd, i, args[i]);
				if (isChanged) {
					change = true;
				}
			}
			loop++;
		}
		if (loop >= 100) {
			System.out.println("too many loop: " + loop);
		}
	}

	public class HillClimbingTask implements Runnable {

		public int id;
		public Random r;
		
		int converge;
		int earlyStop;

		@Override
		public void run() {
			
			converge = options.numHcConverge;
			earlyStop = options.earlyStop;
			
			double goldScore = -Double.MAX_VALUE;
			if (addLoss) {
				goldScore = calcScore(inst);
			}

			DependencyInstance now = new DependencyInstance(inst);
			
			while (!stopped) {
				SemanticFrame[] predFrames = new SemanticFrame[numframes];
				SemanticFrame[] frames = inst.frames;
				
				// sequential sampling
				TIntArrayList[] args = getArgs();
				for (int i = 0; i < numframes; ++i) {
					predFrames[i] = sequentialSampling(frames[i], goldlbids[i], sfd, i, args[i]);
				}
				now.frames = predFrames;
				
				// hill climbing
				hillClimbing(now, args);

				double score = calcScore(now);
				synchronized (pred) {
					++totRuns;
					if (score > bestScore) {
						bestScore = score;

						if (addLoss && unchangedRuns >= earlyStop + options.numHcThreads && bestScore >= goldScore + 1e-6)
							System.out.print("(" + unchangedRuns + ") ");

						unchangedRuns = 0;
                        pred.frames = now.frames;
					} else {
						++unchangedRuns;
						if (unchangedRuns >= converge)
							stopped = true;
						
						if (addLoss && unchangedRuns >= earlyStop && bestScore >= goldScore + 1e-6)
							stopped = true;
					}
					
				}
			}
			
		}

		private SemanticFrame sequentialSampling(SemanticFrame goldFrame, int[] goldlbids, SRLFeatureData sfd,
				 int p, TIntArrayList args) {
			// p is the predicate id, 1 to #frames
			
			SemanticFrame predict = new SemanticFrame(goldFrame);
			predict.arglbids = new int[goldlbids.length];
			for (int i = 0, L = predict.arglbids.length; i < L; ++i)
				predict.arglbids[i] = -1;
			
			// construct the bi-partite graph
			int N = args.size(), M = sfd.L;
			//System.out.println("N: " + N + " M: " + M);
			int T = N + M;
			boolean[] usedRel = new boolean[T];
			double[] score = new double[T];
			
			// sample mappings for argument
			for (int a = 0; a < N; ++a) {
				for (int r = 0; r < T; ++r) {
					if (usedRel[r]) {
						score[r] = -Double.MAX_VALUE;
					}
					else {
						score[r] = r < M ? sfd.getArcScore(goldFrame.predid, args.get(a), r) : nullWeight;
						if (addLoss)
							score[r] += loss(goldlbids[args.get(a)], r < M ? r : -1);
					}
				}
				int sample = samplePoint(score, usedRel);
				usedRel[sample] = true;

				// convert bi-partite graph to semantic frame
				predict.arglbids[args.get(a)] = sample < M ? sample : -1;
			}
			
			// no need to sample mappings for null node
			
			return predict;
		}
		
		private int samplePoint(double[] score, boolean[] used) {
			double sumScore = Double.NEGATIVE_INFINITY;
			for (int i = 0; i < score.length; i++) {
				if (used[i])
					continue;
				sumScore = Utils.logSumExp(sumScore, score[i]);
			}
			double logp = Math.log(r.nextDouble() + 1e-60);
			double cur = Double.NEGATIVE_INFINITY;
			int ret = 0;
			for (; ret < score.length; ret++) {
				if (used[ret])
					continue;
				cur = Utils.logSumExp(cur, score[ret]);
				if (logp + sumScore - 1e-8 < cur)
					break;
			}
			Utils.Assert(ret < score.length && !used[ret]);
			//System.out.println(score.length + " " + ret);
			return ret;
		}

	}
	
	private boolean findOptChange(int[] goldlbids, DependencyInstance inst,
			SRLFeatureData sfd, int p, TIntArrayList args) {
		// p is the predicate id, 1 to #frames
		// find the local optimum match of p by swap
		
		boolean change = false;
		SemanticFrame predict = inst.frames[p];
		
		// re-construct the matching
		int N = args.size(), M = sfd.L;
		int T = N + M;
		int[] match = new int[T];
		int[] invMatch = new int[T];
		for (int i = 0; i < T; ++i)
			invMatch[i] = -1;
		
		int nullIndex = 0;
		for (int i = 0; i < N; ++i) {
			int label = predict.arglbids[args.get(i)];
			if (label >= 0) {
				Utils.Assert(invMatch[label] == -1);
				match[i] = label;
				invMatch[label] = i;
			}
			else {
				Utils.Assert(invMatch[M + nullIndex] == -1);
				match[i] = M + nullIndex;
				invMatch[M + nullIndex] = i;
				nullIndex++;
			}
		}
		
		nullIndex = 0;
		for (int i = N; i < T; ++i) {
			while (nullIndex < T && invMatch[nullIndex] != -1)
				nullIndex++;
			Utils.Assert(nullIndex < T);
			Utils.Assert(invMatch[nullIndex] == -1);
			match[i] = nullIndex;
			invMatch[nullIndex] = i;
			nullIndex++;
		}
		
		for (int i = 0; i < T; ++i)
			Utils.Assert(invMatch[i] != -1);
		
		// climb each argument
		for (int i = 0; i < T; ++i) {
			int oldMatch = match[i];
			//double fullBestScore = sfd.getScore(inst);
			double bestScore = sfd.getPartialScore(inst, p);
			int bestMatch = oldMatch;
			if (addLoss) {
				for (int j = 0; j < N; ++j) {
					bestScore += loss(goldlbids[args.get(j)], predict.arglbids[args.get(j)]);
					//fullBestScore += loss(goldFrame.arglbids[args.get(j)], predict.arglbids[args.get(j)]);
				}
			}
			
			for (int j = 0; j < T; ++j) {
				if (j == oldMatch)
					continue;
				Utils.Assert(match[i] == oldMatch);
				if (swap(predict, N, M, match, invMatch, i, j, args)) {
					double currScore = sfd.getPartialScore(inst, p);
					//double fullCurrScore = sfd.getScore(inst);
					if (addLoss) {
						for (int k = 0; k < N; ++k) {
							currScore += loss(goldlbids[args.get(k)], predict.arglbids[args.get(k)]);
							//fullCurrScore += loss(goldFrame.arglbids[args.get(k)], predict.arglbids[args.get(k)]);
						}
					}
					
					//Utils.Assert(Math.abs(currScore - bestScore - fullCurrScore + fullBestScore) < 1e-6);
					
					if (currScore > bestScore + 1e-6) {
						bestMatch = j;
						bestScore = currScore;
						//fullBestScore = fullCurrScore;
						change = true;
					}
					
					// recover
					swap(predict, N, M, match, invMatch, i, oldMatch, args);
				}
			}
			
			if (bestMatch != oldMatch)
				Utils.Assert(swap(predict, N, M, match, invMatch, i, bestMatch, args));
		}
		
		return change;
	}
	
	private boolean swap(SemanticFrame frame, int N, int M, int[] match, int[] invMatch, int p, int b, TIntArrayList args) {
		// p->a, q->b to p->b, q->a
		int a = match[p];
		int q = invMatch[b];
		
		if (p >= N && q >= N) {
			// both role labels are matched to null arguments
			return false;
		}
		else if (a >= M && b >= M) {
			// both arguments are matched to null labels
			return false;
		}
		
		if (p < N)
			frame.arglbids[args.get(p)] = b < M ? b : -1;
		if (q < N)
			frame.arglbids[args.get(q)] = a < M ? a : -1;
		
		match[p] = b; invMatch[b] = p;
		match[q] = a; invMatch[a] = q;
		
		return true;
	}
	
	private double loss(int gold, int pred) {
		if (gold != pred) {
			if (gold < 0) {
				// wrong arc
				return 1.0;
			}
			else if (pred < 0) {
				return 2.0;
			}
			else {
				return 0.5;
			}
		}
		else {
			return 0.0;
		}
	}
	
	private double calcScore(DependencyInstance now) 
	{
		double score = sfd.getScore(now);
        if (addLoss) 
        	score += sfd.getSRLCost(inst.frames, now.frames);

		return score;
	}
	
}
