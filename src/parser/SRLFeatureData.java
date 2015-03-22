package parser;

import gnu.trove.list.array.TIntArrayList;
import parser.feature.SemanticFeatureFactory;
import utils.DictionarySet.DictionaryTypes;
import utils.FeatureVector;
import utils.Utils;

public class SRLFeatureData {
	
	public DependencyInstance inst;
	public DependencyPipe pipe;
	public Options options;
	public Parameters parameters;
	
	public int F, N, L;
	
	public int numPAs, numPARs;
	public FeatureVector[] pathFvs, wordFvs, parFvs, contextFvs;
	public int[] p2id;
	public double[] parScores;
	public double[][] wpU, wpV, ppW, cpX;
	boolean[] isPruned;				// whether a (h->m) arc is pruned
	
	public double gamma2 = 0.5;
	public int rank2 = 50;
	
	int sibnum;
	int[] sib2id;
	int coparnum;
	int[] copar2id;
	int gpnum;
	int[] gp2id;
	
	FeatureDataItem[] gpPred;			// grandparent-predicate-argument
	
	FeatureDataItem[] consSib;		// consecutive argument
	
	FeatureDataItem[] consPred;	// consecutive predicate

	public SRLFeatureData(DependencyInstance inst, Options options,
			DependencyPipe pipe, Parameters parameters)
	{
		this.inst = inst;
		this.options = options;
		this.parameters = parameters;
		this.pipe = pipe;
		
		F = inst.numframes;
		N = inst.length;
		L = pipe.dictionaries.size(DictionaryTypes.AUGLABEL);
		
		gamma2 = options.gamma2;
		rank2 = options.R2;

		// construct un-pruned arc list. All arcs are kept if there is no pruner.
		initArcPruningMap();
		
		// calculate 1st order feature vectors and scores
		initFirstOrderTable();
		
		// allocate memory for tables of high order features 
		initHighOrderFeatureTables();
	}
	
	public void initArcPruningMap() {
		sib2id = new int[F * N * N];
		copar2id = new int[F * F * N];
		gp2id = new int[F * F * N];
		isPruned = new boolean[F * N];
		
		for (int i = 0; i < isPruned.length; ++i) {
			isPruned[i] = true;
		}
		
		for (int i = 0; i < sib2id.length; ++i)
			sib2id[i] = -1;

		for (int i = 0; i < copar2id.length; ++i)
			copar2id[i] = -1;

		for (int i = 0; i < gp2id.length; ++i)
			gp2id[i] = -1;

		// hard constraint to prune arc
		sibnum = 0;
		coparnum = 0;
		gpnum = 0;
		for (int i = 0; i < F; ++i) {
			int pid = inst.frames[i].predid;
			int[] augs = inst.frames[i].arglbids;
			for (int j = 0; j < N; ++j) {
				if (isValidPredAugPair(pid, j)
						|| augs[j] >= 0) {
					// not pruned or gold arcs
					isPruned[i * N + j] = !isValidPredAugPair(pid, j);
					
					if (options.useSRL2O) {
						// build sib
						for (int k = j + 1; k < N; ++k) 
							if (isValidPredAugPair(pid, k)
								|| augs[k] >= 0) {
								sib2id[(i * N + j) * N + k] = sibnum;
								sibnum++;
							}
						
						// build co-par
						for (int k = i + 1; k < F; ++k) 
							if (SemanticFeatureFactory.isValidPredAugPair(inst, inst.frames[k].predid, j)
									|| inst.frames[k].arglbids[j] >= 0) {
									copar2id[(i * F + k) * N + j] = coparnum;
									coparnum++;
								}
						
						// build gp
						//for (int k = 0; k < F; ++k)
						//	if (SemanticFeatureFactory.isValidPredAugPair(inst, inst.frames[k].predid, pid)
						//			|| inst.frames[k].arglbids[pid] >= 0) {
						//			gp2id[(k * F + i) * N + j] = gpnum;
						//			gpnum++;
						//		}
					}
				}
			}
		}
	}
	
	private void initFirstOrderTable()
	{		
		numPAs = F * N;
		numPARs = F * N * L;
		p2id = new int[N];
		wordFvs = new FeatureVector[N];
		wpU = new double[N][rank2];
		wpV = new double[N][rank2];
		ppW = new double[numPARs][];
		cpX = new double[numPARs][];
		contextFvs = new FeatureVector[numPARs];
		pathFvs = new FeatureVector[numPARs];
		parFvs = new FeatureVector[numPARs];
		parScores = new double[numPARs];
		
		for (int i = 0; i < N; ++i)
			p2id[i] = -1;
		
		for (int i = 0; i < F; ++i) {
			int p = inst.frames[i].predid;
			Utils.Assert(p2id[p] == -1);
			p2id[p] = i;
			Utils.Assert(p2id[p] == inst.predIndex[p]);
			
			if (wordFvs[p] == null) {
				wordFvs[p] = pipe.smnFactory.createWordFeatureVector(inst, p);
				parameters.projectU2(wordFvs[p], wpU[p]);
				parameters.projectV2(wordFvs[p], wpV[p]);
			}
			
			for (int a = 0; a < N; ++a)
				if (!isPruned(p, a)) {
					
					if (wordFvs[a] == null) {
						wordFvs[a] = pipe.smnFactory.createWordFeatureVector(inst, a);
						parameters.projectU2(wordFvs[a], wpU[a]);
						parameters.projectV2(wordFvs[a], wpV[a]);
					}
					
					for (int r = 0; r < L; ++r) {
						int id = i*N*L + a*L + r;
						contextFvs[id] = pipe.smnFactory.createContextFeatureVector(inst, p, a, r);
						cpX[id] = new double[rank2];
						parameters.projectX2(contextFvs[id], cpX[id]);
						pathFvs[id] = pipe.smnFactory.createPathFeatureVector(inst, p, a, r);
						parFvs[id] = pipe.smnFactory.createPredArgLinkFeatures(inst, p, a, r);
						ppW[id] = new double[rank2];
						parameters.projectW2(pathFvs[id], ppW[id]);
						parScores[id] = gamma2 * parameters.dotProduct2(parFvs[id]) +
										(1-gamma2) * parameters.dotProduct2(wpU[p], wpV[a], ppW[id], cpX[id]);
					}
				}
		}
	}

	public void initHighOrderFeatureTables() {
		// init non-first-order feature tables
		
		if (options.useSRL2O) {
			//gpPred = new FeatureDataItem[gpnum * L * L];		// with the label of p->a
			
			consSib = new FeatureDataItem[sibnum * L * L];
			
			consPred = new FeatureDataItem[coparnum * L * L];
		}
		
	}

	public boolean isValidPredAugPair(int p, int a)
	{
		return SemanticFeatureFactory.isValidPredAugPair(inst, p, a);
	}
	
	public boolean isPruned(int p, int a) 
	{
		// p = 1..N
		return isPruned[p2id[p] * N + a];
	}
	
	public double getArcScore(int p, int a, int r)
	{
		// p = 1..N
		int i = p2id[p];
		int id = i * N * L + a * L + r;
		Utils.Assert(i != -1 && parFvs[id] != null);
		return parScores[id];
	}

	public FeatureVector getArcFeatureVector(int p, int a, int r)
	{
		// p = 1..N
		int i = p2id[p];
		int id = i * N * L + a * L + r;
		Utils.Assert(i != -1 && parFvs[id] != null);
		return parFvs[id];
	}

	public double getSibScore(int p, int a1, int a2, int r1, int r2) 
	{
		// a1 < a2, p = 1..N
		// maybe r2 will be included
		int id = sib2id[(p2id[p] * N + a1) * N + a2];
		
		Utils.Assert(id >= 0);
		
		int pos = (id * L + r1) * L + r2;
		if (consSib[pos] == null)
			getSibFeatureVector(p, a1, a2, r1, r2);
		
		return consSib[pos].score;
	}
	
	public FeatureVector getSibFeatureVector(int p, int a1, int a2, int r1, int r2)
	{
		int id = sib2id[(p2id[p] * N + a1) * N + a2];
		
		Utils.Assert(id >= 0);
		
		int pos = (id * L + r1) * L + r2;
		FeatureDataItem item = consSib[pos];
		if (item == null) {
			FeatureVector fv = pipe.smnFactory.addSibFeatures(inst, p, a1, a2, r1, r2);
			double score = parameters.dotProduct2(fv) * gamma2;
			item = new FeatureDataItem(fv, score);
			consSib[pos] = item;			
		}
		return item.fv;
	}
	
	public double getGPScore(int gp,  int p, int a, int r1, int r2) 
	{
		// p, gp are in the index 1..N
		// maybe r2 will be included
		int id = gp2id[(p2id[gp] * F + p2id[p]) * N + a];
		
		Utils.Assert(id >= 0);
		
		int pos = (id * L + r1) * L + r2;
		if (gpPred[pos] == null)
			getGPFeatureVector(gp, p, a, r1, r2);
		
		return gpPred[pos].score;
	}
	
	public FeatureVector getGPFeatureVector(int gp, int p, int a, int r1, int r2)
	{
		int id = gp2id[(p2id[gp] * F + p2id[p]) * N + a];
		
		Utils.Assert(id >= 0);
		
		int pos = (id * L + r1) * L + r2;
		FeatureDataItem item = gpPred[pos];
		if (item == null) {
			FeatureVector fv = pipe.smnFactory.addGPFeatures(inst, gp, p, a, r1, r2);
			double score = parameters.dotProduct2(fv) * gamma2;
			item = new FeatureDataItem(fv, score);
			gpPred[pos] = item;			
		}
		return item.fv;
	}
	
	public double getCoPredScore(int p, int p2, int a, int r1, int r2) 
	{
		// p, p2 are in the index 1..N
		// maybe r2 will be included
		int id = copar2id[(p2id[p] * F + p2id[p2]) * N + a];
		
		Utils.Assert(id >= 0);
		
		int pos = (id * L + r1) * L + r2;
		if (consPred[pos] == null)
			getCoPredFeatureVector(p, p2, a, r1, r2);
		
		return consPred[pos].score;
	}
	
	public FeatureVector getCoPredFeatureVector(int p, int p2, int a, int r1, int r2)
	{
		int id = copar2id[(p2id[p] * F + p2id[p2]) * N + a];
		
		Utils.Assert(id >= 0);
		
		int pos = (id * L + r1) * L + r2;
		FeatureDataItem item = consPred[pos];
		if (item == null) {
			FeatureVector fv = pipe.smnFactory.addCoPredFeatures(inst, p, p2, a, r1, r2);
			double score = parameters.dotProduct2(fv) * gamma2;
			item = new FeatureDataItem(fv, score);
			consPred[pos] = item;			
		}
		return item.fv;
	}
	
	public double getScore(DependencyInstance inst)
	{
		// get score of all predicates
		double score = 0.0;
		SemanticFrame[] frames = inst.frames;
		
		for (int p = 0; p < F; ++p) {
			SemanticFrame frame = frames[p];
			int len = frame.arglbids.length;
			
			TIntArrayList childList = new TIntArrayList();
			for (int c = 0; c < len; ++c)
				if (!isPruned[p * len + c] && frame.arglbids[c] != -1) {
					//Utils.Assert(!isPruned(frame.predid, c));
					childList.add(c);
				}
			
			TIntArrayList gpList = new TIntArrayList();
			for (int gp = 0, L = frames.length; gp < L; ++gp)
				if (!isPruned[gp * len + frame.predid] && frames[gp].arglbids[frame.predid] != -1)
					gpList.add(gp);

			for (int i = 0, L = childList.size(); i < L; ++i) {
				int r = frame.arglbids[childList.get(i)];
				int c = childList.get(i);
				
				// first order
				score += getArcScore(frame.predid, c, r);

				if (options.useSRL2O) {
					// consecutive sib
					if (i < L - 1) {
						score += getSibScore(frame.predid, c, childList.get(i + 1), r, frame.arglbids[childList.get(i + 1)]);
					}
					
					// grandparent, head of p
					//for (int j = 0, L2 = gpList.size(); j < L2; ++j) {
					//	score += getGPScore(frames[gpList.get(j)].predid, frame.predid, c, frames[gpList.get(j)].arglbids[frame.predid], r);
					//}
					
					// co-predicate, next
					if (p < frames.length - 1) {
						SemanticFrame nextFrame = frames[p + 1];
						if (!isPruned[(p + 1) * len + c] && nextFrame.arglbids[c] != -1) {
							score += getCoPredScore(frame.predid, nextFrame.predid, c, r, nextFrame.arglbids[c]);
						}
					}
				}
			}

			if (options.useSRLHO && inst.voice[frame.predid] < 2 && inst.voice[frame.predid] >= 0) {
				FeatureVector fv = pipe.smnFactory.createGlobalFeatures(inst, isPruned, p);
				score += parameters.dotProduct2(fv) * gamma2;
			}
		}
		
		return score;
	}
	
	
	public double getPartialScore(DependencyInstance inst, int p) {
		// get partial score when climbing predicates p
		double score = 0.0;
		
		SemanticFrame[] frames = inst.frames;
		//int[] predIndex = inst.predIndex;
		
		SemanticFrame frame = frames[p];
		int len = frame.arglbids.length;
		
		TIntArrayList childList = new TIntArrayList();
		for (int c = 0; c < len; ++c)
			if (!isPruned[p * len + c] && frame.arglbids[c] != -1) {
				//Utils.Assert(!isPruned(frame.predid, c));
				childList.add(c);
			}
		
		TIntArrayList gpList = new TIntArrayList();
		for (int gp = 0, L = frames.length; gp < L; ++gp)
			if (!isPruned[gp * len + frame.predid] && frames[gp].arglbids[frame.predid] != -1)
				gpList.add(gp);

		for (int i = 0, L = childList.size(); i < L; ++i) {
			int r = frame.arglbids[childList.get(i)];
			int c = childList.get(i);
			
			// first order
			score += getArcScore(frame.predid, c, r);

			if (options.useSRL2O) {
				// consecutive sib
				
				if (i < L - 1) {
					score += getSibScore(frame.predid, c, childList.get(i + 1), r, frame.arglbids[childList.get(i + 1)]);
				}
				
				// grandparent, head of p
				//for (int j = 0, L2 = gpList.size(); j < L2; ++j) {
				//	score += getGPScore(frames[gpList.get(j)].predid, frame.predid, c, frames[gpList.get(j)].arglbids[frame.predid], r);
				//}
				
				// grandparent, child of c
				//if (predIndex[c] >= 0) {
				//	SemanticFrame childFrame = frames[predIndex[c]];
				//	for (int j = 0; j < len; ++j)
				//		if (childFrame.arglbids[j] >= 0) {
				//			Utils.Assert(!isPruned(predIndex[c], j));
				//			score += getGPScore(frame.predid, c, j, r, childFrame.arglbids[j]);
				//		}
				//}
				
				// co-predicate, prev
				if (p > 0) {
					SemanticFrame prevFrame = frames[p - 1];
					if (!isPruned[(p - 1) * len + c] && prevFrame.arglbids[c] != -1) {
						score += getCoPredScore(prevFrame.predid, frame.predid, c, prevFrame.arglbids[c], r);
					}
				}

				// co-predicate, next
				if (p < frames.length - 1) {
					SemanticFrame nextFrame = frames[p + 1];
					if (!isPruned[(p + 1) * len + c] && nextFrame.arglbids[c] != -1) {
						score += getCoPredScore(frame.predid, nextFrame.predid, c, r, nextFrame.arglbids[c]);
					}
				}
				
			}
		}
		
		if (options.useSRLHO && inst.voice[frame.predid] < 2 && inst.voice[frame.predid] >= 0) {
			FeatureVector fv = pipe.smnFactory.createGlobalFeatures(inst, isPruned, p);
			score += parameters.dotProduct2(fv) * gamma2;
		}
		
		return score;
	}

	public FeatureVector getFeatureVector(DependencyInstance inst)
	{
		FeatureVector fv = new FeatureVector(pipe.smnFactory.numLinkFeats);
		
		SemanticFrame[] frames = inst.frames;
		int numframes = inst.numframes;
		int len = inst.length;
		
		Utils.Assert(numframes == F && len == N);

		for (int p = 0; p < numframes; ++p) {
			SemanticFrame frame = frames[p];
			
			TIntArrayList childList = new TIntArrayList();
			for (int c = 0; c < len; ++c)
				if (!isPruned[p * len + c] && frame.arglbids[c] != -1) {
					childList.add(c);
				}
			
			TIntArrayList gpList = new TIntArrayList();
			for (int gp = 0, L = frames.length; gp < L; ++gp)
				if (!isPruned[gp * len + frame.predid] && frames[gp].arglbids[frame.predid] != -1)
					gpList.add(gp);

			for (int i = 0, L = childList.size(); i < L; ++i) {
				int r = frame.arglbids[childList.get(i)];
				int c = childList.get(i);
				
				// first order
				fv.addEntries(getArcFeatureVector(frame.predid, c, r));

				if (options.useSRL2O) {
					// consecutive sib
					if (i < L - 1) {
						fv.addEntries(getSibFeatureVector(frame.predid, c, childList.get(i + 1), r, frame.arglbids[childList.get(i + 1)]));
					}
					
					// grandparent, head of p
					//for (int j = 0, L2 = gpList.size(); j < L2; ++j) {
					//	fv.addEntries(getGPFeatureVector(frames[gpList.get(j)].predid, frame.predid, c, frames[gpList.get(j)].arglbids[frame.predid], r));
					//}
					
					// co-predicate, next
					if (p < frames.length - 1) {
						SemanticFrame nextFrame = frames[p + 1];
						if (!isPruned[(p + 1) * len + c] && nextFrame.arglbids[c] != -1) {
							fv.addEntries(getCoPredFeatureVector(frame.predid, nextFrame.predid, c, r, nextFrame.arglbids[c]));
						}
					}
				}
			}
			
			if (options.useSRLHO && inst.voice[frame.predid] < 2 && inst.voice[frame.predid] >= 0) {
				fv.addEntries(pipe.smnFactory.createGlobalFeatures(inst, isPruned, p));
			}
		}

		return fv;
	}

	public FeatureVector getFeatureDifference(DependencyInstance gold,
				DependencyInstance pred) {
		
		FeatureVector fv = new FeatureVector(pipe.smnFactory.numLinkFeats);

		Utils.Assert(pred.heads == inst.heads);
		Utils.Assert(pred.deplbids == inst.deplbids);
		Utils.Assert(F == gold.numframes && N == gold.length);
		
		/*
		for (int i = 0; i < F; ++i) {
			SemanticFrame frame = gold.frames[i];
			int p = frame.predid;
			for (int a = 0; a < N; ++a) {
				int r = frame.arglbids[a];
				if (r >= 0 && isValidPredAugPair(p, a)) {
					int id = i*N*L + a*L + r;
					fv.addEntries(parFvs[id], 1.0);
				}
			}
		}
		
		for (int i = 0; i < F; ++i) {
			SemanticFrame frame = pred.frames[i];
			int p = frame.predid;
			for (int a = 0; a < N; ++a) {
				int r = frame.arglbids[a];
				if (r >= 0 && isValidPredAugPair(p, a)) {
					int id = i*N*L + a*L + r;
					fv.addEntries(parFvs[id], -1.0);
				}
			}
		}
		*/
		
		SemanticFrame[] temp = pred.frames;
		
		pred.frames = gold.frames;
        fv.addEntries(getFeatureVector(pred), 1.0);

        pred.frames = temp;
        fv.addEntries(getFeatureVector(pred), -1.0);

        return fv;
	}
	
//	public FeatureVector getFeatureDifference(DependencyInstance gold,
//			DependencyInstance pred) {
//		FeatureVector fv = new FeatureVector(pipe.smnFactory.numLinkFeats);
//		
//		Utils.Assert(pred == inst);
//		
//		SemanticFrame[] temp = inst.frames;
//		
//		inst.frames = gold.frames;
//		fv.addEntries(pipe.smnFactory.getFeatureVector(inst), 1.0);
//		
//		inst.frames = temp;
//		fv.addEntries(pipe.smnFactory.getFeatureVector(inst), -1.0);
//		
//		return fv;
//	}
	
	public double getSRLCost(SemanticFrame[] gold, SemanticFrame[] pred)
	{
		Utils.Assert(gold.length == pred.length);
		double dis = 0;
		for (int i = 0, N = gold.length; i < N; ++i) {
			
			Utils.Assert(gold[i].predid == pred[i].predid);
			
			int[] ga = gold[i].arglbids, pa = pred[i].arglbids;
			for (int j = 0, L = ga.length; j < L; ++j) {
				if (isPruned(gold[i].predid, j))
					continue;
				
				// "project" gold links to the prediction
				// remove this link if it is not valid on the predicted tree
				int garg = ga[j];
				if (garg != pa[j]) {
					
					// false positive or false negative
					if (garg < 0) dis += 1.0;
					else if (pa[j] < 0) dis += 2.0;
					// incorrect label
					else dis += 0.5;
				}
			}
		}
		return dis;
	}
	
}
