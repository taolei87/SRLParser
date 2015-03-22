package parser.decoding;

import gnu.trove.list.array.TIntArrayList;
import parser.*;
import utils.Utils;

public class MaxMatchingDecoder extends SRLDecoder {

	/*
	 * Semantic role labeling decoder
	 */
	
	public MaxMatchingDecoder(Options options) {
		this.options = options;
		this.nullWeight = options.nullWeight;
	}
	
	@Override
	public DependencyInstance decode(DependencyInstance inst, SRLFeatureData sfd, boolean addLoss) {
		int len = inst.length;
		int F = inst.numframes;
		SemanticFrame[] frames = inst.frames;
		SemanticFrame[] predFrames = new SemanticFrame[F];
		for (int i = 0; i < F; ++i) {
			
			int p = frames[i].predid;
			
			TIntArrayList augs = new TIntArrayList(5);
			for (int j = 0; j < len; ++j)
				if (!sfd.isPruned(p, j))
					augs.add(j);
			
			predFrames[i] = findMaximumMatching(frames[i], sfd, p, augs, addLoss);
		}

		// may need to re-create an instance here
		inst.frames = predFrames;
		return inst;
	}
	
	public SemanticFrame findMaximumMatching(SemanticFrame frame,
			SRLFeatureData sfd, int p, TIntArrayList args, boolean addLoss)
	{
		//addLoss = false;
		int N = args.size(), M = sfd.L;
		int[] arglbids = new int[frame.arglbids.length];
		if (addLoss) {
			for (int i = 0, L = arglbids.length; i < L; ++i) {
				
				// this is gold link
				arglbids[i] = frame.arglbids[i];
				
				// "project" gold links to the prediction
				// remove this link if it is not valid on the predicted tree
				if (arglbids[i] >= 0 && sfd.isPruned(p, i))
					arglbids[i] = -1;
				
			}
		}
		
		
//		int[] arglbids = frame.arglbids;
//		int[] lb2arg = null;
//		if (addLoss) {
//			lb2arg = new int[M];
//			for (int i = 0; i < M; ++i) lb2arg[i] = -1;
//			for (int i = 0; i < arglbids.length; ++i)
//				if (arglbids[i] >= 0)
//					lb2arg[arglbids[i]] = i;
//		}
		
		int T = N+M;
		double minVal = Math.min(0.0, nullWeight);
		double[] f = new double[T*T];
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < M; ++j) {
				double va = sfd.getArcScore(p, args.get(i), j);
				if (addLoss) {
					// false positive
					if (arglbids[args.get(i)] < 0) va += 1.0;
					
					// incorrect label
					else if (j != arglbids[args.get(i)]) va += 0.5;
					//if (j != arglbids[args.get(i)]) va += 1.0;
				}
				f[i*T+j] = va;
				minVal = minVal > va ? va : minVal;
			}
		
		for (int i = N; i < T; ++i)
			for (int j = 0; j < T; ++j)
				f[i*T+j] = 0;
		
		for (int i = 0; i < N; ++i)
			for (int j = M; j < T; ++j) {		
				// false negative
				f[i*T+j] = nullWeight + ((addLoss && i < N && arglbids[args.get(i)] >= 0) ? 2.0 : 0.0);
			}
		
		if (minVal < 0.0)
			for (int i = 0; i < T; ++i)
				for (int j = 0; j < T; ++j)
					f[i*T+j] -= minVal;
		
		/*
		 * Max-weight matching (Kuhn-Munkres algorithm a.k.a. Edmonds-Karp algorithm)
		 */		
		MatchingInstance minst = new MatchingInstance(T, f);
		minst.run();
		
		SemanticFrame predict = new SemanticFrame(frame);
		//arglbids = new int[arglbids.length];
		predict.arglbids = arglbids;
		for (int i = 0; i < arglbids.length; ++i)
			arglbids[i] = -1;
		
//		for (int i = 0; i < N; ++i)
//			if (minst.find2[i] < M)
//				arglbids[args.get(i)] = minst.find2[i];
		
		for (int i = 0; i < M; ++i)
			if (minst.find[i] < N) {
				int a = args.get(minst.find[i]);
				arglbids[a] = i;
			}
		return predict;
	}

}

class MatchingInstance
{
	final static double eps = 1e-13;
	
	int T;
	double[] x, y, f;
	boolean[] visx, visy;
	int[] find, find2;
	
	boolean verbose = false;
	
	public MatchingInstance(int T, double[] f) 
	{
		this.T = T;
		this.f = f;
		x = new double[T];
		y = new double[T];
		visx = new boolean[T];
		visy = new boolean[T];
		find = new int[T];
		find2 = new int[T];
		for (int i = 0; i < T; ++i) {
			find[i] = -1;
			find2[i] = -1;
			x[i] = Double.NEGATIVE_INFINITY;
			y[i] = 0.0;
			for (int j = 0; j < T; ++j) {
				double v = f[i*T+j];
				//x[i] = x[i] < v ? v : x[i];
				if (v > x[i]) {
					x[i] = v;
					find2[i] = j;
				}
			}
		}
	}
	
	private boolean findPath(int u)
	{
		visx[u] = true;
		for (int v = 0; v < T; ++v) {
			Utils.Assert(x[u]+y[v]+eps >= f[u*T+v]);
			if (visy[v] == false && x[u]+y[v] <= f[u*T+v]+eps) {
				visy[v] = true;
				int w = find[v];
				find[v] = u;
				if (w == -1 || findPath(w)) return true;
				find[v] = w;
			}
		}
		return false;
	}
	
	/*
	 * Max-weight matching (Kuhn-Munkres algorithm a.k.a. Edmonds-Karp algorithm)
	 */		
	public void run()
	{

		for (int k = 0; k < T; ++k) {
			
			// reduce dual object sum x[i] + sum y[i] until an argument path is found.
			for (;;) {
//				if (verbose) {
//					System.out.print(score() + " ");
//				}
				for (int i = 0; i < T; ++i) visx[i] = false;
				for (int i = 0; i < T; ++i) visy[i] = false;
				
				// break if find an argument path (one more matching)
				if (findPath(k)) break;
				
				// reduce dual object sum x[i] + sum y[i]
				double minVal = Double.POSITIVE_INFINITY;
				for (int i = 0; i < T; ++i)
					if (visx[i])
						for (int j = 0; j < T; ++j)
							if (!visy[j]) {
								double va = x[i]+y[j]-f[i*T+j];
								minVal = minVal > va ? va : minVal;
							}
				Utils.Assert(minVal > 0.0 && minVal != Double.POSITIVE_INFINITY);
				for (int i = 0; i < T; ++i)
					if (visx[i]) x[i] -= minVal;
				for (int i = 0; i < T; ++i)
					if (visy[i]) y[i] += minVal;
			}
			
		}
		
//		double sum1 = 0.0, sum2 = 0.0;
//		for (int i = 0; i < T; ++i) {
//			Utils.Assert(find[i] >= 0);
//			sum1 += f[find[i]*T+i];
//			sum2 += x[i];
//			sum2 += y[i];
//		}
//		Utils.Assert(Math.abs(sum1-sum2)<=eps);
//		for (int i = 0; i < T; ++i)
//			for (int j = i+1; j < T; ++j)
//				Utils.Assert(find[i] != find[j]);
//		for (int i = 0; i < T; ++i)
//			for (int j = 0; j < T; ++j)
//				Utils.Assert(x[i]+y[j]+eps >= f[i*T+j]);
//		if (verbose) System.out.println();
	}
	
	public double score()
	{
		double sum = 0;
		for (int i = 0; i < T; ++i)
			sum += x[i] + y[i];
		return sum;
	}
}
