package parser;

import java.io.Serializable;

import parser.feature.SemanticFeatureFactory;
import utils.FeatureVector;
import utils.Utils;

public class Parameters implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	public static final int d = 7;
	
	public transient Options options;
	public final int labelLossType;
	public double synC, smnC, gamma, gammaLabel;
	public int size, sizeL;
	public int rank;
	public int N, M, T, D;
		
	public double[] params, paramsL, backupL, totalL;
	public double[][] U, V, W;
	public transient double[] backup, total;
	public transient double[][] totalU, totalV, totalW;
	public transient double[][] backupU, backupV, backupW;	
	public transient FeatureVector[] dU, dV, dW;
	
	public int rank2;
	public double gamma2;
	public int size2, N2, M2, D2, L2;
	public double[] params2, sqrsum2;
	public transient double[] backup2, total2;
	public double[][] U2, V2, W2, X2;
	public transient double[][] totalU2, totalV2, totalW2, totalX2;
	public transient double[][] backupU2, backupV2, backupW2, backupX2;
	public transient FeatureVector[] dU2, dV2, dW2, dX2;
	
	public Parameters(DependencyPipe pipe, Options options) 
	{
		 //T = pipe.types.length;
        D = d * 2 + 1;
		size = pipe.synFactory.numArcFeats;		
		params = new double[size];
		total = new double[size];
		
		if (options.learnLabel) {
			sizeL = pipe.synFactory.numLabeledArcFeats;
			paramsL = new double[sizeL];
			totalL = new double[sizeL];
		}
		
		this.options = options;
		this.labelLossType = options.labelLossType;
		synC = options.synC;
		smnC = options.smnC;
		gamma = options.gamma;
		gammaLabel = options.gammaLabel;
		rank = options.R;
		
		N = pipe.synFactory.numWordFeats;
		M = N;
		U = new double[rank][N];		
		V = new double[rank][M];
		W = new double[rank][D];
		totalU = new double[rank][N];
		totalV = new double[rank][M];
		totalW = new double[rank][D];
		dU = new FeatureVector[rank];
		dV = new FeatureVector[rank];
		dW = new FeatureVector[rank];
		
		rank2 = options.R2;
		gamma2 = options.gamma2;
		L2 = pipe.smnFactory.numContextFeats;
		D2 = pipe.smnFactory.numPathFeats;
		N2 = pipe.smnFactory.numWordFeats;
		M2 = N2;
		size2 = pipe.smnFactory.numLinkFeats;
		params2 = new double[size2];
		total2 = new double[size2];
		U2 = new double[rank2][N2];		
		V2 = new double[rank2][M2];
		W2 = new double[rank2][D2];
		X2 = new double[rank2][L2];
		totalU2 = new double[rank2][N2];
		totalV2 = new double[rank2][M2];
		totalW2 = new double[rank2][D2];
		totalX2 = new double[rank2][L2];
		dU2 = new FeatureVector[rank2];
		dV2 = new FeatureVector[rank2];
		dW2 = new FeatureVector[rank2];
		dX2 = new FeatureVector[rank2];
		sqrsum2 = new double[size2];

	}
	
	public void randomlyInitUVW() 
	{
		for (int i = 0; i < rank; ++i) {
			U[i] = Utils.getRandomUnitVector(N);
			V[i] = Utils.getRandomUnitVector(M);
			W[i] = Utils.getRandomUnitVector(D);
			totalU[i] = U[i].clone();
			totalV[i] = V[i].clone();
			totalW[i] = W[i].clone();
		}
	}
	
	public void randomlyInitUVWX2() 
	{
		for (int i = 0; i < rank2; ++i) {
			U2[i] = Utils.getRandomUnitVector(N2);
			V2[i] = Utils.getRandomUnitVector(M2);
			W2[i] = Utils.getRandomUnitVector(D2);
			X2[i] = Utils.getRandomUnitVector(L2);
			totalU2[i] = U2[i].clone();
			totalV2[i] = V2[i].clone();
			totalW2[i] = W2[i].clone();
			totalX2[i] = X2[i].clone();
		}
	}
	
	public void averageParametersSyn(int T) 
	{
		backup = params;
		double[] avgParams = new double[size];
		for (int i = 0; i < size; ++i) {
			avgParams[i] = (params[i] * (T+1) - total[i])/T;			
		}		
		params = avgParams;
		
		backupL = paramsL;
		double[] avgParamsL = new double[sizeL];
		for (int i = 0; i < sizeL; ++i) {
			avgParamsL[i] = (paramsL[i] * (T+1) - totalL[i])/T;			
		}		
		paramsL = avgParamsL;


		backupU = U;
		double[][] avgU = new double[rank][N];
		for (int i = 0; i < rank; ++i)
			for (int j = 0; j < N; ++j) {
				avgU[i][j] = (U[i][j] * (T+1) - totalU[i][j])/T;
			}
		U = avgU;
		
		backupV = V;
		double[][] avgV = new double[rank][M];
		for (int i = 0; i < rank; ++i)
			for (int j = 0; j < M; ++j) {
				avgV[i][j] = (V[i][j] * (T+1) - totalV[i][j])/T;
			}
		V = avgV;
		
		backupW = W;
		double[][] avgW = new double[rank][D];
		for (int i = 0; i < rank; ++i)
			for (int j = 0; j < D; ++j) {
				avgW[i][j] = (W[i][j] * (T+1) - totalW[i][j])/T;
			}
		W = avgW;
	}
	
	public void unaverageParametersSyn() 
	{
		params = backup;
		paramsL = backupL;
		U = backupU;
		V = backupV;
		W = backupW;
	}
	
	public void averageParameters(int T) 
	{
		backup = params;
		double[] avgParams = new double[size];
		for (int i = 0; i < size; ++i) {
			avgParams[i] = (params[i] * (T+1) - total[i])/T;			
		}		
		params = avgParams;
		
		backupL = paramsL;
		double[] avgParamsL = new double[sizeL];
		for (int i = 0; i < sizeL; ++i) {
			avgParamsL[i] = (paramsL[i] * (T+1) - totalL[i])/T;			
		}		
		paramsL = avgParamsL;
		
		backupU = U;
		double[][] avgU = new double[rank][N];
		for (int i = 0; i < rank; ++i)
			for (int j = 0; j < N; ++j) {
				avgU[i][j] = (U[i][j] * (T+1) - totalU[i][j])/T;
			}
		U = avgU;
		
		backupV = V;
		double[][] avgV = new double[rank][M];
		for (int i = 0; i < rank; ++i)
			for (int j = 0; j < M; ++j) {
				avgV[i][j] = (V[i][j] * (T+1) - totalV[i][j])/T;
			}
		V = avgV;
		
		backupW = W;
		double[][] avgW = new double[rank][D];
		for (int i = 0; i < rank; ++i)
			for (int j = 0; j < D; ++j) {
				avgW[i][j] = (W[i][j] * (T+1) - totalW[i][j])/T;
			}
		W = avgW;
		
		backup2 = params2;
		double[] avgParams2 = new double[size2];
		for (int i = 0; i < size2; ++i) {
			avgParams2[i] = (params2[i] * (T+1) - total2[i])/T;			
		}		
		params2 = avgParams2;
		
		backupU2 = U2;
		double[][] avgU2 = new double[rank2][N2];
		for (int i = 0; i < rank2; ++i)
			for (int j = 0; j < N2; ++j) {
				avgU2[i][j] = (U2[i][j] * (T+1) - totalU2[i][j])/T;
			}
		U2 = avgU2;
		
		backupV2 = V2;
		double[][] avgV2 = new double[rank2][M2];
		for (int i = 0; i < rank2; ++i)
			for (int j = 0; j < M2; ++j) {
				avgV2[i][j] = (V2[i][j] * (T+1) - totalV2[i][j])/T;
			}
		V2 = avgV2;
		
		backupW2 = W2;
		double[][] avgW2 = new double[rank2][D2];
		for (int i = 0; i < rank2; ++i)
			for (int j = 0; j < D2; ++j) {
				avgW2[i][j] = (W2[i][j] * (T+1) - totalW2[i][j])/T;
			}
		W2 = avgW2;
		
		backupX2 = X2;
		double[][] avgX2 = new double[rank2][L2];
		for (int i = 0; i < rank2; ++i)
			for (int j = 0; j < L2; ++j) {
				avgX2[i][j] = (X2[i][j] * (T+1) - totalX2[i][j])/T;
			}
		X2 = avgX2;
	}
	
	public void unaverageParameters() 
	{
		params = backup;
		paramsL = backupL;
		U = backupU;
		V = backupV;
		W = backupW;
		params2 = backup2;
		U2 = backupU2;
		V2 = backupV2;
		W2 = backupW2;
		X2 = backupX2;
	}
	
	public void clearUVW() 
	{
		U = new double[rank][N];
		V = new double[rank][M];
		W = new double[rank][D];
		totalU = new double[rank][N];
		totalV = new double[rank][M];
		totalW = new double[rank][D];
	}
	
	public void clearTheta() 
	{
		params = new double[size];
		total = new double[size];
		params2 = new double[size2];
		total2 = new double[size2];
	}
	
	public void printUStat() 
	{
		double sum = 0;
		double min = Double.POSITIVE_INFINITY, max = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < rank; ++i) {
			sum += Utils.squaredSum(U[i]);
			min = Math.min(min, Utils.min(U[i]));
			max = Math.max(max, Utils.max(U[i]));
		}
		System.out.printf(" |U|^2: %f min: %f\tmax: %f%n", sum, min, max);
	}
	
	public void printVStat() 
	{
		double sum = 0;
		double min = Double.POSITIVE_INFINITY, max = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < rank; ++i) {
			sum += Utils.squaredSum(V[i]);
			min = Math.min(min, Utils.min(V[i]));
			max = Math.max(max, Utils.max(V[i]));
		}
		System.out.printf(" |V|^2: %f min: %f\tmax: %f%n", sum, min, max);
	}
	
	public void printWStat() 
	{
		double sum = 0;
		double min = Double.POSITIVE_INFINITY, max = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < rank; ++i) {
			sum += Utils.squaredSum(W[i]);
			min = Math.min(min, Utils.min(W[i]));
			max = Math.max(max, Utils.max(W[i]));
		}
		System.out.printf(" |W|^2: %f min: %f\tmax: %f%n", sum, min, max);
	}
	
	public void printU2Stat() 
	{
		double sum = 0;
		double min = Double.POSITIVE_INFINITY, max = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < rank2; ++i) {
			sum += Utils.squaredSum(U2[i]);
			min = Math.min(min, Utils.min(U2[i]));
			max = Math.max(max, Utils.max(U2[i]));
		}
		System.out.printf(" |U|^2: %f min: %f\tmax: %f%n", sum, min, max);
	}
	
	public void printV2Stat() 
	{
		double sum = 0;
		double min = Double.POSITIVE_INFINITY, max = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < rank2; ++i) {
			sum += Utils.squaredSum(V2[i]);
			min = Math.min(min, Utils.min(V2[i]));
			max = Math.max(max, Utils.max(V2[i]));
		}
		System.out.printf(" |V|^2: %f min: %f\tmax: %f%n", sum, min, max);
	}
	
	public void printW2Stat() 
	{
		double sum = 0;
		double min = Double.POSITIVE_INFINITY, max = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < rank2; ++i) {
			sum += Utils.squaredSum(W2[i]);
			min = Math.min(min, Utils.min(W2[i]));
			max = Math.max(max, Utils.max(W2[i]));
		}
		System.out.printf(" |W|^2: %f min: %f\tmax: %f%n", sum, min, max);
	}
	
	public void printX2Stat() 
	{
		double sum = 0;
		double min = Double.POSITIVE_INFINITY, max = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < rank2; ++i) {
			sum += Utils.squaredSum(X2[i]);
			min = Math.min(min, Utils.min(X2[i]));
			max = Math.max(max, Utils.max(X2[i]));
		}
		System.out.printf(" |X|^2: %f min: %f\tmax: %f%n", sum, min, max);
	}
	
	public void printThetaStat() 
	{
		double sum = Utils.squaredSum(params);
		double min = Utils.min(params);
		double max = Utils.max(params);		
		System.out.printf(" |\u03b8_d|^2: %f min: %f\tmax: %f%n", sum, min, max);
		
		sum = Utils.squaredSum(params2);
		min = Utils.min(params2);
		max = Utils.max(params2);
		System.out.printf(" |\u03b8_l|^2: %f min: %f\tmax: %f%n", sum, min, max);
	}
	
	public void projectU(FeatureVector fv, double[] proj) 
	{
		for (int r = 0; r < rank; ++r) 
			proj[r] = fv.dotProduct(U[r]);
	}
	
	public void projectV(FeatureVector fv, double[] proj) 
	{
		for (int r = 0; r < rank; ++r) 
			proj[r] = fv.dotProduct(V[r]);
	}
	
	public double dotProduct(FeatureVector fv)
	{
		return fv.dotProduct(params);
	}
	
	public double dotProductL(FeatureVector fv)
	{
		return fv.dotProduct(paramsL);
	}

	public double dotProduct(double[] proju, double[] projv, int dist)
	{
		double sum = 0;
		int binDist = getBinnedDistance(dist);
		for (int r = 0; r < rank; ++r)
			sum += proju[r] * projv[r] * (W[r][binDist] + W[r][0]);
		return sum;
	}
	
	public void projectU2(FeatureVector fv, double[] proj) 
	{
		for (int r = 0; r < rank2; ++r) 
			proj[r] = fv.dotProduct(U2[r]);
	}
	
	public void projectV2(FeatureVector fv, double[] proj) 
	{
		for (int r = 0; r < rank2; ++r) 
			proj[r] = fv.dotProduct(V2[r]);
	}
	
	public void projectW2(FeatureVector fv, double[] proj) 
	{
		for (int r = 0; r < rank2; ++r) 
			proj[r] = fv.dotProduct(W2[r]);
	}
	
	public void projectX2(FeatureVector fv, double[] proj) 
	{
		for (int r = 0; r < rank2; ++r) 
			proj[r] = fv.dotProduct(X2[r]);
	}
	
	public double dotProduct2(FeatureVector fv)
	{
		return fv.dotProduct(params2);
	}
	
	public double dotProduct2(double[] proju, double[] projv, double[] projw, double[] projx)
	{
		double sum = 0;
		for (int r = 0; r < rank2; ++r)
			sum += proju[r]*projv[r]*projw[r]*projx[r];
		return sum;
	}
	
	public double updateSyn(DependencyInstance gold, DependencyInstance pred,
			LocalFeatureData lfd, GlobalFeatureData gfd,
			int updCnt, int offset)
	{
		int N = gold.length;
    	int[] actDeps = gold.heads;
    	int[] actLabs = gold.deplbids;
    	int[] predDeps = pred.heads;
    	int[] predLabs = pred.deplbids;
    	
    	double Fi = getHammingDis(actDeps, actLabs, predDeps, predLabs);
    	if (Fi == 0) return 0.0;
    	
    	FeatureVector dt = lfd.getFeatureDifference(gold, pred);
    	dt.addEntries(gfd.getFeatureDifference(gold, pred));
    	
        double loss = - dt.dotProduct(params)*gamma +Fi;

        double l2norm = dt.Squaredl2NormUnsafe() * gamma * gamma;        			 
    	
        //int updId = (updCnt + offset) % 3;
        //if ( updId == 1 ) 
        {
        	// update U
        	for (int k = 0; k < rank; ++k) {        		
        		FeatureVector dUk = getdU(k, lfd, actDeps, predDeps);
            	l2norm += dUk.Squaredl2NormUnsafe() * (1-gamma) * (1-gamma);            	
            	loss -= dUk.dotProduct(U[k]) * (1-gamma);
            	dU[k] = dUk;
        	}
        }
        //else if ( updId == 2 ) 
        {
        	// update V
        	for (int k = 0; k < rank; ++k) {
        		FeatureVector dVk = getdV(k, lfd, actDeps, predDeps);
            	l2norm += dVk.Squaredl2NormUnsafe() * (1-gamma) * (1-gamma);
            	//loss -= dVk.dotProduct(V[k]) * (1-gamma);
            	dV[k] = dVk;
        	}        	
        } 
        //else 
        {
        	// update W
        	for (int k = 0; k < rank; ++k) {
        		FeatureVector dWk = getdW(k, lfd, actDeps, predDeps);
            	l2norm += dWk.Squaredl2NormUnsafe() * (1-gamma) * (1-gamma);
            	//loss -= dWk.dotProduct(W[k]) * (1-gamma);
            	dW[k] = dWk;
        	}   
        }
        
        double alpha = loss/l2norm; 
    	alpha = Math.min(synC, alpha);
    	if (alpha > 0) {
    		
    		{
    			// update theta
	    		double coeff = alpha * gamma, coeff2 = coeff * updCnt;
	    		for (int i = 0, K = dt.size(); i < K; ++i) {
		    		int x = dt.x(i);
		    		double z = dt.value(i);
		    		params[x] += coeff * z;
		    		total[x] += coeff2 * z;
	    		}
    		}
    		
    		//if ( updId == 1 ) 
    		{
    			// update U
    			double coeff = alpha * (1-gamma), coeff2 = coeff * updCnt;
            	for (int k = 0; k < rank; ++k) {
            		FeatureVector dUk = dU[k];
            		for (int i = 0, K = dUk.size(); i < K; ++i) {
            			int x = dUk.x(i);
            			double z = dUk.value(i);
            			U[k][x] += coeff * z;
            			totalU[k][x] += coeff2 * z;
            		}
            	}
            	
    		} 
    		//else if ( updId == 2 ) 
    		{
    			// update V
    			double coeff = alpha * (1-gamma), coeff2 = coeff * updCnt;
            	for (int k = 0; k < rank; ++k) {
            		FeatureVector dVk = dV[k];
            		for (int i = 0, K = dVk.size(); i < K; ++i) {
            			int x = dVk.x(i);
            			double z = dVk.value(i);
            			V[k][x] += coeff * z;
            			totalV[k][x] += coeff2 * z;
            		}
            	}            	
    		} 
    		//else 
    		{    		
    			// update W
    			double coeff = alpha * (1-gamma), coeff2 = coeff * updCnt;
            	for (int k = 0; k < rank; ++k) {
            		FeatureVector dWk = dW[k];
            		for (int i = 0, K = dWk.size(); i < K; ++i) {
            			int x = dWk.x(i);
            			double z = dWk.value(i);
            			W[k][x] += coeff * z;
            			totalW[k][x] += coeff2 * z;
            		}
            	}  
    		}
    	}
    	
        return loss;
	}
	
	public double updateLabel(DependencyInstance gold, DependencyInstance pred,
			LocalFeatureData lfd, GlobalFeatureData gfd,
			int updCnt, int offset)
	{
		int N = gold.length;
    	int[] actDeps = gold.heads;
    	int[] actLabs = gold.deplbids;
    	int[] predDeps = pred.heads;
    	int[] predLabs = pred.deplbids;
    	
    	double Fi = getLabelDis(actDeps, actLabs, predDeps, predLabs);
        	
    	FeatureVector dtl = lfd.getLabeledFeatureDifference(gold, pred);
    	double loss = - dtl.dotProduct(paramsL) + Fi;
        double l2norm = dtl.Squaredl2NormUnsafe();
    	
        double alpha = loss/l2norm;
    	alpha = Math.min(synC, alpha);
    	if (alpha > 0) {
    		double coeff = alpha;
    		double coeff2 = coeff * updCnt;
    		for (int i = 0, K = dtl.size(); i < K; ++i) {
	    		int x = dtl.x(i);
	    		double z = dtl.value(i);
	    		paramsL[x] += coeff * z;
	    		totalL[x] += coeff2 * z;
    		}
    	}
    	
    	return loss;
	}
	
	public double updateSmn(DependencyInstance gold, DependencyInstance pred,
			SRLFeatureData sfd, int updCnt, int offset)
	{
    	double Fi = getSRLCost(gold, pred);
    	if (Fi == 0) return 0.0;

    	//SRLFeatureData sfd = new SRLFeatureData(pred, options, lfd.pipe, this);
    	FeatureVector dts = sfd.getFeatureDifference(gold, pred);
    
        double loss = -dts.dotProduct(params2)*gamma2 + Fi;
        double l2norm = dts.Squaredl2NormUnsafe()*gamma2*gamma2;
        
        for (int k = 0; k < rank2; ++k) {
        	FeatureVector dUk = getdU2(k, sfd, gold, pred);
        	FeatureVector dVk = getdV2(k, sfd, gold, pred);
        	FeatureVector dWk = getdW2(k, sfd, gold, pred);
        	FeatureVector dXk = getdX2(k, sfd, gold, pred);
        	l2norm += dUk.Squaredl2NormUnsafe()*(1-gamma2)*(1-gamma2);
        	l2norm += dVk.Squaredl2NormUnsafe()*(1-gamma2)*(1-gamma2);
        	l2norm += dWk.Squaredl2NormUnsafe()*(1-gamma2)*(1-gamma2);
        	l2norm += dXk.Squaredl2NormUnsafe()*(1-gamma2)*(1-gamma2);
        	loss -= dUk.dotProduct(U2[k])*(1-gamma2);
        	dU2[k] = dUk;
        	dV2[k] = dVk;
        	dW2[k] = dWk;
        	dX2[k] = dXk;
        }
    	        
        double alpha = loss/l2norm;
        //System.out.println(alpha);
    	alpha = Math.min(smnC, alpha);
    	if (alpha > 0) {
    		
    		{
	    		double coeff = alpha * gamma2;
	    		double coeff2 = coeff * updCnt;
	    		for (int i = 0, K = dts.size(); i < K; ++i) {
		    		int x = dts.x(i);
		    		double z = dts.value(i);
		    		params2[x] += coeff * z;
		    		total2[x] += coeff2 * z;
	    		}
    		}
    		
    		{
    			double coeff = alpha * (1-gamma2), coeff2 = coeff * updCnt;
            	for (int k = 0; k < rank2; ++k) {
            		FeatureVector dUk = dU2[k];
            		for (int i = 0, K = dUk.size(); i < K; ++i) {
            			int x = dUk.x(i);
            			double z = dUk.value(i);
            			U2[k][x] += coeff * z;
            			totalU2[k][x] += coeff2 * z;
            		}
            	}
    		}
    		
    		{
    			double coeff = alpha * (1-gamma2), coeff2 = coeff * updCnt;
            	for (int k = 0; k < rank2; ++k) {
            		FeatureVector dVk = dV2[k];
            		for (int i = 0, K = dVk.size(); i < K; ++i) {
            			int x = dVk.x(i);
            			double z = dVk.value(i);
            			V2[k][x] += coeff * z;
            			totalV2[k][x] += coeff2 * z;
            		}
            	}
    		}
    		
    		{
    			double coeff = alpha * (1-gamma2), coeff2 = coeff * updCnt;
            	for (int k = 0; k < rank2; ++k) {
            		FeatureVector dWk = dW2[k];
            		for (int i = 0, K = dWk.size(); i < K; ++i) {
            			int x = dWk.x(i);
            			double z = dWk.value(i);
            			W2[k][x] += coeff * z;
            			totalW2[k][x] += coeff2 * z;
            		}
            	}
    		}
    		
    		{
    			double coeff = alpha * (1-gamma2), coeff2 = coeff * updCnt;
            	for (int k = 0; k < rank2; ++k) {
            		FeatureVector dXk = dX2[k];
            		for (int i = 0, K = dXk.size(); i < K; ++i) {
            			int x = dXk.x(i);
            			double z = dXk.value(i);
            			X2[k][x] += coeff * z;
            			totalX2[k][x] += coeff2 * z;
            		}
            	}
    		}
    	}
    	
        return loss;
	}
	
	public double update(DependencyInstance gold, DependencyInstance pred,
			LocalFeatureData lfd, GlobalFeatureData gfd, SRLFeatureData sfd,
			int updCnt, int offset)
	{
		int N = gold.length;
    	int[] actDeps = gold.heads;
    	int[] actLabs = gold.deplbids;
    	int[] predDeps = pred.heads;
    	int[] predLabs = pred.deplbids;
    	
    	double Fi = getHammingDis(actDeps, actLabs, predDeps, predLabs)
    			  + getSRLCost(gold, pred);
    	if (Fi == 0) return 0.0;
    	
    	FeatureVector dt = lfd.getFeatureDifference(gold, pred);
    	dt.addEntries(gfd.getFeatureDifference(gold, pred));

    	FeatureVector dtl = new FeatureVector(size);
    	if (options.learnLabel) {
    		dtl = lfd.getLabeledFeatureDifference(gold, pred);
    	}
    	
    	//SRLFeatureData sfd = new SRLFeatureData(pred, options, lfd.pipe, this);
    	FeatureVector dts = sfd.getFeatureDifference(gold, pred);
        
    	double loss = - dt.dotProduct(params)*gamma - dtl.dotProduct(params)*gammaLabel 
    				- dts.dotProduct(params2)*gamma2 + Fi;

        double l2norm = dt.Squaredl2NormUnsafe() * gamma * gamma 
        			  + dtl.Squaredl2NormUnsafe() * gammaLabel * gammaLabel
        			  + dts.Squaredl2NormUnsafe() * gamma2 * gamma2;
    	
        {
        	// update U
        	for (int k = 0; k < rank; ++k) {        		
        		FeatureVector dUk = getdU(k, lfd, actDeps, predDeps);
            	l2norm += dUk.Squaredl2NormUnsafe() * (1-gamma) * (1-gamma);            	
            	loss -= dUk.dotProduct(U[k]) * (1-gamma);
            	dU[k] = dUk;
        	}
        }

        {
        	// update V
        	for (int k = 0; k < rank; ++k) {
        		FeatureVector dVk = getdV(k, lfd, actDeps, predDeps);
            	l2norm += dVk.Squaredl2NormUnsafe() * (1-gamma) * (1-gamma);
            	//loss -= dVk.dotProduct(V[k]) * (1-gamma);
            	dV[k] = dVk;
        	}        	
        } 

        {
        	// update W
        	for (int k = 0; k < rank; ++k) {
        		FeatureVector dWk = getdW(k, lfd, actDeps, predDeps);
            	l2norm += dWk.Squaredl2NormUnsafe() * (1-gamma) * (1-gamma);
            	//loss -= dWk.dotProduct(W[k]) * (1-gamma);
            	dW[k] = dWk;
        	}   
        }
        
        for (int k = 0; k < rank2; ++k) {
        	FeatureVector dUk = getdU2(k, sfd, gold, pred);
        	FeatureVector dVk = getdV2(k, sfd, gold, pred);
        	FeatureVector dWk = getdW2(k, sfd, gold, pred);
        	l2norm += dUk.Squaredl2NormUnsafe()*(1-gamma2)*(1-gamma2);
        	l2norm += dVk.Squaredl2NormUnsafe()*(1-gamma2)*(1-gamma2);
        	l2norm += dWk.Squaredl2NormUnsafe()*(1-gamma2)*(1-gamma2);
        	loss -= dUk.dotProduct(U2[k])*(1-gamma2);
        	dU2[k] = dUk;
        	dV2[k] = dVk;
        	dW2[k] = dWk;
        }
        
        double alpha = loss/l2norm;
    	alpha = Math.min(synC, alpha);
    	if (alpha > 0) {
    		
    		{
    			// update theta
	    		double coeff = alpha * gamma, coeff2 = coeff * updCnt;
	    		for (int i = 0, K = dt.size(); i < K; ++i) {
		    		int x = dt.x(i);
		    		double z = dt.value(i);
		    		params[x] += coeff * z;
		    		total[x] += coeff2 * z;
	    		}
	    		coeff = alpha * gammaLabel;
	    		coeff2 = coeff * updCnt;
	    		for (int i = 0, K = dtl.size(); i < K; ++i) {
		    		int x = dtl.x(i);
		    		double z = dtl.value(i);
		    		params[x] += coeff * z;
		    		total[x] += coeff2 * z;
	    		}
	    		coeff = alpha * gamma2;
	    		coeff2 = coeff * updCnt;
	    		for (int i = 0, K = dts.size(); i < K; ++i) {
		    		int x = dts.x(i);
		    		double z = dts.value(i);
		    		params2[x] += coeff * z;
		    		total2[x] += coeff2 * z;
	    		}
    		}
    		
    		{
    			// update U
    			double coeff = alpha * (1-gamma), coeff2 = coeff * updCnt;
            	for (int k = 0; k < rank; ++k) {
            		FeatureVector dUk = dU[k];
            		for (int i = 0, K = dUk.size(); i < K; ++i) {
            			int x = dUk.x(i);
            			double z = dUk.value(i);
            			U[k][x] += coeff * z;
            			totalU[k][x] += coeff2 * z;
            		}
            	}
            	
    		} 

    		{
    			// update V
    			double coeff = alpha * (1-gamma), coeff2 = coeff * updCnt;
            	for (int k = 0; k < rank; ++k) {
            		FeatureVector dVk = dV[k];
            		for (int i = 0, K = dVk.size(); i < K; ++i) {
            			int x = dVk.x(i);
            			double z = dVk.value(i);
            			V[k][x] += coeff * z;
            			totalV[k][x] += coeff2 * z;
            		}
            	}            	
    		} 

    		{    		
    			// update W
    			double coeff = alpha * (1-gamma), coeff2 = coeff * updCnt;
            	for (int k = 0; k < rank; ++k) {
            		FeatureVector dWk = dW[k];
            		for (int i = 0, K = dWk.size(); i < K; ++i) {
            			int x = dWk.x(i);
            			double z = dWk.value(i);
            			W[k][x] += coeff * z;
            			totalW[k][x] += coeff2 * z;
            		}
            	}  
    		}
    		
    		{
    			double coeff = alpha * (1-gamma2), coeff2 = coeff * updCnt;
            	for (int k = 0; k < rank2; ++k) {
            		FeatureVector dUk = dU2[k];
            		for (int i = 0, K = dUk.size(); i < K; ++i) {
            			int x = dUk.x(i);
            			double z = dUk.value(i);
            			U2[k][x] += coeff * z;
            			totalU2[k][x] += coeff2 * z;
            		}
            	}
    		}
    		
    		{
    			double coeff = alpha * (1-gamma2), coeff2 = coeff * updCnt;
            	for (int k = 0; k < rank2; ++k) {
            		FeatureVector dVk = dV2[k];
            		for (int i = 0, K = dVk.size(); i < K; ++i) {
            			int x = dVk.x(i);
            			double z = dVk.value(i);
            			V2[k][x] += coeff * z;
            			totalV2[k][x] += coeff2 * z;
            		}
            	}
    		}
    		
    		{
    			double coeff = alpha * (1-gamma2), coeff2 = coeff * updCnt;
            	for (int k = 0; k < rank2; ++k) {
            		FeatureVector dWk = dW2[k];
            		for (int i = 0, K = dWk.size(); i < K; ++i) {
            			int x = dWk.x(i);
            			double z = dWk.value(i);
            			W2[k][x] += coeff * z;
            			totalW2[k][x] += coeff2 * z;
            		}
            	}
    		}
    	}
    	
        return loss;
	}
	
//	public double update(DependencyInstance gold, DependencyInstance pred,
//			LocalFeatureData lfd, GlobalFeatureData gfd,
//			int updCnt, int offset)
//	{
//		int N = gold.length;
//    	int[] actDeps = gold.heads;
//    	int[] actLabs = gold.deplbids;
//    	int[] predDeps = pred.heads;
//    	int[] predLabs = pred.deplbids;
//    	
//    	double Fi = getHammingDis(actDeps, actLabs, predDeps, predLabs) +
//    				getSRLCost(gold, pred);
//    	if (Fi == 0) return 0.0;
//    	
//    	FeatureVector dt = lfd.getFeatureDifference(gold, pred);
//    	dt.addEntries(gfd.getFeatureDifference(gold, pred));
//
//    	//FeatureVector dts = lfd.pipe.smnFactory.getFeatureVector(gold);
//    	//dts.addEntries(lfd.pipe.smnFactory.getFeatureVector(pred), -1.0);
//    	SRLFeatureData sfd = new SRLFeatureData(pred, options, lfd.pipe, this);
//    	FeatureVector dts = sfd.getFeatureDifference(gold, pred);
//    	
//    	FeatureVector dtl = new FeatureVector(size);
//    	if (options.learnLabel) {
//    		dtl = lfd.getLabeledFeatureDifference(gold, pred);
//    	}
//    	
//        double loss = - dt.dotProduct(params)*gamma - dtl.dotProduct(params)*gammaLabel
//        			  - dts.dotProduct(params2) + Fi;
//
//        double l2norm = dt.Squaredl2NormUnsafe() * gamma * gamma 
//        			  + dtl.Squaredl2NormUnsafe() * gammaLabel * gammaLabel
//        			  + dts.Squaredl2NormUnsafe();
//    	
//        int updId = (updCnt + offset) % 3;
//        if ( updId == 1 ) {
//        	// update U
//        	for (int k = 0; k < rank; ++k) {        		
//        		FeatureVector dUk = getdU(k, lfd, actDeps, predDeps);
//            	l2norm += dUk.Squaredl2NormUnsafe() * (1-gamma) * (1-gamma);            	
//            	loss -= dUk.dotProduct(U[k]) * (1-gamma);
//            	dU[k] = dUk;
//        	}
//        } else if ( updId == 2 ) {
//        	// update V
//        	for (int k = 0; k < rank; ++k) {
//        		FeatureVector dVk = getdV(k, lfd, actDeps, predDeps);
//            	l2norm += dVk.Squaredl2NormUnsafe() * (1-gamma) * (1-gamma);
//            	loss -= dVk.dotProduct(V[k]) * (1-gamma);
//            	dV[k] = dVk;
//        	}        	
//        } else {
//        	// update W
//        	for (int k = 0; k < rank; ++k) {
//        		FeatureVector dWk = getdW(k, lfd, actDeps, predDeps);
//            	l2norm += dWk.Squaredl2NormUnsafe() * (1-gamma) * (1-gamma);
//            	loss -= dWk.dotProduct(W[k]) * (1-gamma);
//            	dW[k] = dWk;
//        	}   
//        }
//        
//        double alpha = loss/l2norm;
//    	alpha = Math.min(C, alpha);
//    	if (alpha > 0) {
//    		
//    		{
//    			// update theta
//	    		double coeff = alpha * gamma, coeff2 = coeff * updCnt;
//	    		for (int i = 0, K = dt.size(); i < K; ++i) {
//		    		int x = dt.x(i);
//		    		double z = dt.value(i);
//		    		params[x] += coeff * z;
//		    		total[x] += coeff2 * z;
//	    		}
//	    		coeff = alpha * gammaLabel;
//	    		coeff2 = coeff * updCnt;
//	    		for (int i = 0, K = dtl.size(); i < K; ++i) {
//		    		int x = dtl.x(i);
//		    		double z = dtl.value(i);
//		    		params[x] += coeff * z;
//		    		total[x] += coeff2 * z;
//	    		}
//	    		coeff = alpha;
//	    		coeff2 = coeff * updCnt;
//	    		for (int i = 0, K = dts.size(); i < K; ++i) {
//		    		int x = dts.x(i);
//		    		double z = dts.value(i);
//		    		params2[x] += coeff * z;
//		    		total2[x] += coeff2 * z;
//	    		}
//    		}
//    		
//    		if ( updId == 1 ) {
//    			// update U
//    			double coeff = alpha * (1-gamma), coeff2 = coeff * updCnt;
//            	for (int k = 0; k < rank; ++k) {
//            		FeatureVector dUk = dU[k];
//            		for (int i = 0, K = dUk.size(); i < K; ++i) {
//            			int x = dUk.x(i);
//            			double z = dUk.value(i);
//            			U[k][x] += coeff * z;
//            			totalU[k][x] += coeff2 * z;
//            		}
//            	}
//            	
//    		} else if ( updId == 2 ) {
//    			// update V
//    			double coeff = alpha * (1-gamma), coeff2 = coeff * updCnt;
//            	for (int k = 0; k < rank; ++k) {
//            		FeatureVector dVk = dV[k];
//            		for (int i = 0, K = dVk.size(); i < K; ++i) {
//            			int x = dVk.x(i);
//            			double z = dVk.value(i);
//            			V[k][x] += coeff * z;
//            			totalV[k][x] += coeff2 * z;
//            		}
//            	}            	
//    		} else {
//    			// update W
//    			double coeff = alpha * (1-gamma), coeff2 = coeff * updCnt;
//            	for (int k = 0; k < rank; ++k) {
//            		FeatureVector dWk = dW[k];
//            		for (int i = 0, K = dWk.size(); i < K; ++i) {
//            			int x = dWk.x(i);
//            			double z = dWk.value(i);
//            			W[k][x] += coeff * z;
//            			totalW[k][x] += coeff2 * z;
//            		}
//            	}  
//    		}
//    	}
//    	
//        return loss;
//	}
	
	public void updateTheta(FeatureVector gold, FeatureVector pred, double loss,
			int updCnt) 
	{
		FeatureVector fv = new FeatureVector(size);
		fv.addEntries(gold);
		fv.addEntries(pred, -1.0);
		
		double l2norm = fv.Squaredl2NormUnsafe();
		double alpha = loss/l2norm;
	    alpha = Math.min(synC, alpha);
	    if (alpha > 0) {
			// update theta
    		double coeff = alpha, coeff2 = coeff * updCnt;
    		for (int i = 0, K = fv.size(); i < K; ++i) {
	    		int x = fv.x(i);
	    		double z = fv.value(i);
	    		params[x] += coeff * z;
	    		total[x] += coeff2 * z;
    		}
	    }
	}
	
    private FeatureVector getdU(int k, LocalFeatureData lfd, int[] actDeps, int[] predDeps) 
    {
    	double[][] wpV = lfd.wpV;
    	FeatureVector[] wordFvs = lfd.wordFvs;
    	int L = wordFvs.length;
    	FeatureVector dU = new FeatureVector(N);
    	for (int mod = 1; mod < L; ++mod) {
    		int head  = actDeps[mod];
    		int head2 = predDeps[mod];
    		if (head == head2) continue;
    		int d = getBinnedDistance(head-mod);
    		int d2 = getBinnedDistance(head2-mod);
    		double dotv = wpV[mod][k]; //wordFvs[mod].dotProduct(V[k]);    		
    		dU.addEntries(wordFvs[head], dotv * (W[k][0] + W[k][d]));
    		dU.addEntries(wordFvs[head2], - dotv * (W[k][0] + W[k][d2]));
    	}
    	return dU;
    }
    
    private FeatureVector getdV(int k, LocalFeatureData lfd, int[] actDeps, int[] predDeps)
    {
    	double[][] wpU = lfd.wpU;
    	FeatureVector[] wordFvs = lfd.wordFvs;
    	int L = wordFvs.length;
    	FeatureVector dV = new FeatureVector(M);
    	for (int mod = 1; mod < L; ++mod) {
    		int head  = actDeps[mod];
    		int head2 = predDeps[mod];
    		if (head == head2) continue;
    		int d = getBinnedDistance(head-mod);
    		int d2 = getBinnedDistance(head2-mod);
    		double dotu = wpU[head][k];   //wordFvs[head].dotProduct(U[k]);
    		double dotu2 = wpU[head2][k]; //wordFvs[head2].dotProduct(U[k]);
    		dV.addEntries(wordFvs[mod], dotu  * (W[k][0] + W[k][d])
    									- dotu2 * (W[k][0] + W[k][d2]));    		
    	}
    	return dV;
    }
    
    private FeatureVector getdW(int k, LocalFeatureData lfd, int[] actDeps, int[] predDeps)
    {
    	double[][] wpU = lfd.wpU, wpV = lfd.wpV;
    	FeatureVector[] wordFvs = lfd.wordFvs;
    	int L = wordFvs.length;
    	double[] dW = new double[D];
    	for (int mod = 1; mod < L; ++mod) {
    		int head  = actDeps[mod];
    		int head2 = predDeps[mod];
    		if (head == head2) continue;
    		int d = getBinnedDistance(head-mod);
    		int d2 = getBinnedDistance(head2-mod);
    		double dotu = wpU[head][k];   //wordFvs[head].dotProduct(U[k]);
    		double dotu2 = wpU[head2][k]; //wordFvs[head2].dotProduct(U[k]);
    		double dotv = wpV[mod][k];  //wordFvs[mod].dotProduct(V[k]);
    		dW[0] += (dotu - dotu2) * dotv;
    		dW[d] += dotu * dotv;
    		dW[d2] -= dotu2 * dotv;
    	}
    	
    	FeatureVector dW2 = new FeatureVector(D);
    	for (int i = 0; i < D; ++i)
    		dW2.addEntry(i, dW[i]);
    	return dW2;
    }
    
    private FeatureVector getdU2(int k, SRLFeatureData sfd,
    		DependencyInstance gold, DependencyInstance pred)
    {
    	double[][] wpU = sfd.wpU, wpV = sfd.wpV, ppW = sfd.ppW, cpX = sfd.cpX;
    	FeatureVector dU2 = new FeatureVector(N2);
    	int F = pred.numframes, N = pred.length, L = sfd.L;
		for (int i = 0; i < F; ++i) {
			SemanticFrame frame = gold.frames[i];
			SemanticFrame frame2 = pred.frames[i];
			int p = frame.predid;
			Utils.Assert(frame.predid == frame2.predid);
			Utils.Assert(sfd.gamma2 == gamma2);
			for (int a = 0; a < N; ++a) {
				if (frame.arglbids[a] == frame2.arglbids[a]) continue;
				//boolean isValid = sfd.isValidPredAugPair(p, a);
				boolean isValid = !sfd.isPruned(p, a);
				
				{
					int r = frame.arglbids[a];
					if (isValid && r >= 0) {
						int id = i*N*L + a*L + r;
						double dot = wpV[a][k] * ppW[id][k] * cpX[id][k];
						dU2.addEntries(sfd.wordFvs[p], dot);
					}
				}
				
				{
					int r = frame2.arglbids[a];
					if (isValid && r >= 0) {
						int id = i*N*L + a*L + r;
						double dot = wpV[a][k] * ppW[id][k] * cpX[id][k];
						dU2.addEntries(sfd.wordFvs[p], -dot);
					}
				}
			}
		}
		return dU2;
    }
    
    private FeatureVector getdV2(int k, SRLFeatureData sfd,
    		DependencyInstance gold, DependencyInstance pred)
    {
    	double[][] wpU = sfd.wpU, wpV = sfd.wpV, ppW = sfd.ppW, cpX = sfd.cpX;
    	FeatureVector dV2 = new FeatureVector(M2);
    	int F = pred.numframes, N = pred.length, L = sfd.L;
		for (int i = 0; i < F; ++i) {
			SemanticFrame frame = gold.frames[i];
			SemanticFrame frame2 = pred.frames[i];
			int p = frame.predid;
			for (int a = 0; a < N; ++a) {
				if (frame.arglbids[a] == frame2.arglbids[a]) continue;
				//boolean isValid = sfd.isValidPredAugPair(p, a);
				boolean isValid = !sfd.isPruned(p, a);
				
				{
					int r = frame.arglbids[a];
					if (isValid && r >= 0) {
						int id = i*N*L + a*L + r;
						double dot = wpU[p][k] * ppW[id][k] * cpX[id][k];
						dV2.addEntries(sfd.wordFvs[a], dot);
					}
				}
				
				{
					int r = frame2.arglbids[a];
					if (isValid && r >= 0) {
						int id = i*N*L + a*L + r;
						double dot = wpU[p][k] * ppW[id][k] * cpX[id][k];
						dV2.addEntries(sfd.wordFvs[a], -dot);
					}
				}
			}
		}
		return dV2;
    }
    
    private FeatureVector getdW2(int k, SRLFeatureData sfd,
    		DependencyInstance gold, DependencyInstance pred)
    {
    	double[][] wpU = sfd.wpU, wpV = sfd.wpV, ppW = sfd.ppW, cpX = sfd.cpX;
    	FeatureVector dW2 = new FeatureVector(D2);
    	int F = pred.numframes, N = pred.length, L = sfd.L;
		for (int i = 0; i < F; ++i) {
			SemanticFrame frame = gold.frames[i];
			SemanticFrame frame2 = pred.frames[i];
			int p = frame.predid;
			for (int a = 0; a < N; ++a) {
				if (frame.arglbids[a] == frame2.arglbids[a]) continue;
				//boolean isValid = sfd.isValidPredAugPair(p, a);
				boolean isValid = !sfd.isPruned(p, a);
				
				{
					int r = frame.arglbids[a];
					if (isValid && r >= 0) {
						int id = i*N*L + a*L + r;
						double dot = wpU[p][k] * wpV[a][k] * cpX[id][k];
						dW2.addEntries(sfd.pathFvs[id], dot);
					}
				}
				
				{
					int r = frame2.arglbids[a];
					if (isValid && r >= 0) {
						int id = i*N*L + a*L + r;
						double dot = wpU[p][k] * wpV[a][k] * cpX[id][k];
						dW2.addEntries(sfd.pathFvs[id], -dot);
					}
				}
			}
		}
		return dW2;
    }
    
    private FeatureVector getdX2(int k, SRLFeatureData sfd,
    		DependencyInstance gold, DependencyInstance pred)
    {
    	double[][] wpU = sfd.wpU, wpV = sfd.wpV, ppW = sfd.ppW, cpX = sfd.cpX;
    	FeatureVector dX2 = new FeatureVector(L2);
    	int F = pred.numframes, N = pred.length, L = sfd.L;
		for (int i = 0; i < F; ++i) {
			SemanticFrame frame = gold.frames[i];
			SemanticFrame frame2 = pred.frames[i];
			int p = frame.predid;
			for (int a = 0; a < N; ++a) {
				if (frame.arglbids[a] == frame2.arglbids[a]) continue;
				//boolean isValid = sfd.isValidPredAugPair(p, a);
				boolean isValid = !sfd.isPruned(p, a);
				
				{
					int r = frame.arglbids[a];
					if (isValid && r >= 0) {
						int id = i*N*L + a*L + r;
						double dot = wpU[p][k] * wpV[a][k] * ppW[id][k];
						dX2.addEntries(sfd.contextFvs[id], dot);
					}
				}
				
				{
					int r = frame2.arglbids[a];
					if (isValid && r >= 0) {
						int id = i*N*L + a*L + r;
						double dot = wpU[p][k] * wpV[a][k] * ppW[id][k];
						dX2.addEntries(sfd.contextFvs[id], -dot);
					}
				}
			}
		}
		return dX2;
    }
    
	public double getHammingDis(int[] actDeps, int[] actLabs,
			int[] predDeps, int[] predLabs) 
	{
		double dis = 0;
		for (int i = 1; i < actDeps.length; ++i)
			if (options.learnLabel) {
				if (labelLossType == 0) {
					if (actDeps[i] != predDeps[i]) dis += 0.5;
					if (actLabs[i] != predLabs[i]) dis += 0.5;
				} else if (actDeps[i] != predDeps[i] || actLabs[i] != predLabs[i]) dis += 1;
			} else {
				if (actDeps[i] != predDeps[i]) dis += 1;
			}
		return dis;
    }
	
	public double getSRLCost(DependencyInstance goldinst, DependencyInstance predinst)
	{
		SemanticFrame[] gold = goldinst.frames, pred = predinst.frames;
		Utils.Assert(gold.length == pred.length);
		double dis = 0;
		for (int i = 0, N = gold.length; i < N; ++i) {
			
			Utils.Assert(gold[i].predid == pred[i].predid);
			
			int pid = gold[i].predid;
			int[] ga = gold[i].arglbids, pa = pred[i].arglbids;
			for (int j = 0, L = ga.length; j < L; ++j) {
				
				// "project" gold links to the prediction
				// remove this link if it is not valid on the predicted tree
				int garg = ga[j];
				if (garg >= 0 && !SemanticFeatureFactory.isValidPredAugPair(predinst, pid, j))
					garg = -1;
				
				if (garg != pa[j]) {
					
					// false positive or false negative
					if (garg < 0) dis += 1.0;
					else if (pa[j] < 0) dis += 2.0;
					// incorrect label
					else dis += 0.5;
					//dis+=1.0;
				}
			}
		}
		return dis;
	}
	
    public int getBinnedDistance(int x)
    {
    	int y = x > 0 ? x : -x;
    	int dis = 0;
    	if (y > 10)
    		dis = 7;
    	else if (y > 5)
    		dis = 6;
    	else dis = y;
    	if (dis > d) dis = d;
    	return x > 0 ? dis : dis + d;
    }

    public double getLabelDis(int[] actDeps, int[] actLabs,
			int[] predDeps, int[] predLabs) 
	{
		double dis = 0;
		for (int i = 1; i < actLabs.length; ++i) {
			assert(actDeps[i] == predDeps[i]);
			if (actLabs[i] != predLabs[i]) dis += 1;
		}
		return dis;
    }
}
