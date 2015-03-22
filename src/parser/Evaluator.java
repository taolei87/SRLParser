package parser;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import parser.feature.SemanticFeatureFactory;
import utils.Utils;

public class Evaluator
{
	int uas, las, tot;
	int whole, nsents;
	int corr, totp, totg;
	
	int vis;
	
	boolean learnLabel;
	
	int times, numArgs;
	int[] timeStamps, argFreqCnts, argAppearCnts;
	String[] argLabels;
	
	TIntIntHashMap goldlengthCounts, predlengthCounts;
	int[] corrPL, totPL, corrGL, totGL;
	
	public Evaluator(Options options, DependencyPipe pipe)
	{
		uas = las = tot = 0;
		corr = totp = totg = 0;
		whole = nsents = 0;
		vis = 0;
		learnLabel = options.learnLabel;
		
		numArgs = pipe.smnFactory.numSemanticLabels;
		argLabels = pipe.args;
		timeStamps = new int[numArgs];
		argFreqCnts = new int[numArgs];
		argAppearCnts = new int[numArgs];
		
		goldlengthCounts = new TIntIntHashMap();
		predlengthCounts = new TIntIntHashMap();
		corrPL = new int[15];
		totPL = new int[15];
		corrGL = new int[15];
		totGL = new int[15];
	}
	
	public double FilteringRecall()
	{
		return vis/(totg+1e-20);
	}
	
	public double UAS()
	{
		return uas/(tot+1e-20);
	}
	
	public double LAS()
	{
		return las/(tot+1e-20);
	}
	
	public double CAS()
	{
		return whole/(nsents+1e-20);
	}
	
	public double Precision()
	{
		return corr/(totp+1e-20);
	}
	
	public double Recall()
	{
		return corr/(totg+1e-20);
	}
	
	public double F1()
	{
		double p = Precision();
		double r = Recall();
		return p*r*2/(p+r);
	}
	
	public void dumpArgStats()
	{
		System.out.println();
		for (int i = 0; i < numArgs; ++i) {
			System.out.printf("\t%s  %.2f (%d/%d)%n",
					argLabels[i],
					argAppearCnts[i]/(argFreqCnts[i]+1e-20),
					argAppearCnts[i],
					argFreqCnts[i]);
		}
		System.out.println();
	}
	
	public void add(DependencyInstance gold, DependencyInstance predict, boolean evalWithPunc)
	{
		evaluateDependencies(gold, predict, evalWithPunc);
		evaluateSemanticLabeling(gold, predict);
	}
	
    public void evaluateDependencies(DependencyInstance gold, 
    		DependencyInstance pred, boolean evalWithPunc) 
    {
    	++nsents;
    	int tt = 0, ua = 0, la = 0;
    	for (int i = 1, N = gold.length; i < N; ++i) {

            if (!evalWithPunc)
            	if (gold.forms[i].matches("[-!\"%&'()*,./:;?@\\[\\]_{}、]+")) continue;

            ++tt;
    		if (gold.heads[i] == pred.heads[i]) {
    			++ua;
    			if (learnLabel && gold.deplbids[i] == pred.deplbids[i]) ++la;
    		}
    	
    	}    		
    	
    	tot += tt;
    	uas += ua;
    	las += la;
    	whole += (tt == ua) && (tt == la || !learnLabel) ? 1 : 0;
    }
    
    public void evaluateSemanticLabeling(DependencyInstance gold, DependencyInstance pred)
    {
    	Utils.Assert(gold.numframes == pred.numframes);
    	int n = gold.length;
    	for (int k = 0; k < gold.numframes; ++k) {
    		Utils.Assert(gold.frames[k].predid == pred.frames[k].predid);
    		int p = gold.frames[k].predid;
    		int[] ga = gold.frames[k].arglbids;
    		int[] pa = pred.frames[k].arglbids;
    		for (int i = 0; i < n; ++i) {
    			if (ga[i] >= 0) ++totg;
    			if (pa[i] >= 0) ++totp;
    			if (ga[i] >= 0 && ga[i] == pa[i]) ++corr;
    			if (ga[i] >= 0 && SemanticFeatureFactory.isValidPredAugPair(pred, p, i))
    				++vis;
    			
    			if (ga[i] >= 0 && SemanticFeatureFactory.isValidPredAugPair(gold, p, i)) {
    				int length = SemanticFeatureFactory.getPathLength(gold, p, i);
    				goldlengthCounts.adjustOrPutValue(length, 1, 1);
    				++totGL[length];
    				if (pa[i] == ga[i]) ++corrGL[length];
    			}
    			if (pa[i] >= 0) {
    				int length = SemanticFeatureFactory.getPathLength(pred, p, i);
    				predlengthCounts.adjustOrPutValue(length, 1, 1);
    				++totPL[length];
    				if (pa[i] == ga[i]) ++corrPL[length]; 
    			}
    		}
    		
    		++times;
    		for (int i = 0; i < n; ++i) {
    			int x = ga[i];
    			if (x == -1) continue;
    			if (timeStamps[x] != times) {
    				timeStamps[x] = times;
    				++argAppearCnts[x];
    			}
    			++argFreqCnts[x];
    		}
    	}
    }

    public void dumpPathStats() throws IOException
    {
		BufferedWriter writer = new BufferedWriter(new FileWriter("path.info"));
		for (int length : goldlengthCounts.keys()) {
			writer.write(""+length);
			writer.write("\t");
			writer.write(""+goldlengthCounts.get(length));
			writer.write("\n");
		}
		for (int length : predlengthCounts.keys()) {
			writer.write(""+length);
			writer.write("\t");
			writer.write(""+predlengthCounts.get(length));
			writer.write("\n");
		}
		for (int i = 0;i < corrPL.length; ++i) {
			double prec = corrPL[i]/(totPL[i]+1e-20);
			writer.write(i+ " " + prec + " " + corrPL[i] + " " + totPL[i] + "\n");
		}
		for (int i = 0;i < corrGL.length; ++i) {
			double rec = corrGL[i]/(totGL[i]+1e-20);
			writer.write(i+ " " + rec + " " + corrGL[i] + " " + totGL[i] + "\n");
		}
		writer.close();
    }
}
