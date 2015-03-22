package parser.feature;

import gnu.trove.list.array.TIntArrayList;

import java.io.Serializable;

import parser.DependencyInstance;
import parser.DependencyPipe;
import parser.Options;
import parser.Parameters;
import parser.SemanticFrame;
import parser.SemanticLowRankParam;
import utils.Alphabet;
import utils.DictionarySet;
import utils.FeatureVector;
import utils.Utils;
import static parser.feature.FeatureTemplate.Link.*;
import static parser.feature.FeatureTemplate.Path.*;
import static parser.feature.FeatureTemplate.Word.*;
import static parser.feature.FeatureTemplate.Context.*;

public class SemanticFeatureFactory implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	public static int MAX_DEPTH = 7;
    //public static int MAX_ARG = 6;
		
	public int TOKEN_START = 1;
	public int TOKEN_END = 2;
	public int TOKEN_MID = 3;
	
	public Options options;
	
	public int numSemanticLabels;
	public int tagNumBits, wordNumBits, deplbNumBits, auglbNumBits;
	public Alphabet smnAlphabet, pathcodeAlphabet, argSeqAlphabet;
	
	public int numLinkFeats;
	
	public boolean usePathWW, usePathWP, usePathPP, useWWW;
	public String[] args;
	
	public Alphabet wordAlphabet, pathAlphabet, contextAlphabet;
	public int numWordFeats, numPathFeats, numContextFeats;
	public double[][] wordVectors = null;
	public double[] unknownWv = null;
	
	public SemanticFeatureFactory(Options options, DependencyPipe pipe)
	{
		this.options = options;
		smnAlphabet = new Alphabet();
		pathcodeAlphabet = new Alphabet();
		argSeqAlphabet = new Alphabet();
		numLinkFeats = 0;
		numSemanticLabels = 0;
		args = pipe.args;
		
		// low-rank tensor atomic features
		wordAlphabet = new Alphabet();
		pathAlphabet = new Alphabet();
		contextAlphabet = new Alphabet();
		numWordFeats = 0;
		numPathFeats = 0;
		numContextFeats = 0;
	}
	
	public void closeAlphabets()
	{
		smnAlphabet.stopGrowth();
		pathcodeAlphabet.stopGrowth();
		argSeqAlphabet.stopGrowth();
		
		wordAlphabet.stopGrowth();
		pathAlphabet.stopGrowth();
        contextAlphabet.stopGrowth();
	}
		
	public void initFeatureAlphabets()
	{
		int pathBits = (deplbNumBits+1)*(MAX_DEPTH+1);
		if (pathBits > 64)
			System.out.println("WARNING: path bits exceed 64!");
		
		//useWWW = wordNumBits*3 + numLinkFeatBits + auglbNumBits <= 64;
		
		System.out.println("useWWW: " + useWWW);	

	}
	
	public void initFeatureAlphabets(DependencyInstance inst) 
    {
		
		//STRATEGY 1: only use gold (p,a,r)
		getFeatureVector(inst);
		
		int n = inst.length;
		//int m = inst.numframes;
		for (SemanticFrame frame : inst.frames) {
			int p = frame.predid;
			createWordFeatureVector(inst, p);
			int[] args = frame.arglbids;
			for (int a = 0; a < n; ++a) {
				boolean isValid = isValidPredAugPair(inst, p, a);
				if (args[a] >= 0 && isValid) {
					int r = args[a];
					createContextFeatureVector(inst, p, a, r);
					createPathFeatureVector(inst, p, a, r);
					createWordFeatureVector(inst, a);
				}
			}
		}
		
		
		//STRATEGY 2: use all (p,a,r) pair
//		int n = inst.length;
//		int m = inst.numframes; 
//		for (SemanticFrame frame : inst.frames) {
//			int p = frame.predid;
//			int[] args = frame.arglbids;
//			for (int a = 0; a < n; ++a) 
//				if (isValidPredAugPair(inst, p, a)) {
//					for (int r = 0; r < numSemanticLabels; ++r)
//						createPredAugLinkFeatures(inst, p, a, r);
//				}
//		}
    }
	
	public FeatureVector createContextFeatureVector(DependencyInstance inst,
			int p, int a, int r)
	{
		Utils.Assert(isValidPredAugPair(inst, p, a));
		
		int[] heads = inst.heads;
		int[] depids = inst.deplbids;
		//int[] posids = inst.postagids;
		FeatureVector fv = new FeatureVector(numContextFeats);
		
		long code = 0;
		code = createContextCodeW(CNTFV_BIAS, 0);
		addContextFeature(code, fv);
		
		code = createContextCodeW(CNTFV_LABEL, r);
		addContextFeature(code, fv);
		
		return fv;
	}
	
	public FeatureVector createPathFeatureVector(DependencyInstance inst,
			int p, int a, int r)
	{
		Utils.Assert(isValidPredAugPair(inst, p, a));
		
		int[] heads = inst.heads;
		int[] depids = inst.deplbids;
		//int[] posids = inst.postagids;
		FeatureVector fv = new FeatureVector(numPathFeats);
		
		long code = 0;
		code = createPathCodeW(PATHFV_BIAS, 0);
		addPathFeature(code, fv);
		
		//code = createPathCodeW(PATHFV_LABEL, r);
		//addPathFeature(code, fv);
		
		code = createPathCode(heads, depids, p, a);
		code = createPathCodeW(PATHFV_PATHCODE, code);
		addPathFeature(code, fv);
		
//		code = getPathLength(inst, p, a);
//		code = createPathCodeW(PATHFV_DIS, code);
//		addPathFeature(code, fv);
//		
//		int pa = heads[a];
//		int gpa = pa >= 0 ? heads[pa] : pa;
//		for (int i = 0, x = p; i <= MAX_DEPTH && x >= 0; ++i) {
//			if (x == a) break;
//			if (x == pa) {
//				int labeldir = (depids[a] << 1) | 1;
//				code = createPathCodeW(PATHFV_DEP, labeldir);
//				addPathFeature(code, fv);				
//				break;
//			}
//			int labeldir = depids[x] << 1;
//			code = createPathCodeW(PATHFV_DEP, labeldir);
//			addPathFeature(code, fv);			
//			x = heads[x];
//		}		
		
		return fv;
	}
	
	public FeatureVector createWordFeatureVector(DependencyInstance inst, int i)
	{
		int[] pos = inst.postagids;
        int[] toks = inst.formids;
        int[][] feats = inst.featids;
        
        int w0 = toks[i];
        int l0 = inst.lemmaids == null ? 0 : inst.lemmaids[i];
        
        FeatureVector fv = new FeatureVector(numWordFeats);
    	
    	long code = 0;
        
    	code = createWordCodeP(WORDFV_BIAS, 0);
    	addWordFeature(code, fv);

    	code = createWordCodeW(WORDFV_W0, w0);
    	addWordFeature(code, fv);
    	
    	int Wp = i == 0 ? TOKEN_START : toks[i-1];
    	int Wn = i == inst.length - 1 ? TOKEN_END : toks[i+1];
    		    	
    	code = createWordCodeW(WORDFV_Wp, Wp);
    	addWordFeature(code, fv);
    	
    	code = createWordCodeW(WORDFV_Wn, Wn);
    	addWordFeature(code, fv);

    	
		if (l0 != 0) {
    		code = createWordCodeW(WORDFV_W0, l0);
    		addWordFeature(code, fv);
    		
	    	int Lp = i == 0 ? TOKEN_START : inst.lemmaids[i-1];
	    	int Ln = i == inst.length - 1 ? TOKEN_END : inst.lemmaids[i+1];
	    		    	
	    	code = createWordCodeW(WORDFV_Wp, Lp);
	    	addWordFeature(code, fv);
	    	
	    	code = createWordCodeW(WORDFV_Wn, Ln);
	    	addWordFeature(code, fv);
		}
		
		if (feats[i] != null) {
    		for (int u = 0; u < feats[i].length; ++u) {
    			int f = feats[i][u];
    			
    			code = createWordCodeP(WORDFV_P0, f);
    			addWordFeature(code, fv);
    			
                if (l0 != 0) {
                	code = createWordCodeWP(WORDFV_W0P0, l0, f);
                	addWordFeature(code, fv);
                }
                
            }
		}
			
        int p0 = pos[i];
    	int pLeft = i > 0 ? pos[i-1] : TOKEN_START;
    	int pRight = i < pos.length-1 ? pos[i+1] : TOKEN_END;
    	
    	code = createWordCodeP(WORDFV_P0, p0);
    	addWordFeature(code, fv);
    	code = createWordCodeP(WORDFV_Pp, pLeft);
    	addWordFeature(code, fv);
    	code = createWordCodeP(WORDFV_Pn, pRight);
    	addWordFeature(code, fv);
    	code = createWordCodePP(WORDFV_PpP0, pLeft, p0);
    	addWordFeature(code, fv);
    	code = createWordCodePP(WORDFV_P0Pn, p0, pRight);
    	addWordFeature(code, fv);
    	code = createWordCodePPP(WORDFV_PpP0Pn, pLeft, p0, pRight);
    	addWordFeature(code, fv);
    		    	
		if (l0 != 0) {
    		code = createWordCodeWP(WORDFV_W0P0, l0, p0);
    		addWordFeature(code, fv);
		}
    	    	
    	if (wordVectors != null) {
    		addWordVectorFeatures(inst, i, 0, fv);
    		//addWordVectorFeatures(inst, i, -1, fv);
    		//addWordVectorFeatures(inst, i, 1, fv);	
    	}
    	
    	int voice = inst.voice[i];
    	code = createWordCodeP(WORDFV_VOICE, voice);
    	addWordFeature(code, fv);
    	
    	if (l0 != 0) {
    		code = createWordCodeP(WORDFV_L0VOICE, (l0 << 2) | voice);
    		addWordFeature(code, fv);
    	}
    	
    	return fv;
	}
	
    public void addWordVectorFeatures(DependencyInstance inst, int i, int dis, FeatureVector fv) {
    	
    	int d = getBinnedDistance(dis);
    	double [] v = unknownWv;
    	int pos = i + dis;
    	
    	if (pos >= 0 && pos < inst.length) {
    		int wvid = inst.wordVecIds[pos];
    		if (wvid > 0) v = wordVectors[wvid];
    	}
    	
		//if (v == unknownWv) ++wvMiss; else ++wvHit;
		
		if (v != null) {
			for (int j = 0; j < v.length; ++j) {
				long code = createWordCodeW(WORDFV_EMB, j);
				addWordFeature(code | d, v[j], fv);
			}
		}
    }
    
    public int getBinnedDistance(int x) {
    	int flag = 0;
    	int add = 0;
    	if (x < 0) {
    		x = -x;
    		//flag = 8;
    		add = 7;
    	}
    	if (x > 10)          // x > 10
    		flag |= 0x7;
    	else if (x > 5)		 // x = 6 .. 10
    		flag |= 0x6;
    	else
    		flag |= x;   	 // x = 1 .. 5
    	return flag+add;
    }
	
	public FeatureVector getFeatureVector(DependencyInstance inst)
	{
		FeatureVector fv = new FeatureVector(numLinkFeats);
/*		
		int n = inst.length;
		int m = inst.numframes;
		for (SemanticFrame frame : inst.frames) {
			int p = frame.predid;
			int[] args = frame.arglbids;
			for (int a = 0; a < n; ++a) {
				boolean isValid = isValidPredAugPair(inst, p, a);
				if (args[a] >= 0 && isValid) {
					int r = args[a];
					fv.addEntries(createPredArgLinkFeatures(inst, p, a, r));
				}
			}
		}
*/		

		SemanticFrame[] frames = inst.frames;
		int numframes = inst.numframes;
		int len = inst.length;

		boolean[] isPruned = new boolean[numframes * len];
		for (int i = 0; i < isPruned.length; ++i) 
			isPruned[i] = true;
		for (int p = 0; p < numframes; ++p) {
			for (int j = 0; j < len; ++j) {
				if (isValidPredAugPair(inst, frames[p].predid, j)) {
					isPruned[p * len + j] = false;
				}
			}
		}
		
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
				fv.addEntries(createPredArgLinkFeatures(inst, frame.predid, c, r));

				if (options.useSRL2O) {
					// consecutive sib
					if (i < L - 1) {
						fv.addEntries(addSibFeatures(inst, frame.predid, c, childList.get(i + 1), r, frame.arglbids[childList.get(i + 1)]));
					}
					
					// grandparent, head of p
					//for (int j = 0, L2 = gpList.size(); j < L2; ++j) {
					//	fv.addEntries(addGPFeatures(inst, frames[gpList.get(j)].predid, frame.predid, c, frames[gpList.get(j)].arglbids[frame.predid], r));
					//}
					
					// co-predicate, next
					if (p < frames.length - 1) {
						SemanticFrame nextFrame = frames[p + 1];
						if (!isPruned[(p + 1) * len + c] && nextFrame.arglbids[c] != -1) {
							fv.addEntries(addCoPredFeatures(inst, frame.predid, nextFrame.predid, c, r, nextFrame.arglbids[c]));
						}
					}
				}
			}
			
			if (options.useSRLHO && inst.voice[frame.predid] < 2 && inst.voice[frame.predid] >= 0) {
				fv.addEntries(createGlobalFeatures(inst, isPruned, p));
			}
		}

		return fv;
	}
	
	public FeatureVector createPredArgLinkFeatures(DependencyInstance inst,
			int p, int a, int r)
	{
		FeatureVector fv = new FeatureVector(smnAlphabet.size());
		

		// features based on (Johansson 2009)
		addBasicPredArgLinkFeatures(fv, inst, p, a, r);
		
		// features based on (Xavier et. al. 2013)
		addBasicPredArgLinkFeatures2(fv, inst, p, a, r);
				
		return fv;
		
	}

	/**
	 * Predicate-cross-argument features; 
	 * mostly inspired by (Johansson 2009) in EMNLP;
	 * added some additional features. 
	 * @param fv
	 * @param inst
	 * @param p
	 * @param a
	 * @param r
	 */
	private void addBasicPredArgLinkFeatures(FeatureVector fv,
			DependencyInstance inst, int p, int a, int r) {
		
		int[] forms = inst.formids;
		int[] lemmas = inst.lemmaids;
		int[] postags = inst.postagids;
		int[] heads = inst.heads;
		int[] labels = inst.deplbids;
		//int[][] feats = inst.featids;
		
		int pw = forms[p], aw = forms[a];
		int pp = postags[p], ap = postags[a];
		int r2 = r+1;	// r is zero-based; add 1 so r2 \in [1 .. num_of_labels]
		long path = createPathCode(heads, labels, p, a);
		
		long code = createLinkCodeW(PrW, pw);
		addLinkFeature(code, fv);
		addLinkFeature(code | r2, fv);
		
		code = createLinkCodeW(ArW, aw);
		addLinkFeature(code, fv);
		addLinkFeature(code | r2, fv);
		
		code = createLinkCodeP(PrP, pp);
		addLinkFeature(code, fv);
		addLinkFeature(code | r2, fv);
		
		code = createLinkCodeP(ArP, ap);
		addLinkFeature(code, fv);
		addLinkFeature(code | r2, fv);
		
		code = createLinkCodeWW(PrW_ArW, pw, aw);
		addLinkFeature(code, fv);
		addLinkFeature(code | r2, fv);
		
		code = createLinkCodeWP(PrW_ArP, pw, ap);
		addLinkFeature(code, fv);
		addLinkFeature(code | r2, fv);
		
		code = createLinkCodePW(PrP_ArW, pp, aw);
		addLinkFeature(code, fv);
		addLinkFeature(code | r2, fv);
		
		code = createLinkCodePP(PrP_ArP, pp, ap);
		addLinkFeature(code, fv);
		addLinkFeature(code | r2, fv);
	
		code = createLinkCodePATH(PATH, path);
		addLinkFeature(code, fv);
		addLinkFeature(code | r2, fv);
		
		code = createLinkCodePATHW(PATH_PrW, path, pw);
		addLinkFeature(code, fv);
		addLinkFeature(code | r2, fv);
		
		code = createLinkCodePATHP(PATH_PrP, path, pp);
		addLinkFeature(code, fv);
		addLinkFeature(code | r2, fv);
		
		code = createLinkCodePATHW(PATH_ArW, path, aw);
		addLinkFeature(code, fv);
		addLinkFeature(code | r2, fv);
		
		code = createLinkCodePATHP(PATH_ArP, path, ap);
		addLinkFeature(code, fv);
		addLinkFeature(code | r2, fv);
		
		if (lemmas != null) {
			int pl = lemmas[p], al = lemmas[a];
			
			code = createLinkCodeW(PrW, pl);
			addLinkFeature(code, fv);
			addLinkFeature(code | r2, fv);
			
			code = createLinkCodeW(ArW, al);
			addLinkFeature(code, fv);
			addLinkFeature(code | r2, fv);
			
			code = createLinkCodeWW(PrW_ArW, pl, al);
			addLinkFeature(code, fv);
			addLinkFeature(code | r2, fv);
			
			code = createLinkCodeWP(PrW_ArP, pl, ap);
			addLinkFeature(code, fv);
			addLinkFeature(code | r2, fv);
			
			code = createLinkCodePW(PrP_ArW, pp, al);
			addLinkFeature(code, fv);
			addLinkFeature(code | r2, fv);
			
			code = createLinkCodePATHW(PATH_PrW, path, pl);
			addLinkFeature(code, fv);
			addLinkFeature(code | r2, fv);
			
			code = createLinkCodePATHW(PATH_ArW, path, al);
			addLinkFeature(code, fv);
			addLinkFeature(code | r2, fv);
		}
		
		
		int dir = a < p ? 0 : 1;
		int voice = inst.voice[p];
		Utils.Assert(voice <= 2 && voice >= 0);
		int flag = ((dir << 2) | voice) + 1;
		
		code = createLinkCodePP(VOICE_DIR_PrP, flag, pp);
		addLinkFeature(code, fv);
		addLinkFeature(code | r2, fv);
		
		code = createLinkCodePW(VOICE_DIR_PrW, flag, pw);
		addLinkFeature(code, fv);
		addLinkFeature(code | r2, fv);
		
		code = createLinkCodePPP(VOICE_DIR_PrP_ArP, flag, pp, ap);
		addLinkFeature(code, fv);
		addLinkFeature(code | r2, fv);
		
		if (lemmas != null) {
			code = createLinkCodePW(VOICE_DIR_PrW, flag, lemmas[p]);
			addLinkFeature(code, fv);
			addLinkFeature(code | r2, fv);
		}
		
		int pathlen = (getPathLength(inst, p, a) << 1) | dir;
		code = createLinkCodeP(PATHLEN, pathlen);
		addLinkFeature(code, fv);
		addLinkFeature(code | r2, fv);
		
		code = createLinkCodeP(POSITION, dir);
		addLinkFeature(code, fv);
		addLinkFeature(code | r2, fv);
		
		
//		// more trigram features
//		code = createLinkCodePATHPP(PATH_PrP_ArP, path, pp, ap);
//		addLinkFeature(code, fv);
//		addLinkFeature(code | r2, fv);
//		
//		code = createLinkCodePATHWP(PATH_PrW_ArP, path, pw, ap);
//		addLinkFeature(code, fv);
//		addLinkFeature(code | r2, fv);
//		
//		code = createLinkCodePATHWP(PATH_ArW_PrP, path, aw, pp);
//		addLinkFeature(code, fv);
//		addLinkFeature(code | r2, fv);
//		
//		code = createLinkCodePATHWW(PATH_PrW_ArW, path, pw, aw);
//		addLinkFeature(code, fv);
//		addLinkFeature(code | r2, fv);
//		
//		if (lemmas != null) {
//			int pl = lemmas[p], al = lemmas[a];
//			
//			code = createLinkCodePATHWP(PATH_PrW_ArP, path, pl, ap);
//			addLinkFeature(code, fv);
//			addLinkFeature(code | r2, fv);
//			
//			code = createLinkCodePATHWP(PATH_ArW_PrP, path, al, pp);
//			addLinkFeature(code, fv);
//			addLinkFeature(code | r2, fv);
//			
//			code = createLinkCodePATHWW(PATH_PrW_ArW, path, pl, al);
//			addLinkFeature(code, fv);
//			addLinkFeature(code | r2, fv);
//		}
	}
	
	private void addBasicPredArgLinkFeatures2(FeatureVector fv,
			DependencyInstance inst, int p, int a, int r) {
		
		int[] forms = inst.formids;
		int[] lemmas = inst.lemmaids;
		int[] postags = inst.postagids;
		int[] heads = inst.heads;
		int[] labels = inst.deplbids;
		int[][] feats = inst.featids;
		
		int pw = forms[p], aw = forms[a];
		int pp = postags[p], ap = postags[a];
		int r2 = r+1;	// r is zero-based; add 1 so r2 \in [1 .. num_of_labels]
		long path = createPathCode(heads, labels, p, a);
		
		long code;
		
		if (feats != null && feats[p] != null) {
			for (int i = 0, L = feats[p].length; i < L; ++i) {
				code = createLinkCodePP(PrP_ArP, feats[p][i], ap);
				addLinkFeature(code, fv);
				addLinkFeature(code | r2, fv);
				
				code = createLinkCodePW(PrP_ArW, feats[p][i], aw);
				addLinkFeature(code, fv);
				addLinkFeature(code | r2, fv);
				
				code = createLinkCodePATHP(PATH_PrP, path, feats[p][i]);
				addLinkFeature(code, fv);
				addLinkFeature(code | r2, fv);
				
				if (lemmas != null) {
					code = createLinkCodePW(PrP_ArW, feats[p][i], lemmas[a]);
					addLinkFeature(code, fv);
					addLinkFeature(code | r2, fv);
				}
			}
		}
		
//		if (feats != null && feats[a] != null) {
//			for (int i = 0, L = feats[a].length; i < L; ++i) {
//				code = createLinkCodePP(PrP_ArP, pp, feats[a][i]);
//				addLinkFeature(code, fv);
//				addLinkFeature(code | r2, fv);
//				
//				code = createLinkCodePW(PrW_ArP, pw, feats[a][i]);
//				addLinkFeature(code, fv);
//				addLinkFeature(code | r2, fv);
//				
//				code = createLinkCodePATHP(PATH_ArP, path, feats[a][i]);
//				addLinkFeature(code, fv);
//				addLinkFeature(code | r2, fv);
//			}
//		}
		
	}
	
	/**
	 * consecutive argument features; 
	 * mostly inspired by Martins;
	 * added some additional features with arc labels. 
	 */
	public FeatureVector addSibFeatures(DependencyInstance inst, int p, int a1, int a2, int r1, int r2) {
		// p is the index 1..sentence length
		FeatureVector fv = new FeatureVector(smnAlphabet.size());
		
		int[] forms = inst.formids;
		int[] lemmas = inst.lemmaids;
		int[] postags = inst.postagids;
		
		int pw = forms[p], a1w = forms[a1], a2w = forms[a2];
		if (lemmas != null) {
			pw = lemmas[p]; a1w = lemmas[a1]; a2w = lemmas[a2];
		}
		int pp = postags[p], a1p = postags[a1], a2p = postags[a2];
		int dir = (((a1 < p ? 1 : 0) << 1) | (a2 < p ? 1 : 0)) + 1;
		int r = ((r1 + 1) << auglbNumBits) | (r2 + 1);
		
		long code = 0;
		code = createLinkCodePPP(PP_A1P_A2P, pp, a1p, a2p);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePPPP(PP_A1P_A2P_DIR, pp, a1p, a2p, dir);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePPW(PW_A1P_A2P, a1p, a2p, pw);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePPPW(PW_A1P_A2P_DIR, a1p, a2p, dir, pw);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePPW(PP_A1W_A2P, pp, a2p, a1w);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePPPW(PP_A1W_A2P_DIR, pp, a2p, dir, a1w);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePPW(PP_A1P_A2W, pp, a1p, a2w);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePPPW(PP_A1P_A2W_DIR, pp, a1p, dir, a2w);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePP(A1P_A2P, a1p, a2p);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePPP(A1P_A2P_DIR, a1p, a2p, dir);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePW(A1W_A2P, a2p, a1w);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePPW(A1W_A2P_DIR, a2p, dir, a1w);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePW(A1P_A2W, a1p, a2w);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePPW(A1P_A2W_DIR, a1p, dir, a2w);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		//code = createLinkCodeWW(A1W_A2W, a1w, a2w);
		//addLinkFeature(code, fv);
		//addLinkFeature(code | r, fv);
		
		//code = createLinkCodePWW(A1W_A2W_DIR, dir, a1w, a2w);
		//addLinkFeature(code, fv);
		//addLinkFeature(code | r, fv);
		
		//code = createLinkCodePWW(PW_A1W_A2P, a2p, pw, a1w);
		//addLinkFeature(code, fv);
		//addLinkFeature(code | r, fv);
		
		//code = createLinkCodePWW(PP_A1W_A2W, pp, a1w, a2w);
		//addLinkFeature(code, fv);
		//addLinkFeature(code | r, fv);
		
		//code = createLinkCodePWW(PW_A1P_A2W, a1p, pw, a2w);
		//addLinkFeature(code, fv);
		//addLinkFeature(code | r, fv);
		
		return fv;
	}
	
	/**
	 * grandparent features; 
	 * mostly inspired by Martins;
	 * added some additional features with arc labels. 
	 */
	public FeatureVector addGPFeatures(DependencyInstance inst, int gp, int p, int a, int r1, int r2) {
		// p is the index 1..sentence length
		FeatureVector fv = new FeatureVector(smnAlphabet.size());
		
		int[] forms = inst.formids;
		int[] lemmas = inst.lemmaids;
		int[] postags = inst.postagids;
		
		int gpw = forms[gp], pw = forms[p], aw = forms[a];
		if (lemmas != null) {
			gpw = lemmas[gp]; pw = lemmas[p]; aw = lemmas[a];
		}
		int gpp = postags[gp], pp = postags[p], ap = postags[a];
		int dir = (((p < gp ? 1 : 0) << 1) | (a < p ? 1 : 0)) + 1;
		int r = ((r1 + 1) << auglbNumBits) | (r2 + 1);
		
		long code = 0;
		code = createLinkCodePPP(GPP_PP_AP, gpp, pp, ap);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePPPP(GPP_PP_AP_DIR, gpp, pp, ap, dir);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePPW(GPW_PP_AP, pp, ap, gpw);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePPPW(GPW_PP_AP_DIR, pp, ap, dir, gpw);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePPW(GPP_PW_AP, gpp, ap, pw);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePPPW(GPP_PW_AP_DIR, gpp, ap, dir, pw);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePPW(GPP_PP_AW, gpp, pp, aw);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePPPW(GPP_PP_AW_DIR, gpp, pp, dir, aw);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePP(GPP_AP, gpp, ap);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePPP(GPP_AP_DIR, gpp, ap, dir);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePW(GPW_AP, ap, gpw);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePPW(GPW_AP_DIR, ap, dir, gpw);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePW(GPP_AW, gpp, aw);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePPW(GPP_AW_DIR, gpp, dir, aw);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		//code = createLinkCodeWW(GPW_AW, gpw, aw);
		//addLinkFeature(code, fv);
		//addLinkFeature(code | r, fv);
		
		//code = createLinkCodePWW(GPW_AW_DIR, dir, gpw, aw);
		//addLinkFeature(code, fv);
		//addLinkFeature(code | r, fv);
		
		//code = createLinkCodePWW(GPW_PW_AP, ap, gpw, pw);
		//addLinkFeature(code, fv);
		//addLinkFeature(code | r, fv);
		
		//code = createLinkCodePWW(GPP_PW_AW, gpp, pw, aw);
		//addLinkFeature(code, fv);
		//addLinkFeature(code | r, fv);
		
		//code = createLinkCodePWW(GPW_PP_AW, pp, gpw, aw);
		//addLinkFeature(code, fv);
	
		return fv;
	}
	
	/**
	 * co-predicate features; 
	 * mostly inspired by Martins;
	 * added some additional features with arc labels. 
	 */
	public FeatureVector addCoPredFeatures(DependencyInstance inst, int p, int p2, int a, int r1, int r2) {
		// p is the index 1..sentence length
		FeatureVector fv = new FeatureVector(smnAlphabet.size());
		
		int[] forms = inst.formids;
		int[] lemmas = inst.lemmaids;
		int[] postags = inst.postagids;
		
		int pw = forms[p], p2w = forms[p2], aw = forms[a];
		if (lemmas != null) {
			pw = lemmas[p]; p2w = lemmas[p2]; aw = lemmas[a];
		}
		int pp = postags[p], p2p = postags[p2], ap = postags[a];
		int dir = (((a < p ? 1 : 0) << 1) | (a < p2 ? 1 : 0)) + 1;
		int r = ((r1 + 1) << auglbNumBits) | (r2 + 1);
		
		long code = 0;
		code = createLinkCodePPP(PP_P2P_AP, pp, p2p, ap);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePPPP(PP_P2P_AP_DIR, pp, p2p, ap, dir);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePPW(PW_P2P_AP, p2p, ap, pw);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePPPW(PW_P2P_AP_DIR, p2p, ap, dir, pw);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePPW(PP_P2W_AP, pp, ap, p2w);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePPPW(PP_P2W_AP_DIR, pp, ap, dir, p2w);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePPW(PP_P2P_AW, pp, p2p, aw);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		code = createLinkCodePPPW(PP_P2P_AW_DIR, pp, p2p, dir, aw);
		addLinkFeature(code, fv);
		addLinkFeature(code | r, fv);
		
		//code = createLinkCodePWW(PW_P2W_AP, ap, pw, p2w);
		//addLinkFeature(code, fv);
		//addLinkFeature(code | r, fv);
		
		//code = createLinkCodePWW(PP_P2W_AW, pp, p2w, aw);
		//addLinkFeature(code, fv);
		//addLinkFeature(code | r, fv);
		
		//code = createLinkCodePWW(PW_P2P_AW, p2p, pw, aw);
		//addLinkFeature(code, fv);
		//addLinkFeature(code | r, fv);
	
		return fv;
	}
	
	/**
	 * global features; 
	 * mostly inspired by Totounova;
	 */
	public FeatureVector createGlobalFeatures(DependencyInstance inst, boolean[] isPruned, int p) {
		// p is the index 1..sentence length
		FeatureVector fv = new FeatureVector(smnAlphabet.size());
		
		long code = 0;
		
		SemanticFrame frame = inst.frames[p];
		int pid = frame.predid;
		Utils.Assert(inst.voice[pid] < 2 && inst.voice[pid] >= 0);

		//code = createLinkCodeArgSeq(p, inst, false);
		//addLinkFeature(code, fv);

		//code = createLinkCodeArgSeq(p, inst, true);
		//addLinkFeature(code, fv);
		
		long argseq = createArgSequenceCode(frame.arglbids, isPruned, p, pid);
		code = createLinkCodeARGSEQP(VOICE_ARGSEQ, argseq, inst.voice[pid]);
		addLinkFeature(code, fv);

		int word = inst.formids[pid];
		if (inst.lemmaids != null)
			word = inst.lemmaids[pid];
		
		code = createLinkCodeARGSEQWP(VOICE_LEMMA_ARGSEQ, argseq, word, inst.voice[pid]);
		addLinkFeature(code, fv);
		
		argseq = createArgPosSequenceCode(frame.arglbids, inst.postagids, isPruned, p, pid);
		code = createLinkCodeARGSEQP(VOICE_POSARGSEQ, argseq, inst.voice[pid]);
		addLinkFeature(code, fv);
		
		code = createLinkCodeARGSEQWP(VOICE_LEMMA_POSARGSEQ, argseq, word, inst.voice[pid]);
		addLinkFeature(code, fv);
		
		argseq = createArgPredSequenceCode(frame.arglbids, inst.predIndex, inst.voice, isPruned, p, pid);
		code = createLinkCodeARGSEQP(VOICE_POSARGPREDSEQ, argseq, inst.voice[pid]);
		addLinkFeature(code, fv);
		
		return fv;
	}
	
	public static boolean isValidPredAugPair(DependencyInstance inst, int p, int a)
	{
		int[] heads = inst.heads;
		int pa = heads[a];
		//int gpa = pa >= 0 ? heads[pa] : pa;
		for (int i = 0, x = p; i <= MAX_DEPTH && x >= 0; ++i) {
			if (x == a || x == pa /*|| x == gpa*/) return true;
			x = heads[x];
		}
		return false;
	}
	
	public static String getPathString(DependencyInstance inst, int p, int a) 
	{
		String code = "";
		int[] heads = inst.heads;
		String[] deps = inst.deprels;
		int pa = heads[a];
		//int gpa = pa >= 0 ? heads[pa] : pa;
		for (int i = 0, x = p; i <= MAX_DEPTH && x >= 0; ++i) {
			if (x == a) break;
			if (x == pa) {
				code += "[" + deps[a] + "|d]";
				break;
			}
			code += "[" + deps[x] + "|u]";
			x = heads[x];
		}		
		return code;
	}
	
	public static int getPathLength(DependencyInstance inst, int p, int a) 
	{
		int length = 0;
		int[] heads = inst.heads;
		//String[] deps = inst.deprels;
		int pa = heads[a];
		//int gpa = pa >= 0 ? heads[pa] : pa;
		for (int i = 0, x = p; i <= MAX_DEPTH && x >= 0; ++i) {
			if (x == a) break;
			if (x == pa) {
				++length;
				break;
			}
			++length;
			x = heads[x];
		}		
		return length;
	}
	
	public void addLinkFeature(long code, FeatureVector fv)
	{
		int id = smnAlphabet.lookupIndex(code, numLinkFeats);
		if (id >= 0) {
			fv.addEntry(id, 1.0);
			if (id == numLinkFeats) ++numLinkFeats;
		}
	}
	
	public void addPathFeature(long code, FeatureVector fv)
	{
		int id = pathAlphabet.lookupIndex(code, numPathFeats);
		if (id >= 0) {
			fv.addEntry(id, 1.0);
			if (id == numPathFeats) ++numPathFeats;
		}
	}
	
	public void addContextFeature(long code, FeatureVector fv)
	{
		int id = contextAlphabet.lookupIndex(code, numContextFeats);
		if (id >= 0) {
			fv.addEntry(id, 1.0);
			if (id == numContextFeats) ++numContextFeats;
		}
	}
	
    public void addWordFeature(long code, FeatureVector mat) {
    	int id = wordAlphabet.lookupIndex(code, numWordFeats);
    	if (id >= 0) {
    		mat.addEntry(id, 1.0);
    		if (id == numWordFeats) ++numWordFeats;
    	}
    }
    
    public void addWordFeature(long code, double value, FeatureVector mat) {
    	int id = wordAlphabet.lookupIndex(code, numWordFeats);
    	if (id >= 0) {
    		mat.addEntry(id, value);
    		if (id == numWordFeats) ++numWordFeats;
    	}
    }
    
	public int createPathCode(int[] heads, int[] deplbs, int p , int a)
	{
		long code = 0;
		int pa = heads[a];
		//int gpa = pa >= 0 ? heads[pa] : pa;
		for (int i = 0, x = p; i <= MAX_DEPTH && x >= 0; ++i) {
			if (x == a) break;
			if (x == pa) {
				code = ((code << deplbNumBits) | (deplbs[a]+1)) << 1;
				break;
			}
			//if (x == gpa) {
			//	code = ((((code << deplbNumBits) | (deplbs[a]+1)) << (deplbNumBits+1)) | (deplbs[pa]+1)) << 1;
			//	break;
			//}
			code = (((code << deplbNumBits) | (deplbs[x]+1)) << 1) | 1;
			x = heads[x];
		}
		int indexCode = pathcodeAlphabet.lookupIndex(code) + 1;
		return indexCode;
	}
	
	public int createArgSequenceCode(int[] arglbids, boolean[] isPruned, int p, int pid) {
		long code = 0;
		//String s = "";
		for (int i = 0; i < arglbids.length; ++i) {
			if (i == pid) {
				code = (code << auglbNumBits) | (DictionarySet.ARGINDEX);
				//s += "PRED "; 
			}

			if (arglbids[i] == -1 || arglbids[i] >= DictionarySet.ARGINDEX - 1 || isPruned[p * arglbids.length + i])
				continue;
			
			code = (code << auglbNumBits) | (arglbids[i] + 1);
			//s += args[arglbids[i]] + " ";
		}
		int indexCode = argSeqAlphabet.lookupIndex(code) + 1;
		
		//System.out.println(s + "\t" + indexCode);
		return indexCode;
	}
	
	public int createArgPredSequenceCode(int[] arglbids, int[] predIndex, int[] voice, boolean[] isPruned, int p, int pid) {
		long code = 0;
		for (int i = 0; i < arglbids.length; ++i) {
			if (i == pid) {
				code = (code << auglbNumBits) | (DictionarySet.ARGINDEX);
			}
			else if (predIndex[i] >= 0) {
				Utils.Assert(voice[i] >= 0 && voice[i] <= 2);
				code = (code << auglbNumBits) | (DictionarySet.ARGINDEX + (voice[i] == 2 ? 1 : 2));
			}

			if (arglbids[i] == -1 || arglbids[i] >= DictionarySet.ARGINDEX - 1 || isPruned[p * arglbids.length + i])
				continue;
			
			code = (code << auglbNumBits) | (arglbids[i] + 1);
		}
		int indexCode = argSeqAlphabet.lookupIndex(code) + 1;
		return indexCode;
	}
	
	public int createArgPosSequenceCode(int[] arglbids, int[] posid, boolean[] isPruned, int p, int pid) {
		long code = 0;
		//String s = "";
		for (int i = 0; i < arglbids.length; ++i) {
			if (i == pid) {
				code = (code << auglbNumBits) | (DictionarySet.ARGINDEX);
				code = (code << tagNumBits) | posid[pid];
			}

			if (arglbids[i] == -1 || arglbids[i] >= DictionarySet.ARGINDEX - 1 || isPruned[p * arglbids.length + i])
				continue;
			
			code = (code << auglbNumBits) | (arglbids[i] + 1);
			code = (code << tagNumBits) | posid[i];
		}
		int indexCode = argSeqAlphabet.lookupIndex(code) + 1;
		return indexCode;
	}
	
	public long createContextCodeW(FeatureTemplate.Context temp, long x)
	{
		return (x << numContextFeatBits) | temp.ordinal();
	}
	
	public long createPathCodeW(FeatureTemplate.Path temp, long x)
	{
		return (x << numPathFeatBits) | temp.ordinal();
	}
	
	public long createLinkCodeW(FeatureTemplate.Link temp, long x)
	{
		return ((x << numLinkFeatBits) | temp.ordinal()) << (auglbNumBits * 2);
	}
	
	public long createLinkCodeP(FeatureTemplate.Link temp, long x)
	{
		return ((x << numLinkFeatBits) | temp.ordinal()) << (auglbNumBits * 2);
	}
	
	public long createLinkCodeWW(FeatureTemplate.Link temp, long x, long y)
	{
		return ((((x << wordNumBits) | y) << numLinkFeatBits) | temp.ordinal()) << (auglbNumBits * 2);
	}

	public long createLinkCodeWP(FeatureTemplate.Link temp, long x, long y)
	{
		return ((((x << tagNumBits) | y) << numLinkFeatBits) | temp.ordinal()) << (auglbNumBits * 2);
	}
	
	public long createLinkCodePW(FeatureTemplate.Link temp, long x, long y)
	{
		return ((((x << wordNumBits) | y) << numLinkFeatBits) | temp.ordinal()) << (auglbNumBits * 2);
	}
	
	public long createLinkCodePP(FeatureTemplate.Link temp, long x, long y)
	{
		return ((((x << tagNumBits) | y) << numLinkFeatBits) | temp.ordinal()) << (auglbNumBits * 2);
	}
	
	public long createLinkCodePPP(FeatureTemplate.Link temp, long x, long y, long z)
	{
		return ((((((x << tagNumBits) | y) << tagNumBits) | z ) << numLinkFeatBits) | temp.ordinal()) << (auglbNumBits * 2);
	}
	
	public long createLinkCodePPW(FeatureTemplate.Link temp, long x, long y, long z)
	{
		return ((((((x << tagNumBits) | y) << wordNumBits) | z ) << numLinkFeatBits) | temp.ordinal()) << (auglbNumBits * 2);
	}
	
	public long createLinkCodePWW(FeatureTemplate.Link temp, long x, long y, long z)
	{
		return ((((((x << wordNumBits) | y) << wordNumBits) | z ) << numLinkFeatBits) | temp.ordinal()) << (auglbNumBits * 2);
	}
	
	public long createLinkCodePPPP(FeatureTemplate.Link temp, long x, long y, long z, long w)
	{
		return ((((((((x << tagNumBits) | y) << tagNumBits) | z ) << tagNumBits) | w) << numLinkFeatBits) | temp.ordinal()) << (auglbNumBits * 2);
	}
	
	public long createLinkCodePPPW(FeatureTemplate.Link temp, long x, long y, long z, long w)
	{
		return ((((((((x << tagNumBits) | y) << tagNumBits) | z ) << wordNumBits) | w) << numLinkFeatBits) | temp.ordinal()) << (auglbNumBits * 2);
	}
	
	public long createLinkCodePATH(FeatureTemplate.Link temp, long x)
	{
		return ((x << numLinkFeatBits) | temp.ordinal()) << (auglbNumBits * 2);
	}
	
	public long createLinkCodePATHW(FeatureTemplate.Link temp, long x, long y)
	{
		return ((((x << wordNumBits) | y) << numLinkFeatBits) | temp.ordinal()) << (auglbNumBits * 2);
	}
	
	public long createLinkCodePATHP(FeatureTemplate.Link temp, long x, long y)
	{
		return ((((x << tagNumBits) | y) << numLinkFeatBits) | temp.ordinal()) << (auglbNumBits * 2);
	}
	
	public long createLinkCodePATHPP(FeatureTemplate.Link temp, long x, long y, long z)
	{
		return ((((((x << tagNumBits) | y) << tagNumBits) | z) << numLinkFeatBits) | temp.ordinal()) << (auglbNumBits * 2);
	}
	
	public long createLinkCodePATHWP(FeatureTemplate.Link temp, long x, long y, long z)
	{
		return ((((((x << wordNumBits) | y) << tagNumBits) | z) << numLinkFeatBits) | temp.ordinal()) << (auglbNumBits * 2);
	}
	
	public long createLinkCodePATHWW(FeatureTemplate.Link temp, long x, long y, long z)
	{
		return ((((((x << wordNumBits) | y) << wordNumBits) | z) << numLinkFeatBits) | temp.ordinal()) << (auglbNumBits * 2);
	}
	
	public long createLinkCodeARGSEQP(FeatureTemplate.Link temp, long x, long y)
	{
		return ((((x << tagNumBits) | y) << numLinkFeatBits) | temp.ordinal()) << (auglbNumBits * 2);
	}
	
	public long createLinkCodeARGSEQWP(FeatureTemplate.Link temp, long x, long y, long z)
	{
		return ((((((x << wordNumBits) | y) << tagNumBits) | z) << numLinkFeatBits) | temp.ordinal()) << (auglbNumBits * 2);
	}
	
	public long createWordCodeW(FeatureTemplate.Word temp, long x) {
		return ((x << numWordFeatBits) | temp.ordinal()) << 4;
	}

	public long createWordCodeP(FeatureTemplate.Word temp, long x) {
		return ((x << numWordFeatBits) | temp.ordinal()) << 4;
	}

	public long createWordCodePP(FeatureTemplate.Word temp, long x, long y) {
		return ((((x << tagNumBits) | y) << numWordFeatBits) | temp.ordinal()) << 4;
	}

	public long createWordCodePPP(FeatureTemplate.Word temp, long x, long y, long z) {
		return ((((((x << tagNumBits) | y) << tagNumBits) | z) << numWordFeatBits)
				| temp.ordinal()) << 4;
	}

	public long createWordCodeWP(FeatureTemplate.Word temp, long x, long y) {
		return ((((x << tagNumBits) | y) << numWordFeatBits) | temp.ordinal()) << 4;
	}
	
	public long extractLinkTemplateCode(long code) {
		return (code >> (auglbNumBits*2)) & ((1 << numLinkFeatBits)-1);
	}
	
	public long extractRoleCode(long code) {
		return code & ((1 << auglbNumBits)-1);
	}
	
	public void extractLinkCodeW(long code, int[] parts) {
		code = code >> (auglbNumBits*2) >> numLinkFeatBits;
		parts[0] = (int) (code & ((1 << wordNumBits)-1));
	}
	
	public void extractLinkCodeP(long code, int[] parts) {
		code = code >> (auglbNumBits*2) >> numLinkFeatBits;
		parts[0] = (int) (code & ((1 << tagNumBits)-1));
	}
	
	public void extractLinkCodePATH(long code, int[] parts) {
		code = code >> (auglbNumBits*2) >> numLinkFeatBits;
		parts[0] = (int) code;
	}
	
	public void extractLinkCodePATHW(long code, int[] parts) {
		code = code >> (auglbNumBits*2) >> numLinkFeatBits;
		parts[1] = (int) (code & ((1 << wordNumBits)-1));
		code = code >> wordNumBits;
		parts[0] = (int) code;
	}
	
	public void extractLinkCodePATHP(long code, int[] parts) {
		code = code >> (auglbNumBits*2) >> numLinkFeatBits;
		parts[1] = (int) (code & ((1 << tagNumBits)-1));
		code = code >> tagNumBits;
		parts[0] = (int) code;
	}
	
	public void extractLinkCodeWW(long code, int[] parts) {
		code = code >> (auglbNumBits*2) >> numLinkFeatBits;
		parts[1] = (int) (code & ((1 << wordNumBits)-1));
		code = code >> wordNumBits;
		parts[0] = (int) (code & ((1 << wordNumBits)-1));
	}
	
	public void extractLinkCodeWP(long code, int[] parts) {
		code = code >> (auglbNumBits*2) >> numLinkFeatBits;
		parts[1] = (int) (code & ((1 << tagNumBits)-1));
		code = code >> tagNumBits;
		parts[0] = (int) (code & ((1 << wordNumBits)-1));
	}
	
	public void extractLinkCodePP(long code, int[] parts) {
		code = code >> (auglbNumBits*2) >> numLinkFeatBits;
		parts[1] = (int) (code & ((1 << tagNumBits)-1));
		code = code >> tagNumBits;
		parts[0] = (int) (code & ((1 << tagNumBits)-1));
	}
	
	public void extractLinkCodePW(long code, int[] parts) {
		code = code >> (auglbNumBits*2) >> numLinkFeatBits;
		parts[1] = (int) (code & ((1 << wordNumBits)-1));
		code = code >> wordNumBits;
		parts[0] = (int) (code & ((1 << tagNumBits)-1));
	}
	
	public void fillParameters(SemanticLowRankParam tensor, Parameters params) {
		
		long[] codes = smnAlphabet.toArray();
		
		int[] parts = new int[4];
		
		for (long code : codes) {
			int id = smnAlphabet.lookupIndex(code);
			if (id < 0) continue;
			
			int temp = (int) extractLinkTemplateCode(code);
			int role = (int) extractRoleCode(code);
			
//			BIAS,
//			PrW,
//			PrP,
//			ArW,
//			ArP,
//			PrW_ArW,
//			PrP_ArP,
//			PrW_ArP,
//			PrP_ArW,
//			PATH,
//			PATH_PrW,
//			PATH_PrP,
//			PATH_ArW,
//			PATH_ArP,
			
			long codex = 0, codey = 0, codez = 0, coder = 0;
			int x = 0, y = 0, z = 0, r = 0;
		
			if (temp == PrW.ordinal()) {
				extractLinkCodeW(code, parts);
				codex = createWordCodeW(WORDFV_W0, parts[0]);
				codey = createWordCodeP(WORDFV_BIAS, 0);
				codez = createPathCodeW(PATHFV_BIAS, 0);
			}
			else if (temp == PrP.ordinal()) {
				extractLinkCodeP(code, parts);
				codex = createWordCodeP(WORDFV_P0, parts[0]);
				codey = createWordCodeP(WORDFV_BIAS, 0);
				codez = createPathCodeW(PATHFV_BIAS, 0);
			}
			else if (temp == ArW.ordinal()) {
				extractLinkCodeW(code, parts);
				codey = createWordCodeW(WORDFV_P0, parts[0]);
				codex = createWordCodeP(WORDFV_BIAS, 0);
				codez = createPathCodeW(PATHFV_BIAS, 0);
			}
			else if (temp == ArP.ordinal()) {
				extractLinkCodeP(code, parts);
				codey = createWordCodeP(WORDFV_P0, parts[0]);
				codex = createWordCodeP(WORDFV_BIAS, 0);
				codez = createPathCodeW(PATHFV_BIAS, 0);
			}
			else if (temp == PATH.ordinal()) {
				extractLinkCodePATH(code, parts);
				codex = createWordCodeP(WORDFV_BIAS, 0);
				codey = createWordCodeP(WORDFV_BIAS, 0);
				codez = createPathCodeW(PATHFV_PATHCODE, parts[0]);
			}
			else if (temp == PATH_PrW.ordinal()) {
				extractLinkCodePATHW(code, parts);
				codex = createWordCodeW(WORDFV_W0, parts[1]);
				codey = createWordCodeP(WORDFV_BIAS, 0);
				codez = createPathCodeW(PATHFV_PATHCODE, parts[0]);
			}
			else if (temp == PATH_PrP.ordinal()) {
				extractLinkCodePATHP(code, parts);
				codex = createWordCodeP(WORDFV_P0, parts[1]);
				codey = createWordCodeP(WORDFV_BIAS, 0);
				codez = createPathCodeW(PATHFV_PATHCODE, parts[0]);
			}
			else if (temp == PATH_ArW.ordinal()) {
				extractLinkCodePATHW(code, parts);
				codey = createWordCodeW(WORDFV_W0, parts[1]);
				codex = createWordCodeP(WORDFV_BIAS, 0);
				codez = createPathCodeW(PATHFV_PATHCODE, parts[0]);
			}
			else if (temp == PATH_ArP.ordinal()) {
				extractLinkCodePATHP(code, parts);
				codey = createWordCodeP(WORDFV_P0, parts[1]);
				codex = createWordCodeP(WORDFV_BIAS, 0);
				codez = createPathCodeW(PATHFV_PATHCODE, parts[0]);
			}
			else if (temp == PrW_ArW.ordinal()) {
				extractLinkCodeWW(code, parts);
				codex = createWordCodeW(WORDFV_W0, parts[0]);
				codey = createWordCodeW(WORDFV_W0, parts[1]);
				codez = createPathCodeW(PATHFV_BIAS, 0);
			}
			else if (temp == PrW_ArP.ordinal()) {
				extractLinkCodeWP(code, parts);
				codex = createWordCodeW(WORDFV_W0, parts[0]);
				codey = createWordCodeP(WORDFV_P0, parts[1]);
				codez = createPathCodeW(PATHFV_BIAS, 0);
			}
			else if (temp == PrP_ArP.ordinal()) {
				extractLinkCodePP(code, parts);
				codex = createWordCodeP(WORDFV_P0, parts[0]);
				codey = createWordCodeP(WORDFV_P0, parts[1]);
				codez = createPathCodeW(PATHFV_BIAS, 0);
			}
			else if (temp == PrP_ArW.ordinal()) {
				extractLinkCodePW(code, parts);
				codex = createWordCodeP(WORDFV_P0, parts[0]);
				codey = createWordCodeW(WORDFV_W0, parts[1]);
				codez = createPathCodeW(PATHFV_BIAS, 0);
			}
			
			if (role > 0)
				coder = createContextCodeW(CNTFV_LABEL, role-1);
			else
				coder = createContextCodeW(CNTFV_BIAS, 0);
			
			x = wordAlphabet.lookupIndex(codex);
			y = wordAlphabet.lookupIndex(codey);
			z = pathAlphabet.lookupIndex(codez);
			r = contextAlphabet.lookupIndex(coder);
			if (x >= 0 && y >= 0 && z >= 0 && r >= 0) {
				double value = params.params2[id];
				tensor.add(x,y,z,r,value);
			}

		}
	}
}
