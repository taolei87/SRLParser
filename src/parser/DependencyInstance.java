package parser;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.Serializable;
import java.util.*;
import java.io.*;

import parser.Options.PossibleLang;
import utils.Alphabet;
import utils.Dictionary;
import utils.DictionarySet;

import static utils.DictionarySet.DictionaryTypes.*;

public class DependencyInstance implements Serializable {
	
	public enum SpecialPos {
		C, P, PNX, V, AJ, N, OTHER,
	}
	
	private static final long serialVersionUID = 1L;
	
	public int length;

	// FORM: the forms - usually words, like "thought"
	public String[] forms;

	// LEMMA: the lemmas, or stems, e.g. "think"
	public String[] lemmas;
	
	// COARSE-POS: the coarse part-of-speech tags, e.g."V"
	public String[] cpostags;

	// FINE-POS: the fine-grained part-of-speech tags, e.g."VBD"
	public String[] postags;
	
	// MOST-COARSE-POS: the coarsest part-of-speech tags (about 11 in total)
	public SpecialPos[] specialPos;
	
	// FEATURES: some features associated with the elements separated by "|", e.g. "PAST|3P"
	public String[][] feats;

	// HEAD: the IDs of the heads for each element
	public int[] heads;

	// DEPREL: the dependency relations, e.g. "SUBJ"
	public String[] deprels;
	
	public int[] formids;
	public int[] lemmaids;
	public int[] postagids;
	public int[] cpostagids;
	public int[][] featids;
	public int[] deplbids;
	public int[] wordVecIds;
	
	
	// semantic frames and roles
	public int numframes;
	public SemanticFrame[] frames;
	public int[] predIndex;
	public int[] voice;

    public DependencyInstance() {}
    
    public DependencyInstance(int length) { this.length = length; }
    
    public DependencyInstance(String[] forms) {
    	length = forms.length;
    	this.forms = forms;
    	this.feats = new String[length][];
    	this.deprels = new String[length];
    }
    
    public DependencyInstance(String[] forms, String[] postags, int[] heads) {
    	this.length = forms.length;
    	this.forms = forms;    	
    	this.heads = heads;
	    this.postags = postags;
    }
    
    public DependencyInstance(String[] forms, String[] postags, int[] heads, String[] deprels) {
    	this(forms, postags, heads);
    	this.deprels = deprels;    	
    }

    public DependencyInstance(String[] forms, String[] lemmas, String[] cpostags, String[] postags,
            String[][] feats, int[] heads, String[] deprels, SemanticFrame[] frames, int[] predIndex, int[] voice) {
    	this(forms, postags, heads, deprels);
    	this.lemmas = lemmas;    	
    	this.feats = feats;
    	this.cpostags = cpostags;
    	this.frames = frames;
    	this.predIndex = predIndex;
    	this.voice = voice;
    	this.numframes = frames.length;
    }
    
    public DependencyInstance(DependencyInstance a) {
    	
    	// only store token ids (not strings) to save memory
    	
    	specialPos = a.specialPos;
    	length = a.length;
    	heads = a.heads;
    	formids = a.formids;
    	lemmaids = a.lemmaids;
    	postagids = a.postagids;
    	cpostagids = a.cpostagids;
    	deplbids = a.deplbids;
    	featids = a.featids;
    	wordVecIds = a.wordVecIds;
    	predIndex = a.predIndex;
    	voice = a.voice;
   	
    	numframes = a.numframes;
    	frames = a.frames;
    	if (numframes > 0 && frames[0].arglbs != null) {
    		frames = new SemanticFrame[numframes];
    		for (int i = 0; i < numframes; ++i)
    			frames[i] = new SemanticFrame(a.frames[i]);
    	}
    	
    }
    
    public void setInstIds(DictionarySet dicts, 
    		HashMap<String, String> coarseMap, HashSet<String> conjWord, PossibleLang lang) {
    	    	
    	formids = new int[length];    	
		deplbids = new int[length];
		postagids = new int[length];
		cpostagids = new int[length];
		
    	for (int i = 0; i < length; ++i) {
    		formids[i] = dicts.lookupIndex(WORD, "form="+forms[i]);
			postagids[i] = dicts.lookupIndex(POS, "pos="+postags[i]);
			cpostagids[i] = dicts.lookupIndex(POS, "cpos="+cpostags[i]);
			deplbids[i] = dicts.lookupIndex(DEPLABEL, deprels[i]) - 1;	// zero-based
    	}
    	
    	if (lemmas != null) {
    		lemmaids = new int[length];
    		for (int i = 0; i < length; ++i)
    			lemmaids[i] = dicts.lookupIndex(WORD, "lemma="+lemmas[i]);
    	}

		featids = new int[length][];
		for (int i = 0; i < length; ++i) if (feats[i] != null) {
			featids[i] = new int[feats[i].length];
			for (int j = 0; j < feats[i].length; ++j)
				featids[i][j] = dicts.lookupIndex(POS, "feat="+feats[i][j]);
		}
		
		if (frames != null) {
			for (int i = 0; i < numframes; ++i)
				for (int j = 0; j < length; ++j)
					if (frames[i].arglbs[j] != null)
						frames[i].arglbids[j] = dicts.lookupIndex(AUGLABEL, frames[i].arglbs[j]) - 1; // zero-based
		}
		
		if (dicts.size(WORDVEC) > 0) {
			wordVecIds = new int[length];
			for (int i = 0; i < length; ++i) {
				int wvid = dicts.lookupIndex(WORDVEC, forms[i]);
				if (wvid <= 0) wvid = dicts.lookupIndex(WORDVEC, forms[i].toLowerCase());
				if (wvid > 0) wordVecIds[i] = wvid; else wordVecIds[i] = -1; 
			}
		}
		
		// set special pos
		specialPos = new SpecialPos[length];
		for (int i = 0; i < length; ++i) {
			if (coarseMap.containsKey(postags[i])) {
				String cpos = coarseMap.get(postags[i]);
				if ((cpos.equals("CONJ")
						|| PossibleLang.Japanese == lang) && conjWord.contains(forms[i])) {
					specialPos[i] = SpecialPos.C;
				}
				else if (cpos.equals("ADP"))
					specialPos[i] = SpecialPos.P;
				else if (cpos.equals("."))
					specialPos[i] = SpecialPos.PNX;
				else if (cpos.equals("VERB"))
					specialPos[i] = SpecialPos.V;
				else
					specialPos[i] = SpecialPos.OTHER;
			}
			else {
				//System.out.println("Can't find coarse map: " + postags[i]);
				coarseMap.put(postags[i], "X");
			}
		}
    }
}


