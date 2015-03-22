package parser.io;

import java.io.IOException;
import java.util.ArrayList;

import parser.DependencyInstance;
import parser.SemanticFrame;
import parser.Options;
import utils.Utils;


public class CONLLReader extends DependencyReader {

	public CONLLReader(Options options) {
		this.options = options;
	}
	
	@Override
	public DependencyInstance nextInstance() throws IOException {
		
	    ArrayList<String[]> lstLines = new ArrayList<String[]>();

	    String line = reader.readLine();
	    while (line != null && !line.equals("") && !line.startsWith("*")) {
	    	lstLines.add(line.trim().split("\t"));
	    	line = reader.readLine();
	    }
	    
	    if (lstLines.size() == 0) return null;
	    
	    int length = lstLines.size();
	    String[] forms = new String[length + 1];
	    String[] lemmas = new String[length + 1];
	    String[] cpos = new String[length + 1];
	    String[] pos = new String[length + 1];
	    String[][] feats = new String[length + 1][];
	    String[] deprels = new String[length + 1];
	    int[] heads = new int[length + 1];
	    int[] predIndex = new int[length + 1];
	    int[] voice = new int[length + 1];
	    
	    forms[0] = "<root>";
	    lemmas[0] = "<root-LEMMA>";
	    cpos[0] = "<root-CPOS>";
	    pos[0] = "<root-POS>";
	    deprels[0] = "<no-type>";
	    heads[0] = -1;
	    for (int i = 0; i < predIndex.length; ++i) {
	    	predIndex[i] = -1;
	    	voice[i] = -1;
	    }
	    
	    boolean hasLemma = false;
	    
	    
	    /*
		    0 ID
		    1 FORM
		    2 LEMMA (not used)
		    3 PLEMMA 
		    4 POS (not used)
		    5 PPOS   
		    6 FEAT (not used)
		    7 PFEAT  
		    8 HEAD
		    9 PHEAD
		    10 DEPREL 
		    11 PDEPREL 
		    12 FILLPRED 
		    13 PRED
		    14... APREDn
	   	*/
	    
	    // 11  points  point   point   NNS NNS _   _   8   8   PMOD    PMOD    Y   point.02    _   _   _   _	    
	    // 1   这  这  这  DT  DT  _   _   6   4   DMOD    ADV _   _   _   _   _   _
	    
	    int numframes = 0, cur = 0;
	    for (int i = 1; i < length + 1; ++i) {
	    	String[] parts = lstLines.get(i-1);
	    	if (!parts[12].equals("_")) ++numframes;
	    }
	    
	    SemanticFrame[] frames = new SemanticFrame[numframes];
	    for (int i = 0; i < numframes; ++i)
	    	frames[i] = new SemanticFrame(length+1);
	    
	    for (int i = 1; i < length + 1; ++i) {
	    	String[] parts = lstLines.get(i-1);
	    	forms[i] = normalize(parts[1]);
	    	if (!parts[3].equals("_")) { 
	    		lemmas[i] = normalize(parts[3]);
	    		hasLemma = true;
	    	} //else lemmas[i] = forms[i];
	    	
	    	pos[i] = parts[5];
	    	
	    	// todo: use coarse map or use first few chars of pos[i]
	    	cpos[i] = parts[5];	
	    	
	    	if (!parts[7].equals("_")) feats[i] = parts[7].split("\\|");
	    	heads[i] = Integer.parseInt(parts[8]);
	    	deprels[i] = parts[10];
	    	
	    	if (!parts[12].equals("_")) {
	    		frames[cur].predid = i;
	    		predIndex[i] = cur;
                //System.out.println(parts[13]);
                //try {
	    		//    frames[cur].sense = Integer.parseInt(parts[13].substring(parts[13].length()-2));
                //} catch (Exception e) {
                //    //System.out.printf("warning: bad gold predicate sense=%s%n", parts[13]);
                //    frames[cur].sense = 1;
                //}	
	    		//Utils.Assert(parts[13].charAt(parts[13].length()-3) == '.');
	    		//Utils.Assert(frames[cur].sense > 0);
	    		
	    		++cur;
	    	}
	    	
	    	for (int j = 14; j < parts.length; ++j)
	    		if (!parts[j].equals("_")) {
	    			frames[j-14].arglbs[i] = parts[j];
	    		}
	    	
	    	// passive voice
	    	if (pos[i].startsWith("V")) {
		    	if ((pos[i].equals("VBN") && lemmas[i - 1].equals("be"))
		    			|| (i > 1 && pos[i].equals("VBN") && lemmas[i - 1].equals("not") && lemmas[i - 2].equals("be")))
		    		voice[i] = 1;
		    	else
		    		voice[i] = 0;
	    	}
	    	else {
	    		voice[i] = 2;
	    	}
	    	
	    }
	    if (!hasLemma) lemmas = null;
	    
		return new DependencyInstance(forms, lemmas, cpos, pos, feats, heads, deprels, frames, predIndex, voice);
	}

}

