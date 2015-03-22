package parser.io;

import java.io.IOException;

import parser.DependencyInstance;
import parser.DependencyPipe;
import parser.Options;
import parser.SemanticFrame;

public class CONLLWriter extends DependencyWriter {
	
	
	public CONLLWriter(Options options, DependencyPipe pipe) {
		this.options = options;
		this.labels = pipe.types;
		this.args = pipe.args;
	}
	
	@Override
	public void writeInstance(DependencyInstance inst) throws IOException {
		
		if (first) 
			first = false;
		else
			writer.write("\n");
		
		String[] forms = inst.forms;
		String[] lemmas = inst.lemmas;
		String[] cpos = inst.cpostags;
		String[] pos = inst.postags;
		int[] heads = inst.heads;
		int[] labelids = inst.deplbids;
		SemanticFrame[] frames = inst.frames;
		int numframes = inst.numframes;
		
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
		
		for (int i = 1, N = inst.length; i < N; ++i) {
			writer.write(i + "\t");
			writer.write(forms[i] + "\t");
			writer.write((lemmas != null && lemmas[i] != "" ? inst.lemmas[i] : "_") + "\t");
			writer.write((lemmas != null && lemmas[i] != "" ? inst.lemmas[i] : "_") + "\t");
			writer.write(pos[i] + "\t");
            writer.write(pos[i] + "\t");
			writer.write("_\t");
			writer.write("_\t");
			writer.write(heads[i] + "\t");
			writer.write("_\t");
			writer.write((isLabeled ? labels[labelids[i]] : "_") + "\t_\t");
			
			boolean isPred = false;
			for (int j = 0; isPred == false && j < numframes; ++j)
				if (frames[j].predid == i) isPred = true;
			if (isPred) 
				writer.write("Y\t"+lemmas[i]+".01\t");
			else 
				writer.write("_\t_");
			
			for (int j = 0; j < numframes; ++j)
				if (frames[j].arglbids[i] >= 0)
					writer.write("\t"+args[frames[j].arglbids[i]]);
				else
					writer.write("\t_");
			
			writer.write("\n");
		}
	}

}
