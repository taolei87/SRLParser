package parser;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStreamWriter;
import java.io.Serializable;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import parser.Options.LearningMode;
import parser.decoding.DependencyDecoder;
import parser.decoding.SRLDecoder;
import parser.feature.SemanticFeatureFactory;
import parser.io.DependencyReader;
import parser.io.DependencyWriter;
import parser.pruning.BasicArcPruner;
import parser.sampling.RandomWalkSampler;
import utils.Utils;

public class DependencyParser implements Serializable {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	
	protected Options options;
	protected DependencyPipe pipe;
	protected Parameters parameters;
	
	DependencyParser pruner;
	
	double pruningGoldHits = 0;
	double pruningTotGold = 1e-30;
	double pruningTotUparcs = 0;
	double pruningTotArcs = 1e-30;
	
	transient boolean saved = false;
	
	public static void main(String[] args) 
		throws IOException, ClassNotFoundException, CloneNotSupportedException
	{
		
		Options options = new Options();
		options.processArguments(args);
		//options.printOptions();
		
		DependencyParser pruner = null;
		if (options.train && options.pruning && options.learningMode != LearningMode.Basic) {
			Options prunerOptions = new Options();
			prunerOptions.processArguments(args);
			prunerOptions.maxNumIters = 10;
			
			prunerOptions.learningMode = LearningMode.Basic;
			prunerOptions.pruning = false;
			prunerOptions.test = false;
			prunerOptions.learnLabel = false;
			prunerOptions.gamma = 1.0;
			prunerOptions.gammaLabel = 1.0;
			
			//pruner = new DependencyParser();
			pruner = new BasicArcPruner();
			pruner.options = prunerOptions;
			
			DependencyPipe pipe = new DependencyPipe(prunerOptions);
			pruner.pipe = pipe;
			
			pipe.createAlphabets(prunerOptions.trainFile);
			DependencyInstance[] lstTrain = pipe.createInstances(prunerOptions.trainFile);
			
			Parameters parameters = new Parameters(pipe, prunerOptions);
			pruner.parameters = parameters;
			
			pruner.train(lstTrain);
		}
		
		if (options.train) {
			DependencyParser parser = new DependencyParser();
			parser.options = options;
			options.printOptions();
			
			DependencyPipe pipe = new DependencyPipe(options);
			parser.pipe = pipe;
			
			if (options.pruning) parser.pruner = pruner;
			
			pipe.createAlphabets(options.trainFile);
			DependencyInstance[] lstTrain = pipe.createInstances(options.trainFile);
			
			Parameters parameters = new Parameters(pipe, options);
			parser.parameters = parameters;
			
			parser.train(lstTrain);
			if (!parser.saved) parser.saveModel();
		}
		
		if (options.test) {
			DependencyParser parser = new DependencyParser();
			parser.options = options;
			
			parser.loadModel();
			parser.options.processArguments(args);
			if (!options.train) parser.options.printOptions();
			
			System.out.printf(" Evaluating: %s%n", options.testFile);
			parser.evaluateSet(true, false);
		}
		
	}
	
    public void saveModel() throws IOException 
    {
    	ObjectOutputStream out = new ObjectOutputStream(
    			new GZIPOutputStream(new FileOutputStream(options.modelFile)));
    	out.writeObject(pipe);
    	out.writeObject(parameters);
    	out.writeObject(options);
    	if (options.pruning && options.learningMode != LearningMode.Basic) 
    		out.writeObject(pruner);
    	out.close();
    }
	
    public void loadModel() throws IOException, ClassNotFoundException 
    {
        ObjectInputStream in = new ObjectInputStream(
                new GZIPInputStream(new FileInputStream(options.modelFile)));    
        pipe = (DependencyPipe) in.readObject();
        parameters = (Parameters) in.readObject();
        options = (Options) in.readObject();
        if (options.pruning && options.learningMode != LearningMode.Basic)
        	//pruner = (DependencyParser) in.readObject();
        	pruner = (BasicArcPruner) in.readObject();
        //pipe.options = options;
        parameters.options = options;        
        in.close();
        pipe.closeAlphabets();
    }
    
	public void printPruningStats()
	{
		System.out.printf("  Pruning Recall: %.4f\tEffcy: %.4f%n",
				pruningGoldHits / pruningTotGold,
				pruningTotUparcs / pruningTotArcs);
	}
	
	public void resetPruningStats()
	{
		pruningGoldHits = 0;
		pruningTotGold = 1e-30;
		pruningTotUparcs = 0;
		pruningTotArcs = 1e-30;
	}
	
    public void train(DependencyInstance[] lstTrain) 
    	throws IOException, CloneNotSupportedException 
    {
    	long start = 0, end = 0;
    	
        if (options.R > 0 && options.gamma < 1 && options.initTensorWithPretrain) {

        	Options optionsBak = (Options) options.clone();
        	options.learningMode = LearningMode.Basic;
        	options.R = 0;
        	options.gamma = 1.0;
        	options.gammaLabel = 1.0;
        	options.maxNumIters = options.numPretrainIters;
            options.useHO = false;
            options.useSRL2O = false;
            options.useSRLHO = false;
        	parameters.gamma = 1.0;
        	parameters.gammaLabel = 1.0;
        	parameters.rank = 0;
        	options.gamma2 = 1.0;
        	parameters.gamma2 = 1.0;
        	parameters.rank2 = 1;
    		System.out.println("=============================================");
    		System.out.printf(" Pre-training:%n");
    		System.out.println("=============================================");
    		
    		start = System.currentTimeMillis();

    		System.out.println("Running MIRA ... ");
    		trainIter(lstTrain, false);
    		System.out.println();
    		
    		System.out.println("Init tensor ... ");
    		LowRankParam tensor = new LowRankParam(parameters);
    		pipe.synFactory.fillParameters(tensor, parameters);
    		tensor.decompose(1, parameters);
            if (options.pminit) {
                SemanticLowRankParam tensor2 = new SemanticLowRankParam(parameters);
                pipe.smnFactory.fillParameters(tensor2, parameters);
                tensor2.decompose(parameters);
            }
            System.out.println();
    		end = System.currentTimeMillis();
    		
    		options.learningMode = optionsBak.learningMode;
    		options.R = optionsBak.R;
    		options.gamma = optionsBak.gamma;
    		options.gammaLabel = optionsBak.gammaLabel;
    		options.maxNumIters = optionsBak.maxNumIters;
            options.useHO = optionsBak.useHO;
            options.useSRL2O = optionsBak.useSRL2O;
            options.useSRLHO = optionsBak.useSRLHO;
    		parameters.rank = optionsBak.R;
    		parameters.gamma = optionsBak.gamma;
    		parameters.gammaLabel = optionsBak.gammaLabel;
    		options.gamma2 = optionsBak.gamma2;
    		parameters.gamma2 = options.gamma2;
    		parameters.rank2 = optionsBak.R2;
    		parameters.clearTheta();
            parameters.printUStat();
            parameters.printVStat();
            parameters.printWStat();
            System.out.println();
            System.out.printf("Pre-training took %d ms.%n", end-start);    		
    		System.out.println("=============================================");
    		System.out.println();	    

        } else {
        	parameters.randomlyInitUVW();
        }
        
        if (!options.pminit)
            parameters.randomlyInitUVWX2();

		System.out.println("=============================================");
		System.out.printf(" Training:%n");
		System.out.println("=============================================");
		
		start = System.currentTimeMillis();

		System.out.println("Running MIRA ... ");
		trainIter(lstTrain, true);
		System.out.println();
		
		end = System.currentTimeMillis();
		
		System.out.printf("Training took %d ms.%n", end-start);    		
		System.out.println("=============================================");
		System.out.println();		    	
    }
    
    public void trainIter(DependencyInstance[] lstTrain, boolean evalAndSave) throws IOException
    {

    	DependencyDecoder synDecoder = DependencyDecoder.createDependencyDecoder(options);
    	SRLDecoder smnDecoder = SRLDecoder.createSRLDecoder(options);
    	
    	double bestF1 = -1;
    	int N = lstTrain.length;
    	int printPeriod = 10000 < N ? N/10 : 1000;
    	
    	for (int iIter = 0; iIter < options.maxNumIters; ++iIter) {
    	    
    		if (pruner != null) pruner.resetPruningStats();
    		
            // use this offset to change the udpate ordering of U, V and W
            // when N is a multiple of 3, such that U, V and W get updated
            // on each sentence.
            int offset = (N % 3 == 0) ? iIter : 0;

    		long start = 0;
    		double loss = 0;
    		Evaluator eval = new Evaluator(options, pipe);
    		//int uas = 0, tot = 0;
    		start = System.currentTimeMillis();
                		    		
    		for (int i = 0; i < N; ++i) {
    			
    			if ((i + 1) % printPeriod == 0) {
				System.out.printf("  %d (time=%ds)", (i+1),
					(System.currentTimeMillis()-start)/1000);
    			}
    			//DependencyInstance inst = new DependencyInstance(lstTrain[i]);
    			DependencyInstance inst = lstTrain[i];
    		    int n = inst.length;

    		    LocalFeatureData lfd = new LocalFeatureData(inst, this, true);
    		    GlobalFeatureData gfd = new GlobalFeatureData(lfd);
    		    
    		    // unlabeled syntactic
    		    DependencyInstance synTree = synDecoder.decode(inst, lfd, gfd, true);
        		loss += parameters.updateSyn(inst, synTree, lfd, gfd,
        				iIter * N + i + 1, offset);
    		    
    		    // label
        		Utils.Assert(options.learnLabel);
        		int[] tmpHeads = synTree.heads;
        		synTree.heads = inst.heads;
        		lfd.predictLabels(synTree.heads, synTree.deplbids, true);
        		int la = evaluateLabelCorrect(inst, synTree);
        		if (la != n-1) {
        			loss += parameters.updateLabel(inst, synTree, lfd, gfd,
        					iIter * N + i + 1, offset);
        		}
        		synTree.heads = tmpHeads;

    		    // semantic
        		lfd.predictLabels(synTree.heads, synTree.deplbids, true);
        		SRLFeatureData sfd = new SRLFeatureData(synTree, options, pipe, parameters);
    		    Utils.Assert(synTree.frames == inst.frames);
    		    DependencyInstance predInst = smnDecoder.decode(synTree, sfd, true);
        		loss += parameters.updateSmn(inst, predInst, sfd,
        				iIter * N + i + 1, offset);
        		
        		eval.add(inst, predInst, true);
    		}
    		System.out.printf("%n  Iter %d\tloss=%.4f\tuas=%.4f\tprec=%.4f\trecl=%.4f\t[%ds]%n",
    				iIter+1,
    				loss, eval.UAS(), eval.Precision(), eval.Recall(),
    				(System.currentTimeMillis() - start)/1000);
    		
    		parameters.printU2Stat();
    		parameters.printV2Stat();
    		parameters.printW2Stat();
    		parameters.printX2Stat();
    		parameters.printThetaStat();
    		
    		if (options.learningMode != LearningMode.Basic && options.pruning && pruner != null)
    			pruner.printPruningStats();
    		
    		// evaluate on a development set
    		if (evalAndSave && options.test && ((iIter+1) % 1 == 0 || iIter+1 == options.maxNumIters)) {		
    			System.out.println();
	  			System.out.println("_____________________________________________");
	  			System.out.println();
	  			System.out.printf(" Evaluation: %s%n", options.testFile);
	  			System.out.println(); 
	  			
                if (options.average) 
                	parameters.averageParameters((iIter+1)*N);
	  			double res = evaluateSet(true, false);
	  			if (res > bestF1) {
	  				saveModel();
	  				bestF1 = res;
	  				saved = true;
	  				System.out.printf("  Best Dev F1=%.6f%n", bestF1);
	  			}
	  			if (options.average) 
                	parameters.unaverageParameters();

                System.out.println();
	  			System.out.println("_____________________________________________");
	  			System.out.println();
    		} 
    	}
    	
    	if (evalAndSave && options.average) {
            parameters.averageParameters(options.maxNumIters * N);
    	}

        synDecoder.shutdown();
        smnDecoder.shutdown();
    }
            
    public int evaluateLabelCorrect(DependencyInstance act, DependencyInstance pred) 
    {
    	int nCorrect = 0;
    	for (int i = 1, N = act.length; i < N; ++i) {
    		if (act.heads[i] == pred.heads[i] && act.deplbids[i] == pred.deplbids[i])
    			++nCorrect;
    	}    		  		
    	return nCorrect;
    }
    
    public double evaluateSet(boolean output, boolean evalWithPunc)
    		throws IOException {
    	
    	if (pruner != null) pruner.resetPruningStats();
    	
    	DependencyReader reader = DependencyReader.createDependencyReader(options);
    	reader.startReading(options.testFile);
    	

    	DependencyWriter writer = null;
    	if (output && options.outFile != null) {
    		writer = DependencyWriter.createDependencyWriter(options, pipe);
    		writer.startWriting(options.outFile);
    	}
    	
    	DependencyDecoder synDecoder = DependencyDecoder.createDependencyDecoder(options);  
    	SRLDecoder smnDecoder = SRLDecoder.createSRLDecoder(options);  

    	Evaluator eval = new Evaluator(options, pipe);
    	
		long start = System.currentTimeMillis();
    	
    	DependencyInstance inst = pipe.createInstance(reader);    	
    	while (inst != null) {
    		LocalFeatureData lfd = new LocalFeatureData(inst, this, true);
    		GlobalFeatureData gfd = new GlobalFeatureData(lfd); 
    		    		    		
		    DependencyInstance synTree = synDecoder.decode(inst, lfd, gfd, false);
		    lfd.predictLabels(synTree.heads, synTree.deplbids, false);

		    SRLFeatureData sfd = new SRLFeatureData(synTree, options, pipe, parameters);
		    DependencyInstance predInst = smnDecoder.decode(synTree, sfd, false);

            eval.add(inst, predInst, evalWithPunc);

    		if (writer != null) {
    			inst.heads = predInst.heads;
    			inst.deplbids = predInst.deplbids;
    			inst.numframes = predInst.numframes;
    			inst.frames = predInst.frames;
    			writer.writeInstance(inst);
    		}
    		
    		inst = pipe.createInstance(reader);
    	}
    	
    	reader.close();
    	if (writer != null) writer.close();
    	
    	System.out.printf("  Tokens: %d%n", eval.tot);
    	System.out.printf("  Sentences: %d%n", eval.nsents);
    	System.out.printf("  UAS=%.6f\tLAS=%.6f\tCAS=%.6f%n" + 
    					  "  Pre=%.6f\tRec=%.6f\t F1=%.6f\t[%ds]%n",
    			eval.UAS(), eval.LAS(), eval.CAS(),
    			eval.Precision(), eval.Recall(), eval.F1(),
    			(System.currentTimeMillis() - start)/1000);
    	System.out.printf("  Pred-Aug pair filtering recall: %.2f%n",
    			eval.FilteringRecall()*100);
    	if (options.pruning && options.learningMode != LearningMode.Basic && pruner != null)
    		pruner.printPruningStats();
    	
    	//eval.dumpArgStats();
    	eval.dumpPathStats();
    	
        synDecoder.shutdown();
        smnDecoder.shutdown();

    	//return eval.UAS();
        return eval.F1();
    }
}

