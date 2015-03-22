package parser.decoding;

import parser.*;

public abstract class SRLDecoder {
	Options options;
	double nullWeight = 0.0;
	
	public static SRLDecoder createSRLDecoder(Options options)
	{
		if (!options.useSRL2O && !options.useSRLHO) {
			return new MaxMatchingDecoder(options);
		} else {
			return new SRLHillClimbingDecoder(options);
		}
	}
    
    public void shutdown()
    {
    }

	public abstract DependencyInstance decode(DependencyInstance inst,
						SRLFeatureData sfd,
						boolean addLoss);

}
