package parser;

public class SemanticFrame
{
	// sense number (1-based, 0 means unknown)
	public int sense;
	
	// predicate id in the sentence
	public int predid;
		
	// argument labels and label ids (NULL means not argument)
	public String[] arglbs;
	public int[] arglbids;
	
	public SemanticFrame(int length)
	{
		sense = 0;
		predid = 0;
		arglbs = new String[length];
		arglbids = new int[length];
		for (int i = 0; i < length; ++i) {
			arglbs[i] = null;
			arglbids[i] = -1;
		}
	}
	
	public SemanticFrame(SemanticFrame a)
	{
		// only store ids (not strings to save memory)
		
		sense = a.sense;
		predid = a.predid;
		arglbids = a.arglbids;
		arglbs = null;
		
	}
	
	public int numArgs()
	{
		int cnt = 0;
		for (int i = 0, n = arglbids.length; i < n; ++i)
			if (arglbids[i] >= 0) ++cnt;
		return cnt;
	}
	
	public int numCoreArgs()
	{
		int cnt = 0;
		for (int i = 0, n = arglbids.length; i < n; ++i)
			if (arglbids[i] >= 0 && arglbs[i].length() < 3 && arglbs[i].startsWith("A")) ++cnt;
		return cnt;
	}
}