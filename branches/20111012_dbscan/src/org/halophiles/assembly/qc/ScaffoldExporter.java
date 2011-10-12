package org.halophiles.assembly.qc;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;
import java.lang.String;

public class ScaffoldExporter {

	private File outFile;
	private PrintStream out;
	/** A map for keeping track of how many subsequences have been printed for this contig */
	Map<String,Integer> counts;
	
	public ScaffoldExporter (File file) throws IOException{
		outFile = file;
		out = new PrintStream(outFile);
		counts = new HashMap<String,Integer>();
	}
	
	/**
	 * print a substring from the given sequence with a given header
	 */
	public void export(String hdr, StringBuilder sequence, int left, int right){
		int subSeq = -1;
		// check for long stretches of N
		String seqstr = sequence.substring(left-1, right);
		String[] brokenScaffolds = seqstr.split("[Nn]{200,}");
		for(int i=0; i<brokenScaffolds.length; i++){
			if (counts.containsKey(hdr))
				subSeq = counts.get(hdr)+1;
			else 
				subSeq = 1;
			out.println(">"+hdr+"-"+subSeq);
			out.println(brokenScaffolds[i]);
			counts.put(hdr, subSeq);
		}
	}
	
	/**
	 * print an entire sequence with the given header and sequence
	 */
	public void export(String hdr, StringBuilder sequence){
		out.println(">"+hdr);
		out.println(sequence.toString());
	}
	
	/**
	 * Close the underlying PrintStream
	 */
	public void close(){
		out.close();
	}
	
}
