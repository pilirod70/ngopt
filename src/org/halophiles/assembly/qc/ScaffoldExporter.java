package org.halophiles.assembly.qc;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;

public class ScaffoldExporter {

	private File outFile;
	private PrintStream out;
	Map<String,Integer> counts;
	
	public ScaffoldExporter (File file) throws IOException{
		outFile = file;
		out = new PrintStream(outFile);
		counts = new HashMap<String,Integer>();
	}
	
	public void export(String hdr, StringBuilder sequence, int left, int right){
		int subSeq = -1;
		if (counts.containsKey(hdr))
			subSeq = counts.get(hdr)+1;
		else 
			subSeq = 1;
		out.println(">"+hdr+"-");
		out.println(sequence.substring(left-1, right));
		counts.put(hdr, subSeq);
	}
	
	public void export(String hdr, StringBuilder sequence){
		out.println(">"+hdr);
		out.println(sequence.toString());
	}
	
	public void close(){
		out.close();
	}
	
}
