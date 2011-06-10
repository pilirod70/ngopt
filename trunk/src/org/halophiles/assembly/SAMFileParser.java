package org.halophiles.assembly;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

public class SAMFileParser {
	
	private HashMap<String,Contig> contigs;
	private HashMap<String,ReadPair> reads;
	
	private int rdlen;
	
	public SAMFileParser(String samPath) throws IOException{
		contigs = new HashMap<String, Contig>();
		reads = new HashMap<String, ReadPair>();
		BufferedReader br = null;
		String ctgStr = null;

		File samFile = new File(samPath);
		if (samFile.getName().endsWith(".gz")||samFile.getName().endsWith(".Z")||samFile.getName().endsWith(".z"))
			br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(samFile))));
		else
			br = new BufferedReader(new FileReader(samFile));
		while(nextCharIs(br, '@')){ 
			String[] hdr = br.readLine().split("\t");
			Contig tmp = new Contig(hdr[1].substring(3),Integer.parseInt(hdr[2].substring(3)));
			contigs.put(tmp.name, tmp);
		}
		int rdlen = 0;
		int den = 0;
		while(br.ready()){
			ReadPair tmp = null;
			String[] line = br.readLine().split("\t");
			int left = Integer.parseInt(line[3]);
			boolean rev = isReverse(line[1]);
			int len = 0;
			if (Pattern.matches("[0-9]{2,}M",line[5])){
				len = cigarLength(line[5]);
				rdlen += len;
				den++;
			}
			ctgStr = line[2];//.split("\\|")[0];
			if (len <= 0) continue;			
			if (reads.containsKey(line[0]))
				tmp = reads.get(line[0]);
			else {
				tmp = new ReadPair(line[0]);
				reads.put(line[0], tmp);
			}
			Contig tmpCtg = contigs.get(ctgStr);
			tmpCtg.addReadPair(tmp);
			int val = tmp.addRead(left, rev, tmpCtg, line);
			if (val == 2){
				
			} else if (val == -1)
				System.err.println("ambiguous mapping for read " + line[0]);
		}
	}

	
	/*public SAMFileParser(String samPath, String gclPath) throws IOException{
		this(samPath);
		File gclFile = new File(gclPath);
		BufferedReader br = null;
		
		br = new BufferedReader(new FileReader(gclFile));
		br.readLine();
		Contig tmp = null;
		while(br.ready()){
			String[] line = br.readLine().split("\t");
			String ctg = line[0].split("\\|")[0];
			if (contigs.containsKey(ctg)){
				contigs.get(ctg).setCov(Double.parseDouble(line[1]));
			} else {
				tmp = new Contig(ctg, Integer.parseInt(line[3]), Double.parseDouble(line[1]));
				contigs.put(tmp.name, tmp);
			}
		}
	}*/

	public Iterator<Contig> getContigs(){
		return contigs.values().iterator();
	}
	
	public Iterator<ReadPair> getReadPairs(){
		return reads.values().iterator();
	}
	
	/**
	 * Return the average read length in this SAM file
	 * @return the average read length in the SAM file
	 */
	public int getReadLength() {
		return rdlen;
	}
	
	public static boolean nextCharIs(BufferedReader br, char c) throws IOException{
		if (!br.ready()){ return false; }
		boolean ret = false;
		br.mark(1);
		char b = (char) br.read();
		if (b == c) ret = true; 
		else ret = false;
		br.reset();
		return ret;
	}		
	
	public static boolean isReverse(String flag){
		int iflag = Integer.parseInt(flag);
		if (getBit(4,iflag) == 1) return true;
			else return false;
	}
		
	public static int getBit (int bit, int flag) { 
       int mod = 0;
       int dig = 0;
	   while( flag != 0 && dig <= bit) {  
		   mod = flag % 2;
		   flag = flag / 2;
		   dig++;
	   }
	   return mod;  
	} 
	/**
	 * Return mapping length for CIGAR String if there is a match (i.e. string contains 'M') else return -1
	 * @param cig the CIGAR string to parse
	 * @return the length of the match indicated by the CIGAR string. Return -1 if no match
	 */
	public static int cigarLength(String cig){
		if (cig.contains("M"))
			return Integer.parseInt(cig.substring(0,cig.indexOf('M')));
		else return -1;
	}
	
		
}
