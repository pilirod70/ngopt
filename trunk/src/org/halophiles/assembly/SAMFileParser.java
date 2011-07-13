package org.halophiles.assembly;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.StringTokenizer;
import java.util.zip.GZIPInputStream;

public class SAMFileParser {
	
	private HashMap<String,Contig> contigs;
	private HashMap<String,ReadPair> reads;
	
	private int nreads;
	
	private int npairs;
	
	private int rdlen;
	
	private File samFile;
	
	public SAMFileParser(String samPath) throws IOException{
		contigs = new HashMap<String, Contig>();
		reads = new HashMap<String, ReadPair>();
		BufferedReader br = null;
		String ctgStr = null;
		
		nreads = 0;
		npairs = 0;

		samFile = new File(samPath);
		if (samFile.getName().endsWith(".gz")||samFile.getName().endsWith(".Z")||samFile.getName().endsWith(".z"))
			br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(samFile))));
		else
			br = new BufferedReader(new FileReader(samFile));
		while(nextCharIs(br, '@')){ 
			String[] hdr = br.readLine().split("\t");
			if (!hdr[0].equals("@SQ")) continue;
			String name = null;
			int len = -1;
			for (String s: hdr){
				if (s.startsWith("SN")){
					name = s.substring(3);
				} else if (s.startsWith("LN")){
					len = Integer.parseInt(s.substring(3));
				}
			}
			if (name == null) 
				System.err.println("Found nameless contig in SAM header");
			else if (len == -1) 
				System.err.println("Found contig of unknown length in SAM header");
			Contig tmp = new Contig(name,len);
			contigs.put(tmp.name, tmp);
		}
		if (contigs.size() == 0){
			System.err.println("0 contigs found in SAM header.");
		}
		while(br.ready()){
			ReadPair tmp = null;
			String[] line = br.readLine().split("\t");
			int left = Integer.parseInt(line[3]);
			boolean rev = isReverse(line[1]);
			int len = 0;
			if (line[5].contains("M"))
				len = cigarLength(line[5]);
			ctgStr = line[2];
			if (len <= 0) continue;			
			if (ctgStr.equals("*")) continue;
			if (reads.containsKey(line[0]))
				tmp = reads.get(line[0]);
			else {
				tmp = new ReadPair(line[0]);
				reads.put(line[0], tmp);
			}
			Contig tmpCtg = null;
			if (contigs.containsKey(ctgStr))
				tmpCtg = contigs.get(ctgStr);
			else {
				tmpCtg = new Contig(ctgStr);
				contigs.put(ctgStr, tmpCtg);
			}
			int val = tmp.addRead(left, rev, len, tmpCtg, line);
			if (val == 2){
				npairs++;
			} else if (val == -1){
				System.err.println("ambiguous mapping for read " + line[0]);
			}
			tmpCtg.addReadPair(tmp);
			nreads++;
		}
	}
	
	public String getBasename(){
		if (samFile.getName().endsWith(".sam"))
			return samFile.getName().substring(0,samFile.getName().indexOf(".sam")); 
		else 
			return samFile.getName();
	}

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
	
	public int getNumReads(){
		return nreads;
	}
	
	public int getNumPairs() {	
		return npairs;
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
		StringTokenizer tok = new StringTokenizer(cig, "MIDNSHP", true);
		int totalLen = 0;
		int alignLen = 0;
		while (tok.hasMoreTokens()){
			int len = Integer.parseInt(tok.nextToken());
			char op = tok.nextToken().charAt(0);
			if (op == 'M'){
				alignLen += len;
			}
			if (op != 'I')
				totalLen += len;
		}
		return alignLen;
	}
	
		
}
