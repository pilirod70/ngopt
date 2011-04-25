package org.halophiles.assembly;

import java.util.*;

public class Contig implements Comparable<Contig> {
	private static int CTG_COUNT=0;
	private static int CONCAT_START = 1;
	private int id;
	public String name;
	public int len;
	public double cov;
	private int start;
	private Map<Contig,Integer> counts;
	public int numSelfConnect;
	private Map<String,ReadPair> reads;
	public Contig(String name, int len){
		String[] spl = name.split("\\|");
		this.name = spl[0];
		id = Integer.parseInt(this.name.substring(this.name.indexOf("scaffold")+8, this.name.indexOf('.')));
		this.len = len;
		this.cov = -1;
		numSelfConnect = 0;
		start=CONCAT_START;
		CONCAT_START+=len;
		counts = new HashMap<Contig,Integer>();
		reads = new HashMap<String,ReadPair>();
		++CTG_COUNT;
	}
	public Contig(String name, int len, double cov){
		this(name,len);	
		this.cov = cov;
	}
	public void setCov(double cov){
		this.cov = cov;
	}
	public boolean equals(Contig c){ 
		return this.name.equals(c.name);
	}
	public void addLink(Contig c){
		if (!counts.containsKey(c)) 
			counts.put(c, 1);
		else 
			counts.put(c, counts.get(c)+1);
	}
	public int nLinks(Contig c){
		if (!counts.containsKey(c)) 
			return 0;
		else 
			return counts.get(c);
	}
	public int hashCode(){
		return name.hashCode();
	}
	public String toString(){
		return name;
	}
	public void addEndSpanningPair(){
		numSelfConnect++;
	}
	public int getConcatCoord(int pos){
		return start+pos-1;
	}
	public int compareTo(Contig arg0) {
		return this.name.compareTo(arg0.name);
	}
	public void addReadPair(ReadPair pair){
		
		reads.put(pair.hdr, pair);
	}
	public int getNumLinkedContigs(){
		return counts.size();
	}
	public int getNumLinks(Contig c){
		return counts.get(c);
	}
	public int getId(){
		return id;
	}
}
