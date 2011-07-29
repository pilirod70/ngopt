package org.halophiles.assembly;

import java.util.Map;
import java.util.HashMap;

public class Contig implements Comparable<Contig> {
	private static int CTG_COUNT=0;
	private static int CONCAT_START = 1;
	private ContigTerminal start;
	private ContigTerminal end;
	private int id;
	private int rankId;
	public String name;
	public int len;
	private double cov = 0;
	private int concat_start;
	private Map<Contig,Integer> counts;
	public int numSelfConnect;
	private Map<String,ReadPair> reads;
	private double mappedBasesCount;
	public Contig(String name, int len){
	//	String[] spl = name.split("\\|");
	//	this.name = spl[0];
		this(name);
		this.len = len;
		concat_start=CONCAT_START;
		CONCAT_START+=len;
	}
	
	public Contig(String name){
		//	String[] spl = name.split("\\|");
		//	this.name = spl[0];
			this.name = name;
			if (this.name.startsWith("node")){
				// if contigs came from IDBA : increment by 1 so we don't start at 0. We don't want FISH to shit itself.
				id = Integer.parseInt(this.name.substring(this.name.indexOf("node")+4, this.name.indexOf('_'))) + 1; 			
			} else if (this.name.startsWith("scaffold")){
				// if contigs came from SSPACE
				id = Integer.parseInt(this.name.substring(this.name.indexOf("scaffold")+8, this.name.indexOf('.')));			
			} else {
				id = CTG_COUNT+1;
			}
			this.len= -1;
			this.concat_start = -1;
			this.rankId = id;
			this.cov = -1;
			numSelfConnect = 0;
			counts = new HashMap<Contig,Integer>();
			reads = new HashMap<String,ReadPair>();
			++CTG_COUNT;
			mappedBasesCount = 0.0;
			start = new ContigTerminal(this, ContigTerminal.START);
			end = new ContigTerminal(this,ContigTerminal.END);
		}
	
	public double getCov(){
		if (len == -1) return -1;
		return mappedBasesCount/len;
	}
	public boolean hasCov(){
		return cov != 0;
	}
	public boolean equals(Contig c){
		if (c == null){
			System.out.print("");
		} else if (c.name == null){
			System.out.print("");
		} else if (this.name == null){
			System.out.print("");			
		}
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
		return concat_start+pos-1;
	}
	public int compareTo(Contig arg0) {
		return this.name.compareTo(arg0.name);
	}
	public void addReadPair(ReadPair pair){
		if (reads.containsKey(pair.hdr)){
			mappedBasesCount += SAMFileParser.cigarLength(pair.cig2);
		} else {
			mappedBasesCount += 
				SAMFileParser.cigarLength(pair.cig1);
		}
		reads.put(pair.hdr, pair);
	}
	public void removeReadPair(String readHdr){
		reads.remove(readHdr);
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
	public void setRank(int id){
		this.rankId = id;
	}
	public int getRank(){
		return this.rankId;
	}
	
	public ContigTerminal getStartTerminus(){
		return start;
	}
	
	public ContigTerminal getEndTerminus(){
		return end;
	}
}
