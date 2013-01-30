/**
 * This file is part of the A5 pipeline.
 * (c) 2011, 2012 Andrew Tritt and Aaron Darling
 * This software is licensed under the GPL, v3.0. Please see the file LICENSE for details
 */
package org.halophiles.assembly;

import java.util.List;
import java.util.HashMap;
import java.util.Map;
import java.util.Vector;

import org.halophiles.assembly.qc.MisassemblyBlock;

/**
 * A class to represent any sequence in an assembly. The term contig is not the colloquial contig.
 * i.e. an Object of type Contig can represent a scaffold or a contig 
 * 
 * @author ajtritt
 *
 */
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
	private int mappedReadCount;
	
	private int out;
	private int in;
	
	private MisassemblyBlock leftBlock = null;
	private MisassemblyBlock rightBlock = null;
	
	public Contig(String name, int len){
		this(name);
		this.len = len;
		concat_start=CONCAT_START;
		CONCAT_START+=len;
	}
	
	public Contig(String name){
		//	String[] spl = name.split("\\|");
		//	this.name = spl[0];
			this.name = name;
			id = CTG_COUNT+1;
			this.len= -1;
			this.concat_start = -1;
			this.rankId = id;
			this.cov = -1;
			numSelfConnect = 0;
			counts = new HashMap<Contig,Integer>();
			reads = new HashMap<String,ReadPair>();
			++CTG_COUNT;
			mappedBasesCount = 0.0;
			mappedReadCount = 0;
			start = new ContigTerminal(this, ContigTerminal.START);
			end = new ContigTerminal(this,ContigTerminal.END);
			out = 0;
			in = 0;
	}
	
	public void addOut(){
		out++;
	}
	public int getOut(){
		return out;
	}
	public void addIn(){
		in++;
	}
	public int getIn(){
		return in;
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
	public int numReads(){
		return mappedReadCount;
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
	public void addRead(int len){
		mappedBasesCount += len;
		mappedReadCount++;
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
	
	public MisassemblyBlock get5PrimeBlock(){
		return leftBlock;
	}
	
	public void add5PrimeBlock(MisassemblyBlock block){
		leftBlock = block;
	}
	
	public boolean has5PrimeBlock(){
		return leftBlock != null;
	}
	
	public MisassemblyBlock get3PrimeBlock(){
		return rightBlock;
	}
	
	public void add3PrimeBlock(MisassemblyBlock block){
		rightBlock = block;
	}
	
	public boolean has3PrimeBlock(){
		return rightBlock != null;
	}
}
