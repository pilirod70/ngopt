package org.halophiles.assembly;

import java.text.NumberFormat;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Vector;


public class ReadSet {
	private static NumberFormat NF = NumberFormat.getInstance();
	private static int COUNT = 0;
	
	private Map<String,ReadPair> reads;
	private double sd;
	private double mu;
	private final int ID;
	private boolean mod;
	
	public ReadSet(int id){
		this.ID = id;
		reads = new HashMap<String,ReadPair>();
		mod = false;
	}
	
	public ReadSet(){
		this.ID = COUNT++;
		reads = new HashMap<String,ReadPair>();
		mod = false;
	}
	
	public void add(ReadPair read){
		reads.put(read.hdr,read);
		mod = true;
	}
		
	public void remove(ReadPair read){
		reads.remove(read.hdr);
		mod = true;
	}
	
	public Collection<ReadPair> getReads(){
		return new Vector<ReadPair>(reads.values());
	}
	
	public Collection<String> getReadHdrs() {
		return new Vector<String>(reads.keySet());	
	}
	
	public String toString(){
		update();
		NF.setMaximumFractionDigits(0);
		NF.setGroupingUsed(false);
		return "id="+ID+":mu="+NF.format(mu)+":sd="+NF.format(sd)+":n="+reads.size();
	}
	
	public double mean(){
		update();
		return mu;
	}
	
	public double sd(){
		update();
		return sd;
	}
	
	public int size(){
		return reads.size();
	}
	
	public int getId(){
		return ID;
	}
	
	private void update(){
		if (mod){
			double[] ins = ReadPair.estimateInsertSize(reads.values());
			mu = ins[0];
			sd = ins[1];
			mod = false;
		}
	}
	
	public double p(double ins){
		double den = 2*sd*sd;
		return Math.exp(-Math.pow(mu-ins,2)/den)/Math.sqrt(Math.PI*den);
	}
	
}
