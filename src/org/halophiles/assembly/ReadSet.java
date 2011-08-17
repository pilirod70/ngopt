package org.halophiles.assembly;

import java.text.NumberFormat;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;
import java.util.Vector;


public class ReadSet {
	private static NumberFormat NF = NumberFormat.getInstance();
	
	private Set<ReadPair> reads;
	private double sd;
	private double mu;
	final int ID;
	private boolean mod;
	
	
	public ReadSet(int id){
		this.ID = id;
		reads = new HashSet<ReadPair>();
		mod = false;
	}
	
	public void add(ReadPair read){
		reads.add(read);
		mod = true;
	}
	
	public void addAll(Set<ReadPair> toAdd){
		reads.addAll(toAdd);
		mod = true;
	}
	
	public void remove(ReadPair read){
		reads.remove(read);
		mod = true;
	}
	
	public Collection<ReadPair> getReads(){
		return new Vector<ReadPair>(reads);
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
	
	private void update(){
		if (mod){
			double[] ins = ReadPair.estimateInsertSize(reads);
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