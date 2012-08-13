/**
 * This file is part of the A5 pipeline.
 * (c) 2011, 2012 Andrew Tritt and Aaron Darling
 * This software is licensed under the GPL, v3.0. Please see the file LICENSE for details
 */
package org.halophiles.assembly.qc;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Comparator;
import java.util.Iterator;
import java.util.TreeSet;

import org.halophiles.assembly.Contig;

/**
 * A deprecated class
 * @author ajtritt
 * @deprecated See SpatialCluster.java
 */
public class MatchBuilder {
	// a comparator for sorting reads by position
	private static final Comparator<int[]> COMP = new Comparator<int[]>(){
		public int compare(int[] arg1, int[] arg2) {
			if (arg1[0] == arg2[0])
				return 0;	
			else if (arg1[1] == arg2[1])
				return 0;
			else if (arg1[0] < arg2[0])
				return -1;
			else 
				return 1;
		} 
		
	};
	private Contig ctg1;
	private Contig ctg2;
	private int nPairs;
	
	private TreeSet<int[]> v1;
	
	public MatchBuilder(Contig ctg1, Contig ctg2) throws IOException{
		nPairs = 0;
		if (ctg1==ctg2){
			this.ctg1 = ctg1;
			this.ctg2 = ctg2;
		} else { 
			this.ctg1 = ctg1;
			this.ctg2 = ctg2;
		}
		v1 = new TreeSet<int[]>(COMP);
	}
			
	/**
	 * Assumes pos1 is from ctg1 and pos2 is from ctg2
	 * 
	 * @param pos1 position on <code>ctg1</code> from constructor 
	 * @param pos2 position on <code>ctg2</code> from constructor
	 * @param qual
	 */
	public void addMatch(int pos1, int pos2) {
		int[] ar = {pos1,pos2};
		v1.add(ar);
	}
	
	public void print(File file) throws IOException{
		file.createNewFile();
		PrintStream out = new PrintStream(file);
		Iterator<int[]> it = v1.iterator();
		while(it.hasNext()){
			int[] tmp = it.next();
			out.println(tmp[0]+"\t"+tmp[1]);
		}
		out.close();
	}
	
	public int getNumPairs(){
		return nPairs;
	}
	
	public int[][] getMatches(){
		int[][] ret = new int[v1.size()][];
		v1.toArray(ret);
		return ret;
	}
	
	public Contig getContig1(){
		return ctg1;
	}
	
	public Contig getContig2(){
		return ctg2;
	}
}
