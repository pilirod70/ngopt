package org.halophiles.assembly.dcj;


import java.io.PrintStream;

import java.util.HashMap;
import java.util.StringTokenizer;

import org.halophiles.assembly.stat.LCB;
//import java.util.Vector;
//import java.util.regex.Pattern;


public class Genome {
	
	
	/*
	//private static String WHITESPACE = "[ \t\n\f\r]";
	
	/** A map for holding identifiers */
	public static HashMap<String, Integer> blockIdMap = new HashMap<String, Integer>();
	
	public static int BLOCK_COUNT = 0;
	
	public static int GENOME_COUNT = 1;
	
	private int numChrom;
	
//	private Vector<Adjacency> adjacencies;
	
	private Adjacency[] adj;
	
	/**
	 *  a BLOCK_COUNT x 2 array storing positions of heads and tails in the Adjacency array
	 *  loc[][0] := tail loc[][1] := head 
	 */
	private int[][] loc;
	
	private FastAccessTable fat;
	
//	private Chromosome[] chrom;
	
	
	private int numLinear;
	
	private String name;
	
	private int id;
	
	public Genome (String g, String name){
		id = GENOME_COUNT++;
		this.name = name;
		StringTokenizer tok = new StringTokenizer(g,"$");
		numChrom = tok.countTokens();
/*		chrom = new Chromosome[numChrom];
		int i = 0;
		while(tok.hasMoreTokens()){
			String c = tok.nextToken();
				chrom[i] = new Chromosome(c);
				if (!chrom[i].isCirc() && chrom[i].hasBlocks())
					numLinear++;
			i++;
		}*/
	//	System.out.println("Found " + BLOCK_COUNT + " blocks ");
		
		loc = new int[BLOCK_COUNT][2];
		
		// for N blocks and k linear chromosomes, the number of adjacencies = N + k
		adj = new Adjacency[BLOCK_COUNT+numLinear];
		addAdjacencies();  // we're doing a lot here... lets make a function for this.
		fat = new FastAccessTable(adj, loc);
		
	}
	
/*	public Genome(Chromosome[] chrom){
		id = GENOME_COUNT++;
		this.chrom = chrom;
		
	}*/
	
	public void addLCBs(LCB[] lcbs){
		
	}
	
	private void addAdjacencies(){
	//	Vector<Adjacency> ret = new Vector<Adjacency>();
		int adjIdx = 0;
/*		for (int i = 0; i < chrom.length; i++){
			Block[] blk = chrom[i].getBlocks();
			if (blk.length == 0){
				continue;
			}
		//	System.out.println(chrom[i].toString());
			if (chrom[i].isCirc()){
				addLocations(blk[blk.length-1], blk[0], adjIdx);
			//	ret.add(new Adjacency(blk[blk.length-1].getRightEnd(),blk[0].getLeftEnd()));
				adj[adjIdx] = new Adjacency(blk[blk.length-1].getRightEnd(),blk[0].getLeftEnd());
				adjIdx++;
			} else {
				// add left telomere location
				if (blk[0].isInverted()){
					loc[blockIdMap.get(blk[0].getName())][1] = adjIdx;
				} else {
					loc[blockIdMap.get(blk[0].getName())][0] = adjIdx;
				}
			//	ret.add(new Adjacency(blk[0].getLeftEnd()));
				adj[adjIdx] = new Adjacency(blk[0].getLeftEnd());
				adjIdx++;
			}
			for (int j = 1; j < blk.length; j++){
				addLocations(blk[j-1], blk[j], adjIdx);
			//	ret.add(new Adjacency(blk[j-1].getRightEnd(), blk[j].getLeftEnd()));
				adj[adjIdx] = new Adjacency(blk[j-1].getRightEnd(), blk[j].getLeftEnd());
				adjIdx++;
			}
			if (!chrom[i].isCirc()){
				//  add right telomere location
				if (blk[blk.length-1].isInverted()){
					loc[blockIdMap.get(blk[blk.length-1].getName())][0] = adjIdx;
				} else {
					loc[blockIdMap.get(blk[blk.length-1].getName())][1] = adjIdx;
				}
			//	ret.add(new Adjacency(blk[blk.length-1].getRightEnd()));
				adj[adjIdx] = new Adjacency(blk[blk.length-1].getRightEnd());
				adjIdx++;
			}
		}*/
	//	return ret;
	}
	
	public String getName(){
		return name;
	}

	private void addLocations(Block first, Block second, int idx){
		// a,b,c*$ a,b,-c*$ -a,b,-c*$ 
		// first is c
		// second is a
		
		if(first.isInverted()){  // we'll get c_t
			loc[blockIdMap.get(first.getName())][0] = idx;
		} else {                // we'll get c_h
			loc[blockIdMap.get(first.getName())][1] = idx;
		}
		if (second.isInverted()){  // we'll get a_h
			loc[blockIdMap.get(second.getName())][1] = idx;
		} else {                   // we'll get a_t
			loc[blockIdMap.get(second.getName())][0] = idx;
		}
	}

	public Adjacency[] getAdjacencies(){
		return adj;
	}
	
	public FastAccessTable getFAT(){
		return fat;
	}
	
	public void printDesc(PrintStream out){
		out.println(name + ": " + this.toString());
	}
	
	public String toString(){
		StringBuilder sb = new StringBuilder();
//		for (int i = 0; i < chrom.length; i++){
//			sb.append(chrom[i].toString());
		//	sb.append("$");
//		}
		return sb.toString();
	}
	
	
}
