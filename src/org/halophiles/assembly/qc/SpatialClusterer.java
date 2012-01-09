/**
 * This file is part of the A5 pipeline.
 * (c) 2011, 2012 Andrew Tritt and Aaron Darling
 * This software is licensed under the GPL, v3.0. Please see the file LICENSE for details
 */
package org.halophiles.assembly.qc;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeSet;
import java.util.Vector;

import org.halophiles.assembly.Contig;

public class SpatialClusterer {
	
	static double EPS;
	
	static int MIN_PTS = 10;
	
	public static Comparator<MatchPoint> xSort = new Comparator<MatchPoint>(){
		@Override
		public int compare(MatchPoint arg0, MatchPoint arg1) {
			if (arg0.x() == arg1.x())
				return arg0.y() - arg1.y();
			else 
				return arg0.x() - arg1.x();
		}
	};

	public static Comparator<ReadCluster> CLUST_COMP = new Comparator<ReadCluster>(){
		@Override
		public int compare(ReadCluster arg0, ReadCluster arg1) {
			return arg0.xMax - arg1.xMax;
		}
	};
	
	private Contig ctg1;
	
	private Contig ctg2;

	/**
	 * All points in <code>matrix</code>
	 */
	private TreeSet<MatchPoint> currPoints;
	
	private int numPoints; 
	
	/**
	 * The KClumps resulting from chaining points.
	 * KClumps should be mutually exclusive
	 * 
	 */
	private Set<ReadCluster> kclumps;
	
	private class MatchComparator implements Comparator<Integer>{
		public MatchComparator(int contig, MatchPoint[] matches){
			this.contig = contig;
			this.matches = matches;
		}
		public int compare(Integer a, Integer b){
			if (contig==0)
				return matches[a.intValue()].x()-matches[b.intValue()].x();
			else 
				return matches[a.intValue()].y()-matches[b.intValue()].y();
		}
		int contig;
		MatchPoint[] matches;
	}
	
	public SpatialClusterer(Contig contig1, Contig contig2){
		this.ctg1 = contig1;
		this.ctg2 = contig2;
		currPoints = new TreeSet<MatchPoint>(xSort);
		kclumps = new TreeSet<ReadCluster>(CLUST_COMP);
		numPoints = 0;
	}
	
	public void exportCurrState(File file) throws IOException{
		file.createNewFile();
		PrintStream out = new PrintStream(file);
		Iterator<MatchPoint> it = currPoints.iterator();
		while(it.hasNext()){
			MatchPoint tmp = it.next();
			out.println(tmp.x()+"\t"+Math.abs(tmp.y())+"\t0");
		}
		
		Iterator<ReadCluster> kcIt = kclumps.iterator();
		while(kcIt.hasNext()){
			ReadCluster tmpKc = kcIt.next();
			it = tmpKc.getMatchPoints().iterator();
			while(it.hasNext()){
				MatchPoint tmp = it.next();
				out.println(tmp.x()+"\t"+Math.abs(tmp.y())+"\t"+tmpKc.id);
			}
		}
		out.close();
	}
	
	public int numPoints(){
		return numPoints;
	}
	
	public boolean addMatch(int x, int y){
		numPoints++;
		return currPoints.add(new MatchPoint(x, y));
	}
	
	public void buildKClumps(){
		int win = Math.max(1000,MisassemblyBreaker.MEAN_BLOCK_LEN);
		int[][] grid = new int[ctg1.len/win+(ctg1.len%win==0?0:1)]
		                       [ctg2.len/win+(ctg2.len%win==0?0:1)];
		Iterator<MatchPoint> it = currPoints.iterator();
		while(it.hasNext()){
			MatchPoint tmp = it.next();
			int xIdx = tmp.x()/win+(tmp.x()%win==0?-1:0);
			int yIdx = tmp.y()/win+(tmp.y()%win==0?-1:0); 
			grid[xIdx][yIdx]++;
		}
		double winDouble = win;
		double minDens = Double.POSITIVE_INFINITY;
		
		for (int i = 0; i < grid.length; i++){
			for (int j = 0; j < grid[i].length; j++){
				double tmp = grid[i][j]/(winDouble*winDouble);
				if (tmp > 0 && tmp < minDens)
					minDens = tmp;
			}
		}
		locateNeighbors();
		runDBSCAN();
		if (kclumps.isEmpty())
			return;
		System.out.println("[a5_qc] Found "+kclumps.size()+" initial blocks between contigs "+ctg1.getId()+" and "+ctg2.getId());
		Iterator<ReadCluster> kcIt = kclumps.iterator();
		while(kcIt.hasNext())
			System.out.println("        "+kcIt.next().toString());
	}

	/**
	 * Return a set of mutually exclusive KClumps
	 * @return a set of mutually exclusive KClumps
	 */
	public ReadCluster[] getKClumps() {
		ReadCluster[] ar = new ReadCluster[kclumps.size()];
		kclumps.toArray(ar);
		return ar;
	}
	
	public Contig getContig1(){
		return ctg1;
	}
	
	public Contig getContig2(){
		return ctg2;
	}
	
	private void locateNeighbors(){
		MatchPoint[] matchpoints = new MatchPoint[currPoints.size()];
		Integer[] x_order_int = new Integer[matchpoints.length];
		Iterator<MatchPoint> it = currPoints.iterator();
		int[] x_order = new int[matchpoints.length];
		HashMap<MatchPoint, Integer> xref = new HashMap<MatchPoint, Integer>();
		for(int i=0; i<matchpoints.length; i++){
			matchpoints[i] = it.next();
			x_order_int[i]=i;
		}
		MatchComparator mcx = new MatchComparator(0, matchpoints);
		
		Arrays.sort(x_order_int, mcx);
				
		for(int i=0; i<x_order.length; i++){
			x_order[i] = x_order_int[i];
			xref.put(matchpoints[x_order[i]], i);
		}
		/* FIND NEIGHBORS
		 * 
		 * For each point in <code>matrix</code>, find all other points in <code>matrix</code> whose nieghborhood the point is in
		 * point j is in the neighborhood of point i if
		 * 			j_x - i_x < MAX_INTERPOINT_DIST and
		 * 			j_y - i_y < MAX_INTERPOINT_DIST
		 * 
		 */
		for( int i=0; i<matchpoints.length; i++){
			// where is this point in x?
			int i_in_x = xref.get(matchpoints[i]);
			for(int j_x=i_in_x+1; j_x < x_order.length && 
				matchpoints[x_order[j_x]].x() - matchpoints[i].x() <= EPS; 
																	j_x++)
				if (Math.abs(matchpoints[x_order[j_x]].y() - matchpoints[i].y()) <= EPS)
					matchpoints[i].addNeighbor(matchpoints[x_order[j_x]]);
			for(int j_x=i_in_x-1; j_x >= 0 && 
				matchpoints[i].x() - matchpoints[x_order[j_x]].x() <= EPS; 
																	j_x--)
				if (Math.abs(matchpoints[x_order[j_x]].y() - matchpoints[i].y()) <= EPS)
					matchpoints[i].addNeighbor(matchpoints[x_order[j_x]]);
		}
	}
	
	
	private void runDBSCAN(){
		Iterator<MatchPoint> it = currPoints.iterator();
		Set<MatchPoint> currClust = null;
		MatchPoint tmp = null;
		while(it.hasNext()){
			tmp = it.next();
			if (tmp.isVisited())
				continue;
			if (tmp.size() < MIN_PTS){
				tmp.setNoise();
				tmp.setAssigned();
			} else {
				currClust = new TreeSet<MatchPoint>(xSort);
				expandClusters(tmp, currClust);
				kclumps.add(new ReadCluster(currClust));
			}
			tmp.setVisited();
		}
	}
	
	private void expandClusters(MatchPoint p, Set<MatchPoint> clust){
		clust.add(p);
		p.setAssigned();
		Vector<MatchPoint> neighbors = new Vector<MatchPoint>(p.getNeighbors());
		int i = 0;
		while(i < neighbors.size()){
			MatchPoint tmp = neighbors.get(i);
			if (!tmp.isVisited()) {
				tmp.setVisited();
				if (tmp.size() >= MIN_PTS)
					neighbors.addAll(tmp.getNeighbors());
			} 
			if (!tmp.isAssigned()){
				clust.add(tmp);
				tmp.setAssigned();
			}
			i++;
		}
	}
	
}
