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
/**
 * A data structure for storing match points (i.e. paired-read mapping locations)
 * between two contigs, and additional methods for running the DBSCAN spatial clustering
 * algorithm.
 * 
 * @author Andrew Tritt
 *
 */
public class SpatialClusterer {
	
	/** 
	 * The maximum allowed distance between two MatchPoints 
	 * on either axis. Two points, A and B, are considered neighbors
	 * if and only if |A_x-B_x| <= EPS && |A_y-B_y| <= EPS
	 */
	static double EPS;
	
	/**
	 * The minimum number of neighboring points to any given point
	 * for identifying said point as not being noise.
	 */
	static int MIN_PTS = 10;
	
	/**
	 * A Comparator for sorting MatchPoint in ascending order, first by
	 * x-coordinate (aka Contig 1 coordinate) and second by y-coordinate 
	 * (aka Contig 2 coordinate) to break x-coordinate ties.
	 */
	public static Comparator<MatchPoint> xSort = new Comparator<MatchPoint>(){
		@Override
		public int compare(MatchPoint arg0, MatchPoint arg1) {
			if (arg0.x() == arg1.x())
				return arg0.y() - arg1.y();
			else 
				return arg0.x() - arg1.x();
		}
	};

	/**
	 * A Comparator for sorting ReadClusters in ascending order according
	 * to the maximum (most downstream) mapping location on Contig 1 (aka Contig x)
	 */
	public static Comparator<ReadCluster> CLUST_COMP = new Comparator<ReadCluster>(){
		@Override
		public int compare(ReadCluster arg0, ReadCluster arg1) {
			return arg0.xMax - arg1.xMax;
		}
	};

	/**
	 * The first contig in this SpatialClusterer
	 */
	private Contig ctg1;
	
	/**
	 * The second contig in this SpatialClusterer
	 */
	private Contig ctg2;

	/**
	 * All points in <code>matrix</code>
	 */
	private TreeSet<MatchPoint> currPoints;
	
	/**
	 * The number of points this SpatialClusterer is operating on
	 */
	private int numPoints; 
	
	/**
	 * The ReadPair Clusters resulting from chaining points.
	 * ReadClusters should be mutually exclusive.
	 * 
	 */
	private Set<ReadCluster> readPairClusters;
	
	/**
	 * A Comparator for ordering MatchPoints according to a given
	 * Contig's coordinates.
	 */
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
	
	/**
	 * Create a new SpatialCluster with the given Contigs
	 * @param contig1 the first contig to include
	 * @param contig2 the second contig to include
	 */
	public SpatialClusterer(Contig contig1, Contig contig2){
		this.ctg1 = contig1;
		this.ctg2 = contig2;
		currPoints = new TreeSet<MatchPoint>(xSort);
		readPairClusters = new TreeSet<ReadCluster>(CLUST_COMP);
		numPoints = 0;
	}
	
	/**
	 * Print the current state of this SpatialCluster to the given File.
	 * 
	 * Used mainly for debugging.
	 * 
	 * @param file the File to write the current state to.
	 * @throws IOException if an I/O occurred when trying to write to <code>file</code>
	 */
	void exportCurrState(File file) throws IOException{
		file.createNewFile();
		PrintStream out = new PrintStream(file);
		Iterator<MatchPoint> it = currPoints.iterator();
		while(it.hasNext()){
			MatchPoint tmp = it.next();
			out.println(tmp.x()+"\t"+Math.abs(tmp.y())+"\t0");
		}
		
		Iterator<ReadCluster> kcIt = readPairClusters.iterator();
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
	
	/**
	 * The number of points this SpatialCluster is operating on.
	 * @return number of points this SpatialCluster is operating on
	 */
	public int numPoints(){
		return numPoints;
	}
	
	/**
	 * Add a match to this SpatialClusterer
	 * 
	 * @param x the location of the match on Contig 1
	 * @param y the location of the match on Contig 2
	 * @return true if the point was not already contained in the underlying set of match points, false otherwise
	 */
	public boolean addMatch(int x, boolean xRev, int y, boolean yRev){
		numPoints++;
		return currPoints.add(new MatchPoint(x, xRev, y, yRev));
	}
	
	/**
	 * Build ReadPair clusters: First, identify neighbors of each MatchPoint
	 * using the value of <code>EPS</code>, and then run the DBSCAN algorithm. 
	 */
	public void buildReadPairClusters(){
		locateNeighbors();
		runDBSCAN();
		if (readPairClusters.isEmpty())
			return;
		System.out.println("[a5_qc] Found "+readPairClusters.size()+" initial blocks between contigs "+ctg1.getId()+" and "+ctg2.getId());
		Iterator<ReadCluster> kcIt = readPairClusters.iterator();
		while(kcIt.hasNext())
			System.out.println("        "+kcIt.next().toString());
	}

	/**
	 * Return a set of mutually exclusive ReadClusters
	 * @return a set of mutually exclusive ReadClusters
	 */
	public ReadCluster[] getReadClusters() {
		ReadCluster[] ar = new ReadCluster[readPairClusters.size()];
		readPairClusters.toArray(ar);
		return ar;
	}
	
	/**
	 * Return Contig 1 in this SpatialClusterer
	 * @return Contig 1 in this SpatialClusterer
	 */
	public Contig getContig1(){
		return ctg1;
	}

	/**
	 * Return Contig 2 in this SpatialClusterer
	 * @return Contig 2 in this SpatialClusterer
	 */
	public Contig getContig2(){
		return ctg2;
	}
	
	/**
	 * Locate neighboring map points for each ReadPair
	 */
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
				matchpoints[x_order[j_x]].x() - matchpoints[i].x() <= EPS; j_x++)
				if (Math.abs(matchpoints[x_order[j_x]].y() - matchpoints[i].y()) <= EPS &&
						matchpoints[x_order[j_x]].ori() == matchpoints[i].ori())
					matchpoints[i].addNeighbor(matchpoints[x_order[j_x]]);
			for(int j_x=i_in_x-1; j_x >= 0 && 
				matchpoints[i].x() - matchpoints[x_order[j_x]].x() <= EPS; j_x--)
				if (Math.abs(matchpoints[x_order[j_x]].y() - matchpoints[i].y()) <= EPS && 
						matchpoints[x_order[j_x]].ori() == matchpoints[i].ori())
					matchpoints[i].addNeighbor(matchpoints[x_order[j_x]]);
		}
	}
	
	/**
	 * Run the DBSCAN algorithm
	 */
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
				readPairClusters.add(new ReadCluster(currClust));
			}
			tmp.setVisited();
		}
	}
	
	/**
	 * Add the point <code>p</code> to the cluster and expand the cluster
	 * using the neighbors of MatchPoint <code>p</code> 
	 * 
	 * @param p the MatchPoint to expand upon
	 * @param clust 
	 */
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
