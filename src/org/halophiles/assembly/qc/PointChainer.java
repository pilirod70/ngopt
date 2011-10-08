package org.halophiles.assembly.qc;

import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Stack;
import java.util.Vector;
import java.util.HashMap;
import java.lang.Integer;

import org.halophiles.assembly.Contig;

public class PointChainer {
	
	/**
	 *  maximum residual for adding a MatchPoint to a KClump
	 */ 
	static int MAX_RES;
	
	private Contig ctg1;
	
	private Contig ctg2;

	/**
	 * All points in <code>matrix</code>
	 */
	HashSet<MatchPoint> currPoints;
	
	/**
	 * The KClumps resulting from chaining points.
	 * KClumps should be mutually exclusive
	 * 
	 */
	private HashSet<KClump> kclumps;
	
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
	
	public PointChainer(Contig contig1, Contig contig2){
		this.ctg1 = contig1;
		this.ctg2 = contig2;
		currPoints = new HashSet<MatchPoint>();
		kclumps = new HashSet<KClump>();
	}
	
	public boolean addMatch(int x, int y){
		return currPoints.add(new MatchPoint(x, y));
	}
	
	public void buildKClumps(){
		runAlgo();
		/*
		 * clear the neighborhoods of the remaining points
		 * invert, and re-run the algorithm
		 */
		Iterator<MatchPoint> it = currPoints.iterator();
		while(it.hasNext()){
			MatchPoint tmp = it.next();
			tmp.clearNeighborhood();
			tmp.invert();
		}
		runAlgo();
	}
	

	/**
	 * Return a set of mutually exclusive KClumps
	 * @return a set of mutually exclusive KClumps
	 */
	public KClump[] getKClumps() {
		KClump[] ar = new KClump[kclumps.size()];
		kclumps.toArray(ar);
		return ar;
	}
	
	public Contig getContig1(){
		return ctg1;
	}
	
	public Contig getContig2(){
		return ctg2;
	}
	
	private void runAlgo(){
		MatchPoint[] matchpoints = new MatchPoint[currPoints.size()];
		Integer[] x_order_int = new Integer[matchpoints.length];
		Integer[] y_order_int = new Integer[matchpoints.length];
		Iterator<MatchPoint> it = currPoints.iterator();
		int[] x_order = new int[matchpoints.length];
		int[] y_order = new int[matchpoints.length];	
		HashMap<MatchPoint, Integer> xref = new HashMap<MatchPoint, Integer>();
		HashMap<MatchPoint, Integer> yref = new HashMap<MatchPoint, Integer>();
		for(int i=0; i<matchpoints.length; i++){
			matchpoints[i] = it.next();
			x_order_int[i]=i;
			y_order_int[i]=i;
		}
		MatchComparator mcx = new MatchComparator(0, matchpoints);
		MatchComparator mcy = new MatchComparator(1, matchpoints);
		
		Arrays.sort(x_order_int, mcx);
		Arrays.sort(y_order_int, mcy);
				
		for(int i=0; i<x_order.length; i++){
			x_order[i] = x_order_int[i];
			y_order[i] = y_order_int[i];
			xref.put(matchpoints[x_order[i]], i);
			yref.put(matchpoints[y_order[i]], i);
		}
		/* BUILD NEIGHBORHOODS
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
				matchpoints[x_order[j_x]].x() - matchpoints[i].x() <= MisassemblyBreaker.MAX_INTERPOINT_DIST; 
				j_x++ ){
					int j_in_y = yref.get(matchpoints[x_order[j_x]]);
					if(matchpoints[y_order[j_in_y]].y() > matchpoints[i].y() && 
							matchpoints[y_order[j_in_y]].y() - matchpoints[i].y() <= MisassemblyBreaker.MAX_INTERPOINT_DIST)
						matchpoints[i].addNeighborhood(matchpoints[x_order[j_x]]);
				}
		}
		// run DP algorithm
		for (int i = 0; i < matchpoints.length; i++) {
			matchpoints[i].getScore();
		}
		// first find our Kclumps ( DP traceback )
		Vector<KClump> kclumpSet = new Vector<KClump>();
		for(int i=0; i<matchpoints.length; i++){
			MatchPoint tmp = matchpoints[i];
			if (tmp.pred == null) {
				Stack<MatchPoint> cc = tmp.getCC();
				if (cc.size() < 5)
					continue;
				kclumpSet.add(new KClump(new HashSet<MatchPoint>(cc), MAX_RES));
			}
		}
		// now collect all points that look like they belong to this line.
		Iterator<KClump> kcIt = kclumpSet.iterator();
		while (kcIt.hasNext()) {
			KClump tmpKc = kcIt.next();
			for(int i=0; i<matchpoints.length; i++){
				tmpKc.add(matchpoints[i]);
			}
			currPoints.removeAll(tmpKc.getMatchPoints());
			kclumps.add(tmpKc);
		}
	}
}