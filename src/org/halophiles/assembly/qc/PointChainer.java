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

public class PointChainer {
	
	/**
	 * A Comparator for sorting KClumps by their id.
	 */
	private static Comparator<KClump> COMP = new Comparator<KClump>() {
		@Override
		public int compare(KClump arg0, KClump arg1) {
			return arg1.id - arg0.id;
		}
	};
	
	/**
	 *  maximum residual for adding a MatchPoint to a KClump
	 */ 
	static int MAX_RES;

	/**
	 * The KClumps resulting from chaining points
	 */
	private KClump[] kclumps;
	
	/**
	 * The matrix of points to use for calculating KClumps 
	 */
	MatchPoint[][] matrix;
	
	/**
	 * All points in <code>matrix</code>
	 */
	Collection<MatchPoint> points;
	
	HashMap<MatchPoint, Integer> xref = new HashMap<MatchPoint, Integer>();
	HashMap<MatchPoint, Integer> yref = new HashMap<MatchPoint, Integer>();
	int[][] xy_matches;
	int[] x_order;
	int[] y_order;
	private class MatchComparator implements Comparator{
		public MatchComparator(int contig, int[][] matches){
			this.contig = contig;
			this.matches = matches;
		}
		public int compare(Object a, Object b){
			return matches[((Integer)a).intValue()][contig]-matches[((Integer)b).intValue()][contig];
		}
		int contig;
		int[][] matches;
	}
	
	/**
	 * Constructs a new PointChainer. Builds the matrix of MatchPoints and runs
	 * dynamic programming algorithm. 
	 * 
	 * @param p1 a sorted array of points in contig 1 
	 * @param p2 a sorted array of points in contig 2
	 * @param matches an array of length 2 arrays. each array contains the point in each contig that comprise this match
	 */
	public PointChainer(int[][] matches) {
		xy_matches = matches;
		
		Integer[] x_order_int = new Integer[matches.length];
		Integer[] y_order_int = new Integer[matches.length];
		x_order = new int[matches.length];
		y_order = new int[matches.length];		
		MatchPoint[] matchpoints = new MatchPoint[matches.length];
		for(int i=0; i<matches.length; i++){
			x_order_int[i]=i;
			y_order_int[i]=i;
			matchpoints[i] = new MatchPoint(matches[i][0],matches[i][1]);
		}
		MatchComparator mcx = new MatchComparator(0, matches);
		MatchComparator mcy = new MatchComparator(1, matches);
		
		Arrays.sort(x_order_int, mcx);
		Arrays.sort(y_order_int, mcy);
				
		for(int i=0; i<x_order.length; i++){
			x_order[i] = x_order_int[i];
			y_order[i] = y_order_int[i];
			xref.put(matchpoints[x_order[i]], i);
			yref.put(matchpoints[y_order[i]], i);
		}

		
		
//		matrix = new MatchPoint[p1.length][p2.length];
//		points = buildMatrix(matrix, p1, p2, matches);
		
		buildNeighborHoods(matchpoints);
		getScores(matchpoints);
		// first find our Kclumps
		Vector<KClump> kclumpSet = new Vector<KClump>();
		for(int i=0; i<matchpoints.length; i++){
			MatchPoint tmp = matchpoints[i];
			//tmp.print(System.out);
			if (tmp.pred == null) {
				Stack<MatchPoint> cc = tmp.getCC();
				if (cc.size() < 10)
					continue;
				kclumpSet.add(new KClump(new HashSet<MatchPoint>(cc), MAX_RES));
			}
		}
		// now collect all points that look like they belong to this line.
		Iterator<KClump> kcIt = kclumpSet.iterator();
		while (kcIt.hasNext()) {
			KClump tmpKc = kcIt.next();
			if(tmpKc.id == 11)
				System.out.print("");
			if(tmpKc.id == 9)
				System.out.print("");
			for(int i=0; i<matchpoints.length; i++){
				tmpKc.add(matchpoints[i]);
			}
		}
		
		this.kclumps = new KClump[kclumpSet.size()];
		kclumpSet.toArray(this.kclumps);
		Arrays.sort(this.kclumps, COMP);
	}

	public KClump[] getKClumps() {
		return kclumps;
	}

	/**
	 * For each point in <code>matrix</code>, find all other points in <code>matrix</code> whose nieghborhood the point is in
	 * point j is in the neighborhood of point i if
	 * 			j_x - i_x < MAX_INTERPOINT_DIST and
	 * 			j_y - i_y < MAX_INTERPOINT_DIST
	 * 
	 */
	private void buildNeighborHoods(MatchPoint[] matchpoints){
		for( int i=0; i<matchpoints.length; i++){
			// where is this point in x?
			int i_in_x = xref.get(matchpoints[i]);
			int i_in_y = yref.get(matchpoints[i]);
			for(int j_x=i_in_x+1; j_x < x_order.length && 
				xy_matches[x_order[j_x]][0] - xy_matches[i][0] <= MisassemblyBreaker.MAX_INTERPOINT_DIST; 
				j_x++ ){
					int j_in_y = yref.get(matchpoints[x_order[j_x]]);
					if(xy_matches[y_order[j_in_y]][1] - xy_matches[i][1] <= MisassemblyBreaker.MAX_INTERPOINT_DIST)
						matchpoints[i].addNeighborhood(matchpoints[x_order[j_x]]);
				}
			for(int j_y=i_in_y+1; j_y < y_order.length && 
				xy_matches[y_order[j_y]][1] - xy_matches[i][1] <= MisassemblyBreaker.MAX_INTERPOINT_DIST; 
				j_y++ ){
				int j_in_x = xref.get(matchpoints[y_order[j_y]]);
				if(xy_matches[x_order[j_in_x]][0] - xy_matches[i][0] <= MisassemblyBreaker.MAX_INTERPOINT_DIST)
					matchpoints[i].addNeighborhood(matchpoints[y_order[j_y]]);
			}
		}
	}
/*
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[i].length; j++) {
				if (matrix[i][j] == null)
					continue;
				for (int y = 1; j + y < matrix[i].length
						&& p2[j + y] - p2[j] <= MisassemblyBreaker.MAX_INTERPOINT_DIST; y++) {
					if (matrix[i][j + y] == null)
						continue;
					matrix[i][j + y].addNeighborhood(matrix[i][j]);
				}
				for (int y = 1; j - y >= 0
						&& p2[j] - p2[j - y] <= MisassemblyBreaker.MAX_INTERPOINT_DIST; y++) {
					if (matrix[i][j - y] == null)
						continue;
					matrix[i][j - y].addNeighborhood(matrix[i][j]);
				}

				// stay in i,j's column, and go up
				for (int y = 1; j + y < matrix[i].length && 
						p2[j + y] - p2[j] <= MisassemblyBreaker.MAX_INTERPOINT_DIST; y++) {
					
					if (matrix[i][j + y] == null)
						continue;
					matrix[i][j + y].addNeighborhood(matrix[i][j]);
				}
				// now go down
				for (int y = 1; j - y >= 0 && 
						p2[j] - p2[j - y] <= MisassemblyBreaker.MAX_INTERPOINT_DIST; y++) {
					
					if (matrix[i][j - y] == null)
						continue;
					matrix[i][j - y].addNeighborhood(matrix[i][j]);
				}

				for (int x = 1; i + x < matrix.length
						&& p1[i + x] - p1[i] <= MisassemblyBreaker.MAX_INTERPOINT_DIST; x++) {

					// start at i,j's current row, and add all points within
					// MisassemblyBreaker.MAX_DIST up
					for (int y = 0; j + y < matrix[i + x].length
							&& p2[j + y] - p2[j] <= MisassemblyBreaker.MAX_INTERPOINT_DIST; y++) {
						if (matrix[i + x][j + y] == null)
							continue;
						matrix[i + x][j + y].addNeighborhood(matrix[i][j]);
					}
					// now go down
					for (int y = 0; j - y >= 0
							&& p2[j] - p2[j - y] <= MisassemblyBreaker.MAX_INTERPOINT_DIST; y++) {
						if (matrix[i + x][j - y] == null)
							continue;
						matrix[i + x][j - y].addNeighborhood(matrix[i][j]);
					}

				}

			}
		}
	}
		*/
	/**
	 * A function to compute scores (i.e. run dynamic programming algorithm)
	 * @param matrix the matrix of points to compute scores for
	 */
	private static void getScores(MatchPoint[] matchpoints){
		for (int i = 0; i < matchpoints.length; i++) {
			matchpoints[i].getScore();
		}
	}
}
