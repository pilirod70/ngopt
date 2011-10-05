package org.halophiles.assembly.qc;

import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Stack;
import java.util.Vector;

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
	
	/**
	 * Constructs a new PointChainer. Builds the matrix of MatchPoints and runs
	 * dynamic programming algorithm. 
	 * 
	 * @param p1 a sorted array of points in contig 1 
	 * @param p2 a sorted array of points in contig 2
	 * @param matches an array of length 2 arrays. each array contains the point in each contig that comprise this match
	 */
	public PointChainer(int[] p1, int[] p2, int[][] matches) {
		matrix = new MatchPoint[p1.length][p2.length];
		points = buildMatrix(matrix, p1, p2, matches);
		buildNeighborHoods(matrix, p1, p2);
		getScores(matrix);
		// first find our Kclumps
		Vector<KClump> kclumpSet = new Vector<KClump>();
		Iterator<MatchPoint> it = points.iterator();
		while (it.hasNext()) {
			MatchPoint tmp = it.next();
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
			it = points.iterator();
			while (it.hasNext()) {
				MatchPoint tmp = it.next();
				tmpKc.add(tmp);
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
	 * Fills a matrix with MatchPoints corresponding to matches in <code>matches</code>
	 * 
	 * @param matrix The matrix to fill
	 * @param p1 the points in contig 1
	 * @param p2 the points in contig 2
	 * @param matches matching points between <code>p1</code> and <code>p2</code>
	 * @return the points that are now in <code>matrix</code>
	 */
	private static Collection<MatchPoint> buildMatrix(MatchPoint[][] matrix, int[] p1, int[] p2, int[][] matches) {
		Collection<MatchPoint> points = new Vector<MatchPoint>();
		for (int i = 0; i < matches.length; i++) {
			// Find the index of each point in match i
			int x = Arrays.binarySearch(p1, matches[i][0]);
			int y = Arrays.binarySearch(p2, matches[i][1]);
			if (x < 0) {
				System.err.println("[a5_qc] Can't find " + matches[i][0]
						+ " in map file");
			}
			if (y < 0) {
				System.err.println("[a5_qc] Can't find " + matches[i][1]
						+ " in map file");
			}
			matrix[x][y] = new MatchPoint(matches[i][0], matches[i][1]);
			points.add(matrix[x][y]);
		}
		return points;
	}
	
	/**
	 * For each point in <code>matrix</code>, find all other points in <code>matrix</code> whose nieghborhood the point is in
	 * point j is in the neighborhood of point i if
	 * 			j_x - i_x < MAX_INTERPOINT_DIST and
	 * 			j_y - i_y < MAX_INTERPOINT_DIST
	 * 
	 * @param matrix the matrix of points to build neighborhoods for
	 * @param p1 the points in contig 1
	 * @param p2 the points in contig 2
	 */
	private static void buildNeighborHoods(MatchPoint[][] matrix, int[] p1, int[] p2){
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
	/**
	 * A function to compute scores (i.e. run dynamic programming algorithm)
	 * @param matrix the matrix of points to compute scores for
	 */
	private static void getScores(MatchPoint[][] matrix){
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[i].length - 1; j++) {
				if (matrix[i][j] == null)
					continue;
				matrix[i][j].getScore();
			}
		}
	}
}
