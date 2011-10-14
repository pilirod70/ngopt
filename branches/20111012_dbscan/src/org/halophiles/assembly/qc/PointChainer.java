package org.halophiles.assembly.qc;

import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.Stack;
import java.util.TreeSet;
import java.util.Vector;
import java.util.HashMap;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.lang.Integer;

import org.halophiles.assembly.Contig;
import org.halophiles.tools.SummaryStats;

public class PointChainer {
	
	static double EPS;
	
	static int MIN_PTS = 1;
	
	public static Comparator<MatchPoint> xSort = new Comparator<MatchPoint>(){

		@Override
		public int compare(MatchPoint arg0, MatchPoint arg1) {
			if (arg0.x() == arg1.x())
				return arg0.y() - arg1.y();
			else 
				return arg0.x() - arg1.x();
		}
		
	};

	public static Comparator<KClump> CLUST_COMP = new Comparator<KClump>(){

		@Override
		public int compare(KClump arg0, KClump arg1) {
			return arg0.xMax - arg1.xMax;
		}
		
	};
	
	/**
	 *  maximum residual for adding a MatchPoint to a KClump
	 */ 
	static int MAX_RES;
	
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
	private Set<KClump> kclumps;
	
	private HashSet<MatchPoint> visited;
	private HashSet<MatchPoint> noise;
	private HashSet<MatchPoint> assigned;
	
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
		currPoints = new TreeSet<MatchPoint>(xSort);
		kclumps = new TreeSet<KClump>(CLUST_COMP);
		visited = new HashSet<MatchPoint>();
		noise = new HashSet<MatchPoint>();
		assigned = new HashSet<MatchPoint>();
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
		
		Iterator<KClump> kcIt = kclumps.iterator();
		while(kcIt.hasNext()){
			KClump tmpKc = kcIt.next();
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
	
	public void buildKClumps(File densFile){
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
		PrintStream out = null;
		try {
			out = new PrintStream(densFile);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		for (int i = 0; i < grid.length; i++){
			for (int j = 0; j < grid[i].length; j++){
				double tmp = grid[i][j]/(winDouble*winDouble);
				if (tmp > 0)
					out.println(grid[i][j]);
				if (tmp > 0 && tmp < minDens)
					minDens = tmp;
			}
		}
		out.close();
		
		locateNeighbors();
		runDBSCAN();
		if (kclumps.isEmpty())
			return;
		System.out.println("[a5_qc] Found "+kclumps.size()+" initial blocks between contigs "+ctg1.getId()+" and "+ctg2.getId());
		
		// merge overlapping KClumps
		Vector<KClump> v = new Vector<KClump>(kclumps);
		int i = 0;
		while(i < v.size()-1){
			KClump kcI = v.get(i);
			KClump kcJ = v.get(i+1);
//			if (kcI.xMin < kcJ.xMax && kcJ.xMin < kcI.xMax &&            // x locations overlap
//				kcI.yMin < kcJ.yMax && kcJ.yMin < kcI.yMax &&            // y locations overlap
//				Math.signum(kcI.slope()) == Math.signum(kcJ.slope())) {  // same orientation
			int xGap = Integer.MAX_VALUE;
			int yGap = Integer.MAX_VALUE;
			if (kcI.xMax > kcJ.xMax )
				xGap = kcI.xMin - kcJ.xMax;
			else
				xGap = kcJ.xMin - kcI.xMax;
			if (kcI.yMax > kcJ.yMax )
				yGap = kcI.yMin - kcJ.yMax;
			else
				yGap = kcJ.yMin - kcI.yMax;
			
			
			if (xGap < EPS && yGap < EPS) {
				Set<MatchPoint> merged = new TreeSet<MatchPoint>(xSort);
				merged.addAll(kcI.getMatchPoints());
				merged.addAll(kcJ.getMatchPoints());
				v.remove(i);
				v.set(i, new KClump(merged));
				break;
			} else{
				i++;
				currPoints.removeAll(kcI.getMatchPoints());
			}
		}
		
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
		/* BUILD NEIGHBORHOODS
		 * 
		 * For each point in <code>matrix</code>, find all other points in <code>matrix</code> whose nieghborhood the point is in
		 * point j is in the neighborhood of point i if
		 * 			j_x - i_x < MAX_INTERPOINT_DIST and
		 * 			j_y - i_y < MAX_INTERPOINT_DIST
		 * 
		 */
		double R = EPS*Math.sqrt(2);
		for( int i=0; i<matchpoints.length; i++){
			// where is this point in x?
			int i_in_x = xref.get(matchpoints[i]);
			for(int j_x=i_in_x+1; j_x < x_order.length && 
				euclidean(matchpoints[x_order[j_x]],matchpoints[i]) <= EPS; 
																	j_x++)
					matchpoints[i].addNeighbor(matchpoints[x_order[j_x]]);
			for(int j_x=i_in_x-1; j_x >= 0 && 
				euclidean(matchpoints[x_order[j_x]],matchpoints[i]) <= EPS; 
																	j_x--)
				matchpoints[i].addNeighbor(matchpoints[x_order[j_x]]);
		}
	}
	
	
	private void runDBSCAN(){
		Iterator<MatchPoint> it = currPoints.iterator();
		Set<MatchPoint> currClust = null;
		MatchPoint tmp = null;
		while(it.hasNext()){
			tmp = it.next();
			if (visited.contains(tmp))
				continue;
			if (tmp.getNeighbors().size() < MIN_PTS){
				noise.add(tmp);
				assigned.add(tmp);
				continue;
			} else {
				currClust = new TreeSet<MatchPoint>(xSort);
				expandClusters(tmp, currClust);
				kclumps.add(new KClump(currClust));
			}
			visited.add(tmp);
		}
	}
	
	private void expandClusters(MatchPoint p, Set<MatchPoint> clust){
		clust.add(p);
		assigned.add(p);
		Vector<MatchPoint> neighbors = new Vector<MatchPoint>(p.getNeighbors());
		int i = 0;
		while(i < neighbors.size()){
			MatchPoint tmp = neighbors.get(i);
			if (!visited.contains(tmp)) {
				visited.add(tmp);
				if (tmp.getNeighbors().size() >= MIN_PTS)
					neighbors.addAll(tmp.getNeighbors());
					
			} 
			if (!assigned.contains(tmp)){
				clust.add(tmp);
				assigned.add(tmp);
			}
			i++;
		}
	}
	
	private static double euclidean(MatchPoint p1, MatchPoint p2){
		return Math.sqrt(Math.pow(p1.x()-p2.x(),2)+Math.pow(p1.y()-p2.y(),2));
	}
	
	
}