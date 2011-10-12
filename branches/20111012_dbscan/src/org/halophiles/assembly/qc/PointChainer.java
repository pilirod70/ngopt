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
import java.io.IOException;
import java.io.PrintStream;
import java.lang.Integer;

import org.halophiles.assembly.Contig;
import org.halophiles.tools.SummaryStats;

public class PointChainer {
	
	private static Comparator<MatchPoint> xSort = new Comparator<MatchPoint>(){

		@Override
		public int compare(MatchPoint arg0, MatchPoint arg1) {
			if (arg0.x() == arg1.x())
				return arg0.y() - arg1.y();
			else 
				return arg0.x() - arg1.x();
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
	TreeSet<MatchPoint> currPoints;
	
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
		currPoints = new TreeSet<MatchPoint>(xSort);
		kclumps = new HashSet<KClump>();
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
		
		// merge overlapping KClumps
		KClump[] ar = new KClump[kclumps.size()];
		kclumps.toArray(ar);
		
		boolean[] altered = new boolean[kclumps.size()];
		
		
		
		for (int i = 0; i < ar.length; i++){
			KClump kcI = ar[i];
			if (altered[i])
				continue;
			for (int j = i+1; j < ar.length; j++){
				KClump kcJ = ar[j];
				if (kcI.xMin < kcJ.xMax && kcJ.xMin < kcI.xMax &&            // x locations overlap
					kcI.yMin < kcJ.yMax && kcJ.yMin < kcI.yMax &&            // y locations overlap
					Math.signum(kcI.slope()) == Math.signum(kcJ.slope())) {  // same orientation
					
					altered[i] = true;
					altered[j] = true;
					HashSet<MatchPoint> merged = new HashSet<MatchPoint>();
					merged.addAll(kcI.getMatchPoints());
					merged.addAll(kcJ.getMatchPoints());
					kclumps.add(new KClump(merged, MAX_RES));
					kclumps.remove(kcJ);
					kclumps.remove(kcI);
					break;
				}
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
							matchpoints[y_order[j_in_y]].y() - matchpoints[i].y() <= MisassemblyBreaker.MAX_INTERPOINT_DIST){
						
						/*double slope = slope(matchpoints[i],matchpoints[x_order[j_x]]);
						if (Math.abs(Math.log10(Math.abs(slope))) < 0.25){
							matchpoints[x_order[j_x]].addNeighborhood(matchpoints[i]);							
						}*/
						matchpoints[x_order[j_x]].addNeighborhood(matchpoints[i]);							
						
					} else if(matchpoints[y_order[j_in_y]].y() < matchpoints[i].y() && 
							 matchpoints[i].y() - matchpoints[y_order[j_in_y]].y() <= MisassemblyBreaker.MAX_INTERPOINT_DIST) {
						/*double slope = slope(matchpoints[i],matchpoints[x_order[j_x]]);
						if (Math.abs(Math.log10(Math.abs(slope))) < 0.25){
							matchpoints[x_order[j_x]].addNeighborhood(matchpoints[i]);
						}*/
						matchpoints[x_order[j_x]].addNeighborhood(matchpoints[i]);
					}
						
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
				Set<MatchPoint> cc = tmp.getCC();
				if (cc.size() < 5)
					continue;
				refine(cc);
				KClump tmpKc = new KClump(new HashSet<MatchPoint>(cc), MAX_RES);
				kclumpSet.add(tmpKc);
				currPoints.removeAll(cc);
			}
		}
		// now collect all points that look like they belong to this line.
		Iterator<KClump> kcIt = kclumpSet.iterator();
		while (kcIt.hasNext()) {
			KClump tmpKc = kcIt.next();
			for(int i=0; i<matchpoints.length; i++){
				if (!currPoints.contains(matchpoints[i]) && tmpKc.add(matchpoints[i]))
					currPoints.remove(matchpoints[i]);
				
			}
			kclumps.add(tmpKc);
		}
	}
	
	/**
	 * Remove all points with residual > MAX_RES
	 * @param points the set of points to refine
	 * @return the Set of points that were removed from <code>points</code>
	 */
	private static Set<MatchPoint> refine(Set<MatchPoint> points){
		double[] x = new double[points.size()];
		double[] y = new double[points.size()];

		Iterator<MatchPoint> it = points.iterator();
		int i = 0;
		while(it.hasNext()){
			MatchPoint tmp = it.next();
			x[i] = tmp.x();
			y[i] = tmp.y();
			i++;
		}
		double mu_x = SummaryStats.mean(x);
		double mu_y = SummaryStats.mean(y);
		// compute the linear regression coefficients
		double slope = SummaryStats.covariance(x,mu_x,y,mu_y)/SummaryStats.variance(x,mu_x);
		double intercept = mu_y - slope*mu_x;
		
		it = points.iterator();
		Set<MatchPoint> toRm = new HashSet<MatchPoint>();
		while(it.hasNext()){
			MatchPoint p = it.next();
			if (Math.abs(slope*p.x()+intercept - p.y()) > MAX_RES){
				toRm.add(p);
			}
		}
		
		it = toRm.iterator();
		while(it.hasNext())
			points.remove(it.next());
		
		return toRm;
	}
	
	
	
}