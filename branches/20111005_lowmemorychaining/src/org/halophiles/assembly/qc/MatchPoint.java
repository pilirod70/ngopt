package org.halophiles.assembly.qc;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.Stack;

public class MatchPoint {
	double S = -1;
	private int x;
	private int y;
	MatchPoint pred = null;
	MatchPoint incoming = null;
	Set<MatchPoint> neighborhood;
	Set<MatchPoint> tbEdges;
	/**
	 * Construct a new MatchPoint with the given points x and y
	 * @param x the point in contig 1 
	 * @param y the point in contig 2
	 */
	public MatchPoint(int x, int y){
		this.x = x;
		this.y = y;
		neighborhood = new HashSet<MatchPoint>();
		tbEdges = new HashSet<MatchPoint>();
	}
	/**
	 * Adds a new point to the set of neighborhoods that this MatchPoint is in
	 * 
	 * @param p the MatchPoint who's neighborhood this is in
	 */
	public void addNeighborhood(MatchPoint p){
		neighborhood.add(p);
	}
	
	/**
	 * Computes the score of the maximal chain ending at this MatchPoint
	 * 
	 * @return the score of the maximal chain ending at this MatchPoint
	 */
	public double getScore(){
		if (S != -1) {
			return S;
		} else if (neighborhood.size()==0) {
			S = 0;
			return S;
		} else {			
			// find the predecessor that minimizes the deviation from a slope of 1
			Iterator<MatchPoint> it = neighborhood.iterator();
			double min = Double.POSITIVE_INFINITY;
			double tmpScore = -1;
			// find the minimum value
			while(it.hasNext()){
				MatchPoint tmp = it.next();
				tmpScore = numGaps(tmp,this);
				if (tmpScore < min){
					min = tmpScore;
					pred = tmp;
				}
			}
			pred.incoming = this;
			S = min - 1;
			return S;
		}
	}
	
	/**
	 * Return the x coordinate for this MatchPoint
	 * @return
	 */
	public int x(){
		return x;
	}
	
	/**
	 * Return the y coordinate for this MatchPoint
	 * @return
	 */
	public int y(){
		return y;
	}
	
	/**
	 * Build the connected component starting at this MatchPoint
	 * @return a set of MatchPoints that comprise the connected component that this MatchPoint belongs to.
	 */
	public Stack<MatchPoint> getCC(){
		Stack<MatchPoint> ret = new Stack<MatchPoint>();
		return getCC(ret,this);	
	}
	
	/**
	 * A helper function for getCC()
	 */
	private Stack<MatchPoint> getCC(Stack<MatchPoint> points, MatchPoint p){
		if (p.incoming == null) {
			points.push(p);
			return points;
		} else {
			points.push(p);
			return getCC(points,p.incoming);
		}
	}
	
	/**
	 * Return a String representation of this MatchPoint
	 */
	public String toString(){
		return "("+x+","+y+")";
	}
	
	/**
	 * Computes the number of "gaps" between <code>p1</code> <code>p2</code> using 
	 * the formula given in Haas et al. 2004
	 *  
	 * @param p1
	 * @param p2
	 * @return
	 */
	private static double numGaps(MatchPoint p1, MatchPoint p2){
		double xDiff = (p2.x-p1.x); 
		double yDiff = (p2.y-p1.y);
		return Math.floor((xDiff+yDiff+Math.abs(xDiff-yDiff))
				/ 2*MisassemblyBreaker.MAX_INTERPOINT_DIST + 0.5);
	}

}