package org.halophiles.assembly.qc;

import java.io.PrintStream;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.Stack;

import org.halophiles.tools.SummaryStats;

public class MatchPoint {
	double S = -1;
	private int x;
	private int y;
	private int xRot;
	private int yRot;
	MatchPoint pred = null;
	MatchPoint incoming = null;
	Set<MatchPoint> neighborhood;
	Set<MatchPoint> tbEdges;
	public MatchPoint(int x, int y){
		this.x = x;
		this.y = y;
		xRot = x;
		yRot = y;
		neighborhood = new HashSet<MatchPoint>();
		tbEdges = new HashSet<MatchPoint>();
	}
	public void invert(){
		yRot = -1*yRot;
	}
	/**
	 * this is a member of p's neighborhood
	 * 
	 * @param p the MatchPoint who's neighborhood this is in
	 */
	public void addNeighborhood(MatchPoint p){
		neighborhood.add(p);
	}
	
	public boolean inOwnNbhood(){
		return neighborhood.contains(this);
	}
	
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
	
	public int x(){
		return xRot;
	}
	
	public int y(){
		return yRot;
	}
	
	public static double getCorrelation(Collection<MatchPoint> points, MatchPoint p){
		double[] x = new double[points.size()+1];
		double[] y = new double[points.size()+1];
		Iterator<MatchPoint> it = points.iterator();
		int i = 0;
		while(it.hasNext()){
			MatchPoint tmp = it.next();
			x[i] = tmp.xRot;
			y[i] = tmp.yRot;
			i++;
		}
		return SummaryStats.correlation(x,y);
	}
	
	private static double numGaps(MatchPoint p1, MatchPoint p2){
		double xDiff = (p2.x-p1.x); 
		double yDiff = (p2.y-p1.y);
		return Math.floor((xDiff+yDiff+Math.abs(xDiff-yDiff))
				/ 2*MisassemblyBreaker.MAX_INTERPOINT_DIST + 0.5);
	}
	
	public void print(PrintStream out){
		out.print(this.toString()+" - ");
		if (pred != null)
			out.print(pred.toString()+ " - ");
		Iterator<MatchPoint> it  = neighborhood.iterator();
		if (it.hasNext())
			out.print(it.next().toString());
		while(it.hasNext())
			out.print(","+it.next().toString());
		out.println();
		
	}
	
	public Stack<MatchPoint> getCC(){
		Stack<MatchPoint> ret = new Stack<MatchPoint>();
		return getCC(ret,this,0);	
	}
	
	private Stack<MatchPoint> getCC(Stack<MatchPoint> points, MatchPoint p, int count){
		if (p.incoming == null) {
			points.push(p);
			return points;
		} else {
			points.push(p);
			return getCC(points,p.incoming,++count);
		}
	}
	public String toString(){
		return "("+x+","+y+")";
	}
	public double xDist(MatchPoint p){
		return Math.abs(this.x - p.x);
	}
	public double yDist(MatchPoint p){
		double yDist = Math.abs(this.y - p.y);
		return yDist;
	}
	
	public static double xyDiff(MatchPoint p1, MatchPoint p2) {
		return Math.abs(Math.abs(p1.x-p2.x)-Math.abs(p1.y-p2.y));
	}
	
	public static double manhattan(MatchPoint p1, MatchPoint p2){
		return Math.abs(p1.x-p2.x)+Math.abs(p1.y-p2.y);
	}

	public static double euclidean(MatchPoint p1, MatchPoint p2) {
		return Math.sqrt(Math.pow(p1.x - p2.x, 2) + Math.pow(p1.y - p2.y, 2));
	}

	public static double slopeDev(MatchPoint p1, MatchPoint p2) {
		double xdiff = p2.x - p1.x;
		double ydiff = p2.y - p1.y;
		return Math.abs(1 - Math.abs(ydiff / xdiff));
	}

	public static double slope(MatchPoint p1, MatchPoint p2) {
		double xdiff = p2.x - p1.x;
		double ydiff = p2.y - p1.y;
		return ydiff / xdiff;
	}
}