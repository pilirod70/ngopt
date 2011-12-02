package org.halophiles.assembly.qc;

import java.util.Set;
import java.util.TreeSet;

public class MatchPoint {
	private int x;
	private int y;
	MatchPoint pred = null;
	MatchPoint incoming = null;
	Set<MatchPoint> neighbors;
	
	private boolean visited;
	private boolean noise;
	private boolean assigned;
	
	/**
	 * Construct a new MatchPoint with the given points x and y
	 * @param x the point in contig 1 
	 * @param y the point in contig 2
	 */
	public MatchPoint(int x, int y){
		this.x = x;
		this.y = y;
		visited = false;
		noise = false;
		assigned = false;
	}
	
	public void addNeighbor(MatchPoint p){
		if (neighbors == null)
			neighbors = new TreeSet<MatchPoint>(SpatialClusterer.xSort);
		neighbors.add(p);
	}
	
	public Set<MatchPoint> getNeighbors(){
		return neighbors;
	}
	
	public int size(){
		return (neighbors == null) ? 0 : neighbors.size();
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

	public void setAssigned() {
		this.assigned = true;
	}

	public boolean isAssigned() {
		return assigned;
	}
	
	public void setNoise(){
		noise = true;
	}

	public boolean isNoise(){
		return noise;
	}
	
	public void setVisited(){
		visited = true;
	}
	
	public boolean isVisited(){
		return visited;
	}
}
