/**
 * This file is part of the A5 pipeline.
 * (c) 2011, 2012 Andrew Tritt and Aaron Darling
 * This software is licensed under the GPL, v3.0. Please see the file LICENSE for details
 */
package org.halophiles.assembly.qc;

import java.util.Set;
import java.util.TreeSet;

/**
 * A class for storing paired-read mapping points between two contigs (or within a contig)
 * @author Andrew Tritt
 *
 */
public class MatchPoint {
	public static final int FR = 1;
	public static final int RF = 2;
	public static final int FF = 3;
	public static final int RR = 4;
	
	private int x;
	private int y;
	private int ori;
	MatchPoint pred = null;
	MatchPoint incoming = null;
	Set<MatchPoint> neighbors;
	/**
	 * A boolean to track visit status of this <code>MatchPoint</code> during the DBSCAN algorithm
	 */
	private boolean visited;
	/**
	 * A boolean to mark whether or not this <code>MatchPoint</code> is identified as noise 
	 * during the DBSCAN algorithm
	 */
	private boolean noise;
	/**
	 * A boolean to mark this <code>MatchPoint</code> as having been assigned to a cluster 
	 * yet  
	 */
	private boolean assigned;
	
	/**
	 * Construct a new <code>MatchPoint</code> with the given points x and y
	 * @param x the point in Contig 1 
	 * @param y the point in Contig 2
	 */
	public MatchPoint(int x, boolean xRev, int y, boolean yRev){
		this.x = x;
		this.y = y;
		if (!xRev && yRev)
			ori = FR;
		else if (xRev && !yRev)
			ori = RF;
		else if (xRev && yRev)
			ori = RR;
		else 
			ori = FF;
		visited = false;
		noise = false;
		assigned = false;
	}
	
	/**
	 * Add <code>p</code> to this <code>MatchPoint</code>'s set of neighboring points
	 * @param p the <code>MatchPoint</code> to add
	 */
	public void addNeighbor(MatchPoint p){
		if (neighbors == null)
			neighbors = new TreeSet<MatchPoint>(SpatialClusterer.xSort);
		neighbors.add(p);
	}
	
	/**
	 * Get this <code>MatchPoint</code>'s neighbors
	 * @return this <code>MatchPoint</code>'s neighbors
	 */
	public Set<MatchPoint> getNeighbors(){
		return neighbors;
	}
	
	/**
	 * The number of neighboring points to this <code>MatchPoint</code>
	 * @return
	 */
	public int size(){
		return (neighbors == null) ? 0 : neighbors.size();
	}
	
	/**
	 * Return the x coordinate for this <code>MatchPoint</code>
	 * @return the x coordinate of this <code>MatchPoint</code>
	 */
	public int x(){
		return x;
	}
	
	/**
	 * Return the y coordinate for this <code>MatchPoint</code>
	 * @return the y coordinat of this <code>MatchPoint</code>
	 */
	public int y(){
		return y;
	}
	
	/**
	 * Returns the orientation of this <code>MatchPoint</code>.
	 * The return value can take on one of 4 values:
	 * MatchPoint.FR</br>
	 * MatchPoint.RF</br>
	 * MatchPoint.FF</br>
	 * MatchPoint.RR</br>
	 * @return a flag indicating the orientation of this <code>MatchPoint</code>
	 */
	public int ori(){
		return ori;
	}
	
	/**
	 * Return a String representation of this <code>MatchPoint</code>
	 */
	public String toString(){
		return "("+x+","+y+")";
	}
	
	/**
	 * Set this <code>MatchPoint</code> as assigned to cluster
	 */
	public void setAssigned() {
		this.assigned = true;
	}

	/**
	 * Return true if this <code>MatchPoint</code> is assigned to cluster
	 * @return true if assigned to a cluster, false otherwise
	 */
	public boolean isAssigned() {
		return assigned;
	}
	
	/**
	 * Mark this <code>MatchPoint</code> as being noise (in the context of the DBSCAN algorithm)
	 */
	public void setNoise(){
		noise = true;
	}

	/**
	 * Return true if this <code>MatchPoint</code> has been marked as noise
	 * @return true marked as noise, false otherwise
	 */
	public boolean isNoise(){
		return noise;
	}
	
	/**
	 * Mark this <code>MatchPoint</code> as being visited
	 */
	public void setVisited(){
		visited = true;
	}
	
	/**
	 * Return true if this <code>MatchPoint</code> has been visited yet during the DBSCAN algorithm
	 * @return true if has been visited, false otherwise
	 */
	public boolean isVisited(){
		return visited;
	}
}
