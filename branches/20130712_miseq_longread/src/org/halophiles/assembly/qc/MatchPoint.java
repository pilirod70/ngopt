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
	private int x;
	private int y;
	MatchPoint pred = null;
	MatchPoint incoming = null;
	Set<MatchPoint> neighbors;
	/**
	 * A boolean to track visit status of this MatchPoint during the DBSCAN algorithm
	 */
	private boolean visited;
	/**
	 * A boolean to mark whether or not this MatchPoint is identified as noise 
	 * during the DBSCAN algorithm
	 */
	private boolean noise;
	/**
	 * A boolean to mark this MatchPoint as having been assigned to a cluster 
	 * yet  
	 */
	private boolean assigned;
	
	/**
	 * Construct a new MatchPoint with the given points x and y
	 * @param x the point in Contig 1 
	 * @param y the point in Contig 2
	 */
	public MatchPoint(int x, int y){
		this.x = x;
		this.y = y;
		visited = false;
		noise = false;
		assigned = false;
	}
	
	/**
	 * Add <code>p</code> to this MatchPoint's set of neighboring points
	 * @param p the MatchPoint to add
	 */
	public void addNeighbor(MatchPoint p){
		if (neighbors == null)
			neighbors = new TreeSet<MatchPoint>(SpatialClusterer.xSort);
		neighbors.add(p);
	}
	
	/**
	 * Get this MatchPoint's neighbors
	 * @return this MatchPoint's neighbors
	 */
	public Set<MatchPoint> getNeighbors(){
		return neighbors;
	}
	
	/**
	 * The number of neighboring points to this MatchPoint
	 * @return
	 */
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
	 * Set this MatchPoint as assigned to cluster
	 */
	public void setAssigned() {
		this.assigned = true;
	}

	/**
	 * Return true if this MatchPoint is assigned to cluster
	 * @return true if assigned to a cluster, false otherwise
	 */
	public boolean isAssigned() {
		return assigned;
	}
	
	/**
	 * Mark this MatchPoint as being noise (in the context of the DBSCAN algorithm)
	 */
	public void setNoise(){
		noise = true;
	}

	/**
	 * Return true if this MatchPoint has been marked as noise
	 * @return true marked as noise, false otherwise
	 */
	public boolean isNoise(){
		return noise;
	}
	
	/**
	 * Mark this MatchPoint as being visited
	 */
	public void setVisited(){
		visited = true;
	}
	
	/**
	 * Return true if this MatchPoint has been visited yet during the DBSCAN algorithm
	 * @return true if has been visited, false otherwise
	 */
	public boolean isVisited(){
		return visited;
	}
}
