/**
 * This file is part of the A5 pipeline.
 * (c) 2011, 2012 Andrew Tritt and Aaron Darling
 * This software is licensed under the GPL, v3.0. Please see the file LICENSE for details
 */
package org.halophiles.assembly.qc;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Iterator;
import java.util.Set;

public class ReadCluster {
	private static int COUNT = 0;
	
	static int RDLEN = 50;
	
	/**
	 * The max MatchPoint location on the x Contig (aka Contig 1)
	 */
	int xMax;
	/**
	 * The min MatchPoint location on the x Contig (aka Contig 1)
	 */
	int xMin;
	/**
	 * The orientation of reads on the x Contig (aka Contig 1)
	 */
	boolean xOri;
	
	/**
	 * The max MatchPoint location on the y Contig (aka Contig 2)
	 */
	int yMax;
	/**
	 * The min MatchPoint location on the y Contig (aka Contig 2)
	 */
	int yMin;
	/**
	 * The orientation of reads on the y Contig (aka Contig 2)
	 */
	boolean yOri;
	
	/**
	 * The set of MatchPoints comprising this ReadCluster
	 */
	private Set<MatchPoint> points;
	
	/**
	 * The uniqe id assigned to this ReadCluster
	 */
	final int id;
	
	/**
	 * Create a <code>ReadCluster</code> from the set of given points.
	 * 
	 * @param points the points in this ReadCluster
	 */
	public ReadCluster(Set<MatchPoint> points){
		id = ++COUNT;
		xMax = Integer.MIN_VALUE;
		xMin = Integer.MAX_VALUE;
		yMax = Integer.MIN_VALUE;
		yMin = Integer.MAX_VALUE;
		Iterator<MatchPoint> it = points.iterator();
		
		// Check the first point so we can get orientation information
		MatchPoint tmp = it.next();
		if (tmp.x()+RDLEN > xMax)
			xMax = tmp.x()+RDLEN;
		if (tmp.x() < xMin)
			xMin = tmp.x();
		if (tmp.y()+RDLEN > yMax)
			yMax = tmp.y()+RDLEN;
		if (tmp.y() < yMin)
			yMin = tmp.y();
		int clustOri = tmp.ori();
		switch (clustOri){
			case MatchPoint.FF: xOri = true; yOri = true; break;
			case MatchPoint.RR: xOri = false; yOri = false; break;
			case MatchPoint.FR: xOri = true; yOri = false; break;
			case MatchPoint.RF: xOri = false; yOri = true; break;
			default: throw new IllegalArgumentException("MatchPoint orientation not a recognizable value");
		}
		while(it.hasNext()){
			tmp = it.next();
			if (tmp.x()+RDLEN > xMax)
				xMax = tmp.x()+RDLEN;
			if (tmp.x() < xMin)
				xMin = tmp.x();
			if (tmp.y()+RDLEN > yMax)
				yMax = tmp.y()+RDLEN;
			if (tmp.y() < yMin)
				yMin = tmp.y();			
		
			if (clustOri != tmp.ori())
				throw new IllegalArgumentException("Inconsistent MatchPoint orienations");
		}	
		this.points = points;
	}
	
	
	/**
	 * Returns the number of points in this ReadCluster
	 * @return the number of points in this ReadCluster
	 */
	public int size(){
		return points.size();
	}
	
	/**
	 * Return the points in this ReadCluster
	 * @return the points in this ReadCluster
	 */
	public Set<MatchPoint> getMatchPoints(){
		return points;
	}
	
	/**
	 * Returns a hash code for this ReadCluster. Just returns the unique id assigned to this
	 * ReadCluster.
	 */
	public int hashCode(){
		return id;
	}
	
	/**
	 * Return a String representation of this ReadCluster in the form:</br>
	 * <code>xMin-xMax <-> yMin-yMax</code>
	 */
	public String toString(){
		return xMin+"-"+xMax +" <-> "+yMin+"-"+yMax;
	}
	
	/**
	 * Prints this ReadCluster to the given file
	 * 
	 * @param file
	 * @throws IOException
	 */
	public void print(File file) throws IOException {
		file.createNewFile();
		PrintStream out = new PrintStream(file);
		Iterator<MatchPoint> it = points.iterator();
		while(it.hasNext()){
			MatchPoint tmp = it.next();
			out.println(tmp.x()+"\t"+tmp.y());
		}
		out.close();
	}
	
}
