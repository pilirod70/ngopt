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
	
	int xMax;
	int xMin;
	int yMax;
	int yMin;
	
	
	private Set<MatchPoint> points;
	
	final int id;
	
	/**
	 * Create a KClump from the set of given points.
	 * 
	 * @param points the points in this KClump
	 * @param maxResid the maximum residual for calling a point a member of this KClump
	 */
	public ReadCluster(Set<MatchPoint> points){
		id = ++COUNT;
		xMax = Integer.MIN_VALUE;
		xMin = Integer.MAX_VALUE;
		yMax = Integer.MIN_VALUE;
		yMin = Integer.MAX_VALUE;
		Iterator<MatchPoint> it = points.iterator();
		int i = 0;
		while(it.hasNext()){
			MatchPoint tmp = it.next();
			if (tmp.x()+RDLEN > xMax)
				xMax = tmp.x()+RDLEN;
			if (tmp.x() < xMin)
				xMin = tmp.x();
			if (tmp.y()+RDLEN > yMax)
				yMax = tmp.y()+RDLEN;
			if (tmp.y() < yMin)
				yMin = tmp.y();			
			
			i++;
		}	
		this.points = points;
	}
	
	
	/**
	 * Returns the number of points in this KClump
	 * @return the number of points in this KClump
	 */
	public int size(){
		return points.size();
	}
	
	/**
	 * Return the points in this KClump
	 * @return the points in this KClump
	 */
	public Set<MatchPoint> getMatchPoints(){
		return points;
	}
	
	public int hashCode(){
		return id;
	}
	
	public String toString(){
		return xMin+"-"+xMax +" <-> "+yMin+"-"+yMax;
	}
	
	/**
	 * Prints this KClump to the given file
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
	
	public double density(){
		return points.size()/((double)(xMax-xMin)*(yMax-yMin));
	}
}
