package org.halophiles.assembly.qc;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Iterator;
import java.util.Set;

import org.halophiles.tools.SummaryStats;

public class KClump {
	private static int COUNT = 0;
	
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
	public KClump(Set<MatchPoint> points){
		id = ++COUNT;
		xMax = Integer.MIN_VALUE;
		xMin = Integer.MAX_VALUE;
		yMax = Integer.MIN_VALUE;
		yMin = Integer.MAX_VALUE;
		double[] x = new double[points.size()];
		double[] y = new double[points.size()];
		Iterator<MatchPoint> it = points.iterator();
		int i = 0;
		while(it.hasNext()){
			MatchPoint tmp = it.next();
			x[i] = tmp.x();
			y[i] = tmp.y();
			if (x[i] > xMax)
				xMax = (int) x[i];
			if (x[i] < xMin)
				xMin = (int) x[i];
			if (y[i] > yMax)
				yMax = (int) y[i];
			if (y[i] < yMin)
				yMin = (int) y[i];			
			
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
