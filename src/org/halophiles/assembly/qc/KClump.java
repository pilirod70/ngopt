package org.halophiles.assembly.qc;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Collection;
import java.util.Iterator;
import java.util.Set;

import org.halophiles.tools.SummaryStats;

public class KClump {
	private static int COUNT = 0;
	
	private int maxResid;
	private int xLowerBnd;
	private int xUpperBnd;
	
	int xMax;
	int xMin;
	int yMax;
	int yMin;
	
	private double slope;
	private double intercept;
	
	private Set<MatchPoint> points;
	
	final int id;
	
	/**
	 * Create a KClump from the set of given points.
	 * 
	 * @param points the points in this KClump
	 * @param maxResid the maximum residual for calling a point a member of this KClump
	 */
	public KClump(Set<MatchPoint> points, int maxResid){
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
		double mu_x = SummaryStats.mean(x);
		double mu_y = SummaryStats.mean(y);
		// compute the linear regression coefficients
		slope = SummaryStats.covariance(x,mu_x,y,mu_y)/SummaryStats.variance(x,mu_x);
		intercept = mu_y - slope*mu_x;
		// set upper and lower x limits for this KClump
		xLowerBnd = xMin - MisassemblyBreaker.MAX_INTERPOINT_DIST; 
		xUpperBnd = xMax + MisassemblyBreaker.MAX_INTERPOINT_DIST;
		this.maxResid = maxResid;
		this.points = points;
	}
	
	/**
	 * Determines if MatchPoint p fits in this KClump
	 * @param p the point to fit
	 * @return true if within x limits and residual is small enough
	 */
	public boolean fits(MatchPoint p){
		double yfit = slope*p.x()+intercept;
		return xLowerBnd < p.x() && p.x() < xUpperBnd && Math.abs(yfit - p.y()) < maxResid;
	}
	
	/**
	 * Add MatchPoint p to this KClump if it fits according to the function <code> fit(MatchPoint p) </code>
	 * @param p
	 * @return false if this point is already in this KClump, or if this point does not fit
	 */
	public boolean add(MatchPoint p){
		if (points.contains(p))
			return false;
		if (fits(p)){
			points.add(p);
			if (p.x() > xMax)
				xMax = p.x();
			else if (p.x() < xMin)
				xMin = p.x();
			if (p.y() > yMax)
				yMax = p.y();
			else if (p.y() < yMin)
				yMin = p.y();
			return true;
		} else {
			return false;
		}
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
	
}
