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
	
	private Collection<MatchPoint> points;
	
	final int id;
	
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
			if (i==0){
			//	ctgX = tmp.
			}
			i++;
		}
		double mu_x = SummaryStats.mean(x);
		double mu_y = SummaryStats.mean(y);
		slope = SummaryStats.covariance(x,mu_x,y,mu_y)/SummaryStats.variance(x,mu_x);
		intercept = mu_y - slope*mu_x;
		xLowerBnd = Math.max(xMin, ((int) mu_x) - MisassemblyBreaker.MAX_INTERPOINT_DIST/2); 
		xUpperBnd = Math.min(xMax, ((int) mu_x) + MisassemblyBreaker.MAX_INTERPOINT_DIST/2);
		this.maxResid = maxResid;
		this.points = points;
	}
	
	public boolean fits(MatchPoint p){
		double yfit = slope*p.x()+intercept;
		return xLowerBnd < p.x() && p.x() < xUpperBnd && Math.abs(yfit - p.y()) < maxResid;
	}
	
	public double density(){
		double dx = ((double)points.size())/(xMax-xMin);
		double dy = ((double)points.size())/(yMax-yMin);
		return (dx+dy)/2;
	}
	
	public boolean add(MatchPoint p){
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
			//System.out.println("Adding point "+p.toString()+" to k-clump "+id);
			return true;
		} else {
			return false;
		}
	}
	
	public int size(){
		return points.size();
	}
	
	public Collection<MatchPoint> getMatchPoints(){
		return points;
	}
	
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
