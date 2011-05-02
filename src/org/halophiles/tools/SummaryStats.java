package org.halophiles.tools;

public class SummaryStats {
	
	
	public static double mean(double[] dat) {
		return sum(dat) / ((double)dat.length);
	}
	
	public static double sum(double[] dat) {
		double total = 0;
		for (int i = 0; i < dat.length; i++) {
			total += dat[i];
		}
		return total;
	}
	
	public static int sum(int[] dat) {
		int total = 0;
		for (int i = 0; i < dat.length; i++) {
			total += dat[i];
		}
		return total;
	}
	
	public static double variance(double[] dat, double mean) {
		double total = 0;
		for (int i = 0; i < dat.length; i++) {
			total += Math.pow(dat[i]-mean, 2);
		}
		return total / ((double)dat.length);
	}
	
	public static double variance(double[] dat) {
		return variance(dat,mean(dat));
	}
	
	public static double skew(double[] dat, double mean, double var) {
		double m3 = 0;
		double n = dat.length;
		for (int i = 0; i < dat.length; i++) {
			m3 += Math.pow(dat[i]-mean, 3);
		}
		m3 = m3/n;
		double g1 = m3/(Math.pow(var, 1.5));
		return Math.sqrt(n*(n-1.0))*g1/(n-2);
	}
	
	public static double skew(double[] dat) {
		double mean = mean(dat);
		double var = variance(dat,mean);
		return skew(dat,mean,var);
	}
	
	

}
