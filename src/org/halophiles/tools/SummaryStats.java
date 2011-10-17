package org.halophiles.tools;

import java.util.Comparator;

import jsc.correlation.KendallCorrelation;
import jsc.correlation.SpearmanCorrelation;
import jsc.datastructures.PairedData;

public class SummaryStats {
	
	private static class IntegerSorter implements Comparator<Integer>{
		public IntegerSorter(int[] dat){
			this.dat = dat;
		}
		public int compare(Integer a, Integer b){
			return dat[a.intValue()] - dat[b.intValue()];
		}
		private int[] dat;
	}
	
	private static class DoubleSorter implements Comparator<Double>{
		public DoubleSorter(double[] dat){
			this.dat = dat;
		}
		public int compare(Double a, Double b){
			return (int) (dat[a.intValue()] - dat[b.intValue()]);
		}
		private double[] dat;
	}
	
	public static double mean(double[] dat) {
		return sum(dat) / ((double)dat.length);
	}
	
	public static double mean(Double[] dat) {
		return sum(dat) / ((double)dat.length);
	}
	
	public static double sum(Double[] dat) {
		double total = 0;
		for (int i = 0; i < dat.length; i++) {
			total += dat[i];
		}
		return total;
	}
	
	public static double sum(double[] dat) {
		double total = 0;
		for (int i = 0; i < dat.length; i++) {
			total += dat[i];
		}
		return total;
	}
	
	public static int sum(Integer[] dat) {
		int total = 0;
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
	
	public static int mean(int[] dat){
		return sum(dat)/dat.length;
	}
	
	public static int mean(Integer[] dat){
		return sum(dat)/dat.length;
	}
	
	public static double variance(double[] dat, double mean) {
		double total = 0;
		for (int i = 0; i < dat.length; i++) {
			total += Math.pow(dat[i]-mean, 2);
		}
		return total / ((double)(dat.length-1));
	}
	
	public static double covariance(double[] dat1, double[] dat2){
		double mean1 = mean(dat1);
		double mean2 = mean(dat2);
		return covariance(dat1, mean1, dat2, mean2);
	}
	
	public static double covariance(double[] dat1, double mean1, double[] dat2, double mean2) {
		if (dat1.length != dat2.length)
			throw new IllegalArgumentException("dat1 and dat2 must be the same length");
		double total = 0;
		for (int i = 0; i < dat1.length; i++) {
			total += (dat1[i]-mean1)*(dat2[i]-mean2);
		}
		return total / ((double)(dat1.length-1));
	}
	
	public static double pearson(double[] dat1, double[] dat2){
		return pearson(dat1,mean(dat1),dat2,mean(dat2));
	}
	
	public static double pearson(double[] dat1, double mean1, double[] dat2, double mean2) {
		if (dat1.length != dat2.length)
			throw new IllegalArgumentException("dat1 and dat2 must be the same length");
		return covariance(dat1, mean1, dat2, mean2)/Math.sqrt(variance(dat1,mean1)*variance(dat2,mean2));
	}
	
	
	public static double spearman(double[] dat1, double[] dat2){
		if (dat1.length != dat2.length)
			throw new IllegalArgumentException("dat1 and dat2 must be the same length");
		return new SpearmanCorrelation(new PairedData(dat1,dat2)).getR();
	}
	
	public static double kendall(double[] dat1, double[] dat2){
		if (dat1.length != dat2.length)
			throw new IllegalArgumentException("dat1 and dat2 must be the same length");
		
		if (dat1.length < 2)
			throw new IllegalArgumentException("dat1 must contain at least 2 elements");
		if (dat2.length < 2)
			throw new IllegalArgumentException("dat2 must contain at least 2 elements");
		return new KendallCorrelation(new PairedData(dat1,dat2)).getR();
	}
	
	public static double covariance(int[] dat1, int[] dat2){
		int mean1 = mean(dat1);
		int mean2 = mean(dat2);
		return covariance(dat1, mean1, dat2, mean2);
	}
	
	public static double covariance(int[] dat1, int mean1, int[] dat2, int mean2) {
		if (dat1.length != dat2.length)
			throw new IllegalArgumentException("dat1 and dat2 must be the same length");
	
		double total = 0;
		for (int i = 0; i < dat1.length; i++) {
			total += (dat1[i]-mean1)*(dat2[i]-mean2);
		}
		return total / ((double)(dat1.length-1));
	}

	
	public static double variance(Double[] dat, double mean) {
		double total = 0;
		for (int i = 0; i < dat.length; i++) {
			total += Math.pow(dat[i]-mean, 2);
		}
		return total / ((double)(dat.length-1));
	}
	
	public static double variance(Double[] dat) {
		return variance(dat,mean(dat));
	}
	
	public static double variance(double[] dat) {
		return variance(dat,mean(dat));
	}
	
	public static double skew(Double[] dat, double mean, double var) {
		double m3 = 0;
		double n = dat.length;
		for (int i = 0; i < dat.length; i++) {
			m3 += Math.pow(dat[i]-mean, 3);
		}
		m3 = m3/n;
		double g1 = m3/(Math.pow(var, 1.5));
		return Math.sqrt(n*(n-1.0))*g1/(n-2);
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
	
	public static double skew(Double[] dat) {
		double mean = mean(dat);
		double var = variance(dat,mean);
		return skew(dat,mean,var);
	}
	
	public static double skew(double[] dat) {
		double mean = mean(dat);
		double var = variance(dat,mean);
		return skew(dat,mean,var);
	}
	
	public static double max(double[] dat){
		double max = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < dat.length; i++)
			if (dat[i]>max)
				max = dat[i];
		return max;
		
	}
	
	public static double min(double[] dat){
		double min = Double.POSITIVE_INFINITY;
		for (int i = 0; i < dat.length; i++)
			if (dat[i]<min)
				min = dat[i];
		return min;
		
	}

	public static int max(int[] dat){
		int max = Integer.MIN_VALUE;
		for (int i = 0; i < dat.length; i++)
			if (dat[i]>max)
				max = dat[i];
		return max;
		
	}
	
	public static int min(int[] dat){
		int min = Integer.MAX_VALUE;
		for (int i = 0; i < dat.length; i++)
			if (dat[i]<min)
				min = dat[i];
		return min;
		
	}

}
