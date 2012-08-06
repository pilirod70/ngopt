package org.halophiles.assembly.qc;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Vector;

import org.halophiles.tools.HelperFunctions;

public class MisassemblyRange {
	private String contig;
	private RunningStat[] stats;
	private int[] left, right, junction;
	private int[] minPos;
	MisassemblyRange (String contig, int[] left, int[] right) {
		if (left.length != right.length) throw new IllegalArgumentException("left and right must be the same length");
		this.contig = contig;
		this.left = left;
		this.right = right;
		this.minPos = new int[left.length];
		this.stats = new RunningStat[left.length];
		/*
		 *  [0] 'weakest' position in each range
		 *  [1] score of 'weakest' position
		 *  [2] mean for this range
		 *  [3] stdev for this range
		 *  [4] number of elements for this range
		 */
		for (int i = 0; i < stats.length; i++)
			stats[i] = new RunningStat();
	}
	
	boolean addPos(int pos, double score){
		int idx = Arrays.binarySearch(right, pos);
		if (idx < 0)
			idx = -1 * (idx + 1);
		if (pos < left[idx])
			return false;
		double min = stats[idx].min();
		if (score < min)
			minPos[idx] = pos;
		stats[idx].addVal(score);
		return true;
	}
	
	boolean contains(int pos){
		int idx = Arrays.binarySearch(right, pos);
		if (idx < 0)
			idx = -1 * (idx + 1);
		if (pos < left[idx])
			return false;
		else
			return true;
	}
	
	void printState(PrintStream out) {
		for (int i = 0; i < left.length; i++){
			out.println(left[i]+"\t"+right[i]+"\t"+stats[i].toString());
		}
	}
	
	int[] getJunctions(double maxScore){
		Vector<Integer> scores = new Vector<Integer>();
		for (int i = 0; i < minPos.length; i++){
			if (stats[i].min() < maxScore)
				scores.add(minPos[i]);
		}
		return HelperFunctions.toArray(scores);
	}
	
	private class RunningStat {
		private int n;
		private double oldM, newM, oldS, newS, min, max;
		public RunningStat() {
			n = 0;
			min = Double.POSITIVE_INFINITY;
			max = Double.NEGATIVE_INFINITY;
		}
		void addVal(double x){
			n++;
			if (n == 1){
				oldM = x;
				newM = x;
				oldS = 0;
				newS = 0;
			} else {
				newM = oldM + (x-oldM)/n;
				newS = oldS + (x-oldM)*(x-newM);
				oldM = newM;
				oldS = newS;
			}
			min = x < min ? x : min;
			max = x > max ? x : max;
		}
		int numPts(){
			return n;
		}
		double mean() {
			return (n > 0) ? newM : 0.0;
		}
		double variance() {
			return (n > 1) ? newS/(n-1) : 0.0;
		}
		double sd() {
			return Math.sqrt(variance());
		}
		double min(){
			return min;
		}
		double max(){
			return max;
		}
		
		public String toString() {
			return mean()+"\t"+sd()+"\t"+min+"\t"+max+"\t"+n;
		}
	}	
}