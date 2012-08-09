package org.halophiles.assembly.qc;

import java.io.PrintStream;
import java.util.Vector;

import org.halophiles.assembly.Contig;

public class MisassemblyRange {
	
	private static int COUNT = 0;
	
	private Contig contig;
	private String id;
	private MisassemblyBlock leftBlock, rightBlock;
	private int left, right, minPos;

	
	private int n;
	private double oldM, newM, oldS, newS, min, max;
	
	public MisassemblyRange (Contig contig, MisassemblyBlock leftBlock, MisassemblyBlock rightBlock) {
		this.contig = contig;
		this.leftBlock = leftBlock;
		this.rightBlock = rightBlock;
		this.minPos = -1;
		if (leftBlock.getRight() == rightBlock.getLeft()){
			left = leftBlock.getRight();
			right = left;
		} else if (leftBlock.getRight() > rightBlock.getLeft()) {
			left = rightBlock.getLeft();
			right = leftBlock.getRight();
		} else {
			left = leftBlock.getRight();
			right = rightBlock.getLeft();
		}
		// expand out because our pairs aren't completely perfect.
		left -= 25;
		right += 25;
		
		// attach a unique id to this range
		id = "r"+Integer.toString(COUNT++);
		
		// stats information
		n = 0;
		min = Double.POSITIVE_INFINITY;
		max = Double.NEGATIVE_INFINITY;
		
	}
	
	public MisassemblyBlock getLeftBlock(){
		return leftBlock;
	}
	
	public MisassemblyBlock getRightBlock(){
		return rightBlock;
	}
	
	public boolean addPos(int pos, double score){
		if (!contains(pos))
			return false;
		if (addVal(score)){
			minPos = pos;
		}
		return true;
	}
	
	public boolean contains(int pos){
		return pos >= left && pos <= right;
	}
	
	public void printState(PrintStream out) {
		out.println(mean()+"\t"+sd()+"\t"+min+"\t"+max+"\t"+n);
	}
	
	public int getMinPos(){
		return minPos;
	}
	
	public double getMinScore(){
		return min;
	}
	
	public String toString(){
		return contig+"\t"+left+"\t"+right+"\t"+id;
	}
	
	public String getId(){
		return id;
	}
	
	/**
	 * Add a new value to the running stats, and return whether or not
	 * we've encountered a new minimum value
	 * 
	 * @param x the value to add toe the running stats
	 * @return true if x is the new minimum value, false otherwise
	 */
	private boolean addVal(double x){
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
		max = x > max ? x : max;
		if (x < min ){
			min = x;
			return true;
		} else
			return false;
	}
	private int numPts(){
		return n;
	}
	private double mean() {
		return (n > 0) ? newM : 0.0;
	}
	private double variance() {
		return (n > 1) ? newS/(n-1) : 0.0;
	}
	private double sd() {
		return Math.sqrt(variance());
	}
	private double min(){
		return min;
	}
	private double max(){
		return max;
	}
	
	/**
	 * 
	 * @param ranges a SORTED Vector of MisassemblyRanges
	 * @param pos the position we are looking for
	 * @return
	 */
	public static int binarySearch(Vector<MisassemblyRange> ranges, int pos){
		int max = ranges.size() - 1;
		int min = 0;
		int mid = (max+min)/2;
		while (max > min){
			if (pos < ranges.get(mid).right)
				max = mid;
			else if (pos > ranges.get(mid).right)
				min = mid+1;
			else
				return mid;
			mid = (max+min)/2;
		}
		if (pos >= ranges.get(mid).left)
			return mid;
		else
			return -1;
	}
}