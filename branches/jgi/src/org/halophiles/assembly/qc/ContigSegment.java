package org.halophiles.assembly.qc;

import java.util.HashMap;
import java.util.Map;
import java.util.Vector;

import org.halophiles.assembly.Contig;

public class ContigSegment {
	
	private static int COUNT = 0;
	
	/**
	 * The Contig this segment is on
	 */
	private Contig contig;
	/**
	 * The starting position of this segment in the contig
	 */
	private int start;
	/**
	 * The ending position of this segment in the contig
	 */
	private int end;
	
	/**
	 * The segment that the left end of of this segment connects to
	 */
	private Vector<ContigSegment> left;
	
	/**
	 * The segment that the right end of this segment connects to
	 */
	private Vector<ContigSegment> right;
	
	private Map<ContigSegment,Integer> ori;

	private boolean visited;
	
	private int id;
	
	
	public ContigSegment (Contig c, int start, int end){
		this.contig = c;
		this.start = start;
		this.end = end;
		this.visited = false;
		this.id = COUNT++;
		ori = new HashMap<ContigSegment, Integer>();
	}
	
	public void addLeftConnection(ContigSegment seg, int ori){
		left.add(seg);
		this.ori.put(seg,ori);
	}
	
	public void addRightConnection(ContigSegment seg, int ori){
		right.add(seg);
		this.ori.put(seg,ori);
	}
	
	public boolean inverted(){
		for (int i = 0; i < left.size(); i++) {
			if (ori.get(left.get(i)) != MatchPoint.RR)
				return false;
		}
		for (int i = 0; i < right.size(); i++) {
			if (ori.get(right.get(i)) != MatchPoint.FF)
				return false;
		}
		return true;
	}
	
	public Vector<ContigSegment> getLeftSegments(){
		return left;
	}
	
	public Vector<ContigSegment> getRightSegments(){
		return right;
	}
	
	public int getOri(ContigSegment seg){
		return ori.get(seg);
	}
	
	public String toString(){
		return contig.name+":"+Integer.toString(start)+"-"+Integer.toString(end);
	}
	
	public int hashCode(){
		return id;
	}
	
	public void setVisited(boolean visited){
		this.visited = visited;
	}
	
	public boolean getVisited(){
		return visited;
	}
	
	public int getStart(){
		return start;
	}
	
	public int getEnd(){
		return end;
	}
	
	public boolean isRepeat(){
		return left.size() > 1 || right.size() > 1;
	}
}
