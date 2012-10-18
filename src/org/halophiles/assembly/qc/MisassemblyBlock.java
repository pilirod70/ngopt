package org.halophiles.assembly.qc;

import org.halophiles.assembly.Contig;

/**
 * This class represents a collection of reads with the same orientation that have been 
 * grouped by spatial clustering on read-pairs. 
 * Within each instance of MisassemblyBlock is a member variable representing the block that this block
 * is connected to. 
 * 
 * @author ajtritt
 *
 */
public class MisassemblyBlock {
	
	private static int COUNT = 0;
	
	private Contig contig;
	private int left;
	private int right;
	private boolean rev;
	
	private int id;
	private MisassemblyBlock connection;
	
	public MisassemblyBlock (Contig contig, int left, int right, boolean rev){
		this.contig = contig;
		this.left = left;
		this.right = right;
		this.rev = rev;
		id = COUNT++;
	}
	
	public MisassemblyBlock (Contig contig, int left, int right, boolean rev, int id){
		this.contig = contig;
		this.left = left;
		this.right = right;
		this.rev = rev;
		this.id = id;
	}
	
	public void addConnection(MisassemblyBlock connection){
		this.connection = connection;
	}
	
	public MisassemblyBlock getConnection(){
		return this.connection;
	}
	
	public int getLeft(){
		return left;
	}
	
	public int getRight(){
		return right;
	}
	
	public boolean getRev(){
		return rev;
	}
	
	public int hashCode(){
		return id;
	}
	
	public String toString(){
		return id+","+contig+","+(rev?"-":"+")+","+Integer.toString(left)+","+Integer.toString(right);
	}
	
	public Contig getContig(){
		return contig;
	}
	
	public int getId(){
		return id;
	}
	
	/**
	 * Returns -1, if this block lies at the left end of the contig it is on,
	 * +1 if this block lies on the right end of the contig, and 0 otherwise.
	 * 
	 * @return -1, 0, 1 for left, middle, right end of the contig, respectively
	 */
	public static int getTerminus(MisassemblyBlock block){
		double perc = ((double) block.right - block.left)/block.right;
		if (perc > 0.5)
			return -1;
		perc = ((double) block.right - block.left)/(block.contig.len-block.left+1);
		if (perc > 0.5)
			return 1;
		return 0;
	}
}
