package org.halophiles.assembly.qc;

import org.halophiles.assembly.Contig;

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
}
