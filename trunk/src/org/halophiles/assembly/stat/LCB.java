package org.halophiles.assembly.stat;

import java.io.PrintStream;

public class LCB implements Comparable<LCB> {
	
	private static int COUNT = 0;
	
	private int lcbIdx;
	
	private char[][] seqs;
	
	private int seq1Left;
	private int seq1Right;
	private boolean comp1;
	
	private int seq2Left;
	private int seq2Right;
	private boolean comp2;
	
	
	private boolean linSpec;
	
	private int taxaIdx;
	
	public LCB(char[][] seqs, int s1L, int s1R, boolean s1Comp, int s2L, int s2R, boolean s2Comp){
		seq1Left = s1L; 
		seq1Right = s1R;
		comp1 = s1Comp;
		seq2Left = s2L;
		seq2Right = s2R;
		comp1 = s2Comp;
		this.seqs = seqs;
		lcbIdx = COUNT++;
		linSpec = false;
	}
	
	public LCB(char[][] seqs, int s1L, int s1R, boolean s1Comp){	
		if (seqs[0].length == 0){
			taxaIdx = 1;
			seq2Left = s1L;
			seq2Right = s1R;
			comp2 = s1Comp;
		} else {
			taxaIdx = 0;
			seq1Left = s1L;
			seq1Right = s1R;
			comp1 = s1Comp;
		}
		linSpec = true;
		this.seqs = seqs;
		lcbIdx = COUNT++;
	}
	
	public LCB(char[] seq){
		this.seqs = new char[1][];
		this.seqs[0] = seq;
		linSpec = false;
	}
	
	public boolean isLinSpec(){
		return linSpec;
	}
	
	public char[][] getSeqs(){
		return seqs;
	}
	
	public int getSeqLen(){
		return seqs[taxaIdx].length;
	}
	
	public int getSeq1Len(){
		return seqs[0].length;
	}
	
	public int getSeq2Len(){
		return seqs[1].length;
	}
	
	public void printSum(PrintStream out){
		if (linSpec){
			out.println("Taxa " + (taxaIdx+1) + ": "+ seqs[taxaIdx].length);
		} else {
			out.println("Taxa 1: " + seqs[0].length);
			out.println("Taxa 2: " + seqs[1].length);
		}
	}

	public void setSeq1Left(int seq1Left) {
		this.seq1Left = seq1Left;
	}

	public int getSeq1Left() {
		return seq1Left;
	}

	public void setSeq1Right(int seq1Right) {
		this.seq1Right = seq1Right;
	}

	public int getSeq1Right() {
		return seq1Right;
	}

	public void setComp1(boolean comp1) {
		this.comp1 = comp1;
	}

	public boolean isComp1() {
		return comp1;
	}

	public void setSeq2Left(int seq2Left) {
		this.seq2Left = seq2Left;
	}

	public int getSeq2Left() {
		return seq2Left;
	}

	public void setSeq2Right(int seq2Right) {
		this.seq2Right = seq2Right;
	}

	public int getSeq2Right() {
		return seq2Right;
	}

	public void setComp2(boolean comp2) {
		this.comp2 = comp2;
	}

	public boolean isComp2() {
		return comp2;
	}

	@Override
	public int compareTo(LCB o) {
		return this.lcbIdx - o.lcbIdx;
	}

}
