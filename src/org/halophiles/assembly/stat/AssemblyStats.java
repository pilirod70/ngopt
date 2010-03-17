package org.halophiles.assembly.stat;

import java.util.Iterator;

import java.util.Vector;


public class AssemblyStats {

	private static final int A = 0;
	private static final int C = 1;
	private static final int T = 2;
	private static final int G = 3;
	
	private static final String[] bases = {"A","C","T","G"};
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
			if (args.length != 1){
				System.err.println("Usage: java AssemblyStats <xmfa_file>");
				System.exit(-1);
			}
			
			
			
			int[][] SUBS = new int[4][4];
			
			XMFALoader xmfal = new XMFALoader(args[0]);
			System.out.println("Sequence 1 length: " + xmfal.numCharSeq1());
			System.out.println("Sequence 2 length: " + xmfal.numCharSeq2());
			
			Vector<LCB> lcbs = xmfal.getLCBs();
			
			Iterator<LCB> it = lcbs.iterator();
			
			System.out.print("Counting substitutions... ");
			while(it.hasNext()){
				LCB tmp = it.next();
				//tmp.printSum(System.out);
				if (!tmp.isLinSpec()){
					char[][] seq = tmp.getSeqs();
					for (int i = 0; i < seq[0].length; i++){
						if (seq[0][i] == '-' || seq[1][i] == '-'){
							continue;
						} else {
							SUBS[getIdx(seq[0][i])][getIdx(seq[1][i])]++;
						}
					}
				}
			}
			System.out.print("done!\n");
			StringBuilder sb = new StringBuilder();
			for (int i = 0; i < bases.length; i++){
				sb.append("\t"+bases[i]);
			}
			sb.append("\n");
			for (int i = 0; i < SUBS.length; i++){
				sb.append(bases[i]);
				for (int j = 0; j < SUBS.length; j++){
					sb.append("\t"+SUBS[i][j]);
				}
				sb.append("\n");
			}
			System.out.print(sb.toString());

	}
	

	private static int getIdx(char c){
		switch(c) {
		  case 'a': return A; 
		  case 'A': return A;
		  case 'c': return C;
		  case 'C': return C;
		  case 't': return T;
		  case 'T': return T;
		  case 'g': return G;
	      case 'G': return G;
		  default: { System.err.print(c); return -1;}
		}
	}
	
	

}
