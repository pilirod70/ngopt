package org.halophiles.assembly.stat;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Vector;


public class XMFALoader {
	
	private int[] seqLen;
	
	private Vector<LCB> lcbs;
	
	
	public XMFALoader(String xmfaFile){
		FileReader fr = null;
		lcbs = new Vector<LCB>();
		try {
			fr = new FileReader(xmfaFile);
		} catch (Exception e){
			e.printStackTrace();
			System.exit(-1);
		}
		BufferedReader br = new BufferedReader(fr);
		seqLen = new int[2];
		int numLCB = 0;
		try {
			// first look over the file to see what we're dealing with
			int numTaxa = 0;
			
			while (br.ready()){
				br.mark(1000);
				char c = (char) br.read();
				if (c == '#'){                    // lets check the comments to make sure theres only 2 sequences here
					String line = br.readLine();
					if (line.startsWith("Sequence")){
						//line = line.substring(1);
						int endSeqNum = line.indexOf("F");
						int seqNum = Integer.parseInt(line.substring(8, endSeqNum));
						if (seqNum > 2){
							System.err.println("Too many taxa. For scoring assemblies, you should be doing only a pairwise alignment");
							System.exit(-1);
						} else {
							numTaxa++;
						}
					} else if (line.startsWith("Format")){
						if (!line.split(" ")[1].equalsIgnoreCase("mauve1")){
							System.err.println(xmfaFile + " is in unrecognized Mauve format. Naively proceeding...");
						}
					}
					
				} else {
					br.reset();
					while(br.ready()) {
						System.out.print("Found new LCB. Counting lengths....");
						br.mark(20000000);
						LCBTemplate temp = calcLen(br);
						System.out.println("done!");
						if (temp.seq1Len == 0 && temp.seq2Len == 0){
							break;
						}
						br.reset();
					
						char[][] seqs = new char[2][];
						seqs[0] = new char[temp.seq1Len];
						seqs[1] = new char[temp.seq2Len];
						loadSeqs(br,seqs);
						LCB tmp = temp.createLCB(seqs);
						lcbs.add(tmp);
						tmp.printSum(System.out);
					}
					
				}
				
			}
		
			
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		System.out.println(xmfaFile + " successfully loaded!");
		
	}
	
	private LCBTemplate calcLen(BufferedReader br) throws IOException {
		int[] len = new int[2]; // if -1 in return, there's no sequence for this taxa on this LCB
		int currSeq = 0;
		char c = '+';
		LCBTemplate temp = new LCBTemplate();
		while(c != '=' && br.ready()){
			c = (char) br.read();
			if (c == '>'){
				br.read(); // get rid of that annoying space
				currSeq = temp.addSeq(br.readLine());
				len[currSeq] = 0; 
			} else if (c == '\n') {
				continue;
			} else if (c != '='){
				len[currSeq]++;
			}
		}
		temp.seq1Len = len[0];
		temp.seq2Len = len[1];
		return temp;
	}
	
	private void loadSeqs(BufferedReader br, char[][] seq) throws IOException{
		int currSeq = 0;
		char c = '>';
		int i = 0;
		//System.err.print("First few characters: ");
		while(c != '='){
			c = (char) br.read();
			
			if (c == '>'){
				br.read(); // get rid of that annoying space
				String line = br.readLine();
				currSeq = Integer.parseInt(line.substring(0,line.indexOf(":")))-1; // figure out which sequence we're on
				i = 0; // reset our counter
			} else if (c == '\n') {
				continue;
			} else if (c != '='){
				
				/*if (i < 5){
					System.err.print(c+" ");
				}
				if (i== 5){
					System.err.println();
				}*/
				seq[currSeq][i] = c;
				i++;
			}
		}
	}
	
	public int numCharSeq1(){
		return seqLen[0];
	}
	
	public int numCharSeq2(){
		return seqLen[1];
	}
	
	public Vector<LCB> getLCBs(){
		return lcbs;
	}
	
	class LCBTemplate {
		int seq1Left;
		int seq1Right;
		int seq1Len;
		boolean comp1;
		
		int seq2Left;
		int seq2Right;
		int seq2Len;
		boolean comp2;
		// 0   1        2    3 
		// 2:3374131-3920004 + 454LargeContigs-Reorder.fas
		
		LCBTemplate(){
			seq1Left = -1;
			seq1Right = -1;
			seq1Len = -1;
			seq2Left = -1;
			seq2Right = -1;
			seq2Len = -1;
			comp1 = false;
			comp2 = false;
		}
			
		int addSeq(String header){
			String[] tmp = header.split("[: -]"); 
			int seq = Integer.parseInt(tmp[0])-1;
			if (seq == 0){
				seq1Left = Integer.parseInt(tmp[1]);
				seq1Right = Integer.parseInt(tmp[2]);
				comp1 = tmp[3].equals("-");
			} else if (seq == 1) {
				seq2Left = Integer.parseInt(tmp[1]);
				seq2Right = Integer.parseInt(tmp[2]);
				comp2 = tmp[3].equals("-");
			} else {
				throw new IllegalArgumentException("seq : " + seq);
			}
			return seq;
		}
		
		LCB createLCB(char[][] seqs){
			if (seq1Left < 0) {
				return new LCB(seqs, seq2Left, seq2Right, comp2);
			} else if (seq2Left < 0) {
				return new LCB(seqs, seq1Left, seq1Right, comp1);
			} else {
 				return new LCB(seqs, seq1Left, seq1Right, comp1, seq2Left, seq2Right, comp2);
			}
		}
		
	}
	

}
