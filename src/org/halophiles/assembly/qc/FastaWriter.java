package org.halophiles.assembly.qc;

import java.io.BufferedWriter;
import java.io.IOException;

public class FastaWriter {
	
	
	private BufferedWriter writer;
	
	private int width;
	
	private int currChar;
	
	public FastaWriter(BufferedWriter writer) {
		this(writer,80);
	}
	
	public FastaWriter(BufferedWriter writer, int width) {
		this.writer = writer;
		this.width = width;
		this.currChar = 0;
	}
	
	public void addNewSequence (String header) throws IOException{
		if (currChar % width > 0)
			writer.write("\n>" + header + "\n");
		else
			writer.write(">" + header + "\n");
	}
	
	public void writeSequence (CharSequence sequence) throws IOException {
		for (int i = 0; i < sequence.length(); i++){
			writer.write(sequence.charAt(i));
			currChar++;
			if (currChar % width == 0)
				writer.write('\n');
		}
	}
	
	public void writeSequenceInverted (CharSequence sequence) throws IOException {
		for (int i = sequence.length()-1; i >= 0; i--){
			writer.write(comp(sequence.charAt(i)));
			currChar++;
			if (currChar % width == 0)
				writer.write('\n');
		}
	}

	private static int comp(int b){
		switch (b) {
			case 'A': return 'T';
			case 'T': return 'A';
			case 'G': return 'C';
			case 'C': return 'G';
			case 'a': return 't';
			case 't': return 'a';
			case 'g': return 'c';
			case 'c': return 'g';
			case 'N': return 'N';
			case 'n': return 'n';
			default: return 'N';
		}
	}
}
