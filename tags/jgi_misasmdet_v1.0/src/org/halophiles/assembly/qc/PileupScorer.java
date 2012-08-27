package org.halophiles.assembly.qc;

import java.util.NoSuchElementException;
import java.util.regex.Pattern;

public class PileupScorer {
	
	private static final Pattern MISMATCH_PATTERN = Pattern.compile("atcgATCG");
	private static final Pattern GAP_PATTERN = Pattern.compile("[+,-][0-9]+[ACGTNacgtn]+");
	private static final Pattern CLIP_PATTERN = Pattern.compile("^[^,$]");
	private static final Pattern MATCH_PATTERN = Pattern.compile("[,.]");

	/**
	 * The gain for a match
	 */
	private double match;
	/**
	 * The mismatch penalty 
	 */
	private double mismatch;
	/**
	 * The penalty for opening a gap
	 */
	private double gapOpen;
	/**
	 * The penalty for extending a gap
	 */
	private double gapExtend;
	/**
	 * The penalty for truncating a read 
	 */
	private double clip;
	/**
	 * Score errors if we have non-zero penalties.
	 * i.e. don't bother adding or subtracting if this is a simple scorer 
	 * and our penalties will not have an effect on score.
	 */
	private boolean scoreErrors;
	
	public PileupScorer (){
		this(1,0,0,0,0);
	}
	
	public PileupScorer(double match, double mismatch, double gapOpen, double gapExtend, double clip){
		this.match = match;
		this.mismatch = mismatch;
		this.gapOpen = gapOpen;
		this.gapExtend = gapExtend;
		this.clip = clip;
		scoreErrors = this.mismatch != 0 && this.gapOpen != 0 && this.gapExtend != 0 && this.clip != 0;
	}
	
	public double scorePileup(String bases) {
		PileupBaseParser pbp = new PileupBaseParser(bases);
		String tmpBase = null;
		double score = 0;
		while(pbp.hasNextBase()){
			tmpBase = pbp.nextBase();
			if (MATCH_PATTERN.matcher(tmpBase).matches()){
				score += match;
			} else if (scoreErrors) {
				if (MISMATCH_PATTERN.matcher(tmpBase).matches()) {
					score -= mismatch;
				} else if (GAP_PATTERN.matcher(tmpBase).matches()) {
					score -= gapOpen;
				} else if (tmpBase.equals("*")) {
					score -= gapExtend;
				} else if (CLIP_PATTERN.matcher(tmpBase).matches()){
					score -= clip;
				}
			}
		}
		return score;
	}
	/**
	 * A class for parsing a the base column of a pileup entry.
	 */
	public static class PileupBaseParser {
		private char[] bases;
		private int idx;
		
		public PileupBaseParser(String bases){
			this.bases = bases.toCharArray();
			this.idx = 0;
		}
		
		public boolean hasNextBase(){
			return idx < bases.length;
		}
		
		public String nextBase() throws NoSuchElementException{
			if (idx == bases.length)
				throw new NoSuchElementException("No more bases to parse");
			String ret = "";
			if (isAlnMatch(bases[idx])){ // match to the reference
				ret = ret + bases[idx++];
				if (idx < bases.length && bases[idx] == '$') // check to see if this is the end of a read
					ret = ret + bases[idx++];
			} else if (bases[idx] == '^'){ // read starts here
				ret = ret + bases[idx++] + bases[idx++] + bases[idx++];
			} else if (bases[idx] == '-' || bases[idx] == '+'){
				ret = ret + bases[idx++];
				// extract the length of the indel, so we know how many characters to parse after the length portion
				String lenStr = "";
				do {
					lenStr += bases[idx++];
				} while (bases[idx] >= 48 && bases[idx] <= 57); // numeric characters, i.e. the length of the indel
				ret = ret + lenStr;
				int len = Integer.parseInt(lenStr);
				for (int i = 0; i < len; i++)
					ret = ret + bases[idx++];
			} else if (bases[idx] == '*'){
				ret = ret + bases[idx++];
			}
			if (ret.length() == 0)
				System.out.print("");
			return ret;
		}
		
		private boolean isAlnMatch (char c){
			return c == ',' || c == '.' || 
					c == 'g' || c == 'G' ||
					c == 'c' || c == 'C' ||
					c == 'a' || c == 'A' ||
					c == 't' || c == 'T' ||
					c == 'n' || c == 'N';
			
		}
		
	}
}
