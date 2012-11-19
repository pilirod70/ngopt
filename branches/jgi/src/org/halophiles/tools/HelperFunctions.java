package org.halophiles.tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.io.RandomAccessFile;
import java.text.NumberFormat;
import java.util.StringTokenizer;
import java.util.Vector;

public class HelperFunctions {

	public static String millisToTimeString(long time){
		NumberFormat nf = NumberFormat.getInstance();
		nf.setGroupingUsed(false);
		nf.setMaximumFractionDigits(2);
		long seconds = (time)/1000;
		int min = (int)(seconds/60);
		seconds = (seconds-min*60);
		nf.setMinimumIntegerDigits(2);
		return Integer.toString(min)+":"+nf.format(seconds);
	}

	/**
	 * Return true of the next character in this RandomAccessFile is equal 
	 * to the given character <code>c</code>
	 * @param raf the RandomAccessFile to check
	 * @param c the character value to check against
	 * @return true if the next character in <code>raf</code> is equal to <code>c</code>
	 * @throws IOException if an I/O error occurs while trying to read from the given RandomAccessFIle
	 */
	public static boolean nextCharIs(RandomAccessFile raf, char c) throws IOException{
		char next = (char) raf.read();
		raf.seek(raf.getFilePointer()-1);
		return next == c;
	}
	
	public static boolean nextCharIs(BufferedReader br, char c) throws IOException{
		if (!br.ready()){ return false; }
		boolean ret = false;
		br.mark(1);
		char b = (char) br.read();
		if (b == c) ret = true; 
		else ret = false;
		br.reset();
		return ret;
	}	
	/**
	 * Pad String <code>s</code> with <code>len</code> spaces on the right.
	 */
	public static String pad(String s, int len){
		String ret = new String(s);
		for (int i = 0; i < len-s.length(); i++){
			ret = ret+" ";
		}
		return ret;
	}

	/**
	 * Return mapping length for CIGAR String if there is a match (i.e. string contains 'M')
	 * @param cig the CIGAR string to parse
	 * @return the length of the match indicated by the CIGAR string. Return 0 if no match
	 */
	public static int cigarLength(String cig){
		StringTokenizer tok = new StringTokenizer(cig, "MIDNSHP", true);
		int alignLen = 0;
		while (tok.hasMoreTokens()){
			int len = Integer.parseInt(tok.nextToken());
			char op = tok.nextToken().charAt(0);
			switch (op) {
				case 'M': alignLen += len; break;
				case 'D': alignLen += len; break;
				default: continue;
			}
		}
		return alignLen;
	}

	/**
	 * Return the reference length for CIGAR String if there is a match (i.e. string contains 'M')
	 * @param cig the CIGAR string to parse
	 * @return the length of the match on the reference indicated by the CIGAR string. Return 0 if no match
	 */
	public static int cigarRefLength(String cig){
		StringTokenizer tok = new StringTokenizer(cig, "MIDNSHP", true);
		int alignLen = 0;
		while (tok.hasMoreTokens()){
			int len = Integer.parseInt(tok.nextToken());
			char op = tok.nextToken().charAt(0);
			switch (op) {
				case 'M': alignLen += len; break;
				case 'D': alignLen += len; break;
				default: continue;
			}
		}
		return alignLen;
	}
	
	/**
	 * Returns true of the 4th bit is set in the flag given in the SAM file
	 */
	public static boolean isReverse(int flag){
		if (getBit(4,flag) == 1) return true;
			else return false;
	}
	
	/**
	 * get the bit from the flag given in a SAM file
	 */
	public static int getBit (int bit, int flag) { 
	   int mod = 0;
	   int dig = 0;
	   while( flag != 0 && dig <= bit) {  
		   mod = flag % 2;
		   flag = flag / 2;
		   dig++;
	   }
	   return mod;  
	}
	/**
	 * Strip the pair number off a fastq header
	 */
	public static String trimPairNumber(String s){
		if (s.contains("/")){
			return s.substring(0,s.lastIndexOf("/"));
		} else
			return s;
	}
	
	public static String basename(String path, String suffix){
		int pos = path.indexOf(suffix);
		if (pos < 0)
			return path.substring(path.lastIndexOf('/')+1);
		else 
			return path.substring(path.lastIndexOf('/')+1,pos);
	}
	
	public static String basename(String path){
		return path.substring(path.lastIndexOf('/')+1);
	}
	
	public static String dirname(String path){
		int pos = path.indexOf("/");
		if (pos < 0)
			return ".";
		else 
			return path.substring(0,path.lastIndexOf('/'));
	}
	
	public static boolean isEmpty (File file){
		if (!file.exists())
			return true;
		try {
			BufferedReader br = new BufferedReader(new FileReader(file));
			if (br.ready())
				return false;
			else
				return true;
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.exit(-1);
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return false;
		
	}
	
	public static void logInputs (String mainClass, String[] args) {
		System.out.print("[a5_qc] "+mainClass);
		if (args.length > 0){
			System.out.println(" "+args[0]);
			String space = "";
			for (int i = 0; i < mainClass.length(); i++)
				space += " ";
			for (int i = 1; i < args.length; i++){
				System.out.println("        "+space+" "+args[i]);
			}
		} else {
			System.out.println();
		}
	}
	
	public static int[] toArray(Vector<Integer> v) {
		int[] ret = new int[v.size()];
		for (int i = 0; i < ret.length; i++)
			ret[i] = v.get(i);
		return ret;
	}
	
	public static boolean hasNonZeroSize(String path){
		File file = new File(path);
		if (file.exists())
			if (isEmpty(file))
				return false;
			else
				return true;
		else
			return false;
	}
	
	public static PrintStream openIfClosed(PrintStream out, String path) throws IOException{
		if (out == null){
			File file = new File(path);
			file.createNewFile();
			return new PrintStream(file);
		} else
			return out;
		
			
	}
	
	public static boolean isUnmapped(int flag){
		return getBit(2, flag) == 1;
	}
}
