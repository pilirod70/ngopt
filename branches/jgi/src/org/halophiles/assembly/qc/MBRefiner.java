package org.halophiles.assembly.qc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Vector;

import org.halophiles.tools.HelperFunctions;

public class MBRefiner {

	private static final String SAMTOOLS = "/jgi/tools/bin/samtools";

	public static void main(String[] args) {
		HelperFunctions.logInputs("MBRefiner", args);
		if (args.length == 4) {
			try {
				String bamPath = args[0];
				String bedPath = args[1];
				String ctgPath = args[2];
				String brokenPath = args[3];

				Map<String, MisassemblyRange> ranges = scoreAtBaseLevel(bamPath, bedPath, ctgPath);
				
				Map<String, int[]> junctions = refine(ranges);
				breakContigs(junctions, ctgPath, brokenPath);

			} catch (IOException e) {
				e.printStackTrace();
				System.exit(-1);
			}	
		} else  
		/*if (args.length == 4) {
			try {
				String bedPath = args[0];
				String plpPath = args[1];
				String ctgPath = args[2];
				String brokenPath = args[3];
	
				Map<String, int[]> junctions = refine(bedPath, plpPath);
				breakContigs(junctions, ctgPath, brokenPath);
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(-1);
			}	
		} else*/ {
			System.out.println("Usage: MBRefiner <bam_file> <bed_file> <contig_file> <broken_ctgs_file>");
			System.exit(-1);
		}
	}
	
	public static Map<String,MisassemblyRange> scoreAtBaseLevel(String bedPath, String plpPath) throws IOException {
		Map<String, MisassemblyRange> regions = getRegions(bedPath);
		PileupEvaluator plpEval = new PileupEvaluator(regions, new PileupScorer());
		BufferedReader br = new BufferedReader(new FileReader(new File(plpPath)));
		while(br.ready()){
			plpEval.handleLine(br.readLine());
		}
		return regions;
	}
	
	public static Map<String, MisassemblyRange> scoreAtBaseLevel(String bamPath, String bedPath, String ctgPath) throws IOException {
		Map<String, MisassemblyRange> regions = getRegions(bedPath);
		try {
			runMPileup(bamPath, bedPath, ctgPath, regions);
			// runAndReadMpileup(bamPath, bedPath, ctgPath, regions);
			// runMPileupWApache(bamPath, bedPath, ctgPath, regions);
		} catch (InterruptedException e) {
			System.err.println("Unable to run samtools mpileup. This is what I know:");
			e.printStackTrace();
			System.exit(-1);
		}
		return regions;
	}
	
	public static Map<String, int[]> refine(Map<String, MisassemblyRange> regions) throws IOException{
		Map<String, int[]> ret = new HashMap<String, int[]>();
		Iterator<String> ctgIt = regions.keySet().iterator();
		String tmpCtg;
		MisassemblyRange reg;
		while (ctgIt.hasNext()) {
			tmpCtg = ctgIt.next();
			reg = regions.get(tmpCtg);
			ret.put(tmpCtg, reg.getJunctions(0.90));
			//reg.printState(System.out);
		}
		return ret;
	}

	public static void breakContigs(Map<String, int[]> breaks, String ctgPath, String destPath) throws IOException {
		File brokenFile = new File(destPath);
		brokenFile.createNewFile();
		ScaffoldExporter out = new ScaffoldExporter(brokenFile);
		String tmpCtg;
		BufferedReader br = new BufferedReader(
				new FileReader(new File(ctgPath)));
		StringBuilder sb = null;
		br.read();
		int[] junctions;
		while (br.ready()) {
			tmpCtg = br.readLine();
			// Take only the first part of the contig header to be consistent
			// with bwa
			if (tmpCtg.contains(" "))
				tmpCtg = tmpCtg.substring(0, tmpCtg.indexOf(" "));
			sb = new StringBuilder();
			char c = (char) br.read();
			// read in this sequence
			while (c != '>') {
				if (MisassemblyBreaker.isNuc(c))
					sb.append(c);
				if (!br.ready())
					break;
				c = (char) br.read();
			}
			if (!breaks.containsKey(tmpCtg)) {
				out.export(tmpCtg, sb);
				continue;
			}
			int left = 1;
			int right = 1;
			junctions = breaks.get(tmpCtg);
			for (int i = 0; i < junctions.length; i++) {
				right = junctions[i];
				System.out.println("[a5_qc] Exporting " + tmpCtg + " at "
						+ left + "-" + right);
				out.export(tmpCtg, sb, left, right);
				left = right + 1;
			}
			right = sb.length();
			System.out.println("[a5_qc] Exporting " + tmpCtg + " at " + left+ "-" + right);
			out.export(tmpCtg, sb, left, right);

		}
	}

	/**
	 * Returns a Map of sets of regions. Each set corresponds to a single
	 * contig, indexed by the contig name (a String). <br>
	 * <br>
	 * Each set is stored in a matrix with 4 rows. The first row is the
	 * left-most coordinate of each region, the second row stores the right-most
	 * coordinate of each region, the third stores the current 'weakest'
	 * position, and the fourth stores the score of the current 'weakest'
	 * position. <br>
	 * <br>
	 * The weakest position is the position with the lowest pileup score. For
	 * the definition of this score, see org.gel.halophiles.qc.PileupScorer.
	 */
	public static Map<String, MisassemblyRange> getRegions(String bedPath)
			throws IOException {
		BufferedReader in = new BufferedReader(
				new FileReader(new File(bedPath)));
		String[] dat = null;

		/*
		 * BEGIN: build data structure storing the boundaries of the ranges
		 * we're investigating
		 */
		// Some Map objects to store our ranges indexed, by contig. Vectors will
		// be turned into
		// arrays that can be binary searched easily
		Map<String, Vector<Integer>> left = new HashMap<String, Vector<Integer>>();
		Map<String, Vector<Integer>> right = new HashMap<String, Vector<Integer>>();
		Vector<Integer> vL = null; // temp variable
		Vector<Integer> vR = null; // temp variable
		int l = 0;
		int r = 0;
		int expectedLines = 0;
		while (in.ready()) {
			dat = in.readLine().split("\t");
			// skip any header in the file
			if (dat[0].startsWith("#"))
				continue;
			if (left.containsKey(dat[0])) {
				vL = left.get(dat[0]);
				vR = right.get(dat[0]);
			} else {
				vL = new Vector<Integer>();
				vR = new Vector<Integer>();
				left.put(dat[0], vL);
				right.put(dat[0], vR);
			}
			l = Integer.parseInt(dat[1]);
			r = Integer.parseInt(dat[2]);
			expectedLines += r - l;
			vL.add(l);
			vR.add(r);
		}
		System.out.println("[a5_qc] Expecting to evaluate " + expectedLines + " lines of pileup");
		/*
		 * END: build boundaries
		 * 
		 * START: convert these Vectors into arrays, so they can be queried
		 * using binary search
		 */
		Map<String, MisassemblyRange> lookups = new HashMap<String, MisassemblyRange>();
		Iterator<String> it = left.keySet().iterator();
		String key = null; // temp variable
		int[] arL = null; // temp variable
		int[] arR = null; // temp variable
		while (it.hasNext()) {
			key = it.next();
			arL = toArray(left.get(key));
			arR = toArray(right.get(key));
			lookups.put(key, new MisassemblyRange(key, arL, arR));
		}
		/*
		 * END: convert Vectors to arrays
		 */
		return lookups;

	}

	private static int[] toArray(Vector<Integer> v) {
		int[] ret = new int[v.size()];
		for (int i = 0; i < ret.length; i++)
			ret[i] = v.get(i);
		return ret;
	}

	private static void runMPileup(String bamPath, String bedPath, String ctgPath, Map<String, MisassemblyRange> regions) throws IOException, InterruptedException {
		// System.err.println("Executing command: " + cmd);
		String cmd = SAMTOOLS + " mpileup -d 350 -f " + ctgPath + " -l "+ bedPath + " " + bamPath;
		System.out.println("[a5_qc] Executing: " + cmd);
		Process p = Runtime.getRuntime().exec(cmd);
		LineHandler errLH = new LineHandler() {
			public void handleLine(String line) {
				System.out.println(line);
			}
		};
		PileupEvaluator plpEval = new PileupEvaluator(regions,new PileupScorer());
		StreamReader errRdr = new StreamReader(p.getErrorStream(), errLH);
		StreamReader outRdr = new StreamReader(p.getInputStream(), plpEval);
		errRdr.start();
		outRdr.start();
		int retcode = p.waitFor();
		System.out.println("[a5_qc] mpileup terminated with exit value " + retcode);
		System.out.println("[a5_qc] Evaluated " + plpEval.numLines + " total lines of pileup");
	}

	static class StreamReader extends Thread {
		InputStream in;
		LineHandler lineHandler;
		boolean finished;

		StreamReader(InputStream in, LineHandler lineHandler) {
			this.in = in;
			this.lineHandler = lineHandler;
		}

		public void run() {
			synchronized(this){
				finished=false;
			}
			try {
				InputStreamReader isr = new InputStreamReader(in);
				BufferedReader br = new BufferedReader(isr);
				String line = null;
				while ((line = br.readLine()) != null) {
					lineHandler.handleLine(line);
					//System.err.println("read line");
				}
			} catch (IOException ioe) {
				ioe.printStackTrace();
			} finally {
				synchronized (this) {
					finished=true;
					notifyAll();
				}
			}
		}
		
		public synchronized boolean isFinished(){
			return finished;
		}
	}

	static class PileupEvaluator implements LineHandler {
		Map<String, MisassemblyRange> regions;
		PileupScorer scorer;
		int numLines;

		PileupEvaluator(Map<String, MisassemblyRange> regions, PileupScorer scorer) {
			this.regions = regions;
			this.scorer = scorer;
			this.numLines = 0;
		}

		public void handleLine(String line) {
			String[] plp = line.split("\t");
			MisassemblyRange tmpReg = regions.get(plp[0]);
			int pos = Integer.parseInt(plp[1]);
			double score = scorer.scorePileup(plp[4])/Double.parseDouble(plp[3]);
			tmpReg.addPos(pos, score);
			numLines++;
		}
	}

	static interface LineHandler {
		public void handleLine(String line);
	}
	
	static class MisassemblyRange {
		private String contig;
		private RunningStat[] stats;
		private int[] left, right;
		private int[] minPos;
		MisassemblyRange (String contig, int[] left, int[] right) {
			if (left.length != right.length) throw new IllegalArgumentException("left and right must be the same length");
			this.contig = contig;
			this.left = left;
			this.right = right;
			this.minPos = new int[left.length];
			this.stats = new RunningStat[left.length];
			/*
			 *  [0] 'weakest' position in each range
			 *  [1] score of 'weakest' position
			 *  [2] mean for this range
			 *  [3] stdev for this range
			 *  [4] number of elements for this range
			 */
			for (int i = 0; i < stats.length; i++)
				stats[i] = new RunningStat();
		}
		
		boolean addPos(int pos, double score){
			int idx = Arrays.binarySearch(right, pos);
			if (idx < 0)
				idx = -1 * (idx + 1);
			if (pos < left[idx])
				return false;
			double min = stats[idx].min();
			if (score < min)
				minPos[idx] = pos;
			stats[idx].addVal(score);
			return true;
		}
		
		boolean contains(int pos){
			int idx = Arrays.binarySearch(right, pos);
			if (idx < 0)
				idx = -1 * (idx + 1);
			if (pos < left[idx])
				return false;
			else
				return true;
		}
		
		void printState(PrintStream out) {
			for (int i = 0; i < left.length; i++){
				out.println(left[i]+"\t"+right[i]+"\t"+stats[i].toString());
			}
		}
		
		int[] getJunctions(double maxScore){
			Vector<Integer> scores = new Vector<Integer>();
			for (int i = 0; i < minPos.length; i++){
				if (stats[i].min() < maxScore)
					scores.add(minPos[i]);
			}
			return toArray(scores);
		}
	}
	
	static class RunningStat {
		private int n;
		private double oldM, newM, oldS, newS, min, max;
		public RunningStat() {
			n = 0;
			min = Double.POSITIVE_INFINITY;
			max = Double.NEGATIVE_INFINITY;
		}
		void addVal(double x){
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
			min = x < min ? x : min;
			max = x > max ? x : max;
		}
		int numPts(){
			return n;
		}
		double mean() {
			return (n > 0) ? newM : 0.0;
		}
		double variance() {
			return (n > 1) ? newS/(n-1) : 0.0;
		}
		double sd() {
			return Math.sqrt(variance());
		}
		double min(){
			return min;
		}
		double max(){
			return max;
		}
		
		public String toString() {
			return mean()+"\t"+sd()+"\t"+min+"\t"+max+"\t"+n;
		}
		
	}
	
}