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
import java.util.Set;
import java.util.TreeSet;
import java.util.Vector;

import net.sf.samtools.SAMFileReader;

import org.halophiles.assembly.Contig;
import org.halophiles.tools.HelperFunctions;

public class MBRefiner {

	private static final String SAMTOOLS = "/jgi/tools/bin/samtools";

	public static void main(String[] args) {
		HelperFunctions.logInputs("MBRefiner", args);
		if (args.length == 4) {
			try {
				File bamFile = new File(args[0]);
				File bedFile = new File(args[1]);
				File connectionsFile = new File(args[2]);
				File ctgFile = new File(args[2]);
				File brokenFile = new File(args[3]);
				Map<String, Contig> contigs = MisassemblyBreaker.getContigs(new SAMFileReader(bamFile));
				Map<String, Vector<MisassemblyRegion>> ranges = getRegions(bedFile, connectionsFile, contigs);
				scoreAtBaseLevel(bamFile, bedFile, ctgFile, ranges, contigs);
				Map<String, int[]> junctions = refine(ranges);
				breakContigs(junctions, ctgFile.getAbsolutePath(), brokenFile.getAbsolutePath());

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
	/**
	 * This function is used if we already have pileup file.
	 * 
	 * @param bedPath 
	 * @param plpPath
	 * @return
	 * @throws IOException
	 */
	public static Map<String,Vector<MisassemblyRegion>> scoreAtBaseLevel(File plpFile, Map<String, Vector<MisassemblyRegion>> regions, Map<String,Contig> contigs) throws IOException {
		PileupEvaluator plpEval = new PileupEvaluator(regions, new PileupScorer());
		BufferedReader br = new BufferedReader(new FileReader(plpFile));
		while(br.ready()){
			plpEval.handleLine(br.readLine());
		}
		return regions;
	}
	
	public static void scoreAtBaseLevel(File bamFile, File bedFile, File ctgFile, Map<String, Vector<MisassemblyRegion>> regions,Map<String,Contig> contigs) throws IOException {
		try {
			runMPileup(bamFile.getAbsolutePath(), bedFile.getAbsolutePath(), ctgFile.getAbsolutePath(), regions);
		} catch (InterruptedException e) {
			System.err.println("Unable to run samtools mpileup. This is what I know:");
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public static Map<String, int[]> refine(Map<String, Vector<MisassemblyRegion>> regions) throws IOException{
		Map<String, int[]> ret = new HashMap<String, int[]>();
		Iterator<String> ctgIt = regions.keySet().iterator();
		String tmpCtg;
		MisassemblyRegion range= null;
		while (ctgIt.hasNext()) {
			tmpCtg = ctgIt.next();
			Vector<Integer> breaks = new Vector<Integer>();
			Iterator<MisassemblyRegion> rangeIt = regions.get(tmpCtg).iterator();
			while (rangeIt.hasNext()){
				range = rangeIt.next();
				if (range.getMinScore() < 0.90)
					breaks.add(range.getMinPos());
			}
			ret.put(tmpCtg, HelperFunctions.toArray(breaks));
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
			// Take only the first part of the contig header to be consistent with bwa
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
	 * Parses MisassemblyRanges stored in a bed file and a connections file
	 * 
	 * @param bedFile
	 * @param connectionsFile
	 * @param contigs
	 * @return
	 * @throws IOException
	 */
	public static Map<String, Vector<MisassemblyRegion>> getRegions(File bedFile, File connectionsFile, Map<String,Contig> contigs) throws IOException {
		
		// contig+","+(rev?"-":"+")+","+Integer.toString(left)+","+Integer.toString(right);

		Map<String, Vector<MisassemblyRegion>> ret = new HashMap<String, Vector<MisassemblyRegion>>();
		
		BufferedReader in = new BufferedReader(new FileReader(connectionsFile));
		
		// Stores the misassembly blocks we've seen so far.
		Map<Integer, MisassemblyBlock> blockSet = new HashMap<Integer,MisassemblyBlock>();
		Map<String,MisassemblyBlock[]> blocksByRegion = new HashMap<String, MisassemblyBlock[]>();
		MisassemblyBlock[] tmpBlks = null;

		String[] entry = null;
		String reg = null;
		String line = null;
		while(in.ready()){
			line = in.readLine();
			entry = line.split("\t");
			reg = entry[0];
			tmpBlks = new MisassemblyBlock[2];
			// parse the block on the left side of this region
			tmpBlks[0] = addFlankingBlocks(entry[1], blockSet, contigs);

			// parse the block on the right side of this region
			tmpBlks[1] = addFlankingBlocks(entry[2], blockSet, contigs);
			
			blocksByRegion.put(reg, tmpBlks);
		}
		
		in = new BufferedReader(new FileReader(bedFile));


		Vector<MisassemblyRegion> ranges = null;
		MisassemblyRegion tmpReg = null;
		while (in.ready()){
			entry = in.readLine().split("\t");
			// skip any header in the file
			if (entry[0].startsWith("#"))
				continue;
			tmpBlks = blocksByRegion.get(entry[3]);
			tmpReg = new MisassemblyRegion(contigs.get(entry[0]), tmpBlks[0], tmpBlks[1]);
			if (ret.containsKey(entry[0])) 
				ret.get(entry[0]).add(tmpReg);
			else {
				ranges = new Vector<MisassemblyRegion>();
				ranges.add(tmpReg);
				ret.put(entry[0], ranges);
			}
		}
		return ret;
	}
	/**
	 * Parses a MisassemblyBlock string, and returns the MisasssemblyBlock object represented by it.
	 * @param blockString a String read in from the a connections file
	 * @param blockSet a map to store the MisassemblyBlocks in
	 * @param contigs the contigs on which these MisassemblyBlocks exist
	 * @return the object represented by <code>blockString</code>
	 */
	private static MisassemblyBlock addFlankingBlocks(String blockString, Map<Integer,MisassemblyBlock> blockSet, Map<String,Contig> contigs){
		String[] dat = blockString.split(",");
		int id = -1;
		String ctg = null;
		boolean rev = false;
		int left = -1;
		int right = -2;

		MisassemblyBlock tmpFlank = null;
		MisassemblyBlock tmpConnect = null;
		id = Integer.parseInt(dat[0]);
		if (blockSet.containsKey(id)) {
			tmpFlank = blockSet.get(id);
		} else {
			ctg = dat[1];
			rev = dat[2].equals("-");
			left = Integer.parseInt(dat[3]);
			right = Integer.parseInt(dat[4]);
			tmpFlank = new MisassemblyBlock(contigs.get(ctg), left, right, rev, id);
			blockSet.put(id, tmpFlank);
			
			// now parse the block it is connected to
			id = Integer.parseInt(dat[5]);
			ctg = dat[6];
			rev = dat[7].equals("-");
			left = Integer.parseInt(dat[8]);
			right = Integer.parseInt(dat[9]);
			tmpConnect = new MisassemblyBlock(contigs.get(ctg), left, right, rev, id);
			blockSet.put(id, tmpConnect);
			
			// now set the connections
			tmpFlank.addConnection(tmpConnect);
			tmpConnect.addConnection(tmpFlank);
		}
		return tmpFlank;
	}

	private static void runMPileup(String bamPath, String bedPath, String ctgPath, Map<String, Vector<MisassemblyRegion>> regions) throws IOException, InterruptedException {
		// System.err.println("Executing command: " + cmd);
		String cmd = SAMTOOLS + " mpileup -d 350 -f " + ctgPath + " -l "+ bedPath + " " + bamPath;
		System.out.println("[a5_qc] Executing: " + cmd);
		LineHandler errLH = new LineHandler() {
			public void handleLine(String line) {
				System.out.println(line);
			}
		};
		PileupEvaluator plpEval = new PileupEvaluator(regions,new PileupScorer());
		Process p = Runtime.getRuntime().exec(cmd);
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
			try {
				InputStreamReader isr = new InputStreamReader(in);
				BufferedReader br = new BufferedReader(isr);
				String line = null;
				if (lineHandler instanceof PileupEvaluator)
					System.out.print("");
				while ((line = br.readLine()) != null) {
					lineHandler.handleLine(line);
					//System.err.println("read line");
				}
			} catch (IOException ioe) {
				ioe.printStackTrace();
			}
		}
	}

	static class PileupEvaluator implements LineHandler {
		Map<String, Vector<MisassemblyRegion>> regions;
		PileupScorer scorer;
		int numLines;

		PileupEvaluator(Map<String, Vector<MisassemblyRegion>> regions, PileupScorer scorer) {
			this.regions = regions;
			this.scorer = scorer;
			this.numLines = 0;
		}

		public void handleLine(String line) {
			String[] plp = line.split("\t");
			Vector<MisassemblyRegion> tmpRegs = regions.get(plp[0]);
			if (tmpRegs == null)
				System.out.print("");
			int pos = Integer.parseInt(plp[1]);
			double score = scorer.scorePileup(plp[4])/Double.parseDouble(plp[3]);
			int regIdx = MisassemblyRegion.binarySearch(tmpRegs, pos);
			if (pos > 1099967 && pos < 1100130)
				System.out.print("");
			MisassemblyRegion reg = null;
			if (regIdx >= 0){
				reg = tmpRegs.get(regIdx);
				reg.addPos(pos, score);
			}
			numLines++;
		}
	}

	static interface LineHandler {
		public void handleLine(String line);
	}
	
	
	
}