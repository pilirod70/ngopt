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

import org.halophiles.assembly.Contig;
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
	/**
	 * This function is used if we already have pileup file.
	 * 
	 * @param bedPath 
	 * @param plpPath
	 * @return
	 * @throws IOException
	 */
	public static Map<String,Vector<MisassemblyRange>> scoreAtBaseLevel(File bedFile, File plpFile, File connectionsFile, Map<String,Contig> contigs) throws IOException {
		Map<String, Vector<MisassemblyRange>> regions = getRegions(bedFile, connectionsFile, contigs);
		PileupEvaluator plpEval = new PileupEvaluator(regions, new PileupScorer());
		BufferedReader br = new BufferedReader(new FileReader(plpFile));
		while(br.ready()){
			plpEval.handleLine(br.readLine());
		}
		return regions;
	}
	
	public static Map<String, Vector<MisassemblyRange>> scoreAtBaseLevel(File bamFile, File bedFile, File ctgFile, File connectionsFile, Map<String,Contig> contigs) throws IOException {
		Map<String, Vector<MisassemblyRange>> regions = getRegions(bedFile, connectionsFile, contigs);
		try {
			runMPileup(bamFile.getAbsolutePath(), bedFile.getAbsolutePath(), ctgFile.getAbsolutePath(), regions);
		} catch (InterruptedException e) {
			System.err.println("Unable to run samtools mpileup. This is what I know:");
			e.printStackTrace();
			System.exit(-1);
		}
		return regions;
	}
	
	public static void scoreAtBaseLevel(String bamPath, String bedPath, String ctgPath, Map<String, Vector<MisassemblyRange>> regions) throws IOException {
		try {
			runMPileup(bamPath, bedPath, ctgPath, regions);
		} catch (InterruptedException e) {
			System.err.println("Unable to run samtools mpileup. This is what I know:");
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	
	public static Map<String, int[]> refine(Map<String, Vector<MisassemblyRange>> regions) throws IOException{
		Map<String, int[]> ret = new HashMap<String, int[]>();
		Iterator<String> ctgIt = regions.keySet().iterator();
		String tmpCtg;
		Vector<MisassemblyRange> reg;
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
	 * Parses MisassemblyRanges stored in a bed file and a connections file
	 * 
	 * @param bedFile
	 * @param connectionsFile
	 * @param contigs
	 * @return
	 * @throws IOException
	 */
	public static Map<String, Vector<MisassemblyRange>> getRegions(File bedFile, File connectionsFile, Map<String,Contig> contigs) throws IOException {
		
		// contig+","+(rev?"-":"+")+","+Integer.toString(left)+","+Integer.toString(right);

		Map<String, Vector<MisassemblyRange>> ret = new HashMap<String, Vector<MisassemblyRange>>();
		
		BufferedReader in = new BufferedReader(new FileReader(connectionsFile));
		
		// Stores the misassembly blocks we've seen so far.
		Map<Integer, MisassemblyBlock> blockSet = new HashMap<Integer,MisassemblyBlock>();
		Map<String,MisassemblyBlock[]> blocksByRegion = new HashMap<String, MisassemblyBlock[]>();
		MisassemblyBlock[] tmpBlks = null;

		String[] line = null;
		String reg = null;
		while(in.ready()){
			line = in.readLine().split("\t");
			reg = line[0];
			tmpBlks = new MisassemblyBlock[2];
			// parse the block on the left side of this region
			tmpBlks[0] = addFlankingBlocks(line[1], blockSet, contigs);

			// parse the block on the right side of this region
			tmpBlks[1] = addFlankingBlocks(line[2], blockSet, contigs);
			
			blocksByRegion.put(reg, tmpBlks);
		}
		
		in = new BufferedReader(new FileReader(bedFile));


		Vector<MisassemblyRange> ranges = null;
		
		while (in.ready()){
			line = in.readLine().split("\t");
			// skip any header in the file
			if (line[0].startsWith("#"))
				continue;
			if (ret.containsKey(line[0]))
				ranges = ret.get(line[0]);
			else
				ranges = new Vector<MisassemblyRange>();
			tmpBlks = blocksByRegion.get(line[3]);
			ranges.add(new MisassemblyRange(contigs.get(line[0]), tmpBlks[0], tmpBlks[1]));
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

	private static void runMPileup(String bamPath, String bedPath, String ctgPath, Map<String, Vector<MisassemblyRange>> regions) throws IOException, InterruptedException {
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
		Map<String, Vector<MisassemblyRange>> regions;
		PileupScorer scorer;
		int numLines;

		PileupEvaluator(Map<String, Vector<MisassemblyRange>> regions, PileupScorer scorer) {
			this.regions = regions;
			this.scorer = scorer;
			this.numLines = 0;
		}

		public void handleLine(String line) {
			String[] plp = line.split("\t");
			Vector<MisassemblyRange> tmpRegs = regions.get(plp[0]);
			int pos = Integer.parseInt(plp[1]);
			double score = scorer.scorePileup(plp[4])/Double.parseDouble(plp[3]);
			int regIdx = MisassemblyRange.binarySearch(tmpRegs, pos);
			if (regIdx >= 0)
				tmpRegs.get(regIdx).addPos(pos, score);
			numLines++;
		}
	}

	static interface LineHandler {
		public void handleLine(String line);
	}
	
	
	
}