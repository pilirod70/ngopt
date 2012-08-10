/**

 * This file is part of the A5 pipeline.
 * (c) 2011, 2012 Andrew Tritt and Aaron Darling
 * This software is licensed under the GPL, v3.0. Please see the file LICENSE for details
 */
package org.halophiles.assembly.qc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.RandomAccessFile;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.Vector;
/*
import net.sf.samtools.AbstractBAMFileIndex;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceRecord;
*/
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMSequenceRecord;

import org.halophiles.assembly.Contig;
import org.halophiles.assembly.ReadPair;
import org.halophiles.assembly.ReadSet;
import org.halophiles.tools.HelperFunctions;

public class MisassemblyBreaker {
	
	/**
	 * A Comparator for sorting blocks by right-most position
	 */
	private static Comparator<MisassemblyBlock> BLOCK_COMP = new Comparator<MisassemblyBlock>(){
		public int compare(MisassemblyBlock arg0, MisassemblyBlock arg1) {
			return arg0.getRight() - arg1.getRight();
		}
	};
	
	/**
	 * Used for randomly sampling the SAM file. Read in
	 * HANDFUL_SIZE reads after each jump in the SAM file.
	 */
	private static final int HANDFUL_SIZE = 100;
	
	/**
	 * The approximate average length of each SAM file entry
	 */
	private static final int SAM_LINE_LEN = 215;
	
	/**
	 * The tag used in the SAM file for indicating a non-ambiguous
	 * read mapping location. Use this to ignore reads that were
	 * ambiguously placed.
	 */
	private static final String TAG_KEEP = "XT:A:U";
	
	/**
	 * The Minimum number of points allowed in a Read cluster
	 */
	private static final int MIN_PTS = 3;
	
	/**
	 * The maximum number of ReadPairs to sample for insert size calculations
	 */
	private static int N_ESTREADPAIRS = 100000;
	
	private static NumberFormat NF;
	
	/**
	 * Represents whether there are a significant amount of
	 * innie read pairs
	 */
	private static boolean INWARD = false;
	
	/**
	 * Represents whether there are a significant amount of
	 * outie read pairs
	 */
	private static boolean OUTWARD = false;
	
	/**
	 * Number of inward facing read pairs (innies)
	 */
	private static double NIN = 0;
	
	/**
	 * Number of outward facing read pairs (outies)
	 */
	private static double NOUT = 0;
	
	/** error rate for estimation of maximum values */
	private static double ALPHA = 0.001;
	
	/**
	 * The probability of success for the Geometric distrbution
	 * used for modelling read mapping locations. This
	 * is set to be the minimum non-zero mapping frequency
	 * over windows across the genome.
	 */
	private static double P;
	
	/**
	 * The maximum allowed distance between two mapping points
	 */
	public static int MAX_INTERPOINT_DIST;
	
	/**
	 * The maximum allowed distance between two blocks
	 * for calling a region with a misassembly
	 */
	private static int MAX_INTERBLOCK_DIST;
	
	/**
	 * The expected length of blocks of over-represented
	 * mapping locations.
	 */
	static int MEAN_BLOCK_LEN;
	
	/**
	 * The minimum allowed length of a block for calling
	 * significant over-representation of mapped reads.
	 */
	static int MIN_BLOCK_LEN;
	
	/**
	 * The maximum allowed length of a block for calling
	 * significant over-representation of mapped reads.
	 */
	static int MAX_BLOCK_LEN;
	
	/**
	 * A Collection of data structures for each pair of contig with
	 * paired-read connections for storing each paired-read location
	 * and running the DBSCAN spatial clustering algorithm.
	 */
	private static Collection<SpatialClusterer> matches;
	
	public static void main(String[] args){
		if (args.length != 4 && args.length != 3){
			System.err.println("Usage: java -jar A5qc.jar <sam_file> <contig_file> <output_file>");
			System.exit(-1);
		}
		try{
			HelperFunctions.logInputs("A5qc", args);
			NF = NumberFormat.getInstance();
			NF.setMaximumFractionDigits(0);
			NF.setGroupingUsed(false);
			String samPath = args[0];
			String ctgPath = args[1];
			String brokenCtgPath = args[2];
			
			File ctgFile = new File(ctgPath);
			
			String base = HelperFunctions.dirname(samPath)+"/"+HelperFunctions.basename(samPath,".sam");
			File samFile = new File(samPath);
			File bamFile = new File(base+".bam");
			File bedFile = new File(base+".regions.bed");
			File connectionsFile = new File(base+".connections");
			Map<String,Vector<MisassemblyRegion>> regions = null;
			
			RandomAccessFile raf = new RandomAccessFile(samFile, "r");
			Map<String,Contig> contigs = readContigs(raf);
			
			
			if (!bedFile.exists() || (bedFile.exists() && HelperFunctions.isEmpty(bedFile))) {
				regions = findMisasmRegions(samFile, bedFile, connectionsFile, contigs, getInsertStats(raf, contigs));
			} else {
				regions = MBRefiner.getRegions(bedFile, connectionsFile, contigs);
				System.out.println("[a5_qc] Found a bed file. Assuming I generated this in a past run and proceeding");
			}
			if (regions != null) { // if we created a bedFile, we have regions containing misassemblies
				MBRefiner.scoreAtBaseLevel(bamFile, bedFile, ctgFile, regions, contigs);
				MBRefiner.breakContigs(MBRefiner.refine(regions), ctgPath, brokenCtgPath);
			}
			
		} catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		} catch(Exception e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	/**
	 * Randomly samples reads from the given RandomAccessFile, and computes insert size stats
	 * @param raf the file to sample reads from for computing insert size stats
	 * @param contigs the contigs that the reads stored in <code>raf</code> are mapped to 
	 * @return array of arrays
	 * @throws IOException
	 */
	private static double[][] getInsertStats(RandomAccessFile raf, Map<String, Contig> contigs) throws IOException{
		
		if (contigs.size() == 0){
			System.out.println("[a5_qc] Could not find any contigs in SAM file. File is either not in SAM format, or is missing header.");
			System.exit(-1);
		} else {
			System.out.println("[a5_qc] Found "+contigs.size()+" contigs");
		}
		System.out.println("[a5_qc] Reading in a subset of reads for insert size estimation.");
		long before = System.currentTimeMillis();
		Map<String,ReadPair> reads = readSubsetByChunk(raf, contigs);
		raf.close();
		
		long after = System.currentTimeMillis();
		System.out.println("[a5_qc] Took "+HelperFunctions.millisToTimeString(after-before)+" to read in "+reads.size()+" read pairs.");
		if (reads.size() <= 0) {
			System.err.println("[a5_qc] No paired reads found -- Cannot detect misassemblies -- Exiting.");
			System.exit(-1);
		}
		setOrientation(reads.values());
		if (INWARD && !OUTWARD)
			System.out.println("[a5_qc] Found a substantial amount of innies, but found no outties.");
		else if (!INWARD && OUTWARD)
			System.out.println("[a5_qc] Found a substantial amount of outties, but found no innies.");
		else 
			System.out.println("[a5_qc] Found both innies and outties.");

		return getLibraryStats(reads, 1);
	}
	
	/**
	 * Identifies regions containing misassemblies using the provided SAM file, and writes
	 * these regions to a BED file under the given path. If no regions are found, nothing 
	 * is written to <code>bedPath</code>.
	 * 
	 * @returns true if any regions containing misassemblies are found, false otherwise
	 */
	private static Map<String,Vector<MisassemblyRegion>> findMisasmRegions(File samFile, File bedFile, File connectionsFile, Map<String,Contig> contigs, double[][] insertStats) throws IOException{

		System.out.println("[a5_qc] Reading "+samFile.getPath());
		
		double[][] insRanges = getFilterRanges(insertStats);
		//setMAXBLOCKLEN(clusterStats);
		//loadBAMData(samPath, contigs, ranges);
		loadData(samFile, contigs, insRanges);
		
		String base = HelperFunctions.dirname(samFile.getAbsolutePath())+"/"+HelperFunctions.basename(samFile.getAbsolutePath(), ".bam");
		setParameters(insertStats);
		
		
		printParams();
		
		// collect all of our blocks for each contig
		Iterator<SpatialClusterer> mbIt = matches.iterator();
		Map<String,Vector<MisassemblyBlock>> blocks = new HashMap<String,Vector<MisassemblyBlock>>();
		Vector<MisassemblyBlock> xBlocks = null;
		Vector<MisassemblyBlock> yBlocks = null;
		while(mbIt.hasNext()){
			SpatialClusterer pc = mbIt.next();
			String statePath = base+".matches."+HelperFunctions.basename(pc.getContig1().name)+".v."+HelperFunctions.basename(pc.getContig2().name)+".txt";
			pc.exportCurrState(statePath);
			// get Vector for holding contig X blocks
			xBlocks = new Vector<MisassemblyBlock>();
			// get Vector for holding contig Y blocks
			yBlocks = new Vector<MisassemblyBlock>();
			pc.buildReadPairClusters();
			addBlocks(pc, xBlocks, yBlocks);
			//pc.exportCurrState(new File(matchDir,"match."+pc.getContig1().getId()+"v"+pc.getContig2().getId()+".txt"));
	//		removeTerminalBlocks(pc.getContig1(), xBlocks);
	//		removeTerminalBlocks(pc.getContig2(), yBlocks);
			if (blocks.containsKey(pc.getContig1().name))
				blocks.get(pc.getContig1().name).addAll(xBlocks);
			else
				blocks.put(pc.getContig1().name, xBlocks);
			
			if (blocks.containsKey(pc.getContig2().name))
				blocks.get(pc.getContig2().name).addAll(yBlocks);
			else
				blocks.put(pc.getContig2().name, yBlocks);
			
		}
		Iterator<String> it = blocks.keySet().iterator();
		/*
		while(it.hasNext())
			removeRepeats(blocks.get(it.next()));
			
		it = blocks.keySet().iterator();
		*/
		Vector<String> toRm = new Vector<String>();
		while(it.hasNext()){
			String tmpCtg = it.next();
			Vector<MisassemblyBlock> tmpBlocks = blocks.get(tmpCtg);
			System.out.println("[a5_qc] Found "+tmpBlocks.size()+" blocks on contig " + contigs.get(tmpCtg).getId());
			if (tmpBlocks.isEmpty())
				toRm.add(tmpCtg);
			else 
				Collections.sort(tmpBlocks, BLOCK_COMP);
			Iterator<MisassemblyBlock> blockIt = tmpBlocks.iterator();
			while(blockIt.hasNext()){
				MisassemblyBlock tmp = blockIt.next();
				System.out.println("        "+tmp.toString());
				
			}
		}
		// clear out contigs from the blocks Map if they don't have any blocks.
		if (toRm.size() > 0)
			removeKeys(blocks, toRm);
		
		/*if (blocks.isEmpty()){
			return false;
		}*/
		PrintStream bedOut = null;
		PrintStream connectionsOut = null;
		Map<String,Vector<MisassemblyRegion>> ret = null;
		MisassemblyRegion tmpRange;
		if (!blocks.isEmpty()) {
			Iterator<String> ctgIt = blocks.keySet().iterator();
			String tmpCtg;
			int left = -1;
			int right = -1;
			while(ctgIt.hasNext()){
				tmpCtg = ctgIt.next();
				Vector<MisassemblyBlock> tmpBlks = blocks.get(tmpCtg);
	
				if (tmpBlks.size() < 2)
					continue;
			
				// Collections.sort(tmpBlks,BLOCK_COMP); // This was already sorted abov 
				Vector<MisassemblyRegion> ranges = null;
				for (int i = 1; i < tmpBlks.size(); i++){
					// if they don't face each other, there isn't a misassembly in between this pair of blocks
					if (tmpBlks.get(i-1).getRev() || !tmpBlks.get(i).getRev()) 
						continue;
					left = -1;
					right = -1;
					if (tmpBlks.get(i-1).getRight() > tmpBlks.get(i).getLeft()) {
						left = tmpBlks.get(i).getLeft();
						right = tmpBlks.get(i-1).getRight();
					} else if (tmpBlks.get(i).getLeft()-tmpBlks.get(i-1).getRight() < MAX_INTERBLOCK_DIST) {
						left = tmpBlks.get(i-1).getRight();
						right = tmpBlks.get(i).getLeft();
					}
					
					if (left != -1 && right != -1){
						// if they overlap, split at the midpoint of overlap
						if (bedOut == null) { // don't create the file if we don't need to
							bedFile.createNewFile();
							connectionsFile.createNewFile();
							bedOut = new PrintStream(bedFile);
							connectionsOut = new PrintStream(connectionsFile);
							System.out.println("[a5_qc] Writing regions containing misassemblies to " + bedFile.getAbsolutePath());
							ret = new HashMap<String, Vector<MisassemblyRegion>>();
						}
						if (ranges == null)
							ranges = new Vector<MisassemblyRegion>();
						tmpRange = new MisassemblyRegion(contigs.get(tmpCtg), tmpBlks.get(i-1), tmpBlks.get(i));
						bedOut.println(tmpRange.toString());
						logConnection(connectionsOut, tmpRange.getId(), tmpBlks.get(i-1), tmpBlks.get(i));
						ranges.add(tmpRange);
					}
				}
				if (ranges != null) {
					if (ret == null)
						ret = new HashMap<String, Vector<MisassemblyRegion>>();
					ret.put(tmpCtg, ranges);
				}
			}
		}
		
		if (bedOut != null) {
			bedOut.close();
			connectionsOut.close();
		}
		return ret;
	}
	
	private static void logConnection(PrintStream out, String range, MisassemblyBlock left, MisassemblyBlock right){
		StringBuilder sb = new StringBuilder();
		sb.append(range+"\t");
		sb.append(left.toString()+",");
		sb.append(left.getConnection().toString());
		sb.append("\t");
		sb.append(right.toString()+",");
		sb.append(right.getConnection().toString());
		out.println(sb.toString());
	}
	
	private static boolean resetParam(String in){
		if (!in.contains("="))
			return false;
		String[] pair = in.split("=");
		String key = pair[0];
		String value = pair[1];
		if (key.equals("ALPHA")){
			ALPHA = Double.parseDouble(value);
			return true;
		} else if (key.equals("P")){
			P = Double.parseDouble(value);
			return true;
		} else if (key.equals("MIN_BLOCK_LEN")) {
			MIN_BLOCK_LEN = Integer.parseInt(value);
			return false;
		} else if (key.equals("MEAN_BLOCK_LEN")) {	
			MEAN_BLOCK_LEN = Integer.parseInt(value);
			return false;
		} else if (key.equals("MAX_BLOCK_LEN")) {
			MAX_BLOCK_LEN = Integer.parseInt(value);
			return false;
		} else if (key.equals("MAX_INTERBLOCK_DIST")) {
			MAX_INTERBLOCK_DIST = Integer.parseInt(value);
			return false;
		} else if (key.equals("MAX_INTERPOINT_DIST")) {
			MAX_INTERPOINT_DIST = Integer.parseInt(value);
			return false;
		} else if (key.startsWith("EPS")) {
			SpatialClusterer.EPS = Double.parseDouble(value);
			return false;
		} else if (key.equals("MIN_POINTS")) {
			SpatialClusterer.MIN_PTS = Integer.parseInt(value);
			return false;
		} else if (key.equals("RDLEN")) {
			ReadCluster.RDLEN = Integer.parseInt(value);
			return false;
		} else {
			System.out.println("Unrecognized key: " + key);
			return false;
		}
	}
	
	/**
	 * Identify and add significant blocks resulting from spatial clustering of read pairs.
	 * @param sc the SpatialClusterer to get blocks from
	 * @param xBlocks the Collection of blocks for Contig 1 (aka the x Contig)
	 * @param yBlocks the Collection of blocks for Contig 2 (aka the y Contig)
	 */
	private static void addBlocks(SpatialClusterer sc, Vector<MisassemblyBlock> xBlocks, Vector<MisassemblyBlock> yBlocks){
		ReadCluster[] rdClust = sc.getReadClusters();
		System.out.println("[a5_qc] Found "+rdClust.length+" initial blocks between contigs "+sc.getContig1().getId()+" and "+sc.getContig2().getId());
		for (ReadCluster r: rdClust)
			System.out.println("        "+r.toString());
		int xlen = 0;
		int ylen = 0;
		MisassemblyBlock x = null;
		MisassemblyBlock y = null;
		for (int i = 0; i < rdClust.length; i++){
			xlen = rdClust[i].xMax-rdClust[i].xMin;
			ylen = rdClust[i].yMax-rdClust[i].yMin;
			/*
			y = new int[3];
			if (rdClust[i].yMin < 0){
				y[0] = Math.abs(rdClust[i].yMax);
				y[1] = Math.abs(rdClust[i].yMin);
				System.err.println ("y values are negative in MisassemblyBreaker.addBlocks");
			} else {
				y[0] = rdClust[i].yMin;
				y[1] = rdClust[i].yMax;
			}
			y[2] = rdClust[i].yOri ? 1 : -1; 
			 */
			if (xlen >= MIN_BLOCK_LEN && xlen <= MAX_BLOCK_LEN && ylen >= MIN_BLOCK_LEN && ylen <= MAX_BLOCK_LEN) {
				x = new MisassemblyBlock(sc.getContig1(),rdClust[i].xMin, rdClust[i].xMax, rdClust[i].xRev);
				y = new MisassemblyBlock(sc.getContig2(),rdClust[i].yMin, rdClust[i].yMax, rdClust[i].yRev);
				x.addConnection(y);
				y.addConnection(x);
				xBlocks.add(x);
				yBlocks.add(y);
			} 
		}
	}
	/**
	 * Remove blocks at the end of the given Contig
	 * @param contig the Contig the given blocks correspond to.
	 * @param blocks the Collection of blocks to filter
	 */
	private static void removeTerminalBlocks(Contig contig, Vector<int[]> blocks){
		int i = 0;
		while (i < blocks.size()){
			if (blocks.get(i)[1] < MAX_BLOCK_LEN || 
					contig.len - blocks.get(i)[0] < MAX_BLOCK_LEN)
				blocks.remove(i);
			else
				i++;
		}
	}
	
	/**
	 * Remove segments that overlap by more than 50%
	 * @param blocks
	 */
	private static void removeRepeats(Vector<MisassemblyBlock> blocks){
		Collections.sort(blocks, BLOCK_COMP);
		MisassemblyBlock block1 = null;
		MisassemblyBlock block2 = null;
		int i = 0;
		while (i < blocks.size()-1 && blocks.size() > 1){
			block1 = blocks.get(i);
			block2 = blocks.get(i+1);
			// if they overlap, check to see if they overlap enough to constitute a repeat 
			if (block1.getLeft() < block2.getRight() && block2.getLeft() < block1.getRight()){
				double intersection = block1.getRight() - block2.getLeft();
				double union = block2.getRight() - block1.getLeft();
				if (intersection/union > 0.5){
					blocks.remove(i);
					blocks.remove(i);
				} else
					i++;
			} else
				i++;
		}
		System.out.print("");
	}
	/**
	 * Return true if the given character represents a non-ambiguity code.
	 * @param c
	 * @return
	 */
	public static boolean isNuc(char c){
		switch (c) {
			case 'a': return true;
			case 'c': return true;
			case 'g': return true;
			case 't': return true;
			case 'n': return true;
			case 'A': return true;
			case 'C': return true;
			case 'G': return true;
			case 'T': return true;
			case 'N': return true;
			default : return false;
		}
	}
	
	
	/**
	 * Read data from a SAM file. Filters out pairs with inserts with in the given ranges
	 * 
	 * @param samPath
	 * @param fishDir
	 * @param base
	 * @param ctgs
	 * @param ranges an array of arrays of max and min insert sizes
	 * @throws IOException
	 */
	public static void loadData(File samFile, Map<String,Contig> ctgs, double[][] ranges) throws IOException{
		for (int i = 0; i < ranges.length; i++)
			System.out.println("[a5_qc_load_data] Filtering read pairs with inserts between "+
					NF.format(ranges[i][0])+"-"+NF.format(ranges[i][1]));
		Map<String,SpatialClusterer> clusterers = new HashMap<String,SpatialClusterer>();
		/* Keep track of all clusterers that each contig is associated with.*/
		Map<String,Vector<String>> ctgClusterers = new HashMap<String,Vector<String>>();
		Map<String,Integer> counts = new HashMap<String,Integer>();
		Vector<String> tmpMBs = null;

		FileInputStream fis = new FileInputStream(samFile);
		long start = fis.getChannel().position();
		long len = fis.getChannel().size() - start;
		BufferedReader br = new BufferedReader (new InputStreamReader(fis));
		
		/* begin: build a lookup table for each contig for tallying coverage in windows */
		
		int contigLen = 0;
		// require the window length to be at least 1000 base pairs
		int windowLen = 1000;
		String[] hdr = null;
		String contigName = null;
		Map<String,int[][]> readCounts = new HashMap<String,int[][]>();
		int[][] tmpWin = null;
		while(HelperFunctions.nextCharIs(br, '@')){
			hdr = br.readLine().split("\t");
			if (hdr[0].equals("@SQ")){
				for (int i = 1; i < hdr.length; i++){
					if (hdr[i].startsWith("LN:")){
						contigLen = Integer.parseInt(hdr[i].substring(hdr[i].indexOf("LN:")+3));
						// calculate the number of windows we need for this contig
						int numWindow = contigLen/windowLen;
						if (contigLen % windowLen != 0)
							numWindow++;
						// now create the lookup table: add the upper bounds of each window to the array
						tmpWin = new int[4][numWindow];
						tmpWin[0][numWindow-1] = contigLen;
						tmpWin[1][numWindow-1] = (contigLen - (numWindow-1)*windowLen);
						for (int j = 0; j < numWindow-1; j++){
							tmpWin[0][j] = (j+1)*windowLen;
							tmpWin[1][j] = windowLen;
						}
					} else if (hdr[i].startsWith("SN:")){
						contigName = hdr[i].substring(hdr[i].indexOf("SN:")+3);
					}
				}
				// store our lookup table for this contig
				readCounts.put(contigName,tmpWin);
			}
		}
		/* end: build a lookup table for tallying coverage in windows */
		
		
		String[] line1 = null;
		String[] line2 = null;
		int left1 = 0;
		int left2 = 0;
		String ctgStr = null;
		String tmp = null;
		SpatialClusterer pc = null;
		Contig ctg1 = null;
		Contig ctg2 = null;
		int ctgNameComp = -10;
		System.out.print("[a5_qc_load_data] Reading SAM file...");
		long currPos = start;
		double perc = 0;
		double ten = 1;
		int numKeep = 0;
		int total = 0;
		int index = 0;
		long before = System.currentTimeMillis();
		int rdLen = 0;
		boolean rev1 = false;
		boolean rev2 = false;
		
		while (br.ready()){
			currPos = fis.getChannel().position()-start;
			if (((double)currPos/len)*10 > ten){
				System.out.print(".."+NF.format(10*ten)+"%");
				ten++;
			}
			line1 = br.readLine().split("\t");
			line2 = br.readLine().split("\t");
			line1[0] = HelperFunctions.trimPairNumber(line1[0]);
			line2[0] = HelperFunctions.trimPairNumber(line2[0]);
			while (!line1[0].equals(line2[0]) && br.ready()){
				line1 = line2;
				line2 = br.readLine().split("\t");				
				line2[0] = HelperFunctions.trimPairNumber(line2[0]);
			}
			total++;
			left1 = Integer.parseInt(line1[3]);
			left2 = Integer.parseInt(line2[3]);
			// if pair didn't map, just jump to next
			// line instead of jumping a full step
			if (left1 == 0 || left2 == 0) continue;
			
			ctg1 = ctgs.get(line1[2]);
			ctg2 = ctgs.get(line2[2]);

			if (line1[5].equals("*") || line2[5].equals("*"))
				continue;
			
			int tmpLen = HelperFunctions.cigarLength(line1[5]);
			if (tmpLen > rdLen)
				rdLen = tmpLen;
			tmpLen = HelperFunctions.cigarLength(line2[5]);
			if (tmpLen > rdLen)
				rdLen = tmpLen;
			
			// Remove repetitive reads
			//if (line1[11].equals("XT:A:R")||line2[11].equals("XT:A:R")){
			//	continue;
			//}
			
			rev1 = HelperFunctions.isReverse(line1[1]);
			rev2 = HelperFunctions.isReverse(line2[1]);
			
			/* begin: tally these read positions */
			// tally read 1
			tmpWin = readCounts.get(ctg1.name);
			index = Arrays.binarySearch(tmpWin[0], left1);
			if (index < 0)
				index = -1*(index+1);
			if (rev1)
				tmpWin[3][index]++;
			else 
				tmpWin[2][index]++;
			// tally read 2
			tmpWin = readCounts.get(ctg2.name);
			index = Arrays.binarySearch(tmpWin[0], left1);
			if (index < 0)
				index = -1*(index+1);
			if (rev1)
				tmpWin[3][index]++;
			else 
				tmpWin[2][index]++;
			/* end: tally these read positions */
			
			
			ctgNameComp = line1[2].compareTo(line2[2]);
			
			/*
			 * if this pair spans two different contigs, we need to sort the names
			 * for consistency before we add the match
			 */
			if (ctgNameComp != 0){
				if (ctgNameComp < 0){
					ctgStr = line1[2]+"-"+line2[2];
					if (clusterers.containsKey(ctgStr))
						pc = clusterers.get(ctgStr);
					else {
						pc = new SpatialClusterer(ctg1, ctg2);
						clusterers.put(ctgStr, pc);
					}
					pc.addMatch(left1, rev1, left2, rev2);
				} else {
					ctgStr = line2[2]+"-"+line1[2];
					if (clusterers.containsKey(ctgStr))
						pc = clusterers.get(ctgStr);
					else { 
						pc = new SpatialClusterer(ctg2, ctg1);
						clusterers.put(ctgStr, pc);
					}
					pc.addMatch(left2, rev2, left1, rev1);
				}
				if (counts.containsKey(ctg1.name))
					counts.put(ctg1.name, counts.get(ctg1.name)+2);
				else
					counts.put(ctg1.name, 2);
				if (counts.containsKey(ctg2.name))
					counts.put(ctg2.name, counts.get(ctg2.name)+2);
				else
					counts.put(ctg2.name, 2);
			} else { // same contig, so check to see it's within the given ranges
				int ins = (left2 > left1 ? left2+HelperFunctions.cigarLength(line2[5])-left1 : left1+HelperFunctions.cigarLength(line1[5])-left2);
				if (inRange(ranges,ins)) 
						continue;
				
				ctgStr = line2[2]+"-"+line1[2];
				if (clusterers.containsKey(ctgStr))
					pc = clusterers.get(ctgStr);
				else { 
					pc = new SpatialClusterer(ctg2, ctg1);
					clusterers.put(ctgStr, pc);
				}
				if (left2 < left1) // order for consistency
					pc.addMatch(left2, rev2, left1, rev1);
				else 
					pc.addMatch(left1, rev1, left2, rev2);
				if (counts.containsKey(ctg1.name)){
					counts.put(ctg1.name, counts.get(ctg1.name)+2);
				} else {
					counts.put(ctg1.name, 2);
				}
				
			}
			// add point for contig1, and keep track of which MatchBuilders are associated with contig1
			if (ctgClusterers.containsKey(ctg1.name)){
				tmpMBs = ctgClusterers.get(ctg1.name);
			} else {
				tmpMBs = new Vector<String>();
				ctgClusterers.put(ctg1.name, tmpMBs);
			}
			tmpMBs.add(ctgStr);
			// add point for contig2, and keep track of which MatchBuilders are associated with contig2
			if (ctgClusterers.containsKey(ctg2.name)){
				tmpMBs = ctgClusterers.get(ctg2.name);
			} else {
				tmpMBs = new Vector<String>();
				ctgClusterers.put(ctg2.name, tmpMBs);
			}
			tmpMBs.add(ctgStr);
			numKeep++;
		}
		long after = System.currentTimeMillis();
		System.out.println("..100%... done!... Took "+HelperFunctions.millisToTimeString(after-before));
		perc = (double) numKeep / total * 100;
		System.out.println("[a5_qc_load_data] Keeping "+NF.format(perc)+"% ("+numKeep+"/"+total+") of reads.");
		ReadCluster.RDLEN = rdLen;
		/*
		 * Set LAMDBA, our Poisson rate parameter. We will use this to 
		 * compute key runtime parameters
		 */
		P = Double.POSITIVE_INFINITY;
		int minWindow = -1;
		File covFile = new File (samFile.getParentFile(),HelperFunctions.basename(samFile.getName(),".sam")+".cov");
		covFile.createNewFile();
		
		PrintStream covOut = new PrintStream(covFile);
		Iterator<String> ctgIt = readCounts.keySet().iterator();
		String tmpCtg = null;
		String minWinCtg = null;
		double tmpFreq = 0.0;
		while(ctgIt.hasNext()){
			tmpCtg = ctgIt.next();
			covOut.println("#"+tmpCtg);
			tmpWin = readCounts.get(tmpCtg);
			for (int i = 0; i < tmpWin[2].length; i++){
				if (tmpWin[2][i] != 0){
					tmpFreq = tmpWin[2][i];
					tmpFreq = tmpFreq/tmpWin[1][i];
					if (P > tmpFreq){
						P = tmpFreq;
						minWindow = i;
						minWinCtg = tmpCtg;
					}
				}
				covOut.println(tmpWin[0][i]+"\t"+tmpWin[1][i]+"\t"+tmpWin[2][i]+"\t"+tmpWin[3][i]);
			}
		}
		covOut.close();
		int minWinL;
		int minWinR;
		if (minWindow == 0) {
			minWinL = 1;
			minWinR = readCounts.get(minWinCtg)[0][minWindow];
		} else {
			minWinL = readCounts.get(minWinCtg)[0][minWindow-1];
			minWinR = readCounts.get(minWinCtg)[0][minWindow];
		}
		System.out.println("[a5_load_data] Window with fewest mapped reads: "+minWinCtg+" "+minWinL+" - "+minWinR);
		
		/*
		 * Remove contigs that don't have enough points on them. 
		 */
		ctgIt = ctgClusterers.keySet().iterator();
		Set<String> ctgToRm = new HashSet<String>();
		Set<String> psPairsToRm = new HashSet<String>();
		while(ctgIt.hasNext()){
			tmp = ctgIt.next();
			if (counts.get(tmp) < MIN_PTS) {
				ctgToRm.add(tmp);
				psPairsToRm.addAll(ctgClusterers.get(tmp));
			} 
		}
		removeKeys(ctgs, ctgToRm);
		removeKeys(clusterers, psPairsToRm);
	
		matches = clusterers.values();
	}
	
	/**
	 * Set parameters and print their values to standard out
	 * @param ranges <code>[ mean, sd, n, nSd , min, max ]</code>
	 */
	private static void setParameters(double[][] ranges){
	//	P = 0.025;
	//	NumberFormat nf = NumberFormat.getInstance();
	//	nf.setMaximumFractionDigits(3);
	//	System.out.println("[a5_qc] Setting P to "+nf.format(P));
		int maxSd = 0;
		for (int i = 0; i < ranges.length; i++)
			if (ranges[i][1] > maxSd)
				maxSd = (int) (ranges[i][1]);
		/*
		 * By modelling read mapping points as a Poisson process, the distance
		 * between points is exponential. Take the 1-ALPHA quantile to get
		 * a maximum distance between two points.
		 *
		 * LAMBDA Rate of mapping points (Poisson rate parameter)
		 */	
		MAX_INTERPOINT_DIST = Math.max(ReadCluster.RDLEN, (int) (Math.log(ALPHA)/Math.log(Math.max(1-P,0)))-1);
		//MIN_BLOCK_LEN = 2*MAX_INTERPOINT_DIST;
		MIN_BLOCK_LEN = (int) (1/P)*2;
		for (double[] cluster: ranges){
			if (cluster[0] >= MEAN_BLOCK_LEN){
				MEAN_BLOCK_LEN = (int) (cluster[0]);
				MAX_BLOCK_LEN = (int) (cluster[0]+cluster[1]*6);
				MIN_BLOCK_LEN = (int) Math.max(cluster[0]-cluster[1]*3,ReadCluster.RDLEN);
			}
		}
		/*
		 * apply some extreme value theory to get the minimum of points 
		 * randomly sampled uniformally across an interval of MAX_BLOCK_LEN
		 * 
		 * MAX_INTERBLOCK_DIST = 2*((1-ALPHA)^(1/EXP_POINTS) * MAX_BLOCK_LEN)
		 * 
		 */
		MAX_INTERBLOCK_DIST = (int)(2*(Math.pow(1-ALPHA, 1/(P*MEAN_BLOCK_LEN))*MEAN_BLOCK_LEN-1));
		MAX_INTERBLOCK_DIST = 2*MEAN_BLOCK_LEN;
		SpatialClusterer.MIN_PTS = (int) (P * 2*MAX_INTERPOINT_DIST);
		SpatialClusterer.MIN_PTS = 20;
		SpatialClusterer.EPS = MAX_INTERPOINT_DIST;
	}
	
	/**
	 * Compute the maximum allowed block length.
	 * @param clusterStats an array of arrays with each cluster's statistics (mean, sd, n, numSd)
	 */
	private static void setMAXBLOCKLEN(double[][] clusterStats){
		for (double[] cluster: clusterStats){
			if (cluster[0] > MEAN_BLOCK_LEN){
				MEAN_BLOCK_LEN = (int) (cluster[0]);
				MAX_BLOCK_LEN = (int) (cluster[0]+cluster[1]*6);
			}
		}
	}
	/**
	 * Print the values of major parameters
	 */
	private static void printParams(){
		System.out.println("[a5_qc] parameters:");
		System.out.println("        ALPHA               = " + ALPHA);
		System.out.println("        P                   = " + P);
		System.out.println("        MIN_BLOCK_LEN       = " + MIN_BLOCK_LEN);
		System.out.println("        MEAN_BLOCK_LEN      = " + MEAN_BLOCK_LEN);
		System.out.println("        MAX_BLOCK_LEN       = " + MAX_BLOCK_LEN);
		System.out.println("        MAX_INTERBLOCK_DIST = " + MAX_INTERBLOCK_DIST);
		System.out.println("        MAX_INTERPOINT_DIST = " + MAX_INTERPOINT_DIST);
		System.out.println("        EPSILON             = " + SpatialClusterer.EPS);
		System.out.println("        MIN_POINTS          = " + SpatialClusterer.MIN_PTS);
		System.out.println("        RDLEN               = " + ReadCluster.RDLEN);

	}
	
	/**
	 * Identify if we have innies and/or outies in this set.
	 * @param reads
	 */
	private static void setOrientation(Collection<ReadPair> reads){
		NIN = 0;
		NOUT = 0;
		Iterator<ReadPair> it = reads.iterator();
		ReadPair tmp = null;
		while(it.hasNext()){
			tmp = it.next();
			if (tmp.inward)
				NIN++;
			else if (tmp.outward)
				NOUT++;
		}
		double total = NIN + NOUT;
		if (NIN/total > 0.1)
			INWARD = true;
		if (NOUT/total > 0.1)
			OUTWARD = true;
	}
	
	/**
	 * Return true of <code>ins</code> fits into any of the given ranges
	 * @param ranges
	 * @param ins
	 * @return
	 */
	private static boolean inRange(double[][] ranges, double ins){
		for (int i = 0; i < ranges.length; i++)
			if (ins >= ranges[i][0] && ins <= ranges[i][1])
				return true;
		return false;
	}

	/**
	 * Returns an array of arrays with each array having the following information:	<br><br>
	 * <code>[ mean, sd, n, nSd , min, max ]</code>
	 * 
	 * @return an array or arrays, where individual arrays contain ReadSet stats
	 */
	private static double[][] getLibraryStats(Map<String,ReadPair> reads, int numLibs){
		int maxK = 20;
		double delta = 0.00000005;
		Vector<ReadPair> toFilt = new Vector<ReadPair>();
		Iterator<ReadPair> rpIt = reads.values().iterator();
		ReadPair tmp = null;
		while(rpIt.hasNext()){
			tmp = rpIt.next();
			if (tmp.paired && tmp.ctg1.equals(tmp.ctg2))
				toFilt.add(tmp);
		}
		// estimate initial insert size to determine if we should look for shadow library.
		double[] ins = ReadPair.estimateInsertSize(reads.values());
		int bestModel = 0;
		EMClusterer[] models = new EMClusterer[maxK];
		models[0] = runEM(toFilt, 2, delta);
		double maxL = models[0].likelihood();
		double prevL = maxL;
		int numWorseSteps = 0;
		//int maxWorseSteps = 1;		
		int maxWorseSteps = 0;

		for (int i = 1; i < maxK && numWorseSteps < maxWorseSteps; i++){
			models[i] = runEM(toFilt, i+2, delta);
			if (models[i].likelihood() > maxL){
				bestModel = i;
				maxL = models[i].likelihood();
			}
			if (prevL > models[i].likelihood())
				numWorseSteps++;
			
		}
		/* end: cluster read pairs  */
		System.out.println("[a5_qc] Found "+(bestModel+1)+" clusters.");
		ReadSet[] clusters = new ReadSet[models[bestModel].getClusters().size()];
		models[bestModel].getClusters().toArray(clusters);
		// sort clusters so we can keep the top numLibs underdispersed clusters
		Arrays.sort(clusters,new Comparator<ReadSet>(){
			public int compare(ReadSet x, ReadSet y) {
				double rx = x.sd()/x.mean();
				double ry = y.sd()/y.mean();
				if (rx < ry) 
					return -1;
				else if (rx > ry) 
					return 1;
				else
					return 0;
			}
		});
		
		System.out.println("[a5_qc] Found the following clusters:");	
		
		/* 
		 * begin: find all clusters that look like true signal 
		 *        Use the cluster with the highest standard deviation 
		 *        to compute the maximum residual allowed for fitting
		 *        a point to a KClump (see KClump.fit(MatchPoint)) 
		 */
		Vector<ReadSet> signal = new Vector<ReadSet>();
		for (int i = 0; i < clusters.length; i++){
			NF.setMaximumFractionDigits(0);
			System.out.print("[a5_qc] cluster"+NF.format(clusters[i].getId())+": mu="+HelperFunctions.pad(NF.format(clusters[i].mean()),10)+
					"sd="+HelperFunctions.pad(NF.format(clusters[i].sd()),10)+"n="+HelperFunctions.pad(NF.format(clusters[i].size()),10));
			NF.setMaximumFractionDigits(2);
			double perc = 100*((double)clusters[i].size())/toFilt.size();
			System.out.print("perc="+HelperFunctions.pad(NF.format(perc),10));
			
				
			// if insert size distribution is under dispersed, add all these reads to the signal pile
			// take up to numLibs underdispersed clusters 
			if (clusters[i].sd() <= clusters[i].mean() && i < numLibs){
				signal.add(clusters[i]);
				System.out.println("  (signal)");
			} else {
				System.out.println("  (noise)");
			}
			
		}
		/* end: find all clusters that look like true signal */
		
		Iterator<ReadSet> sigIt = null;
		ReadSet sigSet = null;

		String rmClusters = "";
		sigIt = signal.iterator();
		int nSd = 6;
		double[][] ret = new double[signal.size()][];
		double min = Double.POSITIVE_INFINITY;
		double max = Double.NEGATIVE_INFINITY;
		int i = 0;
		while(sigIt.hasNext()){
			sigSet = sigIt.next();
			rmClusters += " cluster"+sigSet.getId();
			rpIt = sigSet.getReads().iterator();
			// Find the min and max values for this ReadSet
			while(rpIt.hasNext()){
				ReadPair tmpRp = rpIt.next();
				if (tmpRp.getInsert() < min)
					min = tmpRp.getInsert();
				if (tmpRp.getInsert() > max)
					max = tmpRp.getInsert();
			}
			
			nSd = Math.min(((int)sigSet.mean())/((int)sigSet.sd()),6);
			ret[i] = new double[6];
			ret[i][0] = sigSet.mean();
			ret[i][1] = sigSet.sd();
			ret[i][2] = sigSet.size();
			ret[i][3] = nSd;
			ret[i][4] = min;
			ret[i][5] = max;
			removeKeys(reads,sigSet.getReadHdrs());
			i++;
		}
		if (rmClusters.length()>0)
			System.out.println("[a5_qc] Removing "+rmClusters);
		ins = ReadPair.estimateInsertSize(reads.values());
		System.out.println("[a5_qc] Final stats for sample after filtering: mu="+NF.format(ins[0])+" sd="+NF.format(ins[1])+" n="+NF.format(ins[2]));
		
		return ret;
	}
	
	private static double[][] getFilterRanges(double[][] clusterStats){
		double[][] ranges = new double[clusterStats.length][2];
		for (int i = 0; i < clusterStats.length; i++){
			ranges[i][0] = Math.max(1,clusterStats[i][0]-6*clusterStats[i][1]);
			ranges[i][1] = clusterStats[i][0]+6*clusterStats[i][1];
			/*
			if (clusterStats[i][0]>1000){
				/* Use as many standard deviations as possible without going negative
				ranges[i][0] = clusterStats[i][0]-clusterStats[i][3]*clusterStats[i][1];
				ranges[i][1] = clusterStats[i][0]+clusterStats[i][3]*clusterStats[i][1];
				* /
				// Use 6 standard devations to the right, and as many as possible to the left
				ranges[i][0] = Math.max(1,clusterStats[i][0]-6*clusterStats[i][1]);
				ranges[i][1] = clusterStats[i][0]+6*clusterStats[i][1];
				
			} else {
				ranges[i][0] = 1;
				ranges[i][1] = clusterStats[i][0]*2;
			}
			*/
		}
		return ranges;
	}
	
	
	/**
	 * Run the EM-clustering algorithm with the given number of clusters
	 * @param toFilt the ReadPairs to cluster
	 * @param K the number of clusters to partition <code>toFilt</code> into
	 * @param delta the threshold for determining convergence of the algorithm
	 * @return the EMClusterer object used to run the algorithm
	 */
	private static EMClusterer runEM(Collection<ReadPair> toFilt, int K, double delta){
		System.out.print("[a5_qc] EM-clustering insert sizes with K="+K+"... ");
		/* begin: cluster read pairs by insert size so we can efficiently determine the insert size of this library */
		EMClusterer em = new EMClusterer(toFilt, K);
		long before = System.currentTimeMillis();
		int iters = em.iterate(1000, delta);
		long after = System.currentTimeMillis();
		double LK3 = em.likelihood();
		System.out.println("stopping after "+iters+" iterations with delta="+delta+". L = "+LK3+". Took "+(after-before)/1000+" seconds.");
		return em;
	}
	
	/**
	 * A function for removing keys from a Map
	 */
	private static <V> void removeKeys(Map<String,V> reads, Collection<String> torm){
		Iterator<String> it = torm.iterator();
		while(it.hasNext())
			reads.remove(it.next());
	}
	
	/**
	 * Read in contigs from the header of a SAM file
	 */
	private static Map<String,Contig> readContigs(RandomAccessFile raf) throws IOException {
		Map<String,Contig> contigs = new HashMap<String,Contig>();
		String[] hdr = null;
		String name = null;
		while(HelperFunctions.nextCharIs(raf,'@')){
			hdr = raf.readLine().split("\t");
			if (!hdr[0].equals("@SQ"))
				continue;
			int len = -1;
			for (String s: hdr){
				if (s.startsWith("SN")){
					name = s.substring(3);
				} else if (s.startsWith("LN")){
					len = Integer.parseInt(s.substring(3));
				}
			}
			if (name == null) 
				System.err.println("[a5_qc] Found nameless contig in SAM header");
			else if (len == -1) 
				System.err.println("[a5_qc] Found contig of unknown length in SAM header");
			contigs.put(name, new Contig(name,len));
		}
		return contigs;
	}
	
	/**
	 * Reads in a subset of read pairs from a SAM file.
	 * 
	 * This method assumes that the contig header has been read in, and 
	 * the file pointer is at the beginning of the first read's line.
	 * 
	 * @param raf
	 * @throws IOException 
	 */
	private static Map<String,ReadPair> readSubsetByChunk(RandomAccessFile raf, Map<String,Contig> contigs) throws IOException{
		long pos = raf.getFilePointer();
		long len = raf.length();
		
		long step = (len-pos)*HANDFUL_SIZE/N_ESTREADPAIRS;
		int i = 0;
		String[] line1 = null;
		String[] line2 = null;
		int left1 = 0;
		int left2 = 0;
		Map<String,ReadPair> reads = new HashMap<String, ReadPair>();
		ReadPair tmp = null;
		EfficientSAMFileSampler esfs = new EfficientSAMFileSampler(raf, HANDFUL_SIZE*SAM_LINE_LEN, step);
		String[][] lines = null;
		while (i < N_ESTREADPAIRS && esfs.hasNextPair()){
			lines = esfs.nextPair();
			line1 = lines[0];
			line2 = lines[1];
			left1 = Integer.parseInt(line1[3]);
			if (line2.length < 4){
				System.out.print("");
			}
			left2 = Integer.parseInt(line2[3]);
			// if pair didn't map, just jump to next line instead of jumping a full step
			if (left1 == 0 || left2 == 0) continue;
			tmp = new ReadPair(line1[0]);
			boolean rev1 = HelperFunctions.isReverse(line1[1]);
			boolean rev2 = HelperFunctions.isReverse(line2[1]);
			
			if (contigs.get(line1[2]).equals(contigs.get(line2[2])) && rev1 != rev2 &&
					// If this information is present, require this mapping to be unique. 
					(line1.length > 11 && line2.length > 11)&& line1[11].equals(TAG_KEEP) && line2[11].equals(TAG_KEEP)){
				if (left1 < left2){
					if (rev1){
						contigs.get(line1[2]).addOut();
						contigs.get(line2[2]).addOut();
					} else {
						contigs.get(line1[2]).addIn();
						contigs.get(line2[2]).addIn();
					}
				} else {
					if (rev1){
						contigs.get(line1[2]).addIn();	
						contigs.get(line2[2]).addIn();
					} else {
						contigs.get(line1[2]).addOut();
						contigs.get(line2[2]).addOut();						
					}
				}
			}
			
			if (line1[5].equals("*") || line2[5].equals("*"))
				continue;
			
			tmp.addRead(left1, rev1, HelperFunctions.cigarLength(line1[5]), 
					contigs.get(line1[2]), Integer.parseInt(line1[4]), line1[5]);
			tmp.addRead(left2, rev2, HelperFunctions.cigarLength(line2[5]), 
					contigs.get(line2[2]), Integer.parseInt(line2[4]), line2[5]);
			reads.put(line1[0], tmp);
			i++;
		}
		return reads;
	}	
	
	/**
	 * A class for reading and sampling readpairs from a SAM file via random access. 
	 * @author Andrew Tritt
	 *
	 */
	private static class EfficientSAMFileSampler {
		private RandomAccessFile raf;
		private byte[] buf;
		private StringTokenizer tok;
		private int tokLeft;
		private int currChunkSize;
		private long step;
		public EfficientSAMFileSampler(RandomAccessFile file, int bufSize, long step) throws IOException{
			raf = file;
			buf = new byte[bufSize];
			tokLeft = 0;
			currChunkSize = 0;
			this.step = step;
			resetTok();
		}
		
		/**
		 * Returns an array of length 2, with two reads.
		 * @return
		 * @throws IOException 
		 */
		public String[][] nextPair() throws IOException{
			if (raf.getFilePointer() == raf.length())
				throw new IOException("Reached end of file.");
			if (currChunkSize == HANDFUL_SIZE){ 
				// read our handful, so jump to our next step
				raf.seek(raf.getFilePointer()+step);
				resetTok();
			} else if (tokLeft < 3){
				resetTok();
			}
			String[] line1 = tok.nextToken().split("\t");
			String[] line2 = tok.nextToken().split("\t");
			tokLeft -= 2;
			// get to the first read in a pair 
			while (!line1[0].equals(line2[0]) && line1.length < 12 && line2.length < 12 && tok.hasMoreTokens()){
				line1 = line2;
				line2 = tok.nextToken().split("\t");
				tokLeft--;
				if (!tok.hasMoreTokens())
					resetTok();
			}
			String[][] ret = {line1,line2};
			if (line1.length < 11) 
				System.out.print("");
			if ( line2.length < 4)
				System.out.print("");
			currChunkSize++;
			return ret;
		}
		
		public boolean hasNextPair() throws IOException{
			long bytesLeft = raf.length()-raf.getFilePointer();
			boolean ret = tokLeft >= 3 || (bytesLeft >= buf.length);
			if (!ret){
				System.out.print("");
			}
			return ret;
		}
		
		private void resetTok() throws IOException{
			raf.read(buf);
			tok = new StringTokenizer(new String(buf),"\n");
			tok.nextToken(); // make sure we're at the beginning of a line
			tokLeft = tok.countTokens()-1;
		}
	}
	/*
	public static Map<String,ReadPair> sampleBAMFile(SAMFileReader sam, Map<String,Contig> contigs, long seed){	
		Map<String,ReadPair> pairs = new HashMap<String,ReadPair>();
		Random rand = new Random(seed);
		// First we need to count the number of reads in the BAM file
		// so we know the frequency at which to sample pairs
		int count = 0;
		long before = System.currentTimeMillis();
		AbstractBAMFileIndex index = (AbstractBAMFileIndex) sam.getIndex();
		int nRefs = index.getNumberOfReferences();
		for (int i = 0; i < nRefs; i++)
			count += index.getMetaData(i).getAlignedRecordCount();
		long after = System.currentTimeMillis();
		System.out.println("[a5_qc] Found "+ count +" total reads. Took "+ HelperFunctions.millisToTimeString(after-before)+ ".");
		
		//count = count/2; // divide by two to get the number of pairs
		double freq = 2 * 100000/((double) count);
		System.out.println("[a5_qc] Sampling with frequencey "+new BigDecimal(freq).round(new MathContext(3)).doubleValue());
		SAMRecordIterator it = sam.iterator();
		SAMRecord tmpRecord = null;
		ReadPair tmpPair = null;
		while (it.hasNext()){
			tmpRecord = it.next();
			// if the two reads in this fragment don't map to the same contig
			if (tmpRecord.getReferenceName() != tmpRecord.getMateReferenceName())
				continue;
			// we already came across this guys mate
			if (tmpRecord.getAlignmentStart() > tmpRecord.getMateAlignmentStart()) { 
				// and it didn't get sampled, so move on
				if (!pairs.containsKey(tmpRecord.getReadName())) 
					continue;
				else 
					tmpPair = pairs.get(tmpRecord.getReadName());
			} else { // keep with probability=freq
				if (rand.nextDouble() > freq)
					continue;
				else // keeping this fragment
					tmpPair = new ReadPair(tmpRecord.getReadName());
			}
			tmpPair.addRead(tmpRecord.getAlignmentStart(), 
					        tmpRecord.getReadNegativeStrandFlag(), 
					        tmpRecord.getAlignmentEnd() - tmpRecord.getAlignmentStart(), 
					        contigs.get(tmpRecord.getReferenceName()), 
					        tmpRecord.getMappingQuality(), 
					        tmpRecord.getCigarString());
			pairs.put(tmpRecord.getReadName(), tmpPair);
		}
		return pairs;
	}
	*/
	public static Map<String,Contig> getContigs(SAMFileReader sam) {
		SAMFileHeader hdr = sam.getFileHeader();
		Iterator<SAMSequenceRecord> it = hdr.getSequenceDictionary().getSequences().iterator();
		SAMSequenceRecord tmp = null;
		Map<String,Contig> ret = new HashMap<String,Contig>();
		while(it.hasNext()){
			tmp = it.next();
			ret.put(tmp.getSequenceName(),new Contig(tmp.getSequenceName(),tmp.getSequenceLength()));
		}
		return ret;
	}
	/*
	public static void loadBAMData(String bamPath, Map<String,Contig> ctgs, double[][] ranges) throws IOException{
		for (int i = 0; i < ranges.length; i++)
			System.out.println("[a5_qc] Filtering read pairs with inserts between "+ NF.format(ranges[i][0])+"-"+NF.format(ranges[i][1]));
		Map<String,SpatialClusterer> clusterers = new HashMap<String,SpatialClusterer>();
		// Keep track of all clusterers that each contig is associated with.
		Map<String,Vector<String>> ctgClusterers = new HashMap<String,Vector<String>>();
		Map<String,Integer> counts = new HashMap<String,Integer>();
		Vector<String> tmpMBs = null;

		File samFile = new File(bamPath);
		FileInputStream fis = new FileInputStream(samFile);
		long start = fis.getChannel().position();
		long len = fis.getChannel().size() - start;
		SAMFileReader samRdr = new SAMFileReader(fis);
		
		// begin: build a lookup table for each contig for tallying coverage in windows 
		
		// require the window length to be at least 1000 base pairs
		int windowLen = 1000;
		//Map<String,Integer> coordOffset = new HashMap<String,Integer>();
		Map<String,int[][]> readCounts = new HashMap<String,int[][]>();
		int[][] tmpWin = null;
		Iterator<String> ctgIt = ctgs.keySet().iterator();
		Contig tmpCtg = null;
		while(ctgIt.hasNext()){
			tmpCtg = ctgs.get(ctgIt.next());
			// calculate the number of windows we need for this contig
			int numWindow = tmpCtg.len/windowLen;
			if (tmpCtg.len % windowLen != 0)
				numWindow++;
			// now create the lookup table: add the upper bounds of each window to the array
			tmpWin = new int[4][numWindow];
			tmpWin[0][numWindow-1] = tmpCtg.len;
			tmpWin[1][numWindow-1] = (tmpCtg.len - (numWindow-1)*windowLen);
			for (int j = 0; j < numWindow-1; j++){
				tmpWin[0][j] = (j+1)*windowLen;
				tmpWin[1][j] = windowLen;
			}
			// store our lookup table for this contig
			readCounts.put(tmpCtg.name,tmpWin);
		}
		// end: build a lookup table for tallying coverage in windows 
		
		int left1 = 0;
		int left2 = 0;
		String ctgStr = null;
		String tmp = null;
		SpatialClusterer pc = null;
		String ctg1 = null;
		String ctg2 = null;
		int ctgNameComp = -10;
		System.out.print("[a5_qc] Reading SAM file...");
		long currPos = start;
		double perc = 0;
		double ten = 1;
		int numKeep = 0;
		int total = 0;
		long before = System.currentTimeMillis();
		int rdLen = 0;
		boolean rev1 = false;
		boolean rev2 = false;
		
		Iterator<SAMRecord> samIt = samRdr.iterator();
		SAMRecord sam = null;
		
		int ctgComp = -1;
		while (samIt.hasNext()){
			currPos = fis.getChannel().position()-start;
			if (((double)currPos/len)*10 > ten){
				System.out.print(".."+NF.format(10*ten)+"%");
				ten++;
			}
			sam = samIt.next();
			
			if (sam.getMateUnmappedFlag() || sam.getReadUnmappedFlag())
				continue;
			int tmpLen = sam.getAlignmentEnd()-sam.getAlignmentStart();
			if (tmpLen > rdLen)
				rdLen = tmpLen;
			left1 = sam.getAlignmentStart();
			left2 = sam.getMateAlignmentStart();
			ctg1 = sam.getReferenceName();
			ctg2 = sam.getMateReferenceName();
			
			// I love referential equality
			if (ctg1 == ctg2) {
				// Assume we've already seen this fragment before because the file is sorted
				if (left1 > left2)
					continue;
				else {
					ctgComp = 0;
					ctg2 = ctg1;
					total++;
				}
			} else {
				ctgComp = ctg1.compareTo(ctg2);
				//
				// Assuming contigs are sorted lexicographically, if ctg2 is greater than
				// ctg1, we must have come across this fragment already
				//
				if (ctgComp > 0)
					continue;
			}
			
			
			rev1 = sam.getReadNegativeStrandFlag();
			rev2 = sam.getMateNegativeStrandFlag();
			
			// begin: tally these read positions 
			// tally read 1
			tmpWin = readCounts.get(ctg1);
			if (rev1)
				tmpWin[3][left1/windowLen]++;
			else 
				tmpWin[2][left1/windowLen]++;

			// tally read 2
			tmpWin = readCounts.get(ctg2);
			if (rev2)
				tmpWin[3][left2/windowLen]++;
			else 
				tmpWin[2][left2/windowLen]++;
			// end: tally these read positions 
			
			//
			// if this pair spans two different contigs, we need to sort the names
			// for consistency before we add the match. because we are assuming the BAM
			// is sorted by contig lexicographically, ctg1 should always come first
			//
			if (ctgNameComp != 0){
				ctgStr = ctg1+"-"+ctg2;
				if (clusterers.containsKey(ctgStr))
					pc = clusterers.get(ctgStr);
				else {
					pc = new SpatialClusterer(ctgs.get(ctg1), ctgs.get(ctg2));
					clusterers.put(ctgStr, pc);
				}
				pc.addMatch(left1, rev1, left2, rev2);
				
				if (counts.containsKey(ctg1))
					counts.put(ctg1, counts.get(ctg1)+1);
				else
					counts.put(ctg1, 1);
				if (counts.containsKey(ctg2))
					counts.put(ctg2, counts.get(ctg2)+1);
				else
					counts.put(ctg2, 1);
			} else { // same contig, so check to see it's within the given ranges
				if (inRange(ranges,Math.abs(sam.getInferredInsertSize()))) 
						continue;
				ctgStr = ctg1;
				if (clusterers.containsKey(ctgStr))
					pc = clusterers.get(ctgStr);
				else { 
					pc = new SpatialClusterer(ctgs.get(ctg2), ctgs.get(ctg1));
					clusterers.put(ctgStr, pc);
				}
				if (left2 < left1) // order for consistency
					pc.addMatch(left2, rev2, left1, rev1);
				else 
					pc.addMatch(left1, rev1, left2, rev2);
				if (counts.containsKey(ctg1)){
					counts.put(ctg1, counts.get(ctg1)+2);
				} else {
					counts.put(ctg1, 2);
				}
				
			}
			// add point for contig1, and keep track of which MatchBuilders are associated with contig1
			if (ctgClusterers.containsKey(ctg1)){
				tmpMBs = ctgClusterers.get(ctg1);
			} else {
				tmpMBs = new Vector<String>();
				ctgClusterers.put(ctg1, tmpMBs);
			}
			tmpMBs.add(ctgStr);
			// add point for contig2, and keep track of which MatchBuilders are associated with contig2
			if (ctgClusterers.containsKey(ctg2)){
				tmpMBs = ctgClusterers.get(ctg2);
			} else {
				tmpMBs = new Vector<String>();
				ctgClusterers.put(ctg2, tmpMBs);
			}
			tmpMBs.add(ctgStr);
			numKeep++;
		}
		long after = System.currentTimeMillis();
		System.out.println("..100%... done!... Took "+HelperFunctions.millisToTimeString(after-before));
		perc = (double) numKeep / total * 100;
		System.out.println("[a5_qc] Keeping "+NF.format(perc)+"% ("+numKeep+"/"+total+") of pairs.");
		ReadCluster.RDLEN = rdLen;
		//
		// Set LAMDBA, our Poisson rate parameter. We will use this to 
		// compute key runtime parameters
		//
		P = Double.POSITIVE_INFINITY;
		int minWindow = -1;
		File covFile = new File (samFile.getParentFile(),HelperFunctions.basename(samFile.getName(),".sam")+".cov");
		covFile.createNewFile();
		
		PrintStream covOut = new PrintStream(covFile);
		ctgIt = readCounts.keySet().iterator();
		String tmpCtgName = null;
		String minWinCtg = null;
		double tmpFreq = 0.0;
		while(ctgIt.hasNext()){
			tmpCtgName = ctgIt.next();
			covOut.println("#"+tmpCtgName);
			tmpWin = readCounts.get(tmpCtgName);
			for (int i = 0; i < tmpWin[2].length; i++){
				if (tmpWin[2][i] != 0){
					tmpFreq = tmpWin[2][i];
					tmpFreq = tmpFreq/tmpWin[1][i];
					if (P > tmpFreq){
						P = tmpFreq;
						minWindow = i;
						minWinCtg = tmpCtgName;
					}
				}
				covOut.println(tmpWin[0][i]+"\t"+tmpWin[1][i]+"\t"+tmpWin[2][i]+"\t"+tmpWin[3][i]);
			}
		}
		covOut.close();
		int minWinL;
		int minWinR;
		if (minWindow == 0) {
			minWinL = 1;
			minWinR = readCounts.get(minWinCtg)[0][minWindow];
		} else {
			minWinL = readCounts.get(minWinCtg)[0][minWindow-1];
			minWinR = readCounts.get(minWinCtg)[0][minWindow];
		}
		System.out.println("[a5_load_data] Window with fewest mapped reads: "+minWinCtg+" "+minWinL+" - "+minWinR);
		
		//
		// Remove contigs that don't have enough points on them. 
		//
		ctgIt = ctgClusterers.keySet().iterator();
		Set<String> ctgToRm = new HashSet<String>();
		Set<String> psPairsToRm = new HashSet<String>();
		while(ctgIt.hasNext()){
			tmp = ctgIt.next();
			if (counts.get(tmp) < MIN_PTS) {
				ctgToRm.add(tmp);
				psPairsToRm.addAll(ctgClusterers.get(tmp));
			} 
		}
		removeKeys(ctgs, ctgToRm);
		removeKeys(clusterers, psPairsToRm);
	
		matches = clusterers.values();
	}
	*/
}
