/**

 * This file is part of the A5 pipeline.
 * (c) 2011, 2012 Andrew Tritt and Aaron Darling
 * This software is licensed under the GPL, v3.0. Please see the file LICENSE for details
 */
package org.halophiles.assembly.qc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
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

import org.halophiles.assembly.Contig;
import org.halophiles.assembly.ReadPair;
import org.halophiles.assembly.ReadSet;

public class MisassemblyBreaker {
	
	/**
	 * A Comparator for sorting blocks by right-most position
	 */
	private static Comparator<int[]> BLOCK_COMP = new Comparator<int[]>(){
		public int compare(int[] arg0, int[] arg1) {
			return arg0[1] - arg1[1];
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
			
			NF = NumberFormat.getInstance();
			NF.setMaximumFractionDigits(0);
			NF.setGroupingUsed(false);
			String samPath = args[0];
			int numLibs = 1;
			if (args.length == 4){
				numLibs = Integer.parseInt(args[3]);
			}
			System.out.println("[a5_qc] Reading "+samPath);
			File samFile = new File(samPath);
			RandomAccessFile raf = new RandomAccessFile(samFile, "r");
			Map<String,Contig> contigs = readContigs(raf);
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
			System.out.println("[a5_qc] Took "+((after-before)/1000)+" seconds to read in "+reads.size()+" read pairs.");
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
			
				
			double[][] clusterStats = getLibraryStats(reads, numLibs);
			double[][] ranges = getFilterRanges(clusterStats);
			setMAXBLOCKLEN(clusterStats);
			loadData(samPath, contigs, ranges);
			setParameters(clusterStats);
			printParams();
			
			// collect all of our blocks for each contig
			Iterator<SpatialClusterer> mbIt = matches.iterator();
			Map<String,Vector<int[]>> blocks = new HashMap<String,Vector<int[]>>();
			Vector<int[]> xBlocks = null;
			Vector<int[]> yBlocks = null;
			while(mbIt.hasNext()){
				SpatialClusterer pc = mbIt.next();
				// get Vector for holding contig X blocks
				xBlocks = new Vector<int[]>();
				// get Vector for holding contig Y blocks
				yBlocks = new Vector<int[]>();
				pc.buildReadPairClusters();
				addBlocks(pc, xBlocks, yBlocks);
				//pc.exportCurrState(new File(matchDir,"match."+pc.getContig1().getId()+"v"+pc.getContig2().getId()+".txt"));
				removeTerminalBlocks(pc.getContig1(), xBlocks);
				removeTerminalBlocks(pc.getContig2(), yBlocks);
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
				Vector<int[]> tmpBlocks = blocks.get(tmpCtg);
				System.out.println("[a5_qc] Found "+tmpBlocks.size()+" blocks on contig " + contigs.get(tmpCtg).getId());
				if (tmpBlocks.isEmpty())
					toRm.add(tmpCtg);
				else 
					Collections.sort(tmpBlocks, BLOCK_COMP);
				Iterator<int[]> blockIt = tmpBlocks.iterator();
				while(blockIt.hasNext()){
					int[] tmp = blockIt.next();
					System.out.println("        "+tmp[0]+(tmp[2]==1?" -> ":" <- ")+tmp[1]);
					
				}
				
				
			}
			removeKeys(blocks, toRm);
			
			if (blocks.isEmpty()){
				System.out.println("[a5_qc] No blocks were found. Not breaking scaffolds.");
				System.exit(0);
			}
			
			/*
			 *  break on regions of a minimum distance that are flanked by two blocks
			 */
			File brokenScafFile = new File(args[2]);
			brokenScafFile.createNewFile();
			ScaffoldExporter out = new ScaffoldExporter(brokenScafFile); 
			File ctgFile = new File(args[1]);
			BufferedReader br = new BufferedReader(new FileReader(ctgFile));
			File bedFile = new File(samFile.getParentFile(),basename(samFile.getName(),".sam")+".regions.bed");
			bedFile.createNewFile();
			System.out.println("[a5_qc] Writing regions potentially containing misassemblies to "+bedFile.getAbsolutePath());
			PrintStream bedOut = new PrintStream(bedFile);
			// discard the first '>'
			br.read();
			int[][] tmpAr = null;
			StringBuilder sb = null;
			while(br.ready()){
				String tmpCtg = br.readLine(); 
				// Take only the first part of the contig header to be consistent with bwa
				if (tmpCtg.contains(" "))
					tmpCtg = tmpCtg.substring(0,tmpCtg.indexOf(" "));
				sb = new StringBuilder();
				char c = (char) br.read();
				// read in this sequence
				while(c != '>'){
					if (isNuc(c))
						sb.append(c);
					if (!br.ready())
						break;
					c = (char) br.read();
				}
				
				if (!blocks.containsKey(tmpCtg)){
					out.export(tmpCtg, sb);
					continue;
				}
				
				Vector<int[]> tmpBlks = blocks.get(tmpCtg);

				if (tmpBlks.size() < 2){
					out.export(tmpCtg, sb);
					continue;
				}
			
				tmpAr = new int[tmpBlks.size()][];
				tmpBlks.toArray(tmpAr);
				Arrays.sort(tmpAr,BLOCK_COMP);
				int left = 1;
				int right = 1;
				for (int i = 1; i < tmpAr.length; i++){
					// if they don't face each other, there isn't a misassembly in between this pair of blocks
					if (tmpAr[i-1][2] != 1 || tmpAr[i][2] != -1) 
						continue;
					// if they overlap, split at the midpoint of overlap
					if (tmpAr[i-1][1] > tmpAr[i][0]){ 
						bedOut.println(tmpCtg+"\t"+tmpAr[i][0]+"\t"+tmpAr[i-1][1]);
						right = (tmpAr[i-1][1] + tmpAr[i][0])/2;
						System.out.println("[a5_qc] Exporting "+tmpCtg+" at "+left+"-"+right);
						out.export(tmpCtg, sb, left, right);
						left = right+1;
					} else if (tmpAr[i][0]-tmpAr[i-1][1] < MAX_INTERBLOCK_DIST) {
						bedOut.println(tmpCtg+"\t"+tmpAr[i-1][1]+"\t"+tmpAr[i][0]);
						right = tmpAr[i-1][1];
						System.out.println("[a5_qc] Exporting "+tmpCtg+" at "+left+"-"+right);
						out.export(tmpCtg, sb, left, right);
						left = tmpAr[i][0];
					} 
				}
				right = sb.length();
				System.out.println("[a5_qc] Exporting "+tmpCtg+" at "+left+"-"+right);
				out.export(tmpCtg, sb, left, right);
				bedOut.close();
				
				
				
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
	 * Identify and add significant blocks resulting from spatial clustering of read pairs.
	 * @param sc the SpatialClusterer to get blocks from
	 * @param xBlocks the Collection of blocks for Contig 1 (aka the x Contig)
	 * @param yBlocks the Collection of blocks for Contig 2 (aka the y Contig)
	 */
	private static void addBlocks(SpatialClusterer sc, Vector<int[]> xBlocks, Vector<int[]> yBlocks){
		ReadCluster[] rdClust = sc.getReadClusters();
		int xlen = 0;
		int ylen = 0;
		int[] x = null;
		int[] y = null;
		for (int i = 0; i < rdClust.length; i++){
			xlen = rdClust[i].xMax-rdClust[i].xMin;
			ylen = rdClust[i].yMax-rdClust[i].yMin;
			
			x = new int[3];
			x[0] = rdClust[i].xMin;
			x[1] = rdClust[i].xMax;
			x[2] = rdClust[i].xOri ? 1 : -1; 
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
			if (xlen >= MIN_BLOCK_LEN && xlen <= MAX_BLOCK_LEN && ylen >= MIN_BLOCK_LEN && ylen <= MAX_BLOCK_LEN) {
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
	private static void removeRepeats(Vector<int[]> blocks){
		Collections.sort(blocks, BLOCK_COMP);
		int[] block1 = null;
		int[] block2 = null;
		int i = 0;
		while (i < blocks.size()-1 && blocks.size() > 1){
			block1 = blocks.get(i);
			block2 = blocks.get(i+1);
			// if they overlap, check to see if they overlap enough to constitute a repeat 
			if (block1[0] < block2[1] && block2[0] < block1[1]){
				double intersection = block1[1] - block2[0];
				double union = block2[1] - block1[0];
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
	private static boolean isNuc(char c){
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
	public static void loadData(String samPath, Map<String,Contig> ctgs, double[][] ranges) throws IOException{
		for (int i = 0; i < ranges.length; i++)
			System.out.println("[a5_qc] Filtering read pairs with inserts between "+
					NF.format(ranges[i][0])+"-"+NF.format(ranges[i][1]));
		Map<String,SpatialClusterer> clusterers = new HashMap<String,SpatialClusterer>();
		/* Keep track of all clusterers that each contig is associated with.*/
		Map<String,Vector<String>> ctgClusterers = new HashMap<String,Vector<String>>();
		Map<String,Integer> counts = new HashMap<String,Integer>();
		Vector<String> tmpMBs = null;

		File samFile = new File(samPath);
		FileInputStream fis = new FileInputStream(samFile);
		long start = fis.getChannel().position();
		long len = fis.getChannel().size() - start;
		BufferedReader br = new BufferedReader (new InputStreamReader(fis));
		
		
		/* begin: build a lookup table for each contig for tallying coverage in windows */
		
		int contigLen = 0;
		// require the window length to be at least 1000 base pairs
		int windowLen = Math.max(1000, MEAN_BLOCK_LEN);
		String[] hdr = null;
		String contigName = null;
		//Map<String,Integer> coordOffset = new HashMap<String,Integer>();
		Map<String,int[][]> readCounts = new HashMap<String,int[][]>();
		int offset = 0;
		int[][] tmpWin = null;
		while(nextCharIs(br, '@')){
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
						tmpWin = new int[3][numWindow];
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
		System.out.print("[a5_qc] Reading SAM file...");
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
			line1[0] = trimPairNumber(line1[0]);
			line2[0] = trimPairNumber(line2[0]);
			while (!line1[0].equals(line2[0]) && br.ready()){
				line1 = line2;
				line2 = br.readLine().split("\t");				
				line2[0] = trimPairNumber(line2[0]);
			}
			total++;
			left1 = Integer.parseInt(line1[3]);
			left2 = Integer.parseInt(line2[3]);
			// if pair didn't map, just jump to next
			// line instead of jumping a full step
			if (left1 == 0 || left2 == 0) continue;
			
			ctg1 = ctgs.get(line1[2]);
			ctg2 = ctgs.get(line2[2]);

			int tmpLen = cigarLength(line1[5]);
			if (tmpLen > rdLen)
				rdLen = tmpLen;
			tmpLen = cigarLength(line2[5]);
			if (tmpLen > rdLen)
				rdLen = tmpLen;
			
			// Remove repetitive reads
			if (line1[11].equals("XT:A:R")||line2[11].equals("XT:A:R")){
				continue;
			}
			
			rev1 = isReverse(line1[1]);
			rev2 = isReverse(line2[1]);
			
			/* begin: tally these read positions */
			// tally read 1
			tmpWin = readCounts.get(ctg1.name);
			index = Arrays.binarySearch(tmpWin[0], left1);
			if (index < 0)
				index = -1*(index+1);
			tmpWin[2][index]++;
			// tally read 2
			tmpWin = readCounts.get(ctg2.name);
			index = Arrays.binarySearch(tmpWin[0], left1);
			if (index < 0)
				index = -1*(index+1);
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
				int ins = (left2 > left1 ? left2+cigarLength(line2[5])-left1 : left1+cigarLength(line1[5])-left2);
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
		System.out.println("..100%... done!... Took "+(after-before)/1000+" seconds.");
		perc = (double) numKeep / total * 100;
		System.out.println("[a5_qc] Keeping "+NF.format(perc)+"% ("+numKeep+"/"+total+") of reads.");
		ReadCluster.RDLEN = rdLen;
		/*
		 * Set LAMDBA, our Poisson rate parameter. We will use this to 
		 * compute key runtime parameters
		 */
		P = Double.POSITIVE_INFINITY;
		int minWindow = -1;
		File covFile = new File (samFile.getParentFile(),basename(samFile.getName(),".sam")+".cov");
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
				if (tmpWin[2][i] == 0)
					continue;
				tmpFreq = tmpWin[2][i];
				tmpFreq = tmpFreq/tmpWin[1][i];
				if (P > tmpFreq){
					P = tmpFreq;
					minWindow = i;
					minWinCtg = tmpCtg;
				}
				covOut.println(tmpWin[0][i]+"\t"+tmpWin[1][i]+"\t"+tmpWin[2][i]);
			}
		}
		covOut.close();
		System.out.println("[a5_load_data] Window with fewest mapped reads: "+minWinCtg+" "+readCounts.get(minWinCtg)[0][minWindow-1]+" - "+readCounts.get(minWinCtg)[0][minWindow]);
		
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
		System.out.println("        P                   = " + P);
		System.out.println("        MIN_BLOCK_LEN       = " + MIN_BLOCK_LEN);
		System.out.println("        MEAN_BLOCK_LEN      = " + MEAN_BLOCK_LEN);
		System.out.println("        MAX_BLOCK_LEN       = " + MAX_BLOCK_LEN);
		System.out.println("        MAX_INTERBLOCK_DIST = " + MAX_INTERBLOCK_DIST);
		System.out.println("        MAX_INTERPOINT_DIST = " + MAX_INTERPOINT_DIST);
		System.out.println("        EPSILON             = " + SpatialClusterer.EPS);
		System.out.println("        MIN_POINTS          = " + SpatialClusterer.MIN_PTS);
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
			System.out.print("[a5_qc] cluster"+NF.format(clusters[i].getId())+": mu="+pad(NF.format(clusters[i].mean()),10)+
					"sd="+pad(NF.format(clusters[i].sd()),10)+"n="+pad(NF.format(clusters[i].size()),10));
			NF.setMaximumFractionDigits(2);
			double perc = 100*((double)clusters[i].size())/toFilt.size();
			System.out.print("perc="+pad(NF.format(perc),10));
			
				
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
		while(nextCharIs(raf,'@')){
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
			boolean rev1 = isReverse(line1[1]);
			boolean rev2 = isReverse(line2[1]);
			
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
			tmp.addRead(left1, rev1, cigarLength(line1[5]), 
					contigs.get(line1[2]), Integer.parseInt(line1[4]), line1[5]);
			tmp.addRead(left2, rev2, cigarLength(line2[5]), 
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
	/**
	 * Return true of the next character in this RandomAccessFile is equal 
	 * to the given character <code>c</code>
	 * @param raf the RandomAccessFile to check
	 * @param c the character value to check against
	 * @return true if the next character in <code>raf</code> is equal to <code>c</code>
	 * @throws IOException if an I/O error occurs while trying to read from the given RandomAccessFIle
	 */
	private static boolean nextCharIs(RandomAccessFile raf, char c) throws IOException{
		char next = (char) raf.read();
		raf.seek(raf.getFilePointer()-1);
		return next == c;
	}
	
	private static boolean nextCharIs(BufferedReader br, char c) throws IOException{
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
	private static String pad(String s, int len){
		String ret = new String(s);
		for (int i = 0; i < len-s.length(); i++){
			ret = ret+" ";
		}
		return ret;
	}

	/**
	 * Return mapping length for CIGAR String if there is a match (i.e. string contains 'M') else return -1
	 * @param cig the CIGAR string to parse
	 * @return the length of the match indicated by the CIGAR string. Return -1 if no match
	 */
	private static int cigarLength(String cig){
		StringTokenizer tok = new StringTokenizer(cig, "MIDNSHP", true);
		int alignLen = 0;
		while (tok.hasMoreTokens()){
			int len = Integer.parseInt(tok.nextToken());
			char op = tok.nextToken().charAt(0);
			if (op == 'M'){
				alignLen += len;
			}
		}
		return alignLen;
	}

	/**
	 * Returns true of the 4th bit is set in the flag given in the SAM file
	 */
	private static boolean isReverse(String flag){
		int iflag = Integer.parseInt(flag);
		if (getBit(4,iflag) == 1) return true;
			else return false;
	}
	
	/**
	 * get the bit from the flag given in a SAM file
	 */
	private static int getBit (int bit, int flag) { 
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
	private static String trimPairNumber(String s){
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
	
}
