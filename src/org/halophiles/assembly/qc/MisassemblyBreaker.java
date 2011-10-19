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
	
	private static Comparator<int[]> BLOCK_COMP = new Comparator<int[]>(){
		public int compare(int[] arg0, int[] arg1) {
			return arg0[1] - arg1[1];
		}
	};  
	
	private static final int HANDFUL_SIZE = 100;
	private static final int SAM_LINE_LEN = 215;
	private static final String TAG_KEEP = "XT:A:U";
	private static final int MIN_PTS = 3;
	
	
	private static NumberFormat NF;
	private static int N_ESTREADPAIRS = 100000;
	private static boolean INWARD = false;
	private static boolean OUTWARD = false;
	private static double NIN = 0;
	private static double NOUT = 0;
	
	/** error rate for estimation of maximum values */
	private static double ALPHA = 0.001;
	
	private static double LAMBDA;
	public static int MAX_INTERPOINT_DIST;
	private static int MAX_INTERBLOCK_DIST;
	static int MEAN_BLOCK_LEN;
	static int MIN_BLOCK_LEN;
	static int MAX_BLOCK_LEN;
	
	private static Collection<PointChainer> matches;
	
//	private static Map<String,int[]> points;
	
	public static void main(String[] args){
		if (args.length != 5 && args.length != 4){
			System.err.println("Usage: java -jar A5qc.jar <sam_file> <contig_file> <output_file> <num_libs>");
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
			System.out.println("[a5_qc] Found "+contigs.size()+" contigs");
			System.out.println("[a5_qc] Reading in a subset of reads for insert size estimation.");
			long before = System.currentTimeMillis();
			Map<String,ReadPair> reads = readSubsetByChunk(raf, contigs);
			raf.close(); 
			long after = System.currentTimeMillis();
			System.out.println("[a5_qc] Took "+((after-before)/1000)+" seconds to read in "+reads.size()+" read pairs.");
			if (reads.size() <= 0) {
				System.err.println("[a5_qc] No paired reads found. Cannot generate match files for running FISH misassembly detection.");
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
			double[][] ranges = new double[clusterStats.length][2];
			for (int i = 0; i < clusterStats.length; i++){
				if (clusterStats[i][0]>1000){
					ranges[i][0] = clusterStats[i][0]-clusterStats[i][3]*clusterStats[i][1];
					ranges[i][1] = clusterStats[i][0]+clusterStats[i][3]*clusterStats[i][1];
				} else {
					ranges[i][0] = 1;
					ranges[i][1] = clusterStats[i][0]*2;
				}
			}
			setMAXBLOCKLEN(clusterStats);
			loadData(samPath, contigs, ranges);
			setParameters(clusterStats);
			printParams();
			
			// collect all of our blocks for each contig
			Iterator<PointChainer> mbIt = matches.iterator();
			Map<String,Vector<int[]>> blocks = new HashMap<String,Vector<int[]>>();
			Vector<int[]> xBlocks = null;
			Vector<int[]> yBlocks = null;
			while(mbIt.hasNext()){
				PointChainer pc = mbIt.next();
				// get Vector for holding contig X blocks
				xBlocks = new Vector<int[]>();
				// get Vector for holding contig Y blocks
				yBlocks = new Vector<int[]>();
				pc.buildKClumps();
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
			while(it.hasNext())
				removeRepeats(blocks.get(it.next()));
			it = blocks.keySet().iterator();
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
					System.out.println("        "+tmp[0]+" - "+tmp[1]);
					
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
			br.read();
			int[][] tmpAr = null;
			StringBuilder sb = null;
			while(br.ready()){
				String tmpCtg = br.readLine(); 
				sb = new StringBuilder();
				char c = (char) br.read();
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
					if (tmpAr[i-1][1] > tmpAr[i][0]){ // if they overlap, split at the midpoint of overlap
						right = (tmpAr[i-1][1] + tmpAr[i][0])/2;
						System.out.println("[a5_qc] Exporting "+tmpCtg+" at "+left+"-"+right);
						out.export(tmpCtg, sb, left, right);
						left = right+1;
					} else if (tmpAr[i][0]-tmpAr[i-1][1] < MAX_INTERBLOCK_DIST) {
						
						right = tmpAr[i-1][1];
						System.out.println("[a5_qc] Exporting "+tmpCtg+" at "+left+"-"+right);
						out.export(tmpCtg, sb, left, right);
						left = tmpAr[i][0];
					} 
				}
				right = sb.length();
				System.out.println("[a5_qc] Exporting "+tmpCtg+" at "+left+"-"+right);
				out.export(tmpCtg, sb, left, right);
			}
		} catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		} catch(Exception e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	private static void addBlocks(PointChainer pc, Vector<int[]> xBlocks, Vector<int[]> yBlocks) throws IOException {
		KClump[] kclumps = pc.getKClumps();
		int xlen = 0;
		int ylen = 0;
		int[] x = null;
		int[] y = null;
		for (int i = 0; i < kclumps.length; i++){
			xlen = kclumps[i].xMax-kclumps[i].xMin;
			ylen = kclumps[i].yMax-kclumps[i].yMin;
			double xden = ((double)kclumps[i].size())/xlen;
			double yden = ((double)kclumps[i].size())/ylen;
			
			x = new int[2];
			x[0] = kclumps[i].xMin;
			x[1] = kclumps[i].xMax;
			y = new int[2];
			if (kclumps[i].yMin < 0){
				y[0] = Math.abs(kclumps[i].yMax);
				y[1] = Math.abs(kclumps[i].yMin);
			} else {
				y[0] = kclumps[i].yMin;
				y[1] = kclumps[i].yMax;
			}
			if (xlen >= MIN_BLOCK_LEN && xlen <= MAX_BLOCK_LEN && ylen >= MIN_BLOCK_LEN && ylen <= MAX_BLOCK_LEN) {
				xBlocks.add(x);
				yBlocks.add(y);
			} else {
				if (xden >= LAMBDA)
					System.out.print("");
				if (yden >= LAMBDA)
					System.out.print("");
			}
		}
	}
	
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
	}
	
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
		Map<String,PointChainer> chainers = new HashMap<String,PointChainer>();
		Map<String,Vector<String>> ctgMBs = new HashMap<String,Vector<String>>();
		Map<String,Integer> counts = new HashMap<String,Integer>();
		Vector<String> tmpMBs = null;

		File samFile = new File(samPath);
		FileInputStream fis = new FileInputStream(samFile);
		long start = fis.getChannel().position();
		long len = fis.getChannel().size() - start;
		BufferedReader br = new BufferedReader (new InputStreamReader(fis));
		
		
		/* begin: build a lookup table for tallying coverage in windows */
		int genomeLen = 0;
		String[] hdr = null;
		String contigName = null;
		Map<String,Integer> coordOffset = new HashMap<String,Integer>();
		int offset = 0;
		while(nextCharIs(br, '@')){
			hdr = br.readLine().split("\t");
			offset = genomeLen;
			if (hdr[0].equals("@SQ")){
				for (int i = 1; i < hdr.length; i++){
					if (hdr[i].startsWith("LN:")){
						genomeLen += Integer.parseInt(hdr[i].substring(hdr[i].indexOf("LN:")+3));
					} else if (hdr[i].startsWith("SN:")){
						contigName = hdr[i].substring(hdr[i].indexOf("SN:")+3);
					}
				}
				coordOffset.put(contigName, offset);
			}
		}
		// require the window length to be at least 1000 base pairs
		int windowLen = Math.max(1000, MEAN_BLOCK_LEN);
		int[][] readCounts = null;
		int numWindow = genomeLen/windowLen;
		if (genomeLen % windowLen != 0) {
			readCounts = new int[2][numWindow+1];
			readCounts[0][numWindow] = genomeLen;
		} else {
			readCounts = new int[2][numWindow];
		}
		for (int i = 0; i < numWindow; i++){
			readCounts[0][i] = windowLen*(i+1);
			readCounts[1][i] = 0;
		}
		/* end: build a lookup table for tallying coverage in windows */
		
		
		String[] line1 = null;
		String[] line2 = null;
		int left1 = 0;
		int left2 = 0;
		String ctgStr = null;
		String tmp = null;
		PointChainer pc = null;
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
			
			/* begin: tally these read positions */
			offset = coordOffset.get(ctg1.name);
			index = Arrays.binarySearch(readCounts[0], offset+left1);
			if (index < 0)
				index = -1*(index+1);
			readCounts[1][index]++;
			offset = coordOffset.get(ctg2.name);
			index = Arrays.binarySearch(readCounts[0], offset+left2);
			if (index < 0)
				index = -1*(index+1);
			readCounts[1][index]++;
			/* end: tally these read positions */
			
			
			ctgNameComp = line1[2].compareTo(line2[2]);
			
			if (line1[11].equals("XT:A:R")||line2[11].equals("XT:A:R")){
				continue;
			}
			
			/*
			 * if this pair spans two different contigs, we need to sort the names
			 * for consistency before we add the match
			 */
			if (ctgNameComp != 0){
				if (ctgNameComp < 0){
					ctgStr = line1[2]+"-"+line2[2];
					if (chainers.containsKey(ctgStr))
						pc = chainers.get(ctgStr);
					else {
						pc = new PointChainer(ctg1, ctg2);
						chainers.put(ctgStr, pc);
					}
					pc.addMatch(left1, left2);
				} else {
					ctgStr = line2[2]+"-"+line1[2];
					if (chainers.containsKey(ctgStr))
						pc = chainers.get(ctgStr);
					else { 
						pc = new PointChainer(ctg2, ctg1);
						chainers.put(ctgStr, pc);
					}
					pc.addMatch(left2, left1);
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
				if (chainers.containsKey(ctgStr))
					pc = chainers.get(ctgStr);
				else { 
					pc = new PointChainer(ctg2, ctg1);
					chainers.put(ctgStr, pc);
				}
				if (left2 < left1) // order for consistency
					pc.addMatch(left2, left1);
				else 
					pc.addMatch(left1, left2);
				if (counts.containsKey(ctg1.name)){
					counts.put(ctg1.name, counts.get(ctg1.name)+2);
				} else {
					counts.put(ctg1.name, 2);
				}
				
			}
			// add point for contig1, and keep track of which MatchBuilders are associated with contig1
			if (ctgMBs.containsKey(ctg1.name)){
				tmpMBs = ctgMBs.get(ctg1.name);
			} else {
				tmpMBs = new Vector<String>();
				ctgMBs.put(ctg1.name, tmpMBs);
			}
			tmpMBs.add(ctgStr);
			// add point for contig2, and keep track of which MatchBuilders are associated with contig2
			if (ctgMBs.containsKey(ctg2.name)){
				tmpMBs = ctgMBs.get(ctg2.name);
			} else {
				tmpMBs = new Vector<String>();
				ctgMBs.put(ctg2.name, tmpMBs);
			}
			tmpMBs.add(ctgStr);
			numKeep++;
		}
		long after = System.currentTimeMillis();
		System.out.println("..100%... done!... Took "+(after-before)/1000+" seconds.");
		perc = (double) numKeep / total * 100;
		System.out.println("[a5_qc] Keeping "+NF.format(perc)+"% ("+numKeep+"/"+total+") of reads.");
		
		/*
		 * Set LAMDBA, our Poisson rate parameter. We will use this to 
		 * compute key runtime parameters
		 */
		LAMBDA = Double.POSITIVE_INFINITY;
		for (int i = 0; i < readCounts[1].length; i++){
			if (readCounts[1][i] == 0)
				continue;
			if (LAMBDA > readCounts[1][i])
				LAMBDA = readCounts[1][i];
		}
		LAMBDA = LAMBDA/MEAN_BLOCK_LEN;
		
		Iterator<String> ctgIt = ctgMBs.keySet().iterator();
		Set<String> ctgToRm = new HashSet<String>();
		Set<String> psPairsToRm = new HashSet<String>();
		while(ctgIt.hasNext()){
			tmp = ctgIt.next();
			if (counts.get(tmp) < MIN_PTS) {
				ctgToRm.add(tmp);
				psPairsToRm.addAll(ctgMBs.get(tmp));
			} 
		}
		removeKeys(ctgs, ctgToRm);
		removeKeys(chainers, psPairsToRm);
	
		matches = chainers.values();
	}
	
	/**
	 * Set parameters and print their values to standard out
	 * @param ranges <code>[ mean, sd, n, nSd , min, max ]</code>
	 */
	private static void setParameters(double[][] ranges){
		for (int i = 0; i < ranges.length; i++)
//			if (ranges[i][1]*ranges[i][3] > PointChainer.MAX_RES)
//				PointChainer.MAX_RES = (int) (ranges[i][1]*ranges[i][3]);
			if (ranges[i][1] > PointChainer.MAX_RES)
				PointChainer.MAX_RES = (int) (ranges[i][1]);
		setMAXINTERPOINTDIST();
		setMAXINTERBLOCKDIST();
	}
	
	private static void setMAXBLOCKLEN(double[][] ranges){
		for (double[] cluster: ranges){
			if (cluster[0] > MEAN_BLOCK_LEN){
				MEAN_BLOCK_LEN = (int) (cluster[0]);
				MAX_BLOCK_LEN = (int) (cluster[0]+cluster[1]*cluster[3]);
			}
		}
	}
	
	private static void setMAXINTERPOINTDIST(){
		/*
		 * By modelling read mapping points as a Poisson process, the distance
		 * between points is exponential. Take the 1-ALPHA quantile to get
		 * a maximum distance between two points.
		 *
		 * LAMBDA Rate of mapping points (Poisson rate parameter)
		 */	
//		MAX_INTERPOINT_DIST = (int) (Math.log(ALPHA)/(-LAMBDA));
		MAX_INTERPOINT_DIST = Math.max(2*PointChainer.MAX_RES,(int) (Math.log(ALPHA)/(-LAMBDA)));
		PointChainer.EPS = MAX_INTERPOINT_DIST;
		MIN_BLOCK_LEN = 2*MAX_INTERPOINT_DIST;
		//MAX_INTERPOINT_DIST = 600;
	}
	
	private static void setMAXINTERBLOCKDIST(){
		/*
		 * apply some extreme value theory to get the minimum of points 
		 * randomly sampled uniformally across an interval of MAX_BLOCK_LEN
		 * 
		 * MAX_INTERBLOCK_DIST = 2*((1-ALPHA)^(1/EXP_POINTS) * MAX_BLOCK_LEN)
		 * 
		 */
		MAX_INTERBLOCK_DIST = (int)(2*(Math.pow(1-ALPHA, 1/(LAMBDA*MEAN_BLOCK_LEN))*MEAN_BLOCK_LEN));
	}
	
	private static void printParams(){
		System.out.println("[a5_qc] parameters:");
		System.out.println("        LAMBDA              = " + LAMBDA);
		System.out.println("        MIN_BLOCK_LEN       = " + MIN_BLOCK_LEN);
		System.out.println("        MEAN_BLOCK_LEN      = " + MEAN_BLOCK_LEN);
		System.out.println("        MAX_BLOCK_LEN       = " + MAX_BLOCK_LEN);
		System.out.println("        MAX_INTERBLOCK_DIST = " + MAX_INTERBLOCK_DIST);
		System.out.println("        MAX_INTERPOINT_DIST = " + MAX_INTERPOINT_DIST);
		System.out.println("        MAX_RESIDUAL        = " + PointChainer.MAX_RES);
		System.out.println("        EPSILON             = " + PointChainer.EPS);
	}
	
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
		int K = numLibs+1;
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
		System.out.println("[a5_qc] Initial stats for sample: mu="+NF.format(ins[0])+" sd="+NF.format(ins[1])+" n="+NF.format(ins[2]));
		System.out.print("[a5_qc] EM-clustering insert sizes with K="+K+"... ");
		/* begin: cluster read pairs by insert size so we can efficiently determine the insert size of this library */
		EMClusterer em = new EMClusterer(toFilt, K);
		long before = System.currentTimeMillis();
		double delta = 0.00005;
		int iters = em.iterate(1000, delta);
		long after = System.currentTimeMillis();
		System.out.println("stopping after "+iters+" iterations with delta="+delta+". Took "+(after-before)/1000+" seconds.");
		ReadSet[] clusters = new ReadSet[K];
		em.getClusters().toArray(clusters);
		/* end: cluster read pairs  */
		
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
		PointChainer.MAX_RES = 0;
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
			left2 = Integer.parseInt(line2[3]);
			// if pair didn't map, just jump to next line instead of jumping a full step
			if (left1 == 0 || left2 == 0) continue;
			tmp = new ReadPair(line1[0]);
			boolean rev1 = isReverse(line1[1]);
			boolean rev2 = isReverse(line2[1]);
			
			if (contigs.get(line1[2]).equals(contigs.get(line2[2])) && rev1 != rev2 &&
					// require this mapping to be unique. 
					line1[11].equals(TAG_KEEP) && line2[11].equals(TAG_KEEP)){
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
	 * A class for reading a SAM file via random access. 
	 * @author andrew
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
			while (!line1[0].equals(line2[0]) && tok.hasMoreTokens()){
				line1 = line2;
				line2 = tok.nextToken().split("\t");
				tokLeft--;
			}
			String[][] ret = {line1,line2};
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
			tokLeft = tok.countTokens();
		}
	}
	
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
		int totalLen = 0;
		int alignLen = 0;
		while (tok.hasMoreTokens()){
			int len = Integer.parseInt(tok.nextToken());
			char op = tok.nextToken().charAt(0);
			if (op == 'M'){
				alignLen += len;
			}
			if (op != 'I')
				totalLen += len;
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
	
	private static String trimPairNumber(String s){
		if (s.contains("/")){
			return s.substring(0,s.indexOf("/"));
		} else
			return s;
	}
	
}
