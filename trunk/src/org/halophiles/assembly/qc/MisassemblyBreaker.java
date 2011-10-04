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
import java.util.TreeSet;
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
	private static int MEAN_BLOCK_LEN;
	private static int MIN_BLOCK_LEN;
	private static int MAX_BLOCK_LEN;
	
	private static Collection<MatchBuilder> matches;
	
	private static Map<String,int[]> points;
	
	public static void main(String[] args){
		if (args.length != 3 && args.length != 4){
			System.err.println("Usage: java -jar A5qc.jar <sam_file> <contig_file> <output_file> <num_libs>");
			System.exit(-1);
		}
		try{
			
			NF = NumberFormat.getInstance();
			NF.setMaximumFractionDigits(0);
			NF.setGroupingUsed(false);
			String samPath = args[0];
			String base = args[2];
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
			
				
			double[][] clusterStats = getRangesToFilter(reads, numLibs);
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
			loadData(samPath, base, contigs, ranges);
			setParameters(clusterStats);
			
			// collect all of our blocks for each contig
			Iterator<MatchBuilder> mbIt = matches.iterator();
			Map<String,Vector<int[]>> blocks = new HashMap<String,Vector<int[]>>();
			Vector<int[]> xBlocks = null;
			Vector<int[]> yBlocks = null;
			while(mbIt.hasNext()){
				MatchBuilder mb = mbIt.next();
				//mb.print(new File(samFile.getParentFile(),"match."+mb.getContig1().getId()+"v"+mb.getContig2().getId()+".txt"));
				// get Vector for holding contig X blocks
				xBlocks = new Vector<int[]>();
				// get Vector for holding contig Y blocks
				yBlocks = new Vector<int[]>();
				addBlocks(mb, xBlocks, yBlocks ,samFile.getParentFile(),"clump."+mb.getContig1().getId()+"v"+mb.getContig2().getId());
				mergeAndFilterBlocks(mb.getContig1(), xBlocks);
				mergeAndFilterBlocks(mb.getContig2(), yBlocks);
				if (blocks.containsKey(mb.getContig1().name))
					blocks.get(mb.getContig1().name).addAll(xBlocks);
				else
					blocks.put(mb.getContig1().name, xBlocks);
				
				if (blocks.containsKey(mb.getContig2().name))
					blocks.get(mb.getContig2().name).addAll(yBlocks);
				else
					blocks.put(mb.getContig2().name, yBlocks);
				
			}
			
			Iterator<String> it = blocks.keySet().iterator();
			while(it.hasNext())
				removeRepeats(blocks.get(it.next()));
			it = blocks.keySet().iterator();
			Vector<String> toRm = new Vector<String>();
			while(it.hasNext()){
				String tmpCtg = it.next();
				Vector<int[]> tmpBlocks = blocks.get(tmpCtg);
				if (tmpBlocks.isEmpty())
					toRm.add(tmpCtg);
				else 
					Collections.sort(tmpBlocks, BLOCK_COMP);
				
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
				/*System.out.print(tmpCtg+"\t"+tmpBlks.get(0)[0]+"-"+tmpBlks.get(0)[1]);
				for (int i = 1; i < tmpBlks.size(); i++){
					System.out.print(", "+tmpBlks.get(i)[0]+"-"+tmpBlks.get(i)[1]);
				}
				System.out.println(); */
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
						/*if (tmpAr[i-1][1] - left + 1 < MIN_CONTIG_LENGTH){
							left = tmpAr[i][0];
							continue;
						}*/
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
	
	private static void addBlocks(MatchBuilder mb, Vector<int[]> xBlocks, Vector<int[]> yBlocks, File outDir, String pathBase) throws IOException {
		int[] p1 = points.get(mb.getContig1().name);
		int[] p2 = points.get(mb.getContig2().name);
		int[][] matches = mb.getMatches();
		PointChainer pc = new PointChainer(p1, p2, matches);
		KClump[] kclumps = pc.getKClumps();
		int xlen = 0;
		int ylen = 0;
		int[] x = null;
		int[] y = null;
		for (int i = 0; i < kclumps.length; i++){
			xlen = kclumps[i].xMax-kclumps[i].xMin;
			ylen = kclumps[i].yMax-kclumps[i].yMin;
			if (xlen >= MIN_BLOCK_LEN && xlen <= MAX_BLOCK_LEN && ylen >= MIN_BLOCK_LEN && ylen <= MAX_BLOCK_LEN) {
				x=  new int[2];
				x[0] = kclumps[i].xMin;
				x[1] = kclumps[i].xMax;
				y = new int[2];
				y[0] = kclumps[i].yMin;
				y[1] = kclumps[i].yMax;
				xBlocks.add(x);
				yBlocks.add(y);
//				kclumps[i].print(new File(outDir,pathBase+"."+kclumps[i].id+".txt"));
//				System.out.println(kclumps[i].id+"  "+mb.getContig1().name+" "+x[0]+"-"+x[1]+"\t"+
//						mb.getContig2().name+" "+y[0]+"-"+y[1]+" "+kclumps[i].size()+" "+kclumps[i].density()+"     1");
			}
		}

		/*
		 * Now do the same for the the reverse. 
		 */
		int[] p2Inv = new int[p2.length];
		for (int i = 0; i < p2Inv.length; i++){
			p2Inv[i] = -1*p2[p2.length-i-1];
		}
		int[][] matchesInv = new int[matches.length][3];
		for (int i = 0; i < matchesInv.length; i++) {
			matchesInv[i][0] = matches[i][0];
			matchesInv[i][1] = -1*matches[i][1];
		}
		pc = new PointChainer(p1, p2Inv, matchesInv);
		KClump[] kclumpsInv = pc.getKClumps();
		for (int i = 0; i < kclumpsInv.length; i++){
			xlen = kclumpsInv[i].xMax-kclumpsInv[i].xMin;
			ylen = kclumpsInv[i].yMax-kclumpsInv[i].yMin;
			if (xlen >= MIN_BLOCK_LEN && xlen <= MAX_BLOCK_LEN && ylen >= MIN_BLOCK_LEN && ylen <= MAX_BLOCK_LEN) {
				x = new int[2];
				x[0] = kclumpsInv[i].xMin;
				x[1] = kclumpsInv[i].xMax;
				y = new int[2];
				y[0] = -1*kclumpsInv[i].yMax;
				y[1] = -1*kclumpsInv[i].yMin;
				xBlocks.add(x);
				yBlocks.add(y);
//				kclumpsInv[i].print(new File(outDir,pathBase+"."+kclumpsInv[i].id+".txt"));
//				System.out.println(kclumpsInv[i].id+"  "+mb.getContig1().name+" "+x[0]+"-"+x[1]+"\t"+
//							mb.getContig2().name+" "+y[0]+"-"+y[1]+" "+kclumpsInv[i].size()+" "+kclumpsInv[i].density()+"     -1");
			}
		}
	}
	
	private static void mergeAndFilterBlocks(Contig contig, Vector<int[]> blocks){
		int i = 0;
		int len = 0;
		while (i < blocks.size()){
			len = blocks.get(i)[1] - blocks.get(i)[0];
			if (len < MIN_BLOCK_LEN || len > MAX_BLOCK_LEN)
				blocks.remove(i);
			else
				i++;
		}
		Collections.sort(blocks, BLOCK_COMP);
		int[] block1 = null;
		int[] block2 = null;
		i = 1;
		while (i < blocks.size()){
			block1 = blocks.get(i-1);
			block2 = blocks.get(i);
			if (block1[0] < block2[1] && block2[0] < block1[1]){
				block1[0] = Math.min(block1[0], block2[0]);
				block1[1] = Math.max(block1[1], block2[1]);
				blocks.remove(i);
			} else {
				i++;
			}
		}
		
		// remove terminal blocks, we won't use these to break on
		i = 0;
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
					blocks.remove(i+1);
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
	 * 
	 * @param samPath
	 * @param fishDir
	 * @param base
	 * @param ctgs
	 * @param ranges an array of arrays of max and min insert sizes
	 * @throws IOException
	 */
	public static void loadData(String samPath, String base, Map<String,Contig> ctgs, double[][] ranges) throws IOException{
		for (int i = 0; i < ranges.length; i++)
			System.out.println("[a5_qc] Filtering read pairs with inserts between "+
					NF.format(ranges[i][0])+"-"+NF.format(ranges[i][1]));
		Map<String,TreeSet<Integer>> mapPoints = new HashMap<String,TreeSet<Integer>>();
		Map<String,MatchBuilder> matchBldrs = new HashMap<String,MatchBuilder>();
		Map<String,Vector<String>> ctgMBs = new HashMap<String,Vector<String>>();
		TreeSet<Integer> tmpPts = null;
		Vector<String> tmpMBs = null;

		File samFile = new File(samPath);
		FileInputStream fis = new FileInputStream(samFile);
		long start = fis.getChannel().position();
		long len = fis.getChannel().size() - start;
		BufferedReader br = new BufferedReader (new InputStreamReader(fis));
		
		
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
		
		
		String[] line1 = null;
		String[] line2 = null;
		int left1 = 0;
		int left2 = 0;
		String ctgStr = null;
		String tmp = null;
		MatchBuilder mb = null;
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
			while (!line1[0].equals(line2[0]) && br.ready()){
				line1 = line2;
				line2 = br.readLine().split("\t");				
			}
			total++;
			left1 = Integer.parseInt(line1[3]);
			left2 = Integer.parseInt(line2[3]);
			// if pair didn't map, just jump to next
			// line instead of jumping a full step
			if (left1 == 0 || left2 == 0) continue;
			
			ctg1 = ctgs.get(line1[2]);
			ctg2 = ctgs.get(line2[2]);
			
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
			
			
			ctgNameComp = line1[2].compareTo(line2[2]); // order for consistency
			
			if (ctgNameComp < 0){
				ctgStr = line1[2]+"-"+line2[2];
				if (matchBldrs.containsKey(ctgStr))
					mb = matchBldrs.get(ctgStr);
				else {
					mb = new MatchBuilder(ctg1, ctg2);
					matchBldrs.put(ctgStr, mb);
				}
				mb.addMatch(left1, left2);
			} else if (ctgNameComp > 0) {
				ctgStr = line2[2]+"-"+line1[2];
				if (matchBldrs.containsKey(ctgStr))
					mb = matchBldrs.get(ctgStr);
				else { 
					mb = new MatchBuilder(ctg2, ctg1);
					matchBldrs.put(ctgStr, mb);
				}
				mb.addMatch(left2, left1);
			} else { // same contig, so check to see it's within the given ranges
				int ins = (left2 > left1 ? left2+cigarLength(line2[5])-left1 : left1+cigarLength(line1[5])-left2);
				if (inRange(ranges,ins)) 
						continue;
				
				ctgStr = line2[2]+"-"+line1[2];
				if (matchBldrs.containsKey(ctgStr))
					mb = matchBldrs.get(ctgStr);
				else { 
					mb = new MatchBuilder(ctg2, ctg1);
					matchBldrs.put(ctgStr, mb);
				}
				if (left2 < left1) // order for consistency
					mb.addMatch(left2, left1);
				else 
					mb.addMatch(left1, left2);				
			}
			
			// add point for contig1, and keep track of which MatchBuilders are associated with contig1
			if (mapPoints.containsKey(ctg1.name)){
				tmpPts = mapPoints.get(ctg1.name);
				tmpMBs = ctgMBs.get(ctg1.name);
			} else {
				tmpPts = new TreeSet<Integer>();
				mapPoints.put(ctg1.name, tmpPts);
				tmpMBs = new Vector<String>();
				ctgMBs.put(ctg1.name, tmpMBs);
			}
			tmpPts.add(new Integer(left1));
			tmpMBs.add(ctgStr);
			// add point for contig2, and keep track of which MatchBuilders are associated with contig2
			if (mapPoints.containsKey(ctg2.name)){
				tmpPts = mapPoints.get(ctg2.name);
				tmpMBs = ctgMBs.get(ctg2.name);
			} else {
				tmpPts = new TreeSet<Integer>();
				mapPoints.put(ctg2.name, tmpPts);
				tmpMBs = new Vector<String>();
				ctgMBs.put(ctg2.name, tmpMBs);
			}
			tmpMBs.add(ctgStr);
			tmpPts.add(new Integer(left2));
			numKeep++;
		}
		long after = System.currentTimeMillis();
		System.out.println("..100%... done!... Took "+(after-before)/1000+" seconds.");
		perc = (double) numKeep / total * 100;
		System.out.println("[a5_qc] Keeping "+NF.format(perc)+"% ("+numKeep+"/"+total+") of reads.");
		
		LAMBDA = Double.POSITIVE_INFINITY;
		File covFile = new File(samFile.getParentFile(),"coverage.txt");
		covFile.createNewFile();
		PrintStream out = new PrintStream(covFile);
		for (int i = 0; i < readCounts[1].length; i++){
			if (readCounts[1][i] == 0)
				continue;
			if (LAMBDA > readCounts[1][i])
				LAMBDA = readCounts[1][i];
			out.println(readCounts[1][i]);
		}
		out.close();
		LAMBDA = LAMBDA/MEAN_BLOCK_LEN;
		
		Iterator<String> ctgIt = mapPoints.keySet().iterator();
		Set<String> ctgToRm = new HashSet<String>();
		Set<String> psPairsToRm = new HashSet<String>();
		points = new HashMap<String,int[]>();
		while(ctgIt.hasNext()){
			tmp = ctgIt.next();
			tmpPts = mapPoints.get(tmp);
			if (tmpPts.size() < MIN_PTS) {
				ctgToRm.add(tmp);
				psPairsToRm.addAll(ctgMBs.get(tmp));
			} else {
				int[] ar = new int[tmpPts.size()];
				Iterator<Integer> ptsIterator = tmpPts.iterator();
				for (int p = 0; p < ar.length; p++){
					ar[p] = ptsIterator.next();
				}
				points.put(tmp, ar);
			}
		}
		removeKeys(ctgs, ctgToRm);
		removeKeys(mapPoints, ctgToRm);
		removeKeys(matchBldrs, psPairsToRm);
	
		matches = matchBldrs.values();
	}
	
	private static void setParameters(double[][] ranges){
		setMAXINTERPOINTDIST();
		setMAXINTERBLOCKDIST();
		System.out.println("[a5_qc] parameters:");
		System.out.println("        LAMBDA              = " + LAMBDA);
		System.out.println("        MIN_BLOCK_LEN       = " + MIN_BLOCK_LEN);
		System.out.println("        MEAN_BLOCK_LEN      = " + MEAN_BLOCK_LEN);
		System.out.println("        MAX_BLOCK_LEN       = " + MAX_BLOCK_LEN);
		System.out.println("        MAX_INTERBLOCK_DIST = " + MAX_INTERBLOCK_DIST);
		System.out.println("        MAX_INTERPOINT_DIST = " + MAX_INTERPOINT_DIST);
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
		MAX_INTERPOINT_DIST = (int) (Math.log(ALPHA)/(-LAMBDA));
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
	 * @return an array or arrays, where individual arrays contain cluster stats
	 */
	private static double[][] getRangesToFilter(Map<String,ReadPair> reads, int numLibs){
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
		double delta = 0.00005;
		EMClusterer em = new EMClusterer(toFilt, K);
		long before = System.currentTimeMillis();
		int iters = em.iterate(1000, delta);
		long after = System.currentTimeMillis();
		System.out.println("stopping after "+iters+" iterations with delta="+delta+". Took "+(after-before)/1000+" seconds.");
		ReadSet[] clusters = new ReadSet[K];
		em.getClusters().toArray(clusters);
		
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
		
		Vector<ReadSet> signal = new Vector<ReadSet>();
		Vector<ReadSet> noise = new Vector<ReadSet>();
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
				if (clusters[i].sd() > PointChainer.MAX_RES)
					PointChainer.MAX_RES = (int) clusters[i].sd();
			} else {
				noise.add(clusters[i]);
				System.out.println("  (noise)");
			}
			
		}
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
		
		/*
		for (i = 0; i < ret.length; i++){
			double[] tmpAr = {1,2*ret[i][0]};
			if (tmpAr[0] > min)
				tmpAr[0] = min;
			if (tmpAr[1] < max)
				tmpAr[1] = max;
			ret[i] = tmpAr;
			
		}
		*/
		return ret;
	}
	
	private static <V> void removeKeys(Map<String,V> reads, Collection<String> torm){
		Iterator<String> it = torm.iterator();
		while(it.hasNext())
			reads.remove(it.next());
	}
	
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
	
	public static class MatchFile implements Comparable<MatchFile>{
		private Contig ctg1;
		private Contig ctg2;
		private File file;
		
		public MatchFile(Contig ctg1, Contig ctg2, File file){
			this.ctg1 = ctg1;
			this.ctg2 = ctg2;
			this.file = file;
		}
		
		public String toString(){ 
			return ctg1.getRank()+"\t"+ctg2.getRank();
		}
		
		public String getName(){
			return file.getName();
		}
		
		public int compareTo(MatchFile arg0) {
			if (this.ctg1.getRank() - arg0.ctg1.getRank() == 0){
				return this.ctg2.getRank() - arg0.ctg2.getRank();
			} else {
				return this.ctg1.getRank() - arg0.ctg1.getRank();
			}
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

	private static boolean isReverse(String flag){
		int iflag = Integer.parseInt(flag);
		if (getBit(4,iflag) == 1) return true;
			else return false;
	}

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
	
}
