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
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Vector;

import org.halophiles.assembly.Contig;
import org.halophiles.assembly.ReadPair;
import org.halophiles.assembly.ReadSet;

public class FISHInputExporter {
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
	
	public static void main(String[] args){
		if (args.length != 3 && args.length != 4){
			System.err.println("Usage: java -jar GetFishInput.jar <sam_file> <output_base> <output_dir> <num_libs>");
			System.exit(-1);
		}
		try{
			NF = NumberFormat.getInstance();
			NF.setMaximumFractionDigits(0);
			NF.setGroupingUsed(false);
			String samPath = args[0];
			String base = args[1];
			File outdir = new File(args[2]);
			if (!outdir.exists()){
				outdir.mkdirs();
			}
			int numLibs = 1;
			if (args.length == 4){
				numLibs = Integer.parseInt(args[3]);
			}
			System.out.println("[a5_fie] Reading "+samPath);
			System.out.println("[a5_fie] Writing output to "+outdir.getAbsolutePath()+" with basename "+base);
			RandomAccessFile raf = new RandomAccessFile(new File(samPath), "r");
			Map<String,Contig> contigs = readContigs(raf);
			System.out.println("[a5_fie] Found "+contigs.size()+" contigs");
			System.out.println("[a5_fie] Reading in a subset of reads for insert size estimation.");
			long before = System.currentTimeMillis();
			//Map<String,ReadPair> reads = readSubset(raf, contigs);
			Map<String,ReadPair> reads = readSubsetByChunk(raf, contigs);
			raf.close(); 
			long after = System.currentTimeMillis();
			System.out.println("[a5_fie] Took "+((after-before)/1000)+" seconds to read in "+reads.size()+" read pairs.");
			if (reads.size() <= 0) {
				System.err.println("[a5_fie] No paired reads found. Cannot generate match files for running FISH misassembly detection.");
				System.exit(-1);
			}
			setOrientation(reads.values());
			if (INWARD && !OUTWARD)
				System.out.println("[a5_fie] Found a substantial amount of innies, but found no outties.");
			else if (!INWARD && OUTWARD)
				System.out.println("[a5_fie] Found a substantial amount of outties, but found no innies.");
			else 
				System.out.println("[a5_fie] Found both innies and outties.");
			
			double[][] ranges = getRangesToFilter(reads, numLibs);
			
			export2(samPath, outdir, base, contigs, ranges);
			
		} catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		} catch(Exception e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public static void export(Map<String,ReadPair> reads, File fishDir, String base) throws IOException{
		Map<String,Vector<Double>> mapPoints = new HashMap<String,Vector<Double>>();
		Map<String,PrintStreamPair> psPairs = new TreeMap<String,PrintStreamPair>();
		Vector<Double> tmpPts = null;
		
		File controlFile = new File(fishDir,base+".control.txt");
		Map<String,Contig> contigs = new HashMap<String, Contig>();
		PrintStream ctrlOut = new PrintStream(controlFile);
		Iterator<ReadPair> rpIt = reads.values().iterator();
		ReadPair tmpRp = null;
		Contig tmpCtg = null;
		while(rpIt.hasNext()){
			tmpRp = rpIt.next();
			tmpCtg = tmpRp.ctg1;
			boolean rev = tmpRp.rev1;
			int left = tmpRp.pos1;
			if (mapPoints.containsKey(tmpCtg.name)){
				tmpPts = mapPoints.get(tmpCtg.name);
			} else {
				tmpPts = new Vector<Double>();
				mapPoints.put(tmpCtg.name, tmpPts);
			}
			if (rev){
				left = left * -1;
			}
			contigs.put(tmpCtg.name, tmpCtg);
			tmpPts.add(new Double(left));
			
			tmpCtg = tmpRp.ctg2;
			rev = tmpRp.rev2;
			left = tmpRp.pos2;
			if (mapPoints.containsKey(tmpCtg.name)){
				tmpPts = mapPoints.get(tmpCtg.name);
			} else {
				tmpPts = new Vector<Double>();
				mapPoints.put(tmpCtg.name, tmpPts);
			}
			if (rev){
				left = left * -1;
			}
			contigs.put(tmpCtg.name, tmpCtg);
			tmpPts.add(new Double(left));
		}
		Iterator<String> ctgIt = mapPoints.keySet().iterator();
		Set<String> singles = new HashSet<String>();
		String tmpStr = null;
		while(ctgIt.hasNext()){
			tmpStr = ctgIt.next();
			if (mapPoints.get(tmpStr).size() <= 2)
				singles.add(tmpStr);
		}
		removeKeys(contigs, singles);
		removeKeys(mapPoints, singles);
		
		Contig[] ctgRef = new Contig[contigs.values().size()];
		contigs.values().toArray(ctgRef);
		Arrays.sort(ctgRef,new Comparator<Contig>(){
			public int compare(Contig o1, Contig o2) {
				return o1.getId() - o2.getId();
			}
		});
		// Rank contigs so FISH doesn't freak out about non-consecutive ids
		for (int i = 0; i < ctgRef.length; i++){
			ctgRef[i].setRank(i+1);
		}
		
		File ctgLblFile = new File(fishDir,base+".contig_labels.txt");
		ctgLblFile.createNewFile();
		PrintStream ctgLblOut = new PrintStream(ctgLblFile);
		File mapMatchDir = new File(fishDir,base+".dat");
		mapMatchDir.mkdir();
		
		Iterator<String> it = mapPoints.keySet().iterator();
		TreeMap<Integer,String> maps = new TreeMap<Integer, String>();
		while(it.hasNext()){
			String ctg = it.next();
			Vector<Double> tmp = mapPoints.get(ctg);
			Double[] ar = new Double[tmp.size()];
			tmp.toArray(ar);
			File ctgMapFile = new File(mapMatchDir,"map."+contigs.get(ctg).getRank()+".txt");
			ctgLblOut.println(contigs.get(ctg).getRank()+"\t"+contigs.get(ctg).name);
			PrintStream ctgMapOut = new PrintStream(ctgMapFile);
			Arrays.sort(ar, new Comparator<Double>(){
				public int compare(Double arg0, Double arg1) {
					if (Math.abs(arg0) < Math.abs(arg1)) return -1;
					else if (Math.abs(arg0) > Math.abs(arg1)) return 1;
					else return 0;
				}	
			});
			for (Double d: ar){
				//ctgMapOut.println("c"+contigs.get(ctg).getRank()+"p"+Math.abs(d.intValue())+"\t"+(d < 0 ? -1 : 1));
				ctgMapOut.println("c"+contigs.get(ctg).getRank()+"p"+Math.abs(d.intValue())+"\t0");
			}
			
			ctgMapOut.close();
			maps.put(contigs.get(ctg).getRank(),fishDir.getName()+"/"+mapMatchDir.getName()+"/"+ctgMapFile.getName());
		}
		ctgLblOut.close();
		Iterator<Integer> it2 = maps.keySet().iterator();
		ctrlOut.println("-maps");
		while(it2.hasNext()){
			Integer tmpInt = it2.next();
			ctrlOut.println(tmpInt+"\t"+maps.get(tmpInt));
		}
		
		it = reads.keySet().iterator();
		PrintStreamPair tmpPSP = null;
		ctrlOut.println("-matches");
		String ctgStr;
		while(it.hasNext()) {
			ReadPair tmp = reads.get(it.next());
			if (tmp.ctg1 == null || tmp.ctg2 == null) 
				continue;
			ctgStr = tmp.contigString();
			if (psPairs.containsKey(ctgStr))
				tmpPSP = psPairs.get(ctgStr);
			else {
				tmpPSP = new PrintStreamPair(tmp.ctg1,tmp.ctg2,mapMatchDir);
				psPairs.put(ctgStr, tmpPSP);
			}
			tmpPSP.add(tmp);
		}
		
		it = psPairs.keySet().iterator();
		Vector<String> mfToRm = new Vector<String>();
		
		while(it.hasNext()){
			ctgStr = it.next();
			tmpPSP = psPairs.get(ctgStr);
			if (tmpPSP.nPairs <= 2)
				mfToRm.add(ctgStr);
		}
		removeKeys(psPairs, mfToRm);
		Set<MatchFile> matchFiles = new TreeSet<MatchFile>();
		it = psPairs.keySet().iterator();
		while(it.hasNext()) {
			ctgStr = it.next();
			matchFiles.addAll(psPairs.get(ctgStr).getMatchFiles());
			psPairs.get(ctgStr).close();
		}
		Iterator<MatchFile> mfIt = matchFiles.iterator();
		MatchFile tmpMF = null;
		while(mfIt.hasNext()){
			tmpMF = mfIt.next();
			ctrlOut.println(tmpMF.toString()+"\t"+fishDir.getName()+"/"+mapMatchDir.getName()+"/"+tmpMF.getName());
		}
		
		ctrlOut.println("end");
		ctrlOut.close();
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
	public static void export2(String samPath, File fishDir, String base, Map<String,Contig> ctgs, double[][] ranges) throws IOException{
		for (int i = 0; i < ranges.length; i++)
			System.out.println("[a5_fie] Filtering read pairs with inserts between "+
					NF.format(ranges[i][0])+"-"+NF.format(ranges[i][1]));
		Map<String,TreeSet<Double>> mapPoints = new HashMap<String,TreeSet<Double>>();
		Map<String,PrintStreamPair> psPairs = new HashMap<String,PrintStreamPair>();
		Map<String,Vector<String>> ctgPSPairs = new HashMap<String,Vector<String>>();
		TreeSet<Double> tmpPts = null;
		Vector<String> tmpPSPairs = null;
		Comparator<Double> posComp = new Comparator<Double>(){
			public int compare(Double arg0, Double arg1) {
				if (Math.abs(arg0) < Math.abs(arg1)) return -1;
				else if (Math.abs(arg0) > Math.abs(arg1)) return 1;
				else return 0;
			}
		};
//		BufferedReader br = new BufferedReader (new FileReader(new File(samPath)));
		FileInputStream fis = new FileInputStream(new File(samPath));
		long start = fis.getChannel().position();
		long len = fis.getChannel().size() - start;
		BufferedReader br = new BufferedReader (new InputStreamReader(fis));
		
		while(nextCharIs(br, '@')) 
			br.readLine();
		
		String[] line1 = null;
		String[] line2 = null;
		int left1 = 0;
		int left2 = 0;
		String ctgStr = null;
		//Map<String,ReadPair> reads = new HashMap<String, ReadPair>();
		String tmp = null;
		PrintStreamPair psPair = null;
		Contig ctg1 = null;
		Contig ctg2 = null;
		int qual = -1;
		File mapMatchDir = new File(fishDir,base+".dat");
		mapMatchDir.mkdir();
		int ctgNameComp = -10;
		System.out.print("[a5_fie] Reading SAM file...");
		long currPos = start;
		double perc = 0;
		double ten = 1;
		int numKeep = 0;
		int total = 0;
		long before = System.currentTimeMillis();
		while (br.ready()){
			currPos = fis.getChannel().position()-start;
			//System.err.println(perc);
			if (((double)currPos/len)*10 > ten){
				System.out.print(".."+NF.format(10*ten)+"%");
				ten++;
			}
			line1 = br.readLine().split("\t");
			line2 = br.readLine().split("\t");
		//	System.err.println("pos = " + raf.getFilePointer() + " len = "+raf.length());
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
			
			// phred-scaled geometric-mean of posterior probabilites that mapping position is incorrect
			// mod 255 in case the quality score couldn't be computed.	
			qual = ((Integer.parseInt(line1[4]) % 255) + (Integer.parseInt(line2[4]) % 255))/2;
			ctgNameComp = line1[2].compareTo(line2[2]); // order for consistency
			
			if (ctgNameComp < 0){
				if (inRange(ranges,left1, cigarLength(line1[5]), ctg1, isReverse(line1[1]), left2, cigarLength(line2[5]), ctg2, isReverse(line2[1])))
					continue;
				ctgStr = line1[2]+"-"+line2[2];
				if (psPairs.containsKey(ctgStr))
					psPair = psPairs.get(ctgStr);
				else {
					psPair = new PrintStreamPair(ctg1, ctg2, mapMatchDir);
					psPairs.put(ctgStr, psPair);
				}
				psPair.print(left1, left2, qual);
			} else if (ctgNameComp > 0) {
				if (inRange(ranges,left1, cigarLength(line1[5]), ctg1, isReverse(line1[1]), left2, cigarLength(line2[5]), ctg2, isReverse(line2[1])))
					continue;
				ctgStr = line2[2]+"-"+line1[2];
				if (psPairs.containsKey(ctgStr))
					psPair = psPairs.get(ctgStr);
				else { 
					psPair = new PrintStreamPair(ctg2, ctg1, mapMatchDir);
					psPairs.put(ctgStr, psPair);
				}
				psPair.print(left2, left1, qual);
			} else { // same contig, so check to see it's within the given ranges
				if (inRange(ranges,(left2 > left1 ? left2+cigarLength(line2[5])-left1 : left1+cigarLength(line1[5])-left2))) 
						continue;
				
				ctgStr = line2[2]+"-"+line1[2];
				if (psPairs.containsKey(ctgStr))
					psPair = psPairs.get(ctgStr);
				else { 
					psPair = new PrintStreamPair(ctg2, ctg1, mapMatchDir);
					psPairs.put(ctgStr, psPair);
				}
				if (left2 < left1) // order for consistency
					psPair.print(left2, left1, qual);
				else 
					psPair.print(left1, left2, qual);				
			}
			
			// add points for our map files
			if (mapPoints.containsKey(ctg1.name)){
				tmpPts = mapPoints.get(ctg1.name);
				tmpPSPairs = ctgPSPairs.get(ctg1.name);
			} else {
				tmpPts = new TreeSet<Double>();
				mapPoints.put(ctg1.name, tmpPts);
				tmpPSPairs = new Vector<String>();
				ctgPSPairs.put(ctg1.name, tmpPSPairs);
			}
			tmpPts.add(new Double(left1));
			tmpPSPairs.add(ctgStr);
			if (mapPoints.containsKey(ctg2.name)){
				tmpPts = mapPoints.get(ctg2.name);
				tmpPSPairs = ctgPSPairs.get(ctg2.name);
			} else {
				tmpPts = new TreeSet<Double>();
				mapPoints.put(ctg2.name, tmpPts);
				tmpPSPairs = new Vector<String>();
				ctgPSPairs.put(ctg2.name, tmpPSPairs);
			}
			tmpPSPairs.add(ctgStr);
			tmpPts.add(new Double(left2));
			numKeep++;
		}
		long after = System.currentTimeMillis();
		System.out.println("..100%... done!... Took "+(after-before)/1000+" seconds.");
		perc = (double) numKeep / total * 100;
		System.out.println("[a5_fie] Keeping "+NF.format(perc)+"% ("+numKeep+"/"+total+") of reads.");
		Iterator<String> ctgIt = mapPoints.keySet().iterator();
		Set<String> ctgToRm = new HashSet<String>();
		Set<String> psPairsToRm = new HashSet<String>();
		while(ctgIt.hasNext()){
			tmp = ctgIt.next();
			tmpPts = mapPoints.get(tmp);
			if (tmpPts.size() < MIN_PTS) {
				ctgToRm.add(tmp);
				psPairsToRm.addAll(ctgPSPairs.get(tmp));
			}
		}
		//removeKeys(mapPoints, ctgToRm);
		removeKeys(ctgs, ctgToRm);
		removeKeys(mapPoints, ctgToRm);
		removeKeys(psPairs, psPairsToRm);
		Contig[] ctgRef = new Contig[ctgs.values().size()];
		ctgs.values().toArray(ctgRef);
		Arrays.sort(ctgRef,new Comparator<Contig>(){
			public int compare(Contig o1, Contig o2) {
				return o1.getId() - o2.getId();
			}
		});
		// Rank contigs so FISH doesn't freak out about non-consecutive ids
		for (int i = 0; i < ctgRef.length; i++){
			ctgRef[i].setRank(i+1);
		}
		
		File ctgLblFile = new File(fishDir,base+".contig_labels.txt");
		ctgLblFile.createNewFile();
		PrintStream ctgLblOut = new PrintStream(ctgLblFile);
		
		Iterator<String> it = mapPoints.keySet().iterator();
		TreeMap<Integer,String> maps = new TreeMap<Integer, String>();
		while(it.hasNext()){
			String ctg = it.next();
			tmpPts = mapPoints.get(ctg);
			Double[] ar = new Double[tmpPts.size()];
			tmpPts.toArray(ar);
			File ctgMapFile = new File(mapMatchDir,"map."+ctgs.get(ctg).getRank()+".txt");
			ctgLblOut.println(ctgs.get(ctg).getRank()+"\t"+ctgs.get(ctg).name);
			PrintStream ctgMapOut = new PrintStream(ctgMapFile);
			Arrays.sort(ar, posComp);
			for (Double d: ar){
				ctgMapOut.println("c"+ctgs.get(ctg).getRank()+"p"+Math.abs(d.intValue())+"\t0");
			}
			
			ctgMapOut.close();
			maps.put(ctgs.get(ctg).getRank(),fishDir.getName()+"/"+mapMatchDir.getName()+"/"+ctgMapFile.getName());
		}
		ctgLblOut.close();
		File controlFile = new File(fishDir,base+".control.txt");
		PrintStream ctrlOut = new PrintStream(controlFile);
		Iterator<Integer> it2 = maps.keySet().iterator();
		ctrlOut.println("-maps");
		while(it2.hasNext()){
			Integer tmpInt = it2.next();
			ctrlOut.println(tmpInt+"\t"+maps.get(tmpInt));
		}
		
		PrintStreamPair tmpPSP = null;
		ctrlOut.println("-matches");
		it = psPairs.keySet().iterator();
		Vector<String> mfToRm = new Vector<String>();
		
		while(it.hasNext()){
			ctgStr = it.next();
			tmpPSP = psPairs.get(ctgStr);
			if (tmpPSP.nPairs <= 2)
				mfToRm.add(ctgStr);
		}
		removeKeys(psPairs, mfToRm);
		Set<MatchFile> matchFiles = new TreeSet<MatchFile>();
		it = psPairs.keySet().iterator();
		while(it.hasNext()) {
			ctgStr = it.next();
			if (psPairs.get(ctgStr).getNumPairs() < MIN_PTS)
				continue;
			psPairs.get(ctgStr).close();
			matchFiles.addAll(psPairs.get(ctgStr).getMatchFiles());
		}
		Iterator<MatchFile> mfIt = matchFiles.iterator();
		MatchFile tmpMF = null;
		while(mfIt.hasNext()){
			tmpMF = mfIt.next();
			ctrlOut.println(tmpMF.toString()+"\t"+fishDir.getName()+"/"+mapMatchDir.getName()+"/"+tmpMF.getName());
		}
		
		ctrlOut.println("end");
		ctrlOut.close();
		
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
			if (ins >= ranges[i][0] || ins <= ranges[i][1])
				return true;
		return false;
	}
	
	private static boolean inRange(double[][] ranges, int pos1, int len1, Contig ctg1, boolean rev1, int pos2, int len2, Contig ctg2, boolean rev2){
		if (pos1 == 103131 && pos2 == 46024)
			System.out.print("");
		int ins = 0;
		if (rev1 == rev2) { // one of our contigs is invertedins = 0;
			if (pos1 > ctg1.len/2)
				ins += ctg1.len - pos1 + 1;
			else
				ins += pos1 + len1;
			if (pos2 > ctg2.len/2)
				ins+= ctg2.len - pos2 + 1;
			else
				ins += pos2 + len2;
		} else {
			if (INWARD){
				if (rev2)
					ins = (pos2+len2-1) + (ctg1.len-pos1+1);
				else 
					ins = (pos1+len1-1) + (ctg2.len-pos2+1);
			}
			if (OUTWARD) {
				if (rev2)
					ins = (pos1+len1-1) + (ctg2.len-pos2+1);
				else
					ins = (pos2+len2-1) + (ctg1.len-pos1+1);
			}
		}
		for (int i = 0; i < ranges.length; i++)
			if (ins >= ranges[i][0] && ins <= ranges[i][1])
				return true;
		return false;
	}
	

	/**
	 * Returns an array of arrays with each array having the following information:	<br><br>
	 * <code>[ mean, sd, n, nSd ]</code>
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
		System.out.println("[a5_fie] Initial stats for sample: mu="+NF.format(ins[0])+" sd="+NF.format(ins[1])+" n="+NF.format(ins[2]));
		System.out.print("[a5_fie] EM-clustering insert sizes with K="+K+"... ");
		double delta = 0.00005;
		EMClusterer em = new EMClusterer(toFilt, K);
		long before = System.currentTimeMillis();
		int iters = em.iterate(1000, delta);
		long after = System.currentTimeMillis();
		System.out.println("stopping after "+iters+" iterations with delta="+delta+". Took "+(after-before)/1000+" seconds.");
		ReadSet[] clusters = new ReadSet[K];
		em.getClusters().toArray(clusters);
		double[][] allIns = new double[K][4];
		for (int i = 0; i < clusters.length; i++){
			allIns[i][0] = i;
			allIns[i][1] = clusters[i].mean();
			allIns[i][2] = clusters[i].sd();
			allIns[i][3] = ((double)clusters[i].size())/((double)toFilt.size());
		}
		
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
		
		System.out.println("[a5_fie] Found the following clusters:");	
		
		Vector<ReadSet> signal = new Vector<ReadSet>();
		Vector<ReadSet> noise = new Vector<ReadSet>();
		for (int i = 0; i < clusters.length; i++){
			NF.setMaximumFractionDigits(0);
			System.out.print("[a5_fie] cluster"+NF.format(clusters[i].getId())+": mu="+pad(NF.format(allIns[i][1]),10)+
					"sd="+pad(NF.format(allIns[i][2]),10)+"n="+pad(NF.format(clusters[(int)allIns[i][0]].size()),10));
			NF.setMaximumFractionDigits(2);
			System.out.print("perc="+pad(NF.format(100*allIns[i][3]),10));
			
			// if insert size distribution is under dispersed, add all these reads to the signal pile
			// take up to numLibs underdispersed clusters 
			if (clusters[i].sd() <= clusters[i].mean() && i < numLibs){
				signal.add(clusters[i]);
				System.out.println("  (signal)");
			} else {
				noise.add(clusters[i]);
				System.out.println("  (noise)");
			}
			
		}
		Iterator<ReadSet> sigIt = null;
		ReadSet sigSet = null;
	/*	Iterator<ReadSet> noiseIt = noise.iterator();
		ReadSet currSet = null;
		double outInsert = 0.0;
		int numEndSp = 0;
		while(noiseIt.hasNext()){
			// iterate over all noisy reads, and check to see if their outside insert size fits any of the other clusters
			// do this so signal from circular molecules doesn't get mistaken for noise
			currSet = noiseIt.next();
			rpIt = currSet.getReads().iterator();
			while(rpIt.hasNext()){
				tmp = rpIt.next();
				// outside distance
				outInsert = tmp.ctg1.len-tmp.pos2+tmp.pos1+tmp.len1;
				sigIt = signal.iterator();
				while(sigIt.hasNext()){
					sigSet = sigIt.next();
					if (sigSet.p(outInsert) > currSet.p(tmp.getInsert())){
						tmp.setEndSpanning(true);
						currSet.remove(tmp);
						sigSet.add(tmp);
						numEndSp++;
					} 
				}
			}	
		}
		
		System.out.println("[a5_fie] "+numEndSp+" putatively noise reads look like signal if we let the pair span the end of the contig.");*/
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
			NF.setMaximumFractionDigits(0);
			System.out.print("[a5_fie] cluster"+NF.format(sigSet.getId())+": mu="+pad(NF.format(sigSet.mean()),10)+
					"sd="+pad(NF.format(sigSet.sd()),10)+"n="+pad(NF.format(sigSet.size()),10));
			NF.setMaximumFractionDigits(2);
			System.out.println("perc="+pad(NF.format(100*((double)sigSet.size()/(double)toFilt.size())),10)+"nSd="+nSd);
			i++;
		}
		if (rmClusters.length()>0)
			System.out.println("[a5_fie] Removing "+rmClusters);
		ins = ReadPair.estimateInsertSize(reads.values());
		System.out.println("[a5_fie] Final stats for sample after filtering: mu="+NF.format(ins[0])+" sd="+NF.format(ins[1])+" n="+NF.format(ins[2]));
		for (i = 0; i < ret.length; i++){
			//double[] tmpAr = {ret[i][0]-ret[i][1]*ret[i][3],ret[i][0]+ret[i][1]*ret[i][3]};
			double[] tmpAr = {1,2*ret[i][0]};
			if (tmpAr[0] > min)
				tmpAr[0] = min;
			if (tmpAr[1] < max)
				tmpAr[1] = max;
			ret[i] = tmpAr;
			
		}
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
				System.err.println("Found nameless contig in SAM header");
			else if (len == -1) 
				System.err.println("Found contig of unknown length in SAM header");
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
			// if pair didn't map, just jump to next
			// line instead of jumping a full step
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
	
	public static class PrintStreamPair implements Comparable<PrintStreamPair>{
		// a comparator for sorting reads by position
		private static final Comparator<int[]> COMP = new Comparator<int[]>(){
			public int compare(int[] arg1, int[] arg2) {
				if (arg1[0] == arg2[0])
					return 0;	
				else if (arg1[1] == arg2[1])
					return 0;
				else if (arg1[0] < arg2[0])
					return -1;
				else 
					return 1;
			} 
			
		};
		// File for both directions
		File file1;
		File file2;
		File fishDir;
		private Contig ctg1;
		private Contig ctg2;
		private int nPairs;
		
		private TreeSet<int[]> v1;
		private TreeSet<int[]> v2;
		
		public PrintStreamPair(Contig ctg1, Contig ctg2, File fishDir) throws IOException{
			nPairs = 0;
			this.fishDir = fishDir;
			if (ctg1==ctg2){
				this.ctg1 = ctg1;
				this.ctg2 = ctg2;
			} else { 
				this.ctg1 = ctg1;
				this.ctg2 = ctg2;
			}
			v1 = new TreeSet<int[]>(COMP);
			v2 = new TreeSet<int[]>(COMP);
		}
		
		public void add(ReadPair pair){
			print(pair.pos1,pair.pos2,((pair.qual1%255)+(pair.qual2%255))/2);
		}
				
		/**
		 * Assumes pos1 is from ctg1 and pos2 is from ctg2
		 * 
		 * @param pos1 position on <code>ctg1</code> from constructor 
		 * @param pos2 position on <code>ctg2</code> from constructor
		 * @param qual
		 */
		public void print(int pos1, int pos2, int qual) {
			if (pos1 == 111302 || pos2 == 111302){
				System.out.print("");
			}
			
			if (this.ctg1 != this.ctg2){
				int[] ar = {pos1,pos2,qual};
				boolean added = v1.add(ar);
				if (added){
					int[] ar2 = {pos2,pos1,qual};
					added = v2.add(ar2);
					if (!added)
						v1.remove(ar);
					else
						nPairs++;
				}
			} else { 
				int[] ar = {pos1,pos2,qual};
				boolean added = v1.add(ar);
				if (added) 
					nPairs++;
			}
		}
		
		public int getNumPairs(){
			return nPairs;
		}
		
		public void close(){
			if (v1.size() >= MIN_PTS || v2.size() >= MIN_PTS){
				try {
					file1 = new File(fishDir,"match."+ctg1.getRank()+"v"+ctg2.getRank()+".txt");
					file1.createNewFile();
					PrintStream out = new PrintStream(file1);
					Iterator<int[]> it = v1.iterator();
					int[] tmp = null;
					while(it.hasNext()){
						tmp = it.next();					
						out.print("c"+this.ctg1.getRank()+"p"+tmp[0]+"\tc"+this.ctg2.getRank()+"p"+tmp[1]+"\t"+tmp[2]+"\n");
						
					}
					out.close();
					if (ctg1 != ctg2) {
						file2 = new File(fishDir,"match."+ctg2.getRank()+"v"+ctg1.getRank()+".txt");
						file2.createNewFile();
						out = new PrintStream(file2);
						it = v2.iterator();
						tmp = null;
						while(it.hasNext()){
							tmp = it.next();					
							out.print("c"+this.ctg2.getRank()+"p"+tmp[0]+"\tc"+this.ctg1.getRank()+"p"+tmp[1]+"\t"+tmp[2]+"\n");
							
						}
						out.close();
					}
				} catch (IOException e) {
					e.printStackTrace();
					System.exit(-1);
				}
			}
		}

		
		public Set<MatchFile> getMatchFiles(){
			Set<MatchFile> mfiles = new TreeSet<MatchFile>();
			mfiles.add(new MatchFile(ctg1, ctg2, file1));
			if (!ctg1.equals(ctg2))
				mfiles.add(new MatchFile(ctg2, ctg1, file2));
			return mfiles;
		}

		@Override
		public int compareTo(PrintStreamPair arg0) {
			if (this.ctg1.getRank() - arg0.ctg1.getRank() == 0){
				return this.ctg2.getRank() - arg0.ctg2.getRank();
			} else {
				return this.ctg1.getRank() - arg0.ctg1.getRank();
			}
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
