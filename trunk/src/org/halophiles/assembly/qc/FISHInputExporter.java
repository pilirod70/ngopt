package org.halophiles.assembly.qc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
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
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Vector;

import org.halophiles.assembly.Contig;
import org.halophiles.assembly.ReadPair;
import org.halophiles.assembly.ReadSet;
import org.halophiles.assembly.SAMFileParser;

public class FISHInputExporter {
	private static NumberFormat NF;
	private static int N_ESTREADPAIRS = 100000;
	private static int MINQUAL = 13;
	
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
			int numLibs = 1;
			if (args.length == 4){
				numLibs = Integer.parseInt(args[3]);
			}
			System.out.println("[a5_fie] Reading "+samPath);
			System.out.println("[a5_fie] Writing output to "+outdir.getAbsolutePath()+" with basename "+base);
			RandomAccessFile raf = new RandomAccessFile(new File(samPath), "r");
			Map<String,Contig> contigs = readContigs(raf);
			System.out.println("[a5_fie] Found "+contigs.size()+" contigs");
			System.out.println("[a5_fie] Reading in a subset of reads for insert size estimation");
			long before = System.currentTimeMillis();
			Map<String,ReadPair> reads = readSubset(raf, contigs);
			long after = System.currentTimeMillis();
			System.out.println("[a5_fie] Took "+((after-before)/1000)+" seconds.");
			
			double[][] ranges = getFilterRanges(reads, numLibs);
			
			
			/*
			SAMFileParser sfp = new SAMFileParser(samPath,N_ESTREADPAIRS);

			System.out.println("[a5_fie] Found "+sfp.getNumContigs()+" contigs.");
			System.out.println("[a5_fie] Found "+sfp.getNumReads()+" total reads and "+sfp.getNumPairs()+" read-pairs.");
			
			Map<String,ReadPair> reads = new HashMap<String,ReadPair>();
			System.out.println("[a5_fie] Filtering duplicate connections...");
			Set<ReadPair> uniq = new TreeSet<ReadPair>(new Comparator<ReadPair>(){
				public int compare(ReadPair arg0, ReadPair arg1) {
					if (arg0.ctg1.equals(arg1.ctg1)){
						if (arg0.ctg2.equals(arg1.ctg2)){
							if (arg0.pos1 == arg1.pos1){
								if (arg0.pos2 == arg1.pos2){
									return 0;
								} else
									return arg0.pos2 - arg1.pos2;
							} else
								return arg0.pos1 - arg1.pos1;
						} else 
							return arg0.ctg2.compareTo(arg1.ctg2);
					} else 
						return arg0.ctg1.compareTo(arg1.ctg1);
				}
			});
			Iterator<ReadPair> rpIt = sfp.getReadPairs();
			while(rpIt.hasNext()){
				ReadPair tmp = rpIt.next();
				if (!tmp.paired) 
					continue;
				if (tmp.ctg1.equals(tmp.ctg2) && tmp.getInsert() < MINQUAL)
					continue;
				
				if(uniq.add(tmp))
					reads.put(tmp.hdr, tmp);
			}
			*/
			
			if (reads.size() <= 0) {
				System.err.println("[a5_fie] No paired reads found. Cannot generate match files for running FISH misassembly detection.");
				System.exit(-1);
			}
		/*
			File insFile = new File(outdir, base+".ins_size_unfilt.txt");
			insFile.createNewFile();
			PrintStream insOut = new PrintStream(insFile);
			ReadPair.exportInsertSize(reads.values(), insOut);
			insOut.close();

			System.out.println("[a5_fie] Filtering proper connections... ");
			int before = reads.size();
			Collection<ReadPair> filtered = filterProperConnections(reads, numLibs);
			System.out.println("[a5_fie] removed "+(before - reads.size())+" of "+before +" read pairs.");

			insFile = new File(outdir,base+".ins_size_rm.txt");
			insFile.createNewFile();
			insOut = new PrintStream(insFile);
			ReadPair.exportInsertSize(filtered, insOut);
			insOut.close();
		*/
		/*	
			System.out.print("[a5_fie] Filtering tandem connections... ");
			before = reads.size();
			filterTandemConnections(reads);
			System.out.println("removed "+(before - reads.size())+" of "+before +" read pairs.");
		*/	
		/*
			insFile = new File(outdir,base+".ins_size_filt.txt");
			insFile.createNewFile();
			insOut = new PrintStream(insFile);
			ReadPair.exportInsertSize(reads.values(), insOut);
			insOut.close();
		*/
			export(reads,outdir,base);
			
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
			for (Double d: ar)
				ctgMapOut.println("c"+contigs.get(ctg).getRank()+"p"+Math.abs(d.intValue())+"\t"+(d < 0 ? -1 : 1));
			
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
			if (tmpPSP.reads.size() <= 2)
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
	
	public static void export2(String samPath, File fishDir, String base, Map<String,Contig> ctgs, double[][] ranges) throws IOException{
		Map<String,Vector<Double>> mapPoints = new HashMap<String,Vector<Double>>();
		Map<String,PrintStreamPair> psPairs = new TreeMap<String,PrintStreamPair>();
		Vector<Double> tmpPts = null;

		BufferedReader br = new BufferedReader (new FileReader(new File(samPath)));
		while(nextCharIs(br, '@')) 
			br.readLine();
		
		String[] line1 = null;
		String[] line2 = null;
		int left1 = 0;
		int left2 = 0;
		String ctgStr = null;
		Map<String,ReadPair> reads = new HashMap<String, ReadPair>();
		ReadPair tmp = null;
		PrintStreamPair psPair = null;
		Contig ctg1 = null;
		Contig ctg2 = null;
		int qual = -1;
		File mapMatchDir = new File(fishDir,base+".dat");
		int comp = -10;
		int dist = 0;
		while (br.ready()){
			line1 = br.readLine().split("\t");
			line2 = br.readLine().split("\t");
		//	System.err.println("pos = " + raf.getFilePointer() + " len = "+raf.length());
			while (!line1[0].equals(line2[0]) && br.ready()){
				line1 = line2;
				line2 = br.readLine().split("\t");				
			}
			left1 = Integer.parseInt(line1[3]);
			left2 = Integer.parseInt(line2[3]);
			// if pair didn't map, just jump to next
			// line instead of jumping a full step
			if (left1 == 0 || left2 == 0) continue;
			
			ctg1 = ctgs.get(line1[2]);
			ctg2 = ctgs.get(line2[2]);
			
			qual = ((Integer.parseInt(line1[4]) % 255) + (Integer.parseInt(line2[4]) % 255))/2;
			comp = line1[2].compareTo(line2[2]); // order for consistency
			if (comp < 0){
				ctgStr = line1[2]+"-"+line2[2];
				if (psPairs.containsKey(ctgStr))
					psPair = psPairs.get(ctgStr);
				else {
					psPair = new PrintStreamPair(ctg1, ctg2, mapMatchDir);
					psPairs.put(ctgStr, psPair);
				}
				psPair.print(left1, left2, qual);
			} else if (comp > 0) {
				ctgStr = line2[2]+"-"+line1[2];
				if (psPairs.containsKey(ctgStr))
					psPair = psPairs.get(ctgStr);
				else { 
					psPair = new PrintStreamPair(ctg2, ctg1, mapMatchDir);
					psPairs.put(ctgStr, psPair);
				}
				psPair.print(left2, left1, qual);
			} else { // same contig, so check to see it's within the given ranges
				dist = Math.abs(left1 - left2);
				for (int i = 0; i < ranges.length; i++)
					if (dist <= ranges[i][1] && dist >= ranges[i][0]) 
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
			if (mapPoints.containsKey(ctg1.name))
				tmpPts = mapPoints.get(ctg1.name);
			else {
				tmpPts = new Vector<Double>();
				mapPoints.put(ctg1.name, tmpPts);
			}
			tmpPts.add(new Double(left1));
			if (mapPoints.containsKey(ctg2.name))
				tmpPts = mapPoints.get(ctg2.name);
			else {
				tmpPts = new Vector<Double>();
				mapPoints.put(ctg2.name, tmpPts);
			}
			tmpPts.add(new Double(left2));
		}
		
	}
	

	/**
	 * 
	 * @return a collection of ReadPairs that were removed
	 */
	private static double[][] getFilterRanges(Map<String,ReadPair> reads, int numLibs){
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
		System.out.println("[a5_fie] Initial read set stats: mu="+NF.format(ins[0])+" sd="+NF.format(ins[1])+" n="+NF.format(ins[2]));
		System.out.print("[a5_fie] EM-clustering insert sizes with K="+K+"... ");
		
		EMClusterer em = new EMClusterer(toFilt, K);
		double delta = 0.0001;
		int iters = em.iterate(1000, delta);
		System.out.println("stopping after "+iters+" iterations with delta="+delta);
		ReadSet[] clusters = new ReadSet[K];
		em.getClusters().toArray(clusters);
		double[][] allIns = new double[K][4];
		for (int i = 0; i < clusters.length; i++){
			allIns[i][0] = i;
			allIns[i][1] = clusters[i].mean();
			allIns[i][2] = clusters[i].sd();
			allIns[i][3] = ((double)clusters[i].size())/((double)toFilt.size());
		}
		
		Vector<ReadSet> signal = new Vector<ReadSet>();
		Vector<ReadSet> noise = new Vector<ReadSet>();
		for (int i = 0; i < clusters.length; i++){
			NF.setMaximumFractionDigits(0);
			System.out.print("[a5_fie] cluster"+NF.format(clusters[i].getId())+": mu="+pad(NF.format(allIns[i][1]),10)+
					"sd="+pad(NF.format(allIns[i][2]),10)+"n="+pad(NF.format(clusters[(int)allIns[i][0]].size()),10));
			NF.setMaximumFractionDigits(2);
			System.out.print("perc="+pad(NF.format(100*allIns[i][3]),10));
			
			// if insert size distribution is under dispersed, add all these reads to the signal pile
			if (clusters[i].sd() <= clusters[i].mean()){
				// if we have a small and very high variance, but not overdispersed cluster, don't discard 
				if(clusters[i].sd()/clusters[i].mean() >= 0.90 && allIns[i][3] < 0.05){
					noise.add(clusters[i]);
					System.out.println("  (noise)");
				} else {
					signal.add(clusters[i]);
					System.out.println("  (signal)");
				}
			} else {
				noise.add(clusters[i]);
				System.out.println("  (noise)");
			}
			
		}
		Iterator<ReadSet> noiseIt = noise.iterator();
		Iterator<ReadSet> sigIt = null;
		ReadSet currSet = null;
		ReadSet sigSet = null;
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
		
		System.out.println("[a5_fie] "+numEndSp+" putatively noise reads look like signal if we let the pair span the end of the contig.");
		String rmClusters = "";
		sigIt = signal.iterator();
		//Collection<ReadPair> ret = new Vector<ReadPair>();
		int nSd = 6;
		double[][] ret = new double[signal.size()][];
		int i = 0;
		while(sigIt.hasNext()){
			nSd=6;
			sigSet = sigIt.next();
			rmClusters += " cluster"+sigSet.getId();
			//ret.addAll(sigSet.getReads());
			NF.setMaximumFractionDigits(0);
			System.out.print("[a5_fie] cluster"+NF.format(sigSet.getId())+": mu="+pad(NF.format(sigSet.mean()),10)+
					"sd="+pad(NF.format(sigSet.sd()),10)+"n="+pad(NF.format(sigSet.size()),10));
			NF.setMaximumFractionDigits(2);
			System.out.print("perc="+pad(NF.format(100*((double)sigSet.size()/(double)toFilt.size())),10));
			/*	rpIt = sigSet.getReads().iterator();
			while(rpIt.hasNext()){
				toRm.add(rpIt.next().hdr);
			}*/
			while(nSd*sigSet.sd()>sigSet.mean())
				nSd--;
			ret[i] = new double[4];
			ret[i][0] = sigSet.mean();
			ret[i][1] = sigSet.sd();
			ret[i][2] = sigSet.size();
			ret[i][3] = nSd;
			removeKeys(reads,sigSet.getReadHdrs());
		/*
			System.out.println("[a5_fie] Removing reads with inserts between ("+
					NF.format(sigSet.mean()-nSd*sigSet.sd())+","+
					NF.format(sigSet.mean()+nSd*sigSet.sd())+")   nSd="+nSd);
		*/
			// EM-clustering sometimes puts proper connections with noise
		//	ReadPair.filterRange(reads, sigSet.mean()-nSd*sigSet.sd(), sigSet.mean()+nSd*sigSet.sd());
			i++;
		}
		if (rmClusters.length()>0)
			System.out.println("[a5_fie] Removed clusters"+rmClusters);
		ins = ReadPair.estimateInsertSize(reads.values());
		System.out.println("[a5_fie] Final stats (for syntenic read pairs): mu="+NF.format(ins[0])+" sd="+NF.format(ins[1])+" n="+NF.format(ins[2]));
		return ret;
	}
	
	/**
	 * 
	 * @return a collection of ReadPairs that were removed
	 */
	private static Collection<ReadPair> filterProperConnections(Map<String,ReadPair> reads, int numLibs){
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
		System.out.println("[a5_fie] Initial read set stats: mu="+NF.format(ins[0])+" sd="+NF.format(ins[1])+" n="+NF.format(ins[2]));
		System.out.print("[a5_fie] EM-clustering insert sizes with K="+K+"... ");
		
		EMClusterer em = new EMClusterer(toFilt, K);
		double delta = 0.0001;
        int iters = em.iterate(1000, delta);
        System.out.println("stopping after "+iters+" iterations with delta="+delta);
		ReadSet[] clusters = new ReadSet[K];
		em.getClusters().toArray(clusters);
		double[][] allIns = new double[K][4];
		for (int i = 0; i < clusters.length; i++){
			allIns[i][0] = i;
			allIns[i][1] = clusters[i].mean();
			allIns[i][2] = clusters[i].sd();
			allIns[i][3] = ((double)clusters[i].size())/((double)toFilt.size());
		}
		
		Vector<ReadSet> signal = new Vector<ReadSet>();
		Vector<ReadSet> noise = new Vector<ReadSet>();
		for (int i = 0; i < clusters.length; i++){
			NF.setMaximumFractionDigits(0);
			System.out.print("[a5_fie] cluster"+NF.format(clusters[i].getId())+": mu="+pad(NF.format(allIns[i][1]),10)+
					"sd="+pad(NF.format(allIns[i][2]),10)+"n="+pad(NF.format(clusters[(int)allIns[i][0]].size()),10));
			NF.setMaximumFractionDigits(2);
			System.out.print("perc="+pad(NF.format(100*allIns[i][3]),10));
			
			// if insert size distribution is under dispersed, add all these reads to the signal pile
			if (clusters[i].sd() <= clusters[i].mean()){
				// if we have a small and very high variance, but not overdispersed cluster, don't discard 
				if(clusters[i].sd()/clusters[i].mean() >= 0.90 && allIns[i][3] < 0.05){
					noise.add(clusters[i]);
					System.out.println("  (noise)");
				} else {
					signal.add(clusters[i]);
					System.out.println("  (signal)");
				}
			} else {
				noise.add(clusters[i]);
				System.out.println("  (noise)");
			}
			
		}
		Iterator<ReadSet> noiseIt = noise.iterator();
		Iterator<ReadSet> sigIt = null;
		ReadSet currSet = null;
		ReadSet sigSet = null;
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
		
		System.out.println("[a5_fie] Removed "+numEndSp+" putatively noise reads that look like signal if we let the pair span the end of the contig.");
		String rmClusters = "";
		sigIt = signal.iterator();
		Collection<ReadPair> ret = new Vector<ReadPair>();
		int nSd = 6;
		while(sigIt.hasNext()){
			nSd=6;
			sigSet = sigIt.next();
			rmClusters += " "+sigSet.toString();
			ret.addAll(sigSet.getReads());
			NF.setMaximumFractionDigits(0);
			System.out.print("[a5_fie] cluster"+NF.format(sigSet.getId())+": mu="+pad(NF.format(sigSet.mean()),10)+
					"sd="+pad(NF.format(sigSet.sd()),10)+"n="+pad(NF.format(sigSet.size()),10));
			NF.setMaximumFractionDigits(2);
			System.out.print("perc="+pad(NF.format(100*((double)sigSet.size()/(double)toFilt.size())),10));
		/*	rpIt = sigSet.getReads().iterator();
			while(rpIt.hasNext()){
				toRm.add(rpIt.next().hdr);
			}*/
			while(nSd*sigSet.sd()>sigSet.mean())
				nSd--;
			removeKeys(reads,sigSet.getReadHdrs());
			System.out.println("[a5_fie] Removing reads with inserts between ("+
			                        NF.format(sigSet.mean()-nSd*sigSet.sd())+","+
			                        NF.format(sigSet.mean()+nSd*sigSet.sd())+")   nSd="+nSd);
			// EM-clustering sometimes puts proper connections with noise
			ReadPair.filterRange(reads, sigSet.mean()-nSd*sigSet.sd(), sigSet.mean()+nSd*sigSet.sd());
		}
		if (rmClusters.length()>0)
			System.out.println("[a5_fie] Removed clusters"+rmClusters);
		ins = ReadPair.estimateInsertSize(reads.values());
		System.out.println("[a5_fie] Final stats (for syntenic read pairs): mu="+NF.format(ins[0])+" sd="+NF.format(ins[1])+" n="+NF.format(ins[2]));
		return ret;
	}
	
	private static void filterTandemConnections(Map<String,ReadPair> reads){
		Iterator<String> it = reads.keySet().iterator();
		Set<Contig> contigs = new HashSet<Contig>();
		ReadPair pair = null;
		// get a set of the contigs
		while(it.hasNext()){
			pair = reads.get(it.next());
			contigs.add(pair.ctg1);
			contigs.add(pair.ctg2);
		}
		Iterator<Contig> ctgIt = contigs.iterator();
		List<ReadPoint> points = null;
		Set<String> toRm = null;
		ReadPoint rp1 = null;
		ReadPoint rp2 = null;
		Iterator<ReadPoint> ptIt = null;
		Iterator<ReadPair> rpIt = null;
		while(ctgIt.hasNext()){
			// get all pairs where both reads map to this contig
			rpIt = getReadPairs(reads, ctgIt.next()).iterator();
			points = new Vector<ReadPoint>();
			while(rpIt.hasNext()){
				pair = rpIt.next();
				points.add(new ReadPoint(pair.hdr,pair.pos1));
				points.add(new ReadPoint(pair.hdr,pair.pos2));
			}
			if (points.size() < 2)
				continue;
			// sort by position
			Collections.sort(points);
			toRm = new HashSet<String>();
			ptIt = points.iterator();
			rp1 = ptIt.next();
			while(ptIt.hasNext()){
				rp2 = ptIt.next();
				if (rp1.hdr.equals(rp2.hdr))
					toRm.add(rp1.hdr);
				rp1 = rp2;
			}
			removeKeys(reads, toRm);
		}
	}
	
	private static Collection<ReadPair> getReadPairs(Map<String,ReadPair> reads, Contig ctg){
		Iterator<String> it = reads.keySet().iterator();
		Vector<ReadPair> ret = new Vector<ReadPair>();
		ReadPair pair = null;
		while(it.hasNext()){
			pair = reads.get(it.next());
			if (pair.ctg1.equals(ctg) && pair.ctg2.equals(ctg))
				ret.add(pair);
		}
		return ret;
	}
	
	private static <V> void removeKeys(Map<String,V> reads, Collection<String> torm){
		Iterator<String> it = torm.iterator();
		while(it.hasNext())
			reads.remove(it.next());
	}
	
	/**
	 * This method assumes that the contig header has been read in, and 
	 * the file pointer is at the beginning of the first read's line.
	 * 
	 * @param raf
	 * @throws IOException 
	 */
	private static Map<String,ReadPair> readSubset(RandomAccessFile raf, Map<String,Contig> contigs) throws IOException{
		int HANDFUL_SIZE = 1000;
		long pos = raf.getFilePointer();
		long len = raf.length();
		long step = (len-pos)/(N_ESTREADPAIRS*HANDFUL_SIZE);
		int i = 0;
		String[] line1 = null;
		String[] line2 = null;
		int left1 = 0;
		int left2 = 0;
		Map<String,ReadPair> reads = new HashMap<String, ReadPair>();
		ReadPair tmp = null;
		int j = 0;
		while (i < N_ESTREADPAIRS && pos < raf.length()){
			raf.readLine(); // make sure we're at the beginning of a line
			
			line1 = raf.readLine().split("\t");
			line2 = raf.readLine().split("\t");
			// get to the first read in a pair
			while (!line1[0].equals(line2[0]) && raf.getFilePointer() < raf.length()){
				line1 = line2;
				line2 = raf.readLine().split("\t");				
			}
			left1 = Integer.parseInt(line1[3]);
			left2 = Integer.parseInt(line2[3]);
			// if pair didn't map, just jump to next
			// line instead of jumping a full step
			if (left1 == 0 || left2 == 0) continue;
			tmp = new ReadPair(line1[0]);
			tmp.addRead(left1, isReverse(line1[1]), cigarLength(line1[5]), 
					contigs.get(line1[2]), Integer.parseInt(line1[4]), line1[5]);
			tmp.addRead(left2, isReverse(line2[1]), cigarLength(line2[5]), 
					contigs.get(line2[2]), Integer.parseInt(line2[4]), line2[5]);
			reads.put(line1[0], tmp);
			i++;
			// sample HANDFUL_SIZE pairs from this point
			j = 1;
			while (j < HANDFUL_SIZE && raf.getFilePointer() < raf.length()){
				line1 = raf.readLine().split("\t");
				line2 = raf.readLine().split("\t");
				left1 = Integer.parseInt(line1[3]);
				left2 = Integer.parseInt(line2[3]);
				// if pair didn't map, just jump to next
				// line instead of jumping a full step
				if (left1 == 0 || left2 == 0) continue;
				tmp = new ReadPair(line1[0]);
				tmp.addRead(left1, isReverse(line1[1]), cigarLength(line1[5]), 
						contigs.get(line1[2]), Integer.parseInt(line1[4]), line1[5]);
				tmp.addRead(left2, isReverse(line2[1]), cigarLength(line2[5]), 
						contigs.get(line2[2]), Integer.parseInt(line2[4]), line2[5]);
				reads.put(line1[0], tmp);
				i++;
				j++;
			}
			pos += step;
			raf.seek(pos);
		}
		return reads;
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
	
	
	private static class PrintStreamPair implements Comparable<PrintStreamPair>{
		// a comparator for sorting reads by position
		private static final Comparator<ReadPair> COMP = new Comparator<ReadPair>(){
			public int compare(ReadPair arg0, ReadPair arg1) {
				if (arg0.pos1 == arg1.pos1){
					if (arg0.pos2 == arg1.pos2)
						return 0;
					else if (arg0.pos2 < arg1.pos2)
						return -1;
					else
						return 1;
				} else if (arg0.pos1 < arg1.pos1)
					return -1;
				else 
					return 1;
			} 
			
		};
		private static final int MAX_SB_SIZE = 250;
		private File fishDir;
		// a PrintStream and File for both directions
		private PrintStream out1; 
		private PrintStream out2;
		File file1;
		File file2;
		private Contig ctg1;
		private Contig ctg2;
		private TreeSet<ReadPair> reads;
		private StringBuilder sb1;
		private StringBuilder sb2;
		
		public PrintStreamPair(Contig ctg1, Contig ctg2, File outdir) throws IOException{
			fishDir = outdir;
			reads = new TreeSet<ReadPair>(COMP);
			if (ctg1==ctg2){
				file1 = new File(fishDir,"match."+ctg1.getRank()+"v"+ctg2.getRank()+".txt");
				file1.createNewFile();
				this.ctg1 = ctg1;
				this.ctg2 = ctg2;
			} else { 
				file1 = new File(fishDir,"match."+ctg1.getRank()+"v"+ctg2.getRank()+".txt");
				file1.createNewFile();
				file2 = new File(fishDir,"match."+ctg2.getRank()+"v"+ctg1.getRank()+".txt");
				file2.createNewFile();
				this.ctg1 = ctg1;
				this.ctg2 = ctg2;
			}
			sb1 = new StringBuilder();
			sb2 = new StringBuilder();
		}
		
		public void add(ReadPair pair){
			reads.add(pair);
		}
		
		/**
		 * Assumes pos1 is from ctg1 and pos2 is from ctg2
		 * 
		 * @param pos1 position on <code>ctg1</code> from constructor 
		 * @param pos2 position on <code>ctg2</code> from constructor
		 * @param qual
		 */
		public void print(int pos1, int pos2, int qual) {
			// phred-scaled geometric-mean of posterior probabilites that mapping position is incorrect
			// mod 255 in case the quality score couldn't be computed.	
			sb1.append("c"+this.ctg1.getRank()+"p"+pos1+"\tc"+this.ctg2.getRank()+"p"+pos2+"\t"+qual+"\n");
			if (this.ctg1 != this.ctg2)
				sb2.append("c"+this.ctg2.getRank()+"p"+pos2+"\tc"+this.ctg1.getRank()+"p"+pos1+"\t"+qual+"\n"); 
			// flush the buffer 
			if (sb1.length() >= MAX_SB_SIZE || sb2.length() >= MAX_SB_SIZE){
				try {
					BufferedWriter bw = new BufferedWriter(new FileWriter(file1, true));
					bw.write(sb1.toString());
					sb1.delete(0, sb1.length());
					bw.close();
					if (this.ctg1 != this.ctg2) {
						bw = new BufferedWriter(new FileWriter(file2, true));
						bw.write(sb2.toString());
						sb2.delete(0, sb2.length());
						bw.close();
					}
				} catch (IOException e) {
					e.printStackTrace();
					System.exit(-1);
				}
			}
		}
		
		public void close(){
			if (sb1.length() > 0 || sb2.length() > 0){
				try {
					BufferedWriter bw = new BufferedWriter(new FileWriter(file1, true));
					bw.write(sb1.toString());
					sb1.delete(0, sb1.length());
					bw.close();
					if (ctg1 != ctg2) {
						bw = new BufferedWriter(new FileWriter(file2, true));
						bw.write(sb2.toString());
						sb2.delete(0, sb2.length());
						bw.close();
					}
				} catch (IOException e) {
					e.printStackTrace();
					System.exit(-1);
				}
			}
		/*
			if (reads.size() == 0)
				return;
			if (ctg1 == ctg2){
				file1.createNewFile();
				out1 = new PrintStream(file1);
			} else {
				file1.createNewFile();
				file2.createNewFile();
				out1 = new PrintStream(file1);
				out2 = new PrintStream(file2);
			}
			Iterator<ReadPair> it = reads.iterator();
			while(it.hasNext())
				print(it.next());
			
			out1.close();
			if (ctg1 != ctg2)
				out2.close();
		*/
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
	
	static class MatchFile implements Comparable<MatchFile>{
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
	
	static class ReadPoint implements Comparable<ReadPoint>{
		String hdr;
		int pos;
		public ReadPoint(String hdr, int pos){
			this.hdr = hdr;
			this.pos = pos;
		}
		public int compareTo(ReadPoint rp){
			return this.pos - rp.pos;
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
