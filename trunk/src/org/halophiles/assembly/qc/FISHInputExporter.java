package org.halophiles.assembly.qc;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
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
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Vector;

import org.halophiles.assembly.Contig;
import org.halophiles.assembly.ReadPair;
import org.halophiles.assembly.SAMFileParser;
import org.halophiles.tools.SummaryStats;

public class FISHInputExporter {
	private static NumberFormat NF;
	private static int MINQUAL = 13;
	
	public static void main(String[] args){
		if (args.length < 3){
			System.err.println("Usage: java -jar GetFishInput.jar <sam_file> <output_base> <output_dir> <min_map_qual (optional)>");
			System.exit(-1);
		}
		try{
			NF = NumberFormat.getInstance();
			NF.setMaximumFractionDigits(0);
			NF.setGroupingUsed(false);
			if (args.length == 4){
				MINQUAL = Integer.parseInt(args[2]);
			}
			String base = args[1];
			File outdir = new File(args[2]);
			System.out.println("Reading "+args[0]);
			System.out.println("Writing output to "+outdir.getAbsolutePath()+" with basename "+base);
			SAMFileParser sfp = new SAMFileParser(args[0]);

			System.out.println("Found "+sfp.getNumContigs()+" contigs.");
			System.out.println("Found "+sfp.getNumReads()+" total reads and "+sfp.getNumPairs()+" read-pairs.");
			
			Map<String,ReadPair> reads = new HashMap<String,ReadPair>();
			System.out.println("Filtering duplicate connections...");
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
			if (reads.size() <= 0) {
				System.err.println("No paired reads found. Cannot generate match files for running FISH misassembly detection.");
				System.exit(-1);
			}
			File insFile = new File(outdir, base+".ins_size_unfilt.txt");
			insFile.createNewFile();
			PrintStream insOut = new PrintStream(insFile);
			ReadPair.exportInsertSize(reads.values(), insOut);
			insOut.close();
	/*		
			filterShadowLibrary(reads);
			System.out.print("Estimating insert size... ");
			int nSd = 3;
			double[] ins = ReadPair.estimateInsertSize(reads.values());
			System.out.println("mean insert: " + NF.format(ins[0])+"    stdev: " + NF.format(ins[1]));
			if (ins[1] > ins[0]){
				System.out.print("Insert sizes are overdispersed. Discarding upper and lower 10 percent and recomputing mean... ");
				double[] tmpIns = estimateInsertSizeIQR(reads);
				ins = tmpIns;
				System.out.println("mean insert: " + NF.format(ins[0])+"    stdev: " + NF.format(ins[1]));				
			}
			System.out.print("Filtering diagonal self-links. Discarding pairs with inserts between "+NF.format(ins[0]-nSd*ins[1])+" - "+NF.format(ins[0]+nSd*ins[1])+"... ");			
			int before = reads.size();
			filterDiagReads(reads,ins,nSd);
			System.out.println("removed "+(before - reads.size())+" of "+before +" read pairs.");
			
	*/
			System.out.println("Filtering proper connections... ");
			int before = reads.size();
			filteProperConnections(reads);
			System.out.println("removed "+(before - reads.size())+" of "+before +" read pairs.");
			
			System.out.print("Filtering tandem connections... ");
			before = reads.size();
			filterTandemConnections(reads);
			System.out.println("removed "+(before - reads.size())+" of "+before +" read pairs.");
			
			insFile = new File(outdir,base+".ins_size_filt.txt");
			insFile.createNewFile();
			insOut = new PrintStream(insFile);
			ReadPair.exportInsertSize(reads.values(), insOut);
			insOut.close();
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
		Iterator<Contig> ctgIt = contigs.values().iterator();
		Set<Contig> singles = new TreeSet<Contig>();
		while(ctgIt.hasNext()){
			tmpCtg = ctgIt.next();
			if (tmpCtg.numReads() <= 1)
				singles.add(tmpCtg);
		}
		
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
			maps.put(contigs.get(ctg).getRank(),ctgMapFile.getAbsolutePath());
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
			if (tmpPSP.reads.size() <= 1)
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
		while(mfIt.hasNext())
			ctrlOut.println(mfIt.next().toString());
		
		ctrlOut.println("end");
		ctrlOut.close();
	}
	
	

	
	public static double[] estimateInsertSizeIQR(Map<String,ReadPair> reads){
		Iterator<ReadPair> it = reads.values().iterator();
		ReadPair tmp = null;
		Vector<Double> vals = new Vector<Double>();
		while(it.hasNext()){
			tmp = it.next();
			if (tmp.paired && tmp.ctg1.equals(tmp.ctg2))
				vals.add(new Double(tmp.getInsert()));
		}
		Double[] arD = vals.toArray(new Double[vals.size()]);
		Arrays.sort(arD);
		// discard the upper and lower quartiles to prevent any outliers from distorting our estimate
		double[] dat = new double[(arD.length*8)/10];
		int j = arD.length/10;
		for (int i = 0; i < dat.length; i++)
			dat[i] = arD[j++];
		double mean = SummaryStats.mean(dat);
		double stdev = Math.sqrt(SummaryStats.variance(dat,mean));
		double[] ret = {Math.round(mean),Math.round(stdev),dat.length};
		return ret;
	}
	
	/**
	 * 
	 * @param reads
	 * @param ins
	 * @param nSd
	 * @return a reference to the Map that was filtered
	 */
	public static Map<String,ReadPair> filterDiagReads(Map<String,ReadPair> reads, double[] ins, int nSd){
		Iterator<String> it = reads.keySet().iterator();
		Vector<String> rm = new Vector<String>();
		String key = null;
		ReadPair r = null;
		while(it.hasNext()){
			key = it.next();
			r = reads.get(key);
			if (isDiag(r, ins, nSd)){
				rm.add(key);
			}
		}
		removeKeys(reads,rm);
		return reads;
	}
	
	/**
	 * 
	 * @param reads
	 * @param ins
	 * @param nSd
	 * @return a reference to the Map that was filtered
	 */
	public static Map<String,ReadPair> filterTailPairs(Map<String,ReadPair> reads, double alpha){
		ReadPair[] ar = new ReadPair[reads.size()];
		reads.values().toArray(ar);
		Arrays.sort(ar,new Comparator<ReadPair>() {
			public int compare(ReadPair o1, ReadPair o2) {
				return o1.getInsert() - o2.getInsert();
			}
		});
		
		Vector<String> rm = new Vector<String>();
		int l = (int) ((alpha/2) * ar.length);
		int r = (int) ((1-(alpha/2)) * ar.length);
		for (int i = 0; i < l; i++)
			rm.add(ar[i].hdr);
		for (int i = r; i < ar.length; i++)
			rm.add(ar[i].hdr);
		removeKeys(reads,rm);
		return reads;
	}
	/*
	public static void filterShadowLibrary(Map<String,ReadPair> reads){
		Map<String,ReadPair> outies = new HashMap<String, ReadPair>();
		Map<String,ReadPair> inies = new HashMap<String, ReadPair>();
		Iterator<String> it = reads.keySet().iterator();
		String tmp = null;
		ReadPair read = null;
		int nsyn = 0;
		System.out.print("Checking for shadow library... ");
		while(it.hasNext()){
			tmp = it.next();
			read = reads.get(tmp);
			if (!read.paired) continue;
			if (!read.ctg1.equals(read.ctg2)) continue;
			if (!read.inward && !read.outward)
				continue;
			nsyn++;
			if (read.outward){
				outies.put(tmp, reads.get(tmp));
			} else {
				inies.put(tmp, reads.get(tmp));
			}
		}
		double nOut = outies.size();
		double nIn = inies.size();
		System.out.println(inies.size()+" innies "+outies.size()+" outies");
		nOut = nOut/nsyn;
		nIn = nIn/nsyn;
		double[] inDat = new double[inies.size()];
		double[] outDat = new double[outies.size()];
		it = inies.keySet().iterator();
		for (int i = 0; i < inDat.length; i++)
			inDat[i] = reads.get(it.next()).getInsert();
		it = outies.keySet().iterator();
		for (int i = 0; i < outDat.length; i++)
			outDat[i] = reads.get(it.next()).getInsert();
		double[] inIQR = estimateInsertSizeIQR(inies);
		double[] outIQR = estimateInsertSizeIQR(outies);
		double inMean = SummaryStats.mean(inDat);
		double inSd = Math.sqrt(SummaryStats.variance(inDat, inMean));
		double outMean = SummaryStats.mean(outDat);
		double outSd = Math.sqrt(SummaryStats.variance(outDat, outMean));
		System.out.println("innies: mu="+NF.format(inMean)+" sd="+NF.format(inSd));
		System.out.println("outies: mu="+NF.format(outMean)+" sd="+NF.format(outSd));
		System.out.println("IQR:");
		System.out.println("innies: mu="+NF.format(inIQR[0])+" sd="+NF.format(Math.sqrt(inIQR[1])));
		System.out.println("outies: mu="+NF.format(outIQR[0])+" sd="+NF.format(Math.sqrt(outIQR[1])));
		
		if (nOut > 0.1 && nIn > 0.1){ // we have a shadow library
			System.out.print("Found shadow library. Filtering shadow reads... ");
			int nRemoved = -1;
			//double inSd = Math.sqrt(SummaryStats.variance(inDat, inIQR[0]));
			inSd = Math.sqrt(SummaryStats.variance(inDat, inIQR[0]));
			//double outSd = Math.sqrt(SummaryStats.variance(outDat, outIQR[0]));
			outSd = Math.sqrt(SummaryStats.variance(outDat, outIQR[0]));
			double inCv = inSd/inIQR[0];
			double outCv = outSd/outIQR[0];
			if (inIQR[0] > outIQR[0]){
				if (inIQR[0] > 2*outIQR[0] || outCv > inCv) { // remove outies
					filterOffDiagReads(outies, outIQR, 3);
					nRemoved = outies.size();
					removeKeys(reads,outies.keySet());
				}
			} else if (outIQR[0] > inIQR[0]) { 
				if (outIQR[0] > 2*inIQR[0] || inCv > outCv){ // remove innies
					filterOffDiagReads(inies, inIQR, 3);
					nRemoved = inies.size();
					removeKeys(reads,inies.keySet());
				}
			} else {
				System.out.println("Found innie and outie libraries with same insert.");
			}
			System.out.println(" removed "+nRemoved+" reads");
		} else {
			System.out.println(" unable to detect a shadow library.");
		}
			
	}
	*/
	@SuppressWarnings("unchecked")
	private static Collection<ReadPair>[] filteProperConnections(Map<String,ReadPair> reads){
		int K = 2;
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
		if (ins[0] > 1500){
			System.out.print("EM-clustering insert sizes with K=3 (mean ins > 1500bp)... ");
			K = 3;
		} else {
			System.out.print("EM-clustering insert sizes with K=2... ");
		}
		EMClusterer em = new EMClusterer(toFilt, K);
		int iters = em.iterate(1000, 0.0001);
		System.out.println("stopping after "+iters+" iterations.");
		Collection<ReadPair>[] clusters = new Vector[K];
		em.getClusters().toArray(clusters);
		double[][] allIns = new double[K][4];
		for (int i = 0; i < clusters.length; i++){
			allIns[i][0] = i;
			ins = ReadPair.estimateInsertSize(clusters[i]);
			allIns[i][1] = ins[0];
			allIns[i][2] = ins[1];
			allIns[i][3] = ins[2]/toFilt.size();
		}
		Arrays.sort(allIns, new Comparator<double[]>(){
			public int compare(double[] o1, double[] o2) {
				return (int) (o1[1] - o2[1]);
			}
		});
		Vector<String> toRm = new Vector<String>();
		Map<String,ReadPair> map = null;
		ReadPair tmpRp  = null;
		String rmClusters = "";
		double alpha = 0.01;
		for (int i = 0; i < allIns.length; i++){
			System.out.println("cluster"+NF.format(allIns[i][0])+": mu="+pad(NF.format(allIns[i][1]),10)+"sd="+pad(NF.format(allIns[i][2]),10)+"n="+clusters[(int)allIns[i][0]].size());
			// if we have a small and very high variance, but not overdispersed cluster, don't discard 
			if(allIns[i][2]/allIns[i][1] >= 0.90 && allIns[i][3] < 0.05)
				continue;
			// if insert size distribution is under dispersed, filter tails so we don't throw them out
			if (allIns[i][2]/allIns[i][1] <= 1){
				map = new HashMap<String,ReadPair>();
				for (Iterator<ReadPair> it = clusters[(int)allIns[i][0]].iterator(); it.hasNext();){
					tmpRp = it.next();
					map.put(tmpRp.hdr, tmpRp);
				}
				ins[0] = allIns[i][1];
				ins[1] = allIns[i][2];
				filterTailPairs(map,alpha);
				toRm.addAll(map.keySet());
				if (rmClusters.length()>0)
					rmClusters = rmClusters+" "+NF.format(allIns[i][0])+"("+map.keySet().size()+")";
				else 
					rmClusters = NF.format(allIns[i][0])+"("+map.keySet().size()+")";
					
			} 
		}
		if (rmClusters.length()>0)
			System.out.println("Removed but keeping tails from clusters "+rmClusters);
		
		removeKeys(reads, toRm);
		return clusters;
	}
	
	private static String pad(String s, int len){
		String ret = new String(s);
		for (int i = 0; i < len-s.length(); i++){
			ret = ret+" ";
		}
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
	
	
	private static boolean isDiag(ReadPair r, double[] ins, int nSd){
		double dist = Math.abs(r.pos1 - r.pos2 );
		if (r.ctg2 == null){
			System.out.print("");
		}
		return (dist < ins[0]+nSd*ins[1] && dist > ins[0]-nSd*ins[1]) && r.ctg1.equals(r.ctg2);
	}
	
	private static <V> void removeKeys(Map<String,V> reads, Collection<String> torm){
		Iterator<String> it = torm.iterator();
		while(it.hasNext())
			reads.remove(it.next());
	}

	private static class PrintStreamPair implements Comparable<PrintStreamPair>{
		private File fishDir;
		// a PrintStream and File for both directions
		private PrintStream out1; 
		private PrintStream out2;
		File file1;
		File file2;
		private Contig ctg1;
		private Contig ctg2;
		private TreeSet<ReadPair> reads;
		
		public PrintStreamPair(Contig ctg1, Contig ctg2, File outdir) throws IOException{
			fishDir = outdir;
			reads = new TreeSet<ReadPair>(new Comparator<ReadPair>(){
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
				
			});
			if (ctg1==ctg2){
				file1 = new File(fishDir,"match."+ctg1.getRank()+"v"+ctg2.getRank()+".txt");
				//file1.createNewFile();
				this.ctg1 = ctg1;
				this.ctg2 = ctg2;
			} else { 
				file1 = new File(fishDir,"match."+ctg1.getRank()+"v"+ctg2.getRank()+".txt");
				file2 = new File(fishDir,"match."+ctg2.getRank()+"v"+ctg1.getRank()+".txt");
				//file1.createNewFile();
				//file2.createNewFile();
				this.ctg1 = ctg1;
				this.ctg2 = ctg2;
			}
		}
		
		public void add(ReadPair pair){
			reads.add(pair);
		}
		
		private void print(ReadPair pair){
			// phred-scaled geometric-mean of posterior probabilites that mapping position is incorrect
			// mod 255 in case the quality score couldn't be computed.	
			if ((ctg1 == pair.ctg1 && ctg2 == pair.ctg2)) {
				out1.println("c"+ctg1.getRank()+"p"+pair.pos1+"\tc"+ctg2.getRank()+"p"+pair.pos2+"\t"+pair.getQual());
				if (ctg1 != ctg2)
					out2.println("c"+ctg2.getRank()+"p"+pair.pos2+"\tc"+ctg1.getRank()+"p"+pair.pos1+"\t"+pair.getQual());
			} else if((ctg1 == pair.ctg2 && ctg2 == pair.ctg1)){
				out1.println("c"+ctg1.getRank()+"p"+pair.pos2+"\tc"+ctg2.getRank()+"p"+pair.pos1+"\t"+pair.getQual());
				out2.println("c"+ctg2.getRank()+"p"+pair.pos1+"\tc"+ctg1.getRank()+"p"+pair.pos2+"\t"+pair.getQual());
			} else {
				throw new IllegalArgumentException("pair "+pair.hdr+" does not connect contigs "+ctg1.name+" and "+ctg2.name);
			}
		}
		
		public void close() throws IOException{
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
			return ctg1.getRank()+"\t"+ctg2.getRank()+"\t"+file.getAbsolutePath();
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

}
