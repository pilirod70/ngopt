package org.halophiles.assembly.qc;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Vector;

import org.halophiles.assembly.Contig;
import org.halophiles.assembly.ReadPair;
import org.halophiles.assembly.SAMFileParser;
import org.halophiles.tools.SummaryStats;

import weka.estimators.KernelEstimator;

public class FISHInputExporter {
	private static NumberFormat NF;
	
	public static void main(String[] args){
		if (args.length != 2){
			System.err.println("Usage: java -jar GetFishInput.jar <sam_file> <output_base>");
			System.exit(-1);
		}
		try{
			NF = NumberFormat.getInstance();
			NF.setMaximumFractionDigits(0);
			NF.setGroupingUsed(false);
			System.out.println("Reading "+args[0]);
			System.out.println("Writing output to "+System.getProperty("user.dir"));
			SAMFileParser sfp = new SAMFileParser(args[0]);
			String base = args[1];
			File outdir = new File(System.getProperty("user.dir"));
			Map<String,Contig> contigs = new HashMap<String,Contig>();
			Iterator<Contig> ctgIt = sfp.getContigs();
			while(ctgIt.hasNext()){
				Contig tmp = ctgIt.next();
				contigs.put(tmp.name, tmp);
			}
			Map<String,ReadPair> reads = new HashMap<String,ReadPair>();
			Iterator<ReadPair> rpIt = sfp.getReadPairs();
			while(rpIt.hasNext()){
				ReadPair tmp = rpIt.next();
				if (!tmp.paired) 
					continue;
				reads.put(tmp.hdr, tmp);
			}
			
			System.out.println("Found "+contigs.size()+" contigs.");
			System.out.println("Found "+sfp.getNumReads()+" total reads and "+sfp.getNumPairs()+" read-pairs.");
			if (reads.size() <= 0) {
				System.err.println("No paired reads found. Cannot generate match files for running FISH misassembly detection.");
				System.exit(-1);
			}
			
			System.out.print("Estimating insert size... ");
			double[] ins = estimateInsertSizeIQR(reads);
			System.out.println("mean insert: " + NF.format(ins[0])+"    stdev: " + NF.format(ins[1]));
			int nSd = 3;
			File insFile = new File(outdir, base+".ins_size_unfilt.txt");
			insFile.createNewFile();
			PrintStream insOut = new PrintStream(insFile);
			exportInsertSize(reads, insOut);
			insOut.close();
			
			
			if (ins[0] > 1000) { // if we have a mated library, filter out any small insert shadow library
				System.out.print("Looks like this library was mated. Filtering discordant read pairs... ");
				int before = reads.size();
				filterDiscordantPairs(reads);
				System.out.println("removed "+(before - reads.size())+" of "+before +" read pairs.");
				System.out.print("Re-estimating insert size... ");
				ins = estimateInsertSize(reads);
				System.out.println("mean insert: " + NF.format(ins[0])+"    stdev: " + NF.format(ins[1]));
			}
			System.out.print("Filtering diagonal self-links. Discarding pairs with inserts between "+NF.format(ins[0]-nSd*ins[1])+" - "+NF.format(ins[0]+nSd*ins[1])+"... ");			
			int before = reads.size();
			filterDiagReads(reads,ins,nSd);
			System.out.println("removed "+(before - reads.size())+" of "+before +" read pairs.");
			
			
			insFile = new File(outdir,base+".ins_size_filt.txt");
			insFile.createNewFile();
			insOut = new PrintStream(insFile);
			exportInsertSize(reads, insOut);
			insOut.close();
			export(contigs,reads,outdir,base);
			
		} catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		} catch(Exception e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public static void export(Map<String,Contig> contigs, Map<String,ReadPair> reads, File fishDir, String base) throws IOException{
		HashMap<String,Vector<Double>> mapPoints = new HashMap<String,Vector<Double>>();
		HashMap<String,PrintStreamPair> matchesOut = new HashMap<String,PrintStreamPair>();
		Vector<Double> tmpPts = null;
		
		File controlFile = new File(fishDir,base+".control.txt");
		PrintStream ctrlOut = new PrintStream(controlFile);
		Iterator<ReadPair> rpIt = reads.values().iterator();
		while(rpIt.hasNext()){
			ReadPair tmpRp = rpIt.next();
			Contig ctg = tmpRp.ctg1;
			boolean rev = tmpRp.rev1;
			int left = tmpRp.pos1;
			if (mapPoints.containsKey(ctg.name)){
				tmpPts = mapPoints.get(ctg.name);
			} else {
				tmpPts = new Vector<Double>();
				mapPoints.put(ctg.name, tmpPts);
			}
			if (rev){
				left = left * -1;
			}
			tmpPts.add(new Double(left));
			
			ctg = tmpRp.ctg2;
			rev = tmpRp.rev2;
			left = tmpRp.pos2;
			if (mapPoints.containsKey(ctg.name)){
				tmpPts = mapPoints.get(ctg.name);
			} else {
				tmpPts = new Vector<Double>();
				mapPoints.put(ctg.name, tmpPts);
			}
			if (rev){
				left = left * -1;
			}
			tmpPts.add(new Double(left));
		}
		// Rank contigs so FISH doesn't freak out about non-consecutive ids
		Contig[] ctgRef = new Contig[contigs.values().size()];
		contigs.values().toArray(ctgRef);
		Arrays.sort(ctgRef,new Comparator<Contig>(){
			public int compare(Contig o1, Contig o2) {
				return o1.getId() - o2.getId();
			}
		});
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
		ctrlOut.println("-maps");
		while(it.hasNext()){
			String ctg = it.next();
	//		if (contigs.get(ctg).getId() != 1 && contigs.get(ctg).getId() != 2)
	//			continue;
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
//			maps.put(contigs.get(ctg).getRank(),ctgMapFile.getParentFile().getName()+"/"+ctgMapFile.getName());
			maps.put(contigs.get(ctg).getRank(),ctgMapFile.getAbsolutePath());
		}
		ctgLblOut.close();
		Iterator<Integer> it2 = maps.keySet().iterator();
		while(it2.hasNext()){
			Integer tmpInt = it2.next();
			ctrlOut.println(tmpInt+"\t"+maps.get(tmpInt));
		}
		
		it = reads.keySet().iterator();
		PrintStreamPair matchOut = null;
		ctrlOut.println("-matches");
		String ctgStr;
		TreeSet<MatchFile> matchFiles = new TreeSet<MatchFile>();
		while(it.hasNext()) {
			ReadPair tmp = reads.get(it.next());
			/*if (!tmp.paired || !((tmp.ctg1.getId()== 1 && tmp.ctg2.getId()== 2) || 
					(tmp.ctg1.getId()== 2 && tmp.ctg2.getId()== 2) || 
					(tmp.ctg1.getId()== 1 && tmp.ctg2.getId()== 1)))
					continue;*/
			
			if (tmp.ctg1 == null || tmp.ctg2 == null) 
				continue;
			ctgStr = tmp.contigString();
			if (matchesOut.containsKey(ctgStr))
				matchOut = matchesOut.get(ctgStr);
			else {
				matchOut = new PrintStreamPair(tmp.ctg1,tmp.ctg2,mapMatchDir);
				matchOut.addMatchFile(matchFiles);
				matchesOut.put(ctgStr, matchOut);
			}
			matchOut.add(tmp);
		}
		
		Iterator<MatchFile> mfIt = matchFiles.iterator();
		while(mfIt.hasNext())
			ctrlOut.println(mfIt.next().toString());
		ctrlOut.println("end");
		ctrlOut.close();
	//	System.out.println(matchesOut.size()+" PrintStreamPairs for " + contigs.size() + " contigs");
		it = matchesOut.keySet().iterator(); 
		while(it.hasNext())
			matchesOut.get(it.next()).close();
	}
	
	
	public static double[] estimateInsertSize(Map<String,ReadPair> reads){
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
		double[] dat = new double[arD.length];
		for (int i = 0; i < dat.length; i++)
			dat[i] = arD[i];
		double mean = SummaryStats.mean(dat);
		double stdev = Math.sqrt(SummaryStats.variance(dat,mean));
		double[] ret = {Math.round(mean),Math.round(stdev),dat.length};
		return ret;
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
		double[] dat = new double[arD.length/2];
		int j = arD.length/4;
		for (int i = 0; i < dat.length; i++)
			dat[i] = arD[j++];
		double mean = SummaryStats.mean(dat);
		double stdev = Math.sqrt(SummaryStats.variance(dat,mean));
		double[] ret = {Math.round(mean),Math.round(stdev),dat.length};
		return ret;
	}
	
	private static double[] getInsertData(Map<String,ReadPair> reads){
		Iterator<ReadPair> it = reads.values().iterator();
		ReadPair tmp = null;
		Vector<Double> vals = new Vector<Double>();
		while(it.hasNext()){
			tmp = it.next();
			if (tmp.paired && tmp.ctg1.equals(tmp.ctg2))
				vals.add(new Double(tmp.getInsert()));
		}
		Double[] arD = vals.toArray(new Double[vals.size()]);
		double[] ret = new double[arD.length];
		for (int i = 0; i < ret.length; i++)
			ret[i] = arD[i];
		return ret;
	}
	
	public static boolean isDiag(ReadPair r, double[] ins, int nSd){
		double dist = Math.abs(r.pos1 - r.pos2 );
		if (r.ctg2 == null){
			System.out.print("");
		}
		return (dist < ins[0]+nSd*ins[1] && dist > ins[0]-nSd*ins[1]) && r.ctg1.equals(r.ctg2);
	}
	
	
	public static void filterDiagReads(Map<String,ReadPair> reads, double[] ins, int nSd){
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
	}
	
	public static void filterDiscordantPairs(Map<String,ReadPair> reads){
		Vector<String> outies = new Vector<String>();
		Vector<String> inies = new Vector<String>();
		Iterator<String> it = reads.keySet().iterator();
		String tmp = null;
		ReadPair read = null;
		int nsyn = 0;
		while(it.hasNext()){
			tmp = it.next();
			read = reads.get(tmp);
			if (!read.paired) continue;
			if (!read.ctg1.equals(read.ctg2)) continue;
			if (!read.inward && !read.outward)
				continue;
			nsyn++;
			if (read.outward){
				outies.add(tmp);
			} else {
				inies.add(tmp);
			}
		}
		double nOut = outies.size();
		double nIn = inies.size();
		nOut = nOut/nsyn;
		nIn = nIn/nsyn;
		if (nIn > nOut){ // have mostly inies
			if (nIn > 0.5){
				removeKeys(reads,outies);
			}
		} else if (nOut > nIn){ // have mostly outies
			if (nOut > 0.5)
				removeKeys(reads,inies);
		}
	}
	
	private static void removeKeys(Map<String,ReadPair> reads, Vector<String> torm){
		Iterator<String> it = torm.iterator();
		while(it.hasNext())
			reads.remove(it.next());
	}
	
	public static void exportInsertSize(Map<String,ReadPair> reads, PrintStream out) {
		Iterator<ReadPair> it = reads.values().iterator();
		ReadPair tmp = null;
		while(it.hasNext()){
			tmp = it.next();
			if (tmp.paired && tmp.ctg1.equals(tmp.ctg2))
				out.println(tmp.getInsert());
		}
	}
	
	private static class PrintStreamPair{
		private File fishDir;
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
							return arg0.hdr.compareTo(arg1.hdr);
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
				file1.createNewFile();
				this.ctg1 = ctg1;
				this.ctg2 = ctg2;
			} else { 
				file1 = new File(fishDir,"match."+ctg1.getRank()+"v"+ctg2.getRank()+".txt");
				file2 = new File(fishDir,"match."+ctg2.getRank()+"v"+ctg1.getRank()+".txt");
				file1.createNewFile();
				file2.createNewFile();
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
			int qual = ((Integer.parseInt(pair.sam1[4]) % 255)+ (Integer.parseInt(pair.sam2[4]) % 255))/2;
			if ((ctg1 == pair.ctg1 && ctg2 == pair.ctg2)) {
				out1.println("c"+ctg1.getRank()+"p"+pair.pos1+"\tc"+ctg2.getRank()+"p"+pair.pos2+"\t"+qual);
				if (ctg1 != ctg2)
					out2.println("c"+ctg2.getRank()+"p"+pair.pos2+"\tc"+ctg1.getRank()+"p"+pair.pos1+"\t"+qual);
			} else if((ctg1 == pair.ctg2 && ctg2 == pair.ctg1)){
				out1.println("c"+ctg1.getRank()+"p"+pair.pos2+"\tc"+ctg2.getRank()+"p"+pair.pos1+"\t"+qual);
				out2.println("c"+ctg2.getRank()+"p"+pair.pos1+"\tc"+ctg1.getRank()+"p"+pair.pos2+"\t"+qual);
			} else {
				throw new IllegalArgumentException("pair "+pair.hdr+" does not connect contigs "+ctg1.name+" and "+ctg2.name);
			}
		}
		
		public void close() throws FileNotFoundException{
			if (ctg1==ctg2){
				out1 = new PrintStream(file1);
			} else {
				out1 = new PrintStream(file1);
				out2 = new PrintStream(file2);
			}
			Iterator<ReadPair> it = reads.iterator();
			ReadPair pair = null;
			while(it.hasNext())
				print(it.next());
			
			out1.close();
			if (ctg1 != ctg2)
				out2.close();
		}

		
		public void addMatchFile(Set<MatchFile> mfiles){
			mfiles.add(new MatchFile(ctg1, ctg2, file1));
			if (!ctg1.equals(ctg2))
				mfiles.add(new MatchFile(ctg2, ctg1, file2));
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
//			return ctg1.getRank()+"\t"+ctg2.getRank()+"\t"+file.getParentFile().getName()+"/"+file.getName();
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

}
