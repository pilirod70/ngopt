package org.halophiles.assembly.qc;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Vector;

import org.halophiles.assembly.Contig;
import org.halophiles.assembly.ReadPair;
import org.halophiles.assembly.SAMFileParser;
import org.halophiles.tools.SummaryStats;

public class FISHInputExporter {
	
	public static void main(String[] args){
		if (args.length != 2){
			System.err.println("Usage: java FISHInputExporter <sam_file> <output_dir>");
			System.exit(-1);
		}
		try{
			SAMFileParser sfp = new SAMFileParser(args[0]);
			File outdir = new File(args[1]);
			outdir.mkdirs();
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
				reads.put(tmp.hdr, tmp);
			}
			double[] ins = estimateInsertSize(reads);
			filterDiagReads(reads,ins,3);
			export(contigs,reads,outdir);
			
		} catch(Exception e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public static void export(Map<String,Contig> contigs, Map<String,ReadPair> reads, File fishDir) throws IOException{
		HashMap<String,Vector<Double>> mapPoints = new HashMap<String,Vector<Double>>();
		HashMap<String,PrintStreamPair> matchesOut = new HashMap<String,PrintStreamPair>();
		Vector<Double> tmpPts = null;
		
		File controlFile = new File(fishDir,"control.txt");
		PrintStream ctrlOut = new PrintStream(controlFile);
		Iterator<ReadPair> rpIt = reads.values().iterator();
		while(rpIt.hasNext()){
			ReadPair tmpRp = rpIt.next();
			if (!tmpRp.paired) 
				continue;
			
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
		
		File mapMatchDir = new File(fishDir,"dat");
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
			File ctgMapFile = new File(mapMatchDir,"map."+contigs.get(ctg).getId()+".txt");
			PrintStream ctgMapOut = new PrintStream(ctgMapFile);
			
			Arrays.sort(ar, new Comparator<Double>(){
				public int compare(Double arg0, Double arg1) {
					if (Math.abs(arg0) < Math.abs(arg1)) return -1;
					else if (Math.abs(arg0) > Math.abs(arg1)) return 1;
					else return 0;
				}	
			});
			for (Double d: ar){
				ctgMapOut.println("c"+contigs.get(ctg).getId()+"p"+Math.abs(d.intValue())+"\t"+(d < 0 ? -1 : 1));
//				ctgMapOut.println("c"+contigs.get(ctg).getId()+"p"+Math.abs(d.intValue())+"\t"+1);
			}
			ctgMapOut.close();
			//String dir = ctgMapFile.getParentFile().getName();
			//ctrlOut.println(contigs.get(ctg).getId()+"\t"/*+dir+"/"*/+ctgMapFile.getName());
			maps.put(contigs.get(ctg).getId(),ctgMapFile.getParentFile().getName()+"/"+ctgMapFile.getName());
		}
		
		Iterator<Integer> it2 = maps.keySet().iterator();
		while(it2.hasNext()){
			Integer tmpInt = it2.next();
			ctrlOut.println(tmpInt+"\t"+maps.get(tmpInt));
		}
		
//		int rdlen = sfp.getReadLength();
		it = reads.keySet().iterator();
		HashSet<ReadPair> conflicts = new HashSet<ReadPair>();
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
				//matchOut.printControl(ctrlOut);
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
		System.err.println(matchesOut.size()+" PrintStreamPairs for " + contigs.size() + " contigs");
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
				vals.add(new Double(Math.abs(tmp.pos1-tmp.pos2)));
		}
		Double[] arD = vals.toArray(new Double[vals.size()]);
		double[] dat = new double[arD.length];
		for (int i = 0; i < dat.length; i++)
			dat[i] = arD[i];
		double mean = SummaryStats.mean(dat);
		double stdev = Math.sqrt(SummaryStats.variance(dat,mean));
		
		double[] ret = {mean,stdev};
		return ret;
	}
	
	public static boolean isDiag(ReadPair r, double[] ins, int nSd){
		double dist = Math.abs(r.pos1 - r.pos2);
		return (dist < ins[0]+nSd*ins[1] && dist > ins[0]+nSd*ins[1]);
	}
	
	
	public static void filterDiagReads(Map<String,ReadPair> reads, double[] ins, int nSd){
		Iterator<String> it = reads.keySet().iterator();
		String key = null;
		ReadPair r = null;
		while(it.hasNext()){
			key = it.next();
			r = reads.get(key);
			if (isDiag(r, ins, nSd))
				reads.remove(key);
		}
	}
	
	
	private static class PrintStreamPair{
		private static Random r = new Random(123456789);
		private File fishDir;
		private PrintStream out1;
		private PrintStream out2;
		File file1;
		File file2;
		private Contig ctg1;
		private Contig ctg2;
		private Vector<ReadPair> reads;
		
		public PrintStreamPair(Contig ctg1, Contig ctg2, File outdir) throws IOException{
			fishDir = outdir;
			reads = new Vector<ReadPair>();
			if (ctg1==ctg2){
				file1 = new File(fishDir,"match."+ctg1.getId()+"v"+ctg2.getId()+".txt");
				file1.createNewFile();
				this.ctg1 = ctg1;
				this.ctg2 = ctg2;
			} else { 
				file1 = new File(fishDir,"match."+ctg1.getId()+"v"+ctg2.getId()+".txt");
				file2 = new File(fishDir,"match."+ctg2.getId()+"v"+ctg1.getId()+".txt");
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
			if ((ctg1 == pair.ctg1 && ctg2 == pair.ctg2)) {
				
				out1.println("c"+ctg1.getId()+"p"+pair.pos1+"\tc"+ctg2.getId()+"p"+pair.pos2+"\t"+(r.nextInt(100)+201));
				if (ctg1 != ctg2)
					out2.println("c"+ctg2.getId()+"p"+pair.pos2+"\tc"+ctg1.getId()+"p"+pair.pos1+"\t"+(r.nextInt(100)+201));
				else if (ctg1.getId()=='1'){
					System.out.print("");
				}
			} else if((ctg1 == pair.ctg2 && ctg2 == pair.ctg1)){
				out1.println("c"+ctg1.getId()+"p"+pair.pos2+"\tc"+ctg2.getId()+"p"+pair.pos1+"\t"+(r.nextInt(100)+201));
				out2.println("c"+ctg2.getId()+"p"+pair.pos1+"\tc"+ctg1.getId()+"p"+pair.pos2+"\t"+(r.nextInt(100)+201));
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
				print(pair);
			
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
			return ctg1.getId()+"\t"+ctg2.getId()+"\t"+file.getParentFile().getName()+"/"+file.getName();
		}
		
		public int compareTo(MatchFile arg0) {
			if (this.ctg2.getId() - arg0.ctg2.getId() == 0){
				return this.ctg1.getId() - arg0.ctg1.getId();
			} else {
				return this.ctg2.getId() - arg0.ctg2.getId();
			}
		}
		
		
		
	}

}
