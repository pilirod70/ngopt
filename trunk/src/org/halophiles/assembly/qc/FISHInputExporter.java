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
import java.util.Vector;

import org.halophiles.assembly.Contig;
import org.halophiles.assembly.ReadPair;
//import org.halophiles.assembly.scaffold.CountMaxFlow.PrintStreamPair;

public class FISHInputExporter {
	
	public static void export(Map<String,Contig> contigs, Map<String,ReadPair> reads, File fishDir) throws IOException{
		HashMap<String,Vector<Double>> mapPoints = new HashMap<String,Vector<Double>>();
		HashMap<String,PrintStreamPair> matchesOut = new HashMap<String,PrintStreamPair>();
		Vector<Double> tmpPts = null;
		
		File controlFile = new File(fishDir.getParentFile(),"control.txt");
		PrintStream ctrlOut = new PrintStream(controlFile);
		Iterator<ReadPair> rpIt = reads.values().iterator();
		while(rpIt.hasNext()){
			ReadPair tmpRp = rpIt.next();
			if (!reads.containsKey(tmpRp.hdr))
				reads.put(tmpRp.hdr, tmpRp);
			
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
			rev = tmpRp.rev1;
			left = tmpRp.pos1;
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
		
		Iterator<String> it = mapPoints.keySet().iterator();
		ctrlOut.println("-maps");
		while(it.hasNext()){
			String ctg = it.next();
			Vector<Double> tmp = mapPoints.get(ctg);
			Double[] ar = new Double[tmp.size()];
			tmp.toArray(ar);
			File ctgMapFile = new File(fishDir,"map."+contigs.get(ctg).getId()+".txt");
			PrintStream ctgMapOut = new PrintStream(ctgMapFile);
			
			Arrays.sort(ar, new Comparator<Double>(){
				public int compare(Double arg0, Double arg1) {
					if (Math.abs(arg0) < Math.abs(arg1)) return -1;
					 else if (Math.abs(arg0) > Math.abs(arg1)) return 1;
					 else return 0;
				}	
			});
			for (Double d: ar){
				if (d.doubleValue()==904)
					System.out.print("");
				ctgMapOut.println("c"+contigs.get(ctg).getId()+"p"+Math.abs(d.intValue())+"\t"+(d < 0 ? -1 : 1));
			}
			ctgMapOut.close();
			String dir = ctgMapFile.getParentFile().getName();
			ctrlOut.println(contigs.get(ctg).getId()+"\t"+dir+"/"+ctgMapFile.getName());
		}
		
//		int rdlen = sfp.getReadLength();
		it = reads.keySet().iterator();
		HashSet<ReadPair> conflicts = new HashSet<ReadPair>();
		PrintStreamPair matchOut = null;
		ctrlOut.println("-matches");
		String ctgStr;
		while(it.hasNext()) {
			ReadPair tmp = reads.get(it.next());
			if (tmp.ctg1 == null || tmp.ctg2 == null) 
				continue;
			ctgStr = tmp.contigString();
			if (matchesOut.containsKey(ctgStr))
				matchOut = matchesOut.get(ctgStr);
			else {
				matchOut = new PrintStreamPair(tmp.ctg1,tmp.ctg2,fishDir);
				matchOut.printControl(ctrlOut);
				matchesOut.put(ctgStr, matchOut);
			}
			matchOut.print(tmp);
		}
		ctrlOut.println("end");
		ctrlOut.close();
		System.err.println(matchesOut.size()+" PrintStreamPairs for " + contigs.size() + " contigs");
		it = matchesOut.keySet().iterator(); 
		while(it.hasNext())
			matchesOut.get(it.next()).close();
	}
	static class PrintStreamPair {
		private File fishDir;
		private PrintStream out1;
		private PrintStream out2;
		File file1;
		File file2;
		private Contig ctg1;
		private Contig ctg2;
		public PrintStreamPair(Contig ctg1, Contig ctg2, File outdir) throws IOException{
			fishDir = outdir;
			if (ctg1==ctg2){
				file1 = new File(fishDir,"match."+ctg1.getId()+"v"+ctg2.getId()+".txt");
				file1.createNewFile();
				out1 = new PrintStream(file1);
				this.ctg1 = ctg1;
				this.ctg2 = ctg2;
			} else { 
				file1 = new File(fishDir,"match."+ctg1.getId()+"v"+ctg2.getId()+".txt");
				file2 = new File(fishDir,"match."+ctg2.getId()+"v"+ctg1.getId()+".txt");
				file1.createNewFile();
				out1 = new PrintStream(file1);
				file2.createNewFile();
				out2 = new PrintStream(file2);
				this.ctg1 = ctg1;
				this.ctg2 = ctg2;
			}
		}
		public void print(ReadPair pair){
			if ((ctg1 == pair.ctg1 && ctg2 == pair.ctg2)) {
				out1.println("c"+ctg1.getId()+"p"+pair.pos1+"\tc"+ctg2.getId()+"p"+pair.pos2+"\t100");
				if (ctg1 != ctg2)
					out2.println("c"+ctg2.getId()+"p"+pair.pos2+"\tc"+ctg1.getId()+"p"+pair.pos1+"\t100");
				else if (ctg1.getId()=='1'){
					System.out.print("");
				}
			} else if((ctg1 == pair.ctg2 && ctg2 == pair.ctg1)){
				out1.println("c"+ctg1.getId()+"p"+pair.pos2+"\tc"+ctg2.getId()+"p"+pair.pos1+"\t100");
				out2.println("c"+ctg2.getId()+"p"+pair.pos1+"\tc"+ctg1.getId()+"p"+pair.pos2+"\t100");
			} else {
				throw new IllegalArgumentException("pair "+pair.hdr+" does not connect contigs "+ctg1.name+" and "+ctg2.name);
			}
		}
		
		public void close(){
			out1.close();
			if (ctg1 != ctg2)
				out2.close();
		}
		
		public void printControl(PrintStream ctrlOut){
			if (ctg1 == ctg2){
				ctrlOut.println(ctg1.getId()+"\t"+ctg2.getId()+"\t"+fishDir.getName()+"/"+file1.getName());
			} else {
				ctrlOut.println(ctg1.getId()+"\t"+ctg2.getId()+"\t"+fishDir.getName()+"/"+file1.getName());
				ctrlOut.println(ctg2.getId()+"\t"+ctg1.getId()+"\t"+fishDir.getName()+"/"+file2.getName());
			}
		}
	}

}
