/**
 * This file is part of the A5 pipeline.
 * (c) 2011, 2012 Andrew Tritt and Aaron Darling
 * This software is licensed under the GPL, v3.0. Please see the file LICENSE for details
 */
package org.halophiles.assembly.qc;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Random;

import org.halophiles.assembly.Contig;
import org.halophiles.assembly.ReadPair;
import org.halophiles.assembly.SAMFileParser;



public class LGAInputExporter {

	public static void main(String[] args){
		if (args.length != 2 && args.length != 3){
			System.out.println("Usage: LGAInputExporter <output_dir> <sam_file> <gcl_file(optional)>");
			System.exit(-1);
		}
		String outdirPath = args[0];
		String samFilePath = args[1];
		String gclFilePath = null;
		try {
			
			SAMFileParser sfp = null;
			if (args.length == 3){
				gclFilePath = args[2];
				sfp = new SAMFileParser(samFilePath);
			}
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
			File outdir = new File(outdirPath);
			outdir.mkdirs();
			export(contigs, reads, outdir);
			
		}catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
	}
	
	public static void export(Map<String,Contig> contigs, Map<String,ReadPair> reads, File outdir) throws IOException{
		Map<Contig,Map<Contig,PrintStream>> ps = new HashMap<Contig, Map<Contig,PrintStream>>();
		Iterator<ReadPair> it = reads.values().iterator();
		PrintStream out = null;
		File file = null;
		Map<Contig,PrintStream> map = null;
		String fileName = null;
		while(it.hasNext()){
			ReadPair tmp = it.next();
			if (!tmp.paired){
				continue;
			}
			if (ps.containsKey(tmp.ctg1)){
				if (ps.get(tmp.ctg1).containsKey(tmp.ctg2)){
					out = ps.get(tmp.ctg1).get(tmp.ctg2);
				} else {
					if (tmp.ctg1.equals(tmp.ctg2))
						fileName = tmp.ctg2.name+".txt";
					 else
						fileName = tmp.ctg1.name+"-"+tmp.ctg2.name+".txt";
					out = new PrintStream(new File(outdir,fileName));
					ps.get(tmp.ctg1).put(tmp.ctg2,out);
				}
			} else {
				if (tmp.ctg1.equals(tmp.ctg2))
					fileName = tmp.ctg2.name+".txt";
				 else
					fileName = tmp.ctg1.name+"-"+tmp.ctg2.name+".txt";
				out = new PrintStream(new File(outdir,fileName));
				map = new HashMap<Contig,PrintStream>();
				map.put(tmp.ctg2, out);
				ps.put(tmp.ctg1, map);
			}
			out.println(tmp.hdr+"\t"+tmp.ctg1.name+"\t"+tmp.pos1+"\t"+tmp.ctg2.name+"\t"+tmp.pos2);
		}
		
		Iterator<Contig> it1 = ps.keySet().iterator();
		Contig tmp1;
		String rand = "RANDOM";
		int count = 0;
		Random r = new Random(123456789);
		while(it1.hasNext()){
			tmp1 = it1.next();
			map = ps.get(tmp1);
			Iterator<Contig> it2 = map.keySet().iterator();
			Contig tmp2;
			while(it2.hasNext()){
				tmp2 = it2.next();
				out = map.get(tmp2);
				for (int i = 0; i < 50; i++){
					count++;
					int p1 = r.nextInt(tmp1.len)+1;
					int p2 = r.nextInt(tmp2.len)+1;
					out.println(rand+count+"\t"+tmp1.name+"\t"+p1+"\t"+tmp2.name+"\t"+p2);
				}
				out.close();
			}
		}
	}
	
}
