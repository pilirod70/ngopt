package org.halophiles.assembly.qc;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.halophiles.assembly.Contig;
import org.halophiles.assembly.ReadPair;

public class LGAInputExporter {

	public static void export(Map<String,Contig> contigs, Map<String,ReadPair> reads, File outdir) throws IOException{
		Map<Contig,Map<Contig,PrintStream>> ps = new HashMap<Contig, Map<Contig,PrintStream>>();
		Iterator<ReadPair> it = reads.values().iterator();
		PrintStream out = null;
		File file = null;
		Map<Contig,PrintStream> map = null;
		String fileName = null;
		while(it.hasNext()){
			ReadPair tmp = it.next();
			
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
		
		Iterator<Map<Contig,PrintStream>> it1 = ps.values().iterator();
		while(it1.hasNext()){
			map = it1.next();
			Iterator<PrintStream> it2 = map.values().iterator();
			while(it2.hasNext())
				it2.next().close();
		}
	}
	
}
