/**
 * This file is part of the A5 pipeline.
 * (c) 2011, 2012 Andrew Tritt and Aaron Darling
 * This software is licensed under the GPL, v3.0. Please see the file LICENSE for details
 */
package org.halophiles.assembly;

import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.halophiles.assembly.qc.EMClusterer;

public class InsertSizeExporter {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		if (args.length == 0 ){
			System.out.println("Usage: java GetInsertSize.jar <sam_file1> . . . <sam_fileN>");
			System.exit(-1);
		}
		NumberFormat NF = NumberFormat.getInstance();
		NF.setMaximumFractionDigits(0);
		NF.setGroupingUsed(false);
		
		try {
			for (int i = 0; i < args.length; i++){
				File samFile = new File(args[i]);
				System.err.println("[a5_ise] Getting insert stats for " + samFile.getAbsolutePath());
				//String base = samFile.getName().substring(0,samFile.getName().lastIndexOf("."));
				SAMFileParser sfp = new SAMFileParser(samFile.getAbsolutePath());
				ReadPair tmp = null;
				Map<String,ReadPair> reads = new HashMap<String,ReadPair>();
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
				//System.out.println("Filtering duplicate connections...");
				//boolean first = true;
				while(rpIt.hasNext()){
					tmp = rpIt.next();
					/*if (first){
						System.out.println(tmp.hdr)	;
						first =false;
					}*/
					if (!tmp.paired || !(tmp.ctg1.equals(tmp.ctg2))) 
						continue;	
					if(uniq.add(tmp))
						reads.put(tmp.hdr, tmp);
				}
				//ReadPair.exportInsertSize(reads, System.out);
				/*File distFile = new File(new File(args[0]).getParentFile(),"lib3.dist");
				System.out.println(distFile.getAbsolutePath());
				distFile.createNewFile();
				PrintStream distOut = new PrintStream(distFile);
				exportInsertDistances(reads.values(), distOut);*/
				
				double[] insTot = ReadPair.estimateInsertSize(reads.values());
				double[] ret = {insTot[0],insTot[1],insTot[2],ReadPair.getOrientation(reads.values())};
				if (insTot[0]>1500){
					System.err.println("[a5_ise] EM Clustering with K = 3 to remove noisy reads and shadow library");			
					EMClusterer em = new EMClusterer(reads.values(), 3);
					em.iterate(1000,0.0001);
					//System.out.println(" converged in "+numIt+" iterations");
					Collection<ReadSet> emClusters = em.getClusters();
					Iterator<ReadSet> clustIt = emClusters.iterator();
					//System.out.println(emClusters.size()+" clusters");
					ReadSet tmpClust = null;
					double lowV=2000000000;
					ReadSet highC=null;
					while(clustIt.hasNext()){
						tmpClust = clustIt.next();
						System.err.println( "Cluster mean "+tmpClust.mean()+" sd "+tmpClust.sd()+" size "+tmpClust.size() );
						if (tmpClust.mean()>1500 && lowV > tmpClust.sd()){
							lowV = tmpClust.sd();
							highC = tmpClust;
						} 
					}
					
					ret[0] = highC.mean();
					ret[1] = highC.sd();
					ret[2] = highC.size();
					ret[3] = ReadPair.getOrientation(highC.getReads());					
				}
				System.out.println(NF.format(ret[0])+","+NF.format(ret[1])+","+NF.format(ret[2])+","+NF.format(ret[3]));
				System.out.flush();
			}
		} catch (IOException e){
			e.printStackTrace();
			System.exit(-1);
		}

	}

}
