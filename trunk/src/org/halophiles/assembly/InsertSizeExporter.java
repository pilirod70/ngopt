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
				
				double insTot[] = ReadPair.estimateInsertSize(reads.values());
				if (insTot[0]>1500){
					System.err.println("[a5_ise] EM Clustering with K = 2 to remove noisy reads and shadow library");			
					EMClusterer em = new EMClusterer(reads.values(), 2);
					em.iterate(1000,0.0001);
					//System.out.println(" converged in "+numIt+" iterations");
					Collection<ReadSet> emClusters = em.getClusters();
					Iterator<ReadSet> clustIt = emClusters.iterator();
					//System.out.println(emClusters.size()+" clusters");
					ReadSet tmpClust = null;
					while(clustIt.hasNext()){
						tmpClust = clustIt.next();
						if (tmpClust.mean()>1500){
							insTot[0] = tmpClust.mean();
							insTot[1] = tmpClust.sd();
							insTot[2] = tmpClust.size();
							break;
						} 
					}
					
				}
				System.out.println(NF.format(insTot[0])+","+NF.format(insTot[1])+","+NF.format(insTot[2]));
				System.out.flush();
			}
		} catch (IOException e){
			e.printStackTrace();
			System.exit(-1);
		}

	}

}
