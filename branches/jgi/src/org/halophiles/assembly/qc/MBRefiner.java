package org.halophiles.assembly.qc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Vector;

import org.halophiles.assembly.qc.PileupScorer;


public class MBRefiner {
	
	private static final String SAMTOOLS = "/jgi/tools/bin/samtools";
	
	private static Map<String,double[][]> REGIONS;
	
	public static void main(String[] args) {
		if (args.length != 3){
			System.out.println("Usage: MBRefiner <sam_file> <contig_file> <broken_ctgs_file>");
			System.exit(-1);
		}
		try {		
			String samPath = args[0];
			String refPath = args[1];
			String brokenPath = args[2];
			File ctgFile = new File(refPath);
			File samFile = new File(samPath);
			File brokenFile = new File(brokenPath);
			brokenFile.createNewFile();
			samFile = new File(samFile.getAbsolutePath());
			String base = samFile.getParentFile().getAbsolutePath()+"/"+MisassemblyBreaker.basename(samFile.getName(),".sam");
			String bedPath = base + ".regions.bed";
			String bamPath = base + ".bam";
			
			REGIONS = getRegions(bedPath);

			runMPileup (bamPath,bedPath,refPath,REGIONS);
			
			Iterator<String> it = REGIONS.keySet().iterator();
			String tmp = null;
			double[][] reg = null;
			System.out.println("Contig\tLeft\tRight\tMinPosition\tScore");
			while(it.hasNext()) {
				tmp = it.next();
				reg = REGIONS.get(tmp);
				for (int i = 0; i < reg[0].length; i++){
					System.out.println(tmp + "\t" + (int)reg[0][i]+"\t"+(int)reg[1][i] +"\t"+ (int)reg[3][i] +"\t"+ (reg[2][i]));
				}
			}
			System.out.println("Breaking actual contigs now. Writing them to "+brokenFile.getAbsolutePath());
			BufferedReader br = new BufferedReader(new FileReader(ctgFile));
			StringBuilder sb = null;
			ScaffoldExporter out = new ScaffoldExporter(brokenFile);
			br.read();
			while(br.ready()){
				String tmpCtg = br.readLine(); 
				// Take only the first part of the contig header to be consistent with bwa
				if (tmpCtg.contains(" "))
					tmpCtg = tmpCtg.substring(0,tmpCtg.indexOf(" "));
				sb = new StringBuilder();
				char c = (char) br.read();
				// read in this sequence
				while(c != '>'){
					if (MisassemblyBreaker.isNuc(c))
						sb.append(c);
					if (!br.ready())
						break;
					c = (char) br.read();
				}
				if (!REGIONS.containsKey(tmpCtg)){
					out.export(tmpCtg, sb);
					continue;
				}
				int left = 1;
				int right = 1;
				reg = REGIONS.get(tmpCtg);
				for (int i = 0; i < reg[3].length; i++){
					right = (int)reg[3][i];
					System.out.println("[a5_qc] Exporting "+tmpCtg+" at "+left+"-"+right);
					out.export(tmpCtg, sb, left, right);
					left = right+1;
				}
				right = sb.length();
				System.out.println("[a5_qc] Exporting "+tmpCtg+" at "+left+"-"+right);
				out.export(tmpCtg, sb, left, right);
				
			}
			
			
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		} catch (InterruptedException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	
	public static Map<String,double[]> refine(String bamPath, String bedPath, String ctgPath) throws IOException {
		Map<String,double[][]> regions = getRegions(bedPath);
		try {
			runMPileup(bamPath, bedPath, ctgPath, regions);
		} catch (InterruptedException e) {
			System.err.println("Unable to run samtools mpileup. This is what I know:");
			e.printStackTrace();
			System.exit(-1);
		}
		
		return null;
	}
		
	
	/**
	 * Returns a Map of sets of regions. Each set corresponds to a single contig, 
	 * indexed by the contig name (a String). 
	 * <br><br>
	 * Each set is stored in a matrix with 4 rows.
	 * The first row is the left-most coordinate of each region, the second row stores the right-most
	 * coordinate of each region, the third stores the current 'weakest' position, and the fourth
	 * stores the score of the current 'weakest' position. 
	 * <br><br>
	 * The weakest position is the position with
	 * the lowest pileup score. For the definition of this score, see org.gel.halophiles.qc.PileupScorer.
	 */
	public static Map<String,double[][]> getRegions (String bedPath) throws IOException {
		BufferedReader in = new BufferedReader(new FileReader(new File(bedPath)));
		String[] dat = null;
		
		/*
		 * BEGIN: build data structure storing the boundaries of the ranges we're investigating
		 */
		// Some Map objects to store our ranges indexed, by contig. Vectors will be turned into
		// arrays that can be binary searched easily
		Map<String,Vector<Integer> > left = new HashMap<String, Vector<Integer> >();
		Map<String,Vector<Integer> > right = new HashMap<String, Vector<Integer> >();		
		Vector<Integer> vL = null; // temp variable
		Vector<Integer> vR = null; // temp variable
		while (in.ready()){
			dat = in.readLine().split("\t");
			// skip any header in the file
			if (dat[0].startsWith("#")) 
				continue;
			if (left.containsKey(dat[0])){
				vL = left.get(dat[0]);
				vR = right.get(dat[0]);
			} else {
				vL = new Vector<Integer>();
				vR = new Vector<Integer>();
				left.put(dat[0], vL);
				right.put(dat[0], vR);
			}
			vL.add(Integer.parseInt(dat[1]));
			vR.add(Integer.parseInt(dat[2]));
		}
		/*
		 * END: build boundaries
		 * 
		 * START: convert these Vectors into arrays, so they can be queried using binary search
		 */
		Map<String,double[][]> lookups = new HashMap<String,double[][]>();
		Iterator<String> it = left.keySet().iterator();
		String key = null; // temp variable
		double[] arL = null; // temp variable
		double[] arR = null; // temp variable
		double[][] tempLookup = null; // temp variable
		while(it.hasNext()){
			key = it.next();
			vL = left.get(key);
			vR = right.get(key);
			arL = toArray(vL);
			arR = toArray(vR);
			tempLookup = new double[4][];
			tempLookup[0] = arL; // left coordinates of ranges
			tempLookup[1] = arR; // right coordinates of ranges
			tempLookup[2] = new double[arL.length];  // 'weakest' position in each range
			for (int i = 0; i < tempLookup[2].length; i++)
				tempLookup[2][i] = Double.POSITIVE_INFINITY;
			tempLookup[3] = new double[arL.length];  // score of 'weakest' position
			lookups.put(key, tempLookup);
		}
		/*
		 * END: convert Vectors to arrays
		 */
		return lookups;
		
	}
	
	private static double[] toArray(Vector<Integer> v){
		double[] ret = new double[v.size()];
		for (int i = 0; i < ret.length; i++)
			ret[i] = v.get(i);
		return ret;
	}

	private static void runMPileup (String bamPath, String bedPath, String ctgPath, Map<String,double[][]> regions) throws IOException, InterruptedException {
		//System.err.println("Executing command: " + cmd);
		String cmd = SAMTOOLS + " mpileup -d 350 -f " + ctgPath + " -l " + bedPath + " " + bamPath;
		Process p = Runtime.getRuntime().exec(cmd);
		LineHandler errLH = new LineHandler(){
			public void handleLine(String line) {
				System.err.println(line);
			}
		};
		StreamReader errRdr = new StreamReader(p.getErrorStream(),System.err, errLH);
		StreamReader outRdr = new StreamReader(p.getInputStream(),System.out,new PileupEvaluator(regions, new PileupScorer(0,0,0,0)));
		errRdr.start();
		outRdr.start();
		p.waitFor();
	}
	
	static class StreamReader extends Thread {
		InputStream in;
		PrintStream ps;
		LineHandler lineHandler;
		StreamReader (InputStream in, PrintStream redirect, LineHandler lineHandler){
			this.in = in;
			this.ps = redirect;
			this.lineHandler = lineHandler;
		}
		
		public void run()
	    {
	        try
	        {
	            InputStreamReader isr = new InputStreamReader(in);
	            BufferedReader br = new BufferedReader(isr);
	            String line=null;
	            while ( (line = br.readLine()) != null)
	            {
	            	lineHandler.handleLine(line);
	            }
	        } 
	    	catch (IOException ioe)
	    	{	
	    	    ioe.printStackTrace();
	    	}
	    }
	}
	
	static class PileupEvaluator implements LineHandler {
		Map<String,double[][]> regions;
		PileupScorer scorer;
		PileupEvaluator (Map<String,double[][]> regions, PileupScorer scorer){
			this.regions = regions;
			this.scorer = scorer;
		}

		public void handleLine(String line) {
			String[] plp = line.split("\t");
			//System.out.println("Lookup at regions for contig "+plp[0]);
        	double[][] tmpReg = regions.get(plp[0]);
			int pos = Integer.parseInt(plp[1]);
			int idx = Arrays.binarySearch(tmpReg[1], pos);
			if (idx < 0) 
				idx = -1*(idx+1);
			double score = scorer.scorePileup(plp[4]);
			if (score < tmpReg[2][idx]){
				tmpReg[2][idx] = score;
				tmpReg[3][idx] = pos;
			}
		}
	}
	
	static interface LineHandler {
		public void handleLine(String line);
	}
}