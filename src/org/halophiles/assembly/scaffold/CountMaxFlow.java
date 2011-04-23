package org.halophiles.assembly.scaffold;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.Vector;
import java.util.zip.GZIPInputStream;

import jsc.datastructures.PairedData;
import jsc.regression.PearsonCorrelation;

import org.jgrapht.EdgeFactory;
import org.jgrapht.Graph;
import org.jgrapht.VertexFactory;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.alg.EdmondsKarpMaximumFlow;
import org.jgrapht.ext.ComponentAttributeProvider;
import org.jgrapht.ext.DOTExporter;
import org.jgrapht.ext.EdgeNameProvider;
import org.jgrapht.ext.IntegerNameProvider;
import org.jgrapht.ext.StringNameProvider;
import org.jgrapht.ext.VertexNameProvider;
import org.jgrapht.graph.AsUndirectedGraph;
import org.jgrapht.graph.ClassBasedEdgeFactory;
import org.jgrapht.graph.ClassBasedVertexFactory;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.DirectedWeightedMultigraph;


public class CountMaxFlow {
	private static int MIN_LINK = 1;
	private static int NSAMCOL = 11;
	private static HashMap<String,Contig> contigs;
	private static HashMap<String,ReadPair> reads;
	private static File outdir; 
	private static HashMap<String,Vector<Double>> mapPoints;
	private static File fishDir;
	private static HashMap<String,PrintStreamPair> matchesOut;
	
	public static void main(String[] args){
		if (/*args.length % 3 != 0 ||*/args.length == 0){
			System.err.println("Usage: java CountMaxFlow <outdir> <gcl_contig_file> <sam1> <ins1> <err1> .... <samN> <insN> <errN>");
			System.exit(-1);
		}
		
		VertexFactory<Contig> vf = new ClassBasedVertexFactory<Contig>(Contig.class);
		EdgeFactory<Contig,DefaultWeightedEdge> ef = new ClassBasedEdgeFactory<Contig, DefaultWeightedEdge>(DefaultWeightedEdge.class);
		DirectedWeightedMultigraph<Contig,DefaultWeightedEdge> dg = new DirectedWeightedMultigraph<Contig, DefaultWeightedEdge>(ef);
		
		
		contigs = new HashMap<String, Contig>();
		reads = new HashMap<String, ReadPair>();
		outdir = new File(args[0]);
		fishDir = new File(outdir,"fish_dat");
		
		matchesOut = new HashMap<String,PrintStreamPair>();
		mapPoints = new HashMap<String,Vector<Double>>();
		
		File gclFile = new File(args[1]);
		BufferedReader br = null;
		
		try {
			if (!fishDir.exists()) fishDir.mkdirs();
			br = new BufferedReader(new FileReader(gclFile));
			br.readLine();
			while(br.ready()){
				String[] line = br.readLine().split("\t");
				Contig tmp = vf.createVertex();
				tmp.setInfo(line[0], Integer.parseInt(line[3]), Double.parseDouble(line[1]));
				contigs.put(tmp.name, tmp);
			}
		} catch (Exception e){
			e.printStackTrace();
			System.exit(-1);
		}
		for (int i = 2; i < args.length;i+=3){
			File maxFlowFile = new File(outdir,"max_flow.txt");
			File edgeWeightsFile = new File(outdir,"edge_weights.txt");
			File cnctTypeFile = new File(outdir,"connection_counts.txt");
			File dotFile = new File(outdir,"normalized_graph.dot");
			File insDistFile = new File(outdir,"insert_distribution.txt");
			File conflictSAMFile = new File(outdir,"conflicting_pairs.sam");
			File conflictTDFile = new File(outdir,"conflicting_pairs.txt");
			File correlationsFile = new File(outdir,"correlation.txt");
			File controlFile = new File(outdir,"control.txt");
			try {
				maxFlowFile.createNewFile();
				edgeWeightsFile.createNewFile();
				cnctTypeFile.createNewFile();
				dotFile.createNewFile();
				insDistFile.createNewFile();
				conflictSAMFile.createNewFile();
				controlFile.createNewFile();
				PrintStream mfOut = new PrintStream(maxFlowFile);
				PrintStream ewOut = new PrintStream(edgeWeightsFile);
				PrintStream ctOut = new PrintStream(cnctTypeFile);
				PrintStream dotOut = new PrintStream(dotFile);
				PrintStream idOut = new PrintStream(insDistFile);
				PrintStream cftSAMOut = new PrintStream(conflictSAMFile);
				PrintStream cftTDOut = new PrintStream(conflictTDFile);
				PrintStream corrOut = new PrintStream(correlationsFile);
				PrintStream ctrlOut = new PrintStream(controlFile);

				String ctgStr = null;
				br = null;
				File samFile = new File(args[i]);
				int ins = Integer.parseInt(args[i+1]);
				double err = Double.parseDouble(args[i+2]);
				if (samFile.getName().endsWith(".gz")||samFile.getName().endsWith(".Z")||samFile.getName().endsWith(".z"))
					br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(samFile))));
				else
					br = new BufferedReader(new FileReader(samFile));
				while(nextCharIs(br, '@')) 
					br.readLine();
				int rdlen = 0;
				int den = 0;
				Vector<Double> tmpPts = null;
				while(br.ready()){
					ReadPair tmp = null;
					String[] line = br.readLine().split("\t");
					int left = Integer.parseInt(line[3]);
					boolean rev = isReverse(line[1]);
					int len = cigarLength(line[5]);
					rdlen += len;
					den++;
					ctgStr = line[2].split("\\|")[0];
					if (len <= 0) continue;
					Contig tmpCtg = contigs.get(ctgStr);
					if (reads.containsKey(line[0]))
						tmp = reads.get(line[0]);
					else {
						tmp = new ReadPair(line[0]);
						reads.put(line[0], tmp);
					}
					int val = tmp.addRead(left, rev, tmpCtg, line);
					if (mapPoints.containsKey(tmpCtg.name)){
						tmpPts = mapPoints.get(tmpCtg.name);
					} else {
						tmpPts = new Vector<Double>();
						mapPoints.put(tmpCtg.name, tmpPts);
					}
					if (rev)
						left = left * -1;
					if (tmp.paired)
						tmpPts.add(new Double(left));
					if (val == -1)
						System.err.println("ambiguous mapping for read " + line[0]);
				}
				
				Iterator<String> it = mapPoints.keySet().iterator();
				ctrlOut.println("-maps");
				while(it.hasNext()){
					String ctg = it.next();
					Vector<Double> tmp = mapPoints.get(ctg);
					Double[] ar = new Double[tmp.size()];
					tmp.toArray(ar);
					File ctgMapFile = new File(fishDir,"map."+ctg.charAt(8)+".txt");
					PrintStream ctgMapOut = new PrintStream(ctgMapFile);
					
					Arrays.sort(ar, new Comparator<Double>(){
						public int compare(Double arg0, Double arg1) {
							if (Math.abs(arg0) < Math.abs(arg1)) return -1;
							 else if (Math.abs(arg0) > Math.abs(arg1)) return 1;
							 else return 0;
						}	
					});
					for (Double d: ar)
						ctgMapOut.println("c"+ctg.charAt(8)+"p"+Math.abs(d.intValue())+"\t"+(d < 0 ? -1 : 1));
					ctgMapOut.close();
					String dir = ctgMapFile.getParentFile().getName();
					ctrlOut.println(ctg.charAt(8)+"\t"+dir+"/"+ctgMapFile.getName());
				}
				
				rdlen = rdlen/den;
				it = reads.keySet().iterator();
				HashSet<ReadPair> conflicts = new HashSet<ReadPair>();
				PrintStreamPair matchOut = null;
				ctrlOut.println("-matches");
				while(it.hasNext()) {
					ReadPair tmp = reads.get(it.next());
					if (tmp.ctg1 == null || tmp.ctg2 == null) 
						continue;
					if (tmp.ctg1.equals(tmp.ctg2)) {
						if (spansTerminal(tmp,ins,err,rdlen)){
							tmp.ctg1.addEndSpanningPair();
						} else {
							idOut.println(tmp.hdr+"\t"+tmp.ctg1.getConcatCoord(tmp.pos1)
									      +"\t"+tmp.ctg2.getConcatCoord(tmp.pos2)+"\t"+
									      tmp.ctg2.getConcatCoord(((tmp.pos1+tmp.pos2)/2))
									      +"\t"+getInsertSize(tmp, rdlen));
						}
					} else { 
						if (isTerminal(tmp.pos1, tmp.pos1+rdlen-1, ins, err, tmp.ctg1) && 
							isTerminal(tmp.pos2, tmp.pos2+rdlen-1, ins, err, tmp.ctg2)){
							tmp.ctg1.addLink(tmp.ctg2);
					 		tmp.ctg2.addLink(tmp.ctg1);
						} else {
							conflicts.add(tmp);
						}
					}
					ctgStr = tmp.contigString();
					if (matchesOut.containsKey(ctgStr))
						matchOut = matchesOut.get(ctgStr);
					else {
						matchOut = new PrintStreamPair(tmp.ctg1,tmp.ctg2);
						matchOut.printControl(ctrlOut);
						matchesOut.put(ctgStr, matchOut);
					}
					matchOut.print(tmp);
				}
				ctrlOut.println("end");
				ctrlOut.close();
				
				it = matchesOut.keySet().iterator(); 
				while(it.hasNext())
					matchesOut.get(it.next()).close();
				idOut.close();
				it = contigs.keySet().iterator();
				while (it.hasNext()){
					Contig tmp = contigs.get(it.next());
					cftSAMOut.println("@SQ\tSN:"+tmp.name+"\tLN:"+tmp.len);
				}
				Iterator<ReadPair> rpIt = conflicts.iterator();
				Vector<double[]> dat = null;
				HashMap<String,Vector<double[]>> datSets = new HashMap<String,Vector<double[]>>();
				
				while(rpIt.hasNext()){
					ReadPair tmp = rpIt.next();
					cftSAMOut.print(tmp.sam1[0]);
					for (int sI = 1; sI < NSAMCOL; sI++)
						cftSAMOut.print("\t"+tmp.sam1[sI]);
					cftSAMOut.println();
					cftSAMOut.print(tmp.sam2[0]);
					for (int sI = 1; sI < NSAMCOL; sI++)
						cftSAMOut.print("\t"+tmp.sam2[sI]);
					cftSAMOut.println();
					cftTDOut.println(tmp.hdr+"\t"+tmp.ctg1.name+"\t"+tmp.pos1+"\t"+tmp.ctg2.name+"\t"+tmp.pos2);
					String pair = tmp.ctg1.name+":"+tmp.ctg2.name;
					if (datSets.containsKey(pair))
						dat = datSets.get(pair);
					else {
						dat = new Vector<double[]>();
						datSets.put(pair, dat);
					}
					double[] att = {tmp.pos1,tmp.pos2};
					dat.add(att);
				}
				cftTDOut.close();
				cftSAMOut.close();
				it = datSets.keySet().iterator();
				while(it.hasNext()){
					String ctgKey = it.next();
					String[] ctg = ctgKey.split(":");		
					dat = datSets.get(ctgKey);
					double[][] ar = new double[dat.size()][];
					dat.toArray(ar);
					double[] x = new double[ar.length];
					double[] y = new double[ar.length];
					if (ar.length <= 2) 
						continue;
					for (int pI = 0; pI < ar.length; pI++){
						x[pI] = ar[pI][0];
						y[pI] = ar[pI][1];
					}
					PairedData pdat = new PairedData(x, y);
					try {
						PearsonCorrelation lm = new PearsonCorrelation(pdat);
						corrOut.println(ctg[0]+"\t"+ctg[1]+"\t"+lm.getB()+"\t"+lm.getR()+"\t"+lm.getN());
					} catch (IllegalArgumentException iae){
						iae.printStackTrace();
						System.exit(-1);
					}
					
				}
				corrOut.close();
				
				
				
				
			
				Contig[] ctgRef = new Contig[contigs.values().size()];
				contigs.values().toArray(ctgRef);
				for (int cI = 0; cI < ctgRef.length; cI++){
					ctOut.println(ctgRef[cI].name+"\t"+ctgRef[cI].counts.size()+"\t"+ctgRef[cI].numSelfConnect);
					for (int cJ = cI+1; cJ < ctgRef.length; cJ++){
						int nlink = ctgRef[cI].nLinks(ctgRef[cJ]);
						if (nlink < MIN_LINK) continue;
						if (!dg.containsVertex(ctgRef[cI]))	dg.addVertex(ctgRef[cI]);
						if (!dg.containsVertex(ctgRef[cJ])) dg.addVertex(ctgRef[cJ]);
						double L = Math.min(Math.min(ins*(1.0+err), ctgRef[cI].len), Math.min(ins*(1.0+err), ctgRef[cJ].len));
						double Cov = (ctgRef[cI].cov+ctgRef[cJ].cov)/2;
						double normed = nlink*rdlen/(Cov*L);
						ewOut.println(normed+"\t"+nlink);
						DefaultWeightedEdge adj = dg.addEdge(ctgRef[cI], ctgRef[cJ]);
						dg.setEdgeWeight(adj, normed);
						adj = dg.addEdge(ctgRef[cJ],ctgRef[cI]);
						dg.setEdgeWeight(adj, normed);
					}
				}
				ewOut.close();
				ctOut.close();
				EdmondsKarpMaximumFlow<Contig, DefaultWeightedEdge> ekmf = new EdmondsKarpMaximumFlow<Contig, DefaultWeightedEdge>(dg);
				ConnectivityInspector<Contig, DefaultWeightedEdge> ci = new ConnectivityInspector<Contig, DefaultWeightedEdge>(dg);
				Iterator<Set<Contig>> ccIt = ci.connectedSets().iterator();
				int ccCount = 0;
				while(ccIt.hasNext()){
					Set<Contig> cc = ccIt.next();
					ccCount++;
					ctgRef = new Contig[cc.size()];
					cc.toArray(ctgRef);
					for (int cI = 0; cI < ctgRef.length; cI++){
						if (!dg.containsVertex(ctgRef[cI])) continue;
						for (int cJ = cI+1; cJ < ctgRef.length; cJ++){
							if (!dg.containsVertex(ctgRef[cJ])) continue;
							ekmf.calculateMaximumFlow(ctgRef[cI], ctgRef[cJ]);
							double maxFlow = ekmf.getMaximumFlowValue();
							if (maxFlow > 0.0)
								mfOut.println(ccCount+"\t"+ctgRef[cI]+"\t"+ctgRef[cJ]+"\t"+maxFlow);
						}
					}
				}
				AsUndirectedGraph<Contig,DefaultWeightedEdge> udg = new AsUndirectedGraph<Contig,DefaultWeightedEdge>(dg);
				
				Iterator<DefaultWeightedEdge> edgeIt = udg.edgeSet().iterator();
				Set<DefaultWeightedEdge> disc = new HashSet<DefaultWeightedEdge>(); 
				while(edgeIt.hasNext()){
					DefaultWeightedEdge e = edgeIt.next();
					if (disc.contains(e)) 
						continue;
					Contig src = udg.getEdgeSource(e);
					Contig tgt = udg.getEdgeTarget(e);
					if (udg.containsEdge(tgt, src)) 
						disc.add(udg.getEdge(tgt, src));
				}
				edgeIt = disc.iterator();
				while(edgeIt.hasNext())
					udg.removeEdge(edgeIt.next());
				
				System.err.println(udg.edgeSet().size()+" edges");
				VertexNameProvider<Contig> inp = new IntegerNameProvider<Contig>();
				VertexNameProvider<Contig> snp = new StringNameProvider<Contig>();
				EdgeNameProvider<DefaultWeightedEdge> enp = new DWENameProvider(udg);
				ComponentAttributeProvider<DefaultWeightedEdge> eap = new DWEAttributeProvider(udg);
				ComponentAttributeProvider<Contig> vap = new ContigAttributeProvider();
				DOTExporter<Contig,DefaultWeightedEdge> dotexp = new DOTExporter<Contig,DefaultWeightedEdge>(inp,snp,enp,vap,eap);
				dotexp.export(new PrintWriter(dotOut), udg);
				dotOut.close();
				mfOut.close();
			} catch (Exception e){
				e.printStackTrace();
				System.exit(-1);
			}
		}
	}
	
	static class PrintStreamPair {
		private PrintStream out1;
		private PrintStream out2;
		final File file1;
		final File file2;
		private Contig ctg1;
		private Contig ctg2;
		public PrintStreamPair(Contig ctg1, Contig ctg2) throws IOException{
			file1 = new File(fishDir,"match."+ctg1.name.charAt(8)+"v"+ctg2.name.charAt(8)+".txt");
			file2 = new File(fishDir,"match."+ctg2.name.charAt(8)+"v"+ctg1.name.charAt(8)+".txt");
			file1.createNewFile();
			out1 = new PrintStream(file1);
			file2.createNewFile();
			out2 = new PrintStream(file2);
			this.ctg1 = ctg1;
			this.ctg2 = ctg2;
		}
		public void print(ReadPair pair){
			if ((ctg1 == pair.ctg1 && ctg2 == pair.ctg2)) {
				out1.println("c"+ctg1.name.charAt(8)+"p"+pair.pos1+"\tc"+ctg2.name.charAt(8)+pair.pos2+"\t100");
				out2.println("c"+ctg2.name.charAt(8)+"p"+pair.pos2+"\tc"+ctg1.name.charAt(8)+pair.pos1+"\t100");
			} else if((ctg1 == pair.ctg2 && ctg2 == pair.ctg1)){
				out1.println("c"+ctg1.name.charAt(8)+"p"+pair.pos2+"\tc"+ctg2.name.charAt(8)+pair.pos1+"\t100");
				out2.println("c"+ctg2.name.charAt(8)+"p"+pair.pos1+"\tc"+ctg1.name.charAt(8)+pair.pos2+"\t100");
			} else {
				throw new IllegalArgumentException("pair doesn't connect contigs");
			}
		}
		
		public void close(){
			out1.close();
			out2.close();
		}
		
		public void printControl(PrintStream ctrlOut){
			if (ctg1 == ctg2){
				ctrlOut.println(ctg1.name.charAt(8)+"\t"+ctg2.name.charAt(8)+"\t"+fishDir.getName()+"/"+file1.getName());
			} else {
				ctrlOut.println(ctg1.name.charAt(8)+"\t"+ctg2.name.charAt(8)+"\t"+fishDir.getName()+"/"+file1.getName());
				ctrlOut.println(ctg2.name.charAt(8)+"\t"+ctg1.name.charAt(8)+"\t"+fishDir.getName()+"/"+file2.getName());
			}
		}
	}
	
	static class ReadPair{
		public final String hdr;
		public int pos1 = 0;
		public boolean rev1 = false;
		public Contig ctg1;
		public String[] sam1;
		public int pos2 = 0;
		public boolean rev2 = false;
		public Contig ctg2;
		public String[] sam2;
		public boolean paired;
		public ReadPair(String hdr){
			this.hdr=hdr;
		}
		public String contigString(){
			return ctg1.name+"-"+ctg2.name;
		}
		public int addRead(int pos, boolean rev, Contig ctg, String[] sam){
			if (pos1 == 0){
				pos1 = pos;
				rev1 = rev;
				ctg1 = ctg;
				sam1 = sam;
				paired = false;
				return 1;
			} else if (pos2 == 0){
				if (ctg1.name.compareTo(ctg.name)<0){ // alphabetize for consistency
					pos2 = pos;
					rev2 = rev;
					ctg2 = ctg;
					sam2 = sam;
				} else {
					pos2 = pos1;
					rev2 = rev1;
					ctg2 = ctg1;
					sam2 = sam1;
					pos1 = pos;
					rev1 = rev;
					ctg1 = ctg;
					sam1 = sam;
				}
				paired = true;
				return 2;
			} else {
				return -1;
			}
		}
	}
	public static class Contig implements Comparable<Contig> {
		private static int CONCAT_START = 1;
		public String name;
		public int len;
		public double cov;
		private int start;
		public Map<Contig,Integer> counts;
		public int numSelfConnect;
		public void setInfo(String name, int len, double cov){
			String[] spl = name.split("\\|");
			this.name = spl[0];
			this.len = len;
			this.cov = cov;
			numSelfConnect = 0;
			start=CONCAT_START;
			CONCAT_START+=len;
			counts = new HashMap<Contig,Integer>();
		}
		public boolean equals(Contig c){ 
			return this.name.equals(c.name);
		}
		public void addLink(Contig c){
			if (!counts.containsKey(c)) counts.put(c, 1);
			else counts.put(c, counts.get(c)+1);
		}
		public int nLinks(Contig c){
			if (!counts.containsKey(c)) return 0;
			return counts.get(c);
		}
		public int hashCode(){
			return name.hashCode();
		}
		public String toString(){
			return name;
		}
		public void addEndSpanningPair(){
			numSelfConnect++;
		}
		public int getConcatCoord(int pos){
		/*	if (pos < start) 
				return -1;
			else*/
				return start+pos-1;
		}
		@Override
		public int compareTo(Contig arg0) {
			return this.name.compareTo(arg0.name);
		}
	}
	
	public static class DWENameProvider implements EdgeNameProvider<DefaultWeightedEdge>{
		
		private static NumberFormat nf; 
		//private DirectedWeightedMultigraph<Contig,DefaultWeightedEdge> g;
		private Graph<Contig,DefaultWeightedEdge> g;
		public DWENameProvider(Graph<Contig,DefaultWeightedEdge> dg){
			this.g = dg;
			if (nf == null){
				nf =  NumberFormat.getInstance();
				nf.setMaximumFractionDigits(6);
			}
		}
		@Override
		public String getEdgeName(DefaultWeightedEdge arg0) {
			return nf.format(g.getEdgeWeight(arg0));
			//return Double.toString(g.getEdgeWeight(arg0));
		}
	}
	
	public static class DWEAttributeProvider implements ComponentAttributeProvider<DefaultWeightedEdge>{
		private Graph<Contig,DefaultWeightedEdge> g;
		public DWEAttributeProvider(Graph<Contig, DefaultWeightedEdge> g){
			this.g = g;
		}
		@Override
		public Map<String, String> getComponentAttributes(DefaultWeightedEdge arg0) {
			Map<String,String> map = new HashMap<String,String>();
			map.put("label", g.getEdgeSource(arg0).toString()+":"+g.getEdgeTarget(arg0).toString());
			map.put("weight", Double.toString(g.getEdgeWeight(arg0)));
			map.put("nlinks", Integer.toString(g.getEdgeSource(arg0).counts.get(g.getEdgeTarget(arg0))));
			return map;
		}
	}
	
	public static class ContigAttributeProvider implements ComponentAttributeProvider<Contig>{

		@Override
		public Map<String, String> getComponentAttributes(Contig arg0) {
			Map<String,String> map = new HashMap<String,String>();
			map.put("label", arg0.name.substring(arg0.name.indexOf('|')+1));
			map.put("length", Integer.toString(arg0.len));
			map.put("cov", Double.toString(arg0.cov));
			return map;
		}
		
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
	
	private static boolean isTerminal(int left, int right, int ins, double err, Contig ctg){
		if (ctg.len < (ins)*(1.0+err)) return true;
		if (left < (ins)*(1.0+err)) {
			return true;
		} else if (ctg.len - (ins)*(1.0+err) < right){
			return true;
		}
		return false;
	}
	
	/**
	 * Return mapping length for CIGAR String if there is a match (i.e. string contains 'M') else return -1
	 * @param cig the CIGAR string to parse
	 * @return the length of the match indicated by the CIGAR string. Return -1 if no match
	 */
	private static int cigarLength(String cig){
		if (cig.contains("M"))
			return Integer.parseInt(cig.substring(0,cig.indexOf('M')));
		else return -1;
	}
	
	private static double normLinkWeight(Contig src, Contig sink, int ins, double err, int rdlen){
		double nLink = src.counts.get(sink);
		double L = Math.min(ins*(1.0+err), src.len);
		return (rdlen*nLink)/(src.cov*L);
	}
	
	private static boolean spansTerminal(ReadPair pm, int ins, double err, int rdlen){
		if (!isTerminal(pm.pos1, pm.pos1+rdlen-1, ins, err, pm.ctg1) || 
			!isTerminal(pm.pos2, pm.pos2+rdlen-1, ins, err, pm.ctg2))
			return false;
		int outLen;
		int inLen;
		if (pm.pos1 < pm.pos2){
			outLen = (pm.pos1+rdlen) + (pm.ctg1.len-pm.pos2+1);
			inLen = pm.pos2 - pm.pos1 + rdlen+1;
		} else {
			outLen = (pm.pos2+rdlen) + (pm.ctg1.len-pm.pos1+1);
			inLen = pm.pos1 - pm.pos2 + rdlen+1;
		}
		boolean outLenOk = outLen > (1-err)*ins && outLen < (1+err)*ins;
		boolean inLenOk = inLen > (1-err)*ins && inLen < (1+err)*ins;
		if (outLenOk && inLenOk) 
			return (Math.abs(outLen-ins) < Math.abs(inLen-ins));
		else if (outLenOk) 
			return true;
		else 
			return false;
	}
	
	private static int getInsertSize(ReadPair pm, int rdlen){
		if (!pm.ctg1.equals(pm.ctg2)) 
			return -1;
		return (pm.pos1 < pm.pos2) ? (pm.pos2+rdlen-pm.pos1) : (pm.pos1+rdlen-pm.pos2);
	}
	
	private static double getDistance(ReadPair p1, ReadPair p2){
		double dx = p1.pos1 - p2.pos1;
		double dy = p1.pos2 - p2.pos2;
		return Math.sqrt(Math.pow(dx, 2)+Math.pow(dy,2));
	}
	
}
