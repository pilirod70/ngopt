package org.halophiles.assembly.scaffold;

import java.io.BufferedReader;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
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
import java.util.TreeSet;

import org.jgrapht.EdgeFactory;
import org.jgrapht.Graph;
import org.jgrapht.VertexFactory;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.alg.EdmondsKarpMaximumFlow;
import org.jgrapht.ext.ComponentAttributeProvider;
import org.jgrapht.ext.DOTExporter;
import org.jgrapht.ext.EdgeNameProvider;
import org.jgrapht.ext.IntegerEdgeNameProvider;
import org.jgrapht.ext.IntegerNameProvider;
import org.jgrapht.ext.StringEdgeNameProvider;
import org.jgrapht.ext.StringNameProvider;
import org.jgrapht.ext.VertexNameProvider;
import org.jgrapht.graph.AsUndirectedGraph;
import org.jgrapht.graph.ClassBasedEdgeFactory;
import org.jgrapht.graph.ClassBasedVertexFactory;
import org.jgrapht.graph.DirectedWeightedMultigraph;
import org.jgrapht.graph.DefaultWeightedEdge;

//import cern.jet.random.Normal;
//import cern.jet.random.engine.MersenneTwister64;
public class CountMaxFlow {
	//private static MersenneTwister64 RAND = new MersenneTwister64(946235897);
	private static int MIN_LINK = 1;
	private static HashMap<String,Contig> contigs;
	private static HashMap<String,ReadPair> reads;
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
		File gclFile = new File(args[1]);
		BufferedReader br = null;
		File outdir = new File(args[0]);
		try {
			if (!outdir.exists()) outdir.mkdirs();
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
			try {
				if (!maxFlowFile.exists()) maxFlowFile.createNewFile();
				if (!edgeWeightsFile.exists()) edgeWeightsFile.createNewFile();
				if (!cnctTypeFile.exists()) cnctTypeFile.createNewFile();
				PrintStream mfOut = new PrintStream(maxFlowFile);
				PrintStream ewOut = new PrintStream(edgeWeightsFile);
				PrintStream ctOut = new PrintStream(cnctTypeFile);
				PrintStream dotOut = new PrintStream(dotFile);
				br = null;
				File samFile = new File(args[i]);
				int ins = Integer.parseInt(args[i+1]);
				double err = Double.parseDouble(args[i+2]);
				br = new BufferedReader(new FileReader(samFile));
				while(nextCharIs(br, '@')) 
					br.readLine();
				ReadPair tmp = null;
				int rdlen = 0;
				int den = 0;
				while(br.ready()){
					String[] line = br.readLine().split("\t");
					int left = Integer.parseInt(line[3]);
					boolean rev = isReverse(line[1]);
					int len = cigarLength(line[5]);
					rdlen += len;
					den++;
					if (len < 0) continue;
					Contig tmpCtg = contigs.get(line[2]);
					if (!isTerminal(left,left+len-1,ins,err,tmpCtg)) continue; 
					if (reads.containsKey(line[0]))
						tmp = reads.get(line[0]);
					else {
						tmp = new ReadPair(line[0]);
						reads.put(line[0], tmp);
					}
					int val = tmp.addRead(left, rev, tmpCtg);
					if (val == -1)
						System.err.println("ambiguous mapping for read " + line[0]);
				}
				rdlen = rdlen/den;
				Iterator<String> it = reads.keySet().iterator();
				while(it.hasNext()) {
					tmp = reads.get(it.next());
					if (tmp.ctg1 == null || tmp.ctg2 == null) 
						continue;
					if (tmp.ctg1.equals(tmp.ctg2)) {
						if (spansTerminal(tmp,ins,err,rdlen)){
							tmp.ctg1.addEndSpanningPair();
						}
					}
					tmp.ctg1.addLink(tmp.ctg2);
					tmp.ctg2.addLink(tmp.ctg1);
				}
				Contig[] ctgRef = new Contig[contigs.values().size()];
				contigs.values().toArray(ctgRef);
				for (int cI = 0; cI < ctgRef.length; cI++){
					ctOut.println(ctgRef[cI].name+"\t"+ctgRef[cI].counts.size()+"\t"+ctgRef[cI].numSelfConnect);
					for (int cJ = cI+1; cJ < ctgRef.length; cJ++){
						int nlink = ctgRef[cI].nLinks(ctgRef[cJ]);
						if (nlink < MIN_LINK) continue;
						if (!dg.containsVertex(ctgRef[cI]))	dg.addVertex(ctgRef[cI]);
						if (!dg.containsVertex(ctgRef[cJ])) dg.addVertex(ctgRef[cJ]);
						
						double normed1 = normLinkWeight(ctgRef[cI], ctgRef[cJ], ins, err, rdlen);
						double normed2 = normLinkWeight(ctgRef[cJ], ctgRef[cI], ins, err, rdlen);
						//double normed = (normed1+normed2)/2;
						//double normed = Math.min(normed1,normed2);
						//Math.min(ins*(1.0+err), src.len);
						double L = Math.min(Math.min(ins*(1.0+err), ctgRef[cI].len), Math.min(ins*(1.0+err), ctgRef[cJ].len));
						double Cov = (ctgRef[cI].cov+ctgRef[cJ].cov)/2;
						double normed = nlink*rdlen/(Cov*L);
						ewOut.println(normed+"\t"+nlink);
						//if (normed < 0.1) continue;
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
		//		System.out.println(ci.connectedSets().size());
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
	
	static class ReadPair{
		public final String hdr;
		public int pos1 = 0;
		public boolean rev1 = false;
		public Contig ctg1;
		public int pos2 = 0;
		public boolean rev2 = false;
		public Contig ctg2;
		public ReadPair(String hdr){
			this.hdr=hdr;
		}
		public int addRead(int pos, boolean rev, Contig ctg){
			if (pos1 == 0){
				pos1 = pos;
				rev1 = rev;
				ctg1 = ctg;
				return 1;
			} else if (pos2 == 0){
				pos2 = pos;
				rev2 = rev;
				ctg2 = ctg;
				return 2;
			} else {
				return -1;
			}
		}
	}
	public static class Contig implements Comparable<Contig> {
		public String name;
		public int len;
		public double cov;
		public Map<Contig,Integer> counts;
		public int numSelfConnect;
		public void setInfo(String name, int len, double cov){
			this.name = name;
			this.len = len;
			this.cov = cov;
			numSelfConnect = 0;
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
		@Override
		public int compareTo(Contig arg0) {
			return this.name.compareTo(arg0.name);
		}
	}
	
	/*public static class DWEComparator implements Comparator<DefaultWeightedEdge>{

		private Graph<Contig,DefaultWeightedEdge> g;
		
		public DWEComparator(Graph<Contig,DefaultWeightedEdge> g){
			this.g = g;
		}
		@Override
		public int compare(DefaultWeightedEdge arg0, DefaultWeightedEdge arg1) {
			if (g.getEdgeSource(arg0).equals(g.getEdgeSource(arg1)))
				return g.getEdgeTarget(arg0).compareTo(g.getEdgeTarget(arg1));
			else 
				 return g.getEdgeSource(arg0).compareTo(g.getEdgeSource(arg1));
		}
		
	}*/
	
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
	
}
