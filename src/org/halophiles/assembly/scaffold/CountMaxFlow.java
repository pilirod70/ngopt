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

import org.halophiles.assembly.Contig;
import org.halophiles.assembly.ContigTerminal;
import org.halophiles.assembly.ReadPair;
import org.halophiles.assembly.SAMFileParser;
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
	private static final int START_TERM = 1;
	private static final int END_TERM = 3;
	private static final int NO_TERM = -1;
	
	private static int MIN_LINK = 1;
	private static int NSAMCOL = 11;
	private static HashMap<String,Contig> contigs;
	private static HashMap<String,ReadPair> reads;
	private static File outdir; 
	private static HashMap<String,Vector<Double>> mapPoints;
	
	public static void main(String[] args){
		if (/*args.length % 3 != 0 ||*/args.length == 0){
			System.err.println("Usage: java CountMaxFlow <outdir> <sam1> <ins1> <err1> .... <samN> <insN> <errN>");
			System.exit(-1);
		}
		
		VertexFactory<Contig> vf = new ClassBasedVertexFactory<Contig>(Contig.class);
		EdgeFactory<ContigTerminal,DefaultWeightedEdge> ef = new ClassBasedEdgeFactory<ContigTerminal, DefaultWeightedEdge>(DefaultWeightedEdge.class);
		DirectedWeightedMultigraph<ContigTerminal,DefaultWeightedEdge> dg = new DirectedWeightedMultigraph<ContigTerminal, DefaultWeightedEdge>(ef);
		
		
		contigs = new HashMap<String, Contig>();
		reads = new HashMap<String, ReadPair>();
		outdir = new File(args[0]);
		
		mapPoints = new HashMap<String,Vector<Double>>();
		
	//	File gclFile = new File(args[1]);
		BufferedReader br = null;
		
		try {
/*			br = new BufferedReader(new FileReader(gclFile));
			br.readLine();
			while(br.ready()){
				String[] line = br.readLine().split("\t");
				Contig tmp = vf.createVertex();
				tmp.setInfo(line[0], Integer.parseInt(line[3]), Double.parseDouble(line[1]));
				contigs.put(tmp.name, tmp);
			}*/
		} catch (Exception e){
			e.printStackTrace();
			System.exit(-1);
		}
		for (int i = 1; i < args.length;i+=3){
			// create all the output files we need.
			File maxFlowFile = new File(outdir,"max_flow.txt");
			File edgeWeightsFile = new File(outdir,"edge_weights.txt");
			File cnctTypeFile = new File(outdir,"connection_counts.txt");
			File dotFile = new File(outdir,"normalized_graph.dot");
			File insDistFile = new File(outdir,"insert_distribution.txt");
			File conflictSAMFile = new File(outdir,"conflicting_pairs.sam");
			File conflictTDFile = new File(outdir,"conflicting_pairs.txt");
			try {
				maxFlowFile.createNewFile();
				edgeWeightsFile.createNewFile();
				cnctTypeFile.createNewFile();
				dotFile.createNewFile();
				insDistFile.createNewFile();
				conflictSAMFile.createNewFile();
				PrintStream mfOut = new PrintStream(maxFlowFile);
				PrintStream ewOut = new PrintStream(edgeWeightsFile);
				PrintStream ctOut = new PrintStream(cnctTypeFile);
				PrintStream dotOut = new PrintStream(dotFile);
				PrintStream idOut = new PrintStream(insDistFile);
				PrintStream cftSAMOut = new PrintStream(conflictSAMFile);
				PrintStream cftTDOut = new PrintStream(conflictTDFile);

				br = null;
				// read in SAM file
				SAMFileParser sfp = new SAMFileParser(args[i]);
				int ins = Integer.parseInt(args[i+1]);
				double err = Double.parseDouble(args[i+2]);
				
				// build up a map for contigs
				Iterator<Contig> ctgIt = sfp.getContigs();
				while(ctgIt.hasNext()){
					Contig tmpCtg = ctgIt.next();
					if (!contigs.containsKey(tmpCtg.name))
						contigs.put(tmpCtg.name, tmpCtg);
				}
				// ...and again for reads
				Iterator<ReadPair> rpIt = sfp.getReadPairs();
				while(rpIt.hasNext()){
					ReadPair tmpRp = rpIt.next();
					if (!reads.containsKey(tmpRp.hdr))
						reads.put(tmpRp.hdr, tmpRp);
				}
				
				int rdlen = sfp.getReadLength();
				Iterator<String> it = reads.keySet().iterator();
				HashSet<ReadPair> conflicts = new HashSet<ReadPair>();
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
						int term1 = isTerminal(tmp.pos1, tmp.pos1+rdlen-1, ins, err, tmp.ctg1);
						int term2 = isTerminal(tmp.pos2, tmp.pos2+rdlen-1, ins, err, tmp.ctg2);
						if (term1 != 0 && term2 != 0){
							tmp.ctg1.addLink(tmp.ctg2);
					 		tmp.ctg2.addLink(tmp.ctg1);
					 		ContigTerminal from = null;
					 		ContigTerminal to = null;
					 		int dist = 0;
					 		if (term1 == START_TERM) {
					 			from = tmp.ctg1.getStartTerminus(); 
					 			dist += tmp.pos1+rdlen-1;
					 		} else if (term1 == END_TERM) {
					 			from = tmp.ctg1.getEndTerminus(); 
					 			dist += tmp.ctg1.len - tmp.pos1 + 1;
					 		} 
					 		if (term2 == START_TERM) {
					 			to = tmp.ctg2.getStartTerminus(); 
					 			dist += tmp.pos2+rdlen-1;
					 		} else if (term2 == END_TERM ) { 
					 			to = tmp.ctg2.getEndTerminus();
					 			dist += tmp.ctg2.len - tmp.pos2 + 1;					 		
					 		}
					 			
					 		if (from != null && to != null){
					 			from.addLink(to,dist);
					 			to.addLink(from,dist);
					 		}
					 		
						} else {
							conflicts.add(tmp);
						}
					}
				}
			
				idOut.close();
				it = contigs.keySet().iterator();
				while (it.hasNext()){
					Contig tmp = contigs.get(it.next());
					cftSAMOut.println("@SQ\tSN:"+tmp.name+"\tLN:"+tmp.len);
				}
				rpIt = conflicts.iterator();
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
				
				
				Contig[] ctgRef = new Contig[contigs.values().size()];
				contigs.values().toArray(ctgRef);
				
				
				for (int cI = 0; cI < ctgRef.length; cI++){
					ContigTerminal start = ctgRef[cI].getStartTerminus();
					ContigTerminal end = ctgRef[cI].getEndTerminus();
					dg.addVertex(start);
					dg.addVertex(end);
					DefaultWeightedEdge adj = dg.addEdge(start, end);
					dg.setEdgeWeight(adj, ctgRef[cI].len);
					adj = dg.addEdge(end, start);
					dg.setEdgeWeight(adj, ctgRef[cI].len);
					
				}
				for (int cI = 0; cI < ctgRef.length; cI++){
					ctOut.println(ctgRef[cI].name+"\t"+ctgRef[cI].getNumLinkedContigs()+"\t"+ctgRef[cI].numSelfConnect);
					for (int cJ = cI+1; cJ < ctgRef.length; cJ++){
						int nlink = ctgRef[cI].nLinks(ctgRef[cJ]);
						if (nlink < MIN_LINK) 
							continue; 
						else 
							System.out.print("");
			//			if (!dg.containsVertex(ctgRef[cI]))	dg.addVertex(ctgRef[cI]);
			//			if (!dg.containsVertex(ctgRef[cJ])) dg.addVertex(ctgRef[cJ]);
			//			double L = Math.min(Math.min(ins*(1.0+err), ctgRef[cI].len), Math.min(ins*(1.0+err), ctgRef[cJ].len));
			//			double Cov = (ctgRef[cI].getCov()+ctgRef[cJ].getCov())/2;
			//			double normed = nlink*rdlen/(Cov*L);
			//			ewOut.println(normed+"\t"+nlink);
						ContigTerminal[] ctg1 = {ctgRef[cI].getStartTerminus(),ctgRef[cI].getEndTerminus()};
						ContigTerminal[] ctg2 = {ctgRef[cJ].getStartTerminus(),ctgRef[cJ].getEndTerminus()};
						for (ContigTerminal c1: ctg1) {
							for (ContigTerminal c2: ctg2) {
								if (c1.getDistance(c2) == -1 && c2.getDistance(c1) == -1) {
									continue;
								} else if (c1.getDistance(c2) == -1 || c2.getDistance(c1) == -1){
									System.err.println("unequal distances between symmetric edges between two nodes.");
									System.exit(-1);
								}
								DefaultWeightedEdge adj = dg.addEdge(c1, c2);
								dg.setEdgeWeight(adj, c1.getDistance(c2));
								adj = dg.addEdge(c2, c1);
								dg.setEdgeWeight(adj, c2.getDistance(c1));
							}
						}
//						DefaultWeightedEdge adj = dg.addEdge(ctgRef[cI], ctgRef[cJ]);
//						dg.setEdgeWeight(adj, normed);
//						adj = dg.addEdge(ctgRef[cJ],ctgRef[cI]);
//						dg.setEdgeWeight(adj, normed);
					}
				}
				ewOut.close();
				ctOut.close();
				EdmondsKarpMaximumFlow<ContigTerminal, DefaultWeightedEdge> ekmf = new EdmondsKarpMaximumFlow<ContigTerminal, DefaultWeightedEdge>(dg);
				ConnectivityInspector<ContigTerminal, DefaultWeightedEdge> ci = new ConnectivityInspector<ContigTerminal, DefaultWeightedEdge>(dg);
				Iterator<Set<ContigTerminal>> ccIt = ci.connectedSets().iterator();
				int ccCount = 0;
				ContigTerminal[] termRef = null;
				while(ccIt.hasNext()){
					Set<ContigTerminal> cc = ccIt.next();
					ccCount++;
					termRef = new ContigTerminal[cc.size()];
					cc.toArray(termRef);
					for (int cI = 0; cI < termRef.length; cI++){
						if (!dg.containsVertex(termRef[cI])) continue;
						if (!termRef[cI].isStart()) continue;
						
						for (int cJ = cI+1; cJ < ctgRef.length; cJ++){
							if (!dg.containsVertex(termRef[cJ])) continue;
							if (termRef[cJ].isStart()) continue;
							
							ekmf.calculateMaximumFlow(termRef[cI], termRef[cJ]);
							double maxFlow = ekmf.getMaximumFlowValue();
							if (maxFlow > 0.0)
								mfOut.println(ccCount+"\t"+termRef[cI]+"\t"+termRef[cJ]+"\t"+maxFlow);
						}
					}
				}
				AsUndirectedGraph<ContigTerminal,DefaultWeightedEdge> udg = new AsUndirectedGraph<ContigTerminal,DefaultWeightedEdge>(dg);
				
				Iterator<DefaultWeightedEdge> edgeIt = udg.edgeSet().iterator();
				Set<DefaultWeightedEdge> disc = new HashSet<DefaultWeightedEdge>(); 
				while(edgeIt.hasNext()){
					DefaultWeightedEdge e = edgeIt.next();
					if (disc.contains(e)) 
						continue;
					ContigTerminal src = udg.getEdgeSource(e);
					ContigTerminal tgt = udg.getEdgeTarget(e);
					if (udg.containsEdge(tgt, src)) 
						disc.add(udg.getEdge(tgt, src));
				}
				edgeIt = disc.iterator();
				while(edgeIt.hasNext())
					udg.removeEdge(edgeIt.next());
				
				System.err.println(udg.edgeSet().size()+" edges");
				VertexNameProvider<ContigTerminal> inp = new IntegerNameProvider<ContigTerminal>();
				VertexNameProvider<ContigTerminal> snp = new StringNameProvider<ContigTerminal>();
				EdgeNameProvider<DefaultWeightedEdge> enp = new DWENameProvider(udg);
				ComponentAttributeProvider<DefaultWeightedEdge> eap = new DWEAttributeProvider(udg);
				ComponentAttributeProvider<ContigTerminal> vap = new ContigAttributeProvider();
				DOTExporter<ContigTerminal,DefaultWeightedEdge> dotexp = new DOTExporter<ContigTerminal,DefaultWeightedEdge>(inp,snp,enp,vap,eap);
				dotexp.export(new PrintWriter(dotOut), udg);
				dotOut.close();
				mfOut.close();
			} catch (Exception e){
				e.printStackTrace();
				System.exit(-1);
			}
		}
	}
	
	
	
	public static class DWENameProvider implements EdgeNameProvider<DefaultWeightedEdge>{
		
		private static NumberFormat nf; 
		//private DirectedWeightedMultigraph<Contig,DefaultWeightedEdge> g;
		private Graph<ContigTerminal,DefaultWeightedEdge> g;
		public DWENameProvider(Graph<ContigTerminal,DefaultWeightedEdge> dg){
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
		private Graph<ContigTerminal,DefaultWeightedEdge> g;
		public DWEAttributeProvider(Graph<ContigTerminal, DefaultWeightedEdge> g){
			this.g = g;
		}
		@Override
		public Map<String, String> getComponentAttributes(DefaultWeightedEdge arg0) {
			Map<String,String> map = new HashMap<String,String>();
			map.put("label", g.getEdgeSource(arg0).toString()+":"+g.getEdgeTarget(arg0).toString());
			map.put("weight", Double.toString(g.getEdgeWeight(arg0)));
			map.put("nlinks", Integer.toString(g.getEdgeSource(arg0).nlinks(g.getEdgeTarget(arg0))));
			return map;
		}
	}
	
	public static class ContigAttributeProvider implements ComponentAttributeProvider<ContigTerminal>{

		@Override
		public Map<String, String> getComponentAttributes(ContigTerminal arg0) {
			Map<String,String> map = new HashMap<String,String>();
			map.put("label", arg0.getName());
			//map.put("length", Integer.toString(arg0.len));
			//map.put("cov", Double.toString(arg0.getCov()));
			return map;
		}
		
	}
	
	private static int isTerminal(int left, int right, int ins, double err, Contig ctg){
		if (ctg.len < (ins)*(1.0+err)) return START_TERM;
		if (left < (ins)*(1.0+err)) {
			return START_TERM;
		} else if (ctg.len - (ins)*(1.0+err) < right){
			return END_TERM;
		}
		return NO_TERM;
	}
	
	private static boolean spansTerminal(ReadPair pm, int ins, double err, int rdlen){
		if (isTerminal(pm.pos1, pm.pos1+rdlen-1, ins, err, pm.ctg1) == 0 || 
			isTerminal(pm.pos2, pm.pos2+rdlen-1, ins, err, pm.ctg2) == 0)
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
}
