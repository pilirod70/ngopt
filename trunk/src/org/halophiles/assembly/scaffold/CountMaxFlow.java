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
	private static int MIN_LINK = 1;
	private static int NSAMCOL = 11;
	private static HashMap<String,Contig> contigs;
	private static HashMap<String,ReadPair> reads;
	private static File outdir; 
	private static HashMap<String,Vector<Double>> mapPoints;
	
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
		for (int i = 2; i < args.length;i+=3){
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
				SAMFileParser sfp = new SAMFileParser(args[i], args[1]);
				int ins = Integer.parseInt(args[i+1]);
				double err = Double.parseDouble(args[i+2]);
				
				Iterator<Contig> ctgIt = sfp.getContigs();
				while(ctgIt.hasNext()){
					Contig tmpCtg = ctgIt.next();
					if (!contigs.containsKey(tmpCtg.name))
						contigs.put(tmpCtg.name, tmpCtg);
				}
				
				Iterator<ReadPair> rpIt = sfp.getReadPairs();
				Vector<Double> tmpPts = null;
				while(rpIt.hasNext()){
					ReadPair tmpRp = rpIt.next();
					if (!reads.containsKey(tmpRp.hdr))
						reads.put(tmpRp.hdr, tmpRp);
				}
				
		/*		Iterator<String> it = mapPoints.keySet().iterator();
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
				}*/
				
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
						if (isTerminal(tmp.pos1, tmp.pos1+rdlen-1, ins, err, tmp.ctg1) && 
							isTerminal(tmp.pos2, tmp.pos2+rdlen-1, ins, err, tmp.ctg2)){
							tmp.ctg1.addLink(tmp.ctg2);
					 		tmp.ctg2.addLink(tmp.ctg1);
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
					ctOut.println(ctgRef[cI].name+"\t"+ctgRef[cI].getNumLinkedContigs()+"\t"+ctgRef[cI].numSelfConnect);
					for (int cJ = cI+1; cJ < ctgRef.length; cJ++){
						int nlink = ctgRef[cI].nLinks(ctgRef[cJ]);
						if (nlink < MIN_LINK) 
							continue; 
						else 
							System.out.print("");
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
			map.put("nlinks", Integer.toString(g.getEdgeSource(arg0).getNumLinks(g.getEdgeTarget(arg0))));
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
	
	private static boolean isTerminal(int left, int right, int ins, double err, Contig ctg){
		if (ctg.len < (ins)*(1.0+err)) return true;
		if (left < (ins)*(1.0+err)) {
			return true;
		} else if (ctg.len - (ins)*(1.0+err) < right){
			return true;
		}
		return false;
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
}
