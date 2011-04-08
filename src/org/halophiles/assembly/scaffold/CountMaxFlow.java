package org.halophiles.assembly.scaffold;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.jgrapht.EdgeFactory;
import org.jgrapht.VertexFactory;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.alg.EdmondsKarpMaximumFlow;
import org.jgrapht.graph.ClassBasedEdgeFactory;
import org.jgrapht.graph.ClassBasedVertexFactory;
import org.jgrapht.graph.DirectedWeightedMultigraph;
import org.jgrapht.graph.DefaultWeightedEdge;

public class CountMaxFlow {

	private static HashMap<String,Contig> contigs;
	private static HashMap<String,ReadPair> reads;
	public static void main(String[] args){
		if (args.length % 3 != 0 || args.length == 0){
			System.err.println("Usage: java CountMaxFlow <sam1> <ins1> <err1> .... <samN> <insN> <errN>");
			System.exit(-1);
		}
		VertexFactory<Contig> vf = new ClassBasedVertexFactory<Contig>(Contig.class);
		EdgeFactory<Contig,DefaultWeightedEdge> ef = new ClassBasedEdgeFactory<Contig, DefaultWeightedEdge>(DefaultWeightedEdge.class);
		DirectedWeightedMultigraph<Contig,DefaultWeightedEdge> dg = new DirectedWeightedMultigraph<Contig, DefaultWeightedEdge>(ef);
		contigs = new HashMap<String, Contig>();
		reads = new HashMap<String, ReadPair>();
		for (int i = 0; i < args.length;i+=3){
			File samFile = new File(args[i]);
			double ins = Double.parseDouble(args[i+1]);
			double err = Double.parseDouble(args[i+2]);
			BufferedReader br = null;
			try {
				br = new BufferedReader(new FileReader(samFile));
			
				while (nextCharIs(br,'@')){
					String[] line = br.readLine().split("\t");
					Contig tmp = vf.createVertex();
					tmp.setInfo(line[1].substring(3), Integer.parseInt(line[2].substring(3)));
					contigs.put(tmp.name, tmp);
				}
				ReadPair tmp = null;
				while(br.ready()){
					String[] line = br.readLine().split("\t");
					if (reads.containsKey(line[0]))
						tmp = reads.get(line[0]);
					else {
						tmp = new ReadPair(line[0]);
						reads.put(line[0], tmp);
					}
					int val = tmp.addRead(Integer.parseInt(line[3]), isReverse(line[1]), contigs.get(line[2]));
					if (val == -1)
						System.err.println("ambiguous mapping for read " + line[0]);
				}
				Iterator<String> it = reads.keySet().iterator();
				while(it.hasNext()) {
					tmp = reads.get(it.next());
					if (tmp.ctg1.equals(tmp.ctg2)) continue;
					tmp.ctg1.addLink(tmp.ctg2);
					tmp.ctg2.addLink(tmp.ctg1);
				}
				Contig[] ctgRef = new Contig[contigs.values().size()];
				contigs.values().toArray(ctgRef);
				
				for (int cI = 0; cI < ctgRef.length; cI++){
					for (int cJ = cI+1; cJ < ctgRef.length; cJ++){
						int nlink = ctgRef[cI].nLinks(ctgRef[cJ]);
						if (nlink <= 0) continue;
						if (!dg.containsVertex(ctgRef[cI]))	dg.addVertex(ctgRef[cI]);
						if (!dg.containsVertex(ctgRef[cJ])) dg.addVertex(ctgRef[cJ]);
						
						DefaultWeightedEdge adj = dg.addEdge(ctgRef[cI], ctgRef[cJ]);
						//dg.addEdge(ctgRef[cI], ctgRef[cJ]);
						dg.setEdgeWeight(adj, nlink);
					}
				}
				
				EdmondsKarpMaximumFlow<Contig, DefaultWeightedEdge> ekmf = new EdmondsKarpMaximumFlow<Contig, DefaultWeightedEdge>(dg);
				ConnectivityInspector ci = new ConnectivityInspector<Contig, DefaultWeightedEdge>(dg);
				Iterator<Set<Contig>> ccIt = ci.connectedSets().iterator();
				System.out.println(ci.connectedSets().size());
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
							if (maxFlow <= 0.0) continue;
							System.out.println(ccCount +"\t" + ctgRef[cI] + "\t" + ctgRef[cJ] + "\t" + maxFlow);
						}
					}
				}
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
	
	public static class Adjacency{
		public Contig ctg1;
		public char end1;
		public Contig ctg2;
		public char end2;
		public Adjacency(Contig c1, char e1, Contig c2, char e2){
			ctg1 = c1;
			end1 = e1;
			ctg2 = c2;
			end2 = e2;
		}
	}
	
	public static class Contig {
		public String name;
		public int len;
		public Map<Contig,Integer> counts;
		public void setInfo(String name, int len){
			this.name = name;
			this.len = len;
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
			if (!counts.containsKey(c)) return -1;
			return counts.get(c);
		}
		public int hashCode(){
			return name.hashCode();
		}
		public String toString(){
			return name;
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
	
	private static int getBit ( int bit, int flag ) { 
        int mod = 0;  
        int dig = 0;
        while( flag != 0 && dig <= bit) {  
            mod = flag % 2;
            flag = flag / 2;
            dig++;
        }
        return mod;  
    }  
	
}
