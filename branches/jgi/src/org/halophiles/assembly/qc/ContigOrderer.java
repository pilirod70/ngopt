package org.halophiles.assembly.qc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.Vector;

import org.halophiles.assembly.Contig;

public class ContigOrderer {

	public static final double MIN_PLP_SCORE = 0.9;
	
	private static String GAP_SEQUENCE = "N"; 
	
	private static HashSet<ContigSegment> allSegs = new HashSet<ContigSegment>();
	
	public static void exportNewContigs(Set<SortedSet<ContigSegment>> conComps, String contigFastaPath, String outputPath) throws IOException{
		Iterator<SortedSet<ContigSegment>> it = conComps.iterator();
		SortedSet<ContigSegment> cc = null;
		ContigSegment seg = null;
		Map<String,StringBuilder> sequence = getFastaSequence(contigFastaPath);
		File fixedContigsFile = new File(outputPath);
		fixedContigsFile.createNewFile();
		BufferedWriter bw = new BufferedWriter(new FileWriter(fixedContigsFile));
		FastaWriter fastaWriter = new FastaWriter(bw);
		NumberFormat nf = NumberFormat.getInstance();
		nf.setMinimumIntegerDigits((int)Math.ceil(Math.log10(conComps.size())));
		nf.setMaximumIntegerDigits(4);
		nf.setGroupingUsed(false);
		int numSeqs = 0;
		String header = "";
		boolean inverted = false;
		ContigSegment nextSeg = null;
		while(it.hasNext()){
			inverted = false;
			cc = it.next();
			
			Iterator<ContigSegment> ccIt = cc.iterator();
			while (ccIt.hasNext()){
				seg = ccIt.next();
				if (seg.get5PrimeSegments().isEmpty())
					break;
			}
			
			//seg = cc.first();
			Set<ContigSegment> visited = new HashSet<ContigSegment>();
			header = "scaffold"+nf.format(++numSeqs);
			fastaWriter.addNewSequence(header);
			System.out.println(header);
			while (!visited.contains(seg)) {
				System.out.println("\t"+seg.toString()+(inverted ? " (-)": " (+)"));
				allSegs.remove(seg);
				if (visited.size() > 0 && GAP_SEQUENCE.length() > 0)
					fastaWriter.writeSequence(GAP_SEQUENCE);
				visited.add(seg);
				if (inverted) {
					fastaWriter.writeSequenceInverted(sequence.get(seg.getContig().name).subSequence(seg.getStart()-1, seg.getEnd()));					
					if (seg.get5PrimeSegments().isEmpty())
						continue;
					nextSeg = seg.get5PrimeSegments().firstElement();
					if (seg.getOri(nextSeg) == MatchPoint.RR) 
						inverted = !inverted;
				} else { 
					//StringBuilder seq = sequence.get(seg.getContig().name);
					fastaWriter.writeSequence(sequence.get(seg.getContig().name).subSequence(seg.getStart()-1, seg.getEnd()));
					if (seg.get3PrimeSegments().isEmpty())
						continue;
					nextSeg = seg.get3PrimeSegments().firstElement();
					if (seg.getOri(nextSeg) == MatchPoint.FF) 
						inverted = !inverted;
				}
				seg = nextSeg;
				
			}
		}
		if (allSegs.size() > 0) {
			System.out.println("[a5_qc] "+allSegs.size() +" EXTRA SEGMENTS. Looks like there's a bug somewheres.");
		}
		Iterator<ContigSegment> segIt = allSegs.iterator();
		while(segIt.hasNext()){
			seg = segIt.next();
			System.out.println(seg.toString());
		}
		bw.close();
	}
	
	// TODO: make sure we print sequence that wasn't included during misassembly detection, e.g. contigs that had no misassemblies.
	/**
	 * 
	 * @return
	 * @throws IOException 
	 */
	private static Map<String,StringBuilder> getFastaSequence(String ctgPath) throws IOException {
		BufferedReader br = new BufferedReader(
				new FileReader(new File(ctgPath)));
		StringBuilder sb = null;
		br.read();
		String tmpCtg;
		Map<String,StringBuilder> sequences = new HashMap<String,StringBuilder>();
		while (br.ready()) {
			tmpCtg = br.readLine();
			// Take only the first part of the contig header to be consistent with bwa
			if (tmpCtg.contains(" "))
				tmpCtg = tmpCtg.substring(0, tmpCtg.indexOf(" "));
			sb = new StringBuilder();
			sequences.put(tmpCtg, sb);
			char c = (char) br.read();
			// read in this sequence
			while (c != '>') {
				if (MisassemblyBreaker.isNuc(c))
					sb.append(c);
				if (!br.ready())
					break;
				c = (char) br.read();
			}
		}	
		return sequences;
	}
	
	/**
	 * Builds a graph of ContigSegments from the given MisassemblyRegions. This does so by  
	 * using the underlying MisassemblyBlocks that were used to call the MisassemblyRegions.
	 * 
	 * @param regionMap the MisassemblyRegions for each contig
	 * @param contigs the Contigs that this MisassemblyRanges exits
	 * @return a Set of connected components, where components are ContigSegments
	 */
	public static Set<SortedSet<ContigSegment>> buildGraph(Map<String,Vector<MisassemblyRegion>> regionMap, Map<String,Contig> contigs){
		Iterator<String> ctgIt = regionMap.keySet().iterator();
		Vector<MisassemblyRegion> regions = null;
		String ctgKey = null;
		Contig tmpCtg;
		// A Map to associate ContigSegments with blocks
		Map<MisassemblyBlock, ContigSegment> blockMap = new HashMap<MisassemblyBlock, ContigSegment>();
		Map<ContigSegment, Vector<MisassemblyBlock>> segMap = new HashMap<ContigSegment, Vector<MisassemblyBlock>>();
		Set<SortedSet<ContigSegment>> connectedComponents = new HashSet<SortedSet<ContigSegment>>();
		
		ContigSegment tmpSeg = null;
		MisassemblyRegion tmpRegion = null;
		int start = 1;
		int end = -1;
		MisassemblyBlock blockL = null;
		MisassemblyBlock blockR = null;
		Vector<MisassemblyBlock> tmpBlocks = null;
		/*
		 * For each contig, slice the contig up into segments delimited by misassembly junctions
		 */
		while (ctgIt.hasNext()) {
			ctgKey = ctgIt.next();
			start = 1;
			tmpCtg = contigs.get(ctgKey); 
			regions = regionMap.get(ctgKey);
			blockL = tmpCtg.get5PrimeBlock();
			for (int i = 0; i < regions.size(); i++){
				tmpRegion = regions.get(i);
				if (tmpRegion.getMinScore() < MIN_PLP_SCORE){
					// get the end of this ContigSegment, and the block at that end
					end = tmpRegion.getMinPos();
					blockR = tmpRegion.getForwardBlock();
					// create our ContigSegment object and this segments list of blocks
					if (end < start)
						System.out.print("");
					tmpSeg = new ContigSegment(tmpCtg, start, end);
					allSegs.add(tmpSeg);
					tmpBlocks = new Vector<MisassemblyBlock>();
					segMap.put(tmpSeg, tmpBlocks);
					// add the mapping between this segment and its right block
					if (blockR != null){
						tmpBlocks.add(blockR);
						blockMap.put(blockR, tmpSeg);
					} 
					/*
					 * If a block exists, add a mapping between this segment
					 * and its left block. It should only be non-existent if
					 * the contig is not connected to anything.
					 */
					if (blockL != null){
						tmpBlocks.add(blockL);
						blockMap.put(blockL, tmpSeg);
					}
					// now get the next left block
					blockL = tmpRegion.getReverseBlock();
					start = end+1;                                                                                                    
				}
			}
			end = tmpCtg.len;
			if (end < start)
				System.out.print("");
			tmpSeg = new ContigSegment(tmpCtg, start, end);
			allSegs.add(tmpSeg);
			tmpBlocks = new Vector<MisassemblyBlock>();
			segMap.put(tmpSeg, tmpBlocks);
			end = tmpCtg.len;
			blockR = tmpCtg.get3PrimeBlock();
			// add the mapping between this segment and its right block
			if (blockR != null ){
				tmpBlocks.add(blockR);
				blockMap.put(blockR, tmpSeg);
			} 
			// add theB mapping between this segment and its left block
			if (blockL != null){
				tmpBlocks.add(blockL);
				blockMap.put(blockL, tmpSeg);
			}
		}
		
		Iterator<ContigSegment> segIt = segMap.keySet().iterator();
		/*
		while (segIt.hasNext())
			System.out.println(segIt.next().toString());	
		segIt = segMap.keySet().iterator();
		*/
		
		
		Iterator<MisassemblyBlock> blockIt = null;
		MisassemblyBlock tmpBlock = null;
		ContigSegment tmpCnct = null;
		/*
		 * for each ContigSegment, get its associated blocks, and add
		 * the ContigSegments that are connected via these blocks 
		 */
		while (segIt.hasNext()){
			tmpSeg = segIt.next();
			if (!segMap.containsKey(tmpSeg))
				continue;
			Vector<MisassemblyBlock> vect = segMap.get(tmpSeg);
			vect = getLeftBlocks(tmpSeg, vect); 
			blockIt = vect.iterator();
			while(blockIt.hasNext()){
				tmpBlock = blockIt.next();
				if (!blockMap.containsKey(tmpBlock.getConnection()))
					continue;
				tmpCnct = blockMap.get(tmpBlock.getConnection());
				tmpSeg.add5PrimeConnection(tmpCnct, getOri(tmpBlock.getRev(),tmpBlock.getConnection().getRev()));
			}
			
			vect = segMap.get(tmpSeg);
			vect = getRightBlocks(tmpSeg, vect); 
			blockIt = vect.iterator();
			while(blockIt.hasNext()){
				tmpBlock = blockIt.next();
				if (!blockMap.containsKey(tmpBlock.getConnection()))
					continue;
				tmpCnct = blockMap.get(tmpBlock.getConnection());
				tmpSeg.add3PrimeConnection(tmpCnct, getOri(tmpBlock.getRev(),tmpBlock.getConnection().getRev()));
			}
		}
		
		/*
		 * BEGIN: Build connected components
		 */
		segIt = segMap.keySet().iterator();
		while (segIt.hasNext()) {
			tmpSeg = segIt.next();
			if (tmpSeg.getVisited())
				continue;
			SortedSet<ContigSegment> cc = new TreeSet<ContigSegment>();
			buildConnectedComponent(tmpSeg, cc);
			connectedComponents.add(cc);
		}
		/*
		 * END: Build connected components
		 */
		
		// clear the visit marks
		segIt = segMap.keySet().iterator();
		while(segIt.hasNext())
			segIt.next().setVisited(false);
		
		return connectedComponents;
	}
	
	private static int getOri(boolean rev1, boolean rev2){
		if (rev1)
			if (rev2)
				return MatchPoint.RR;
			else
				return MatchPoint.RF;
		else
			if (rev2)
				return MatchPoint.FR;
			else
				return MatchPoint.FF;
	}
	
	/**
	 * Recursive algorithm for building connected components 
	 * @param seg the segment that needs to be added to this connected component
	 * @param cc the set of components that are connected to build upon
	 */
	private static void buildConnectedComponent(ContigSegment seg, SortedSet<ContigSegment> cc){ 
		if (seg.getVisited())
			return;
		seg.setVisited(true);
		cc.add(seg);
		if (seg.get5PrimeSegments().size() > 0){
			Iterator<ContigSegment> it = seg.get5PrimeSegments().iterator();
			while (it.hasNext())
				buildConnectedComponent(it.next(),cc);
		}
		if (seg.get3PrimeSegments().size() > 0){
			Iterator<ContigSegment> it = seg.get3PrimeSegments().iterator();
			while (it.hasNext())
				buildConnectedComponent(it.next(),cc);
		}
	}
	
	private static Vector<MisassemblyBlock> getLeftBlocks(ContigSegment seg, Vector<MisassemblyBlock> blocks) {
		Vector<MisassemblyBlock> ret = new Vector<MisassemblyBlock>();
		MisassemblyBlock tmp = null;
		for (int i = 0; i < blocks.size(); i++){
			tmp = blocks.get(i);
			// if its closer to the left, call this a left block
			if (tmp.getLeft() - seg.getStart() < seg.getEnd() - tmp.getRight())
				ret.add(tmp);
		}
		return ret;
	}
	
	private static Vector<MisassemblyBlock> getRightBlocks(ContigSegment seg, Vector<MisassemblyBlock> blocks) {
		Vector<MisassemblyBlock> ret = new Vector<MisassemblyBlock>();
		MisassemblyBlock tmp = null;
		for (int i = 0; i < blocks.size(); i++){
			tmp = blocks.get(i);
			// if its closer to the right, call this a right block
			if (tmp.getLeft() - seg.getStart() > seg.getEnd() - tmp.getRight())
				ret.add(tmp);
		}
		return ret;
	}
	
	public static void setGapSequence(int numberOfNs){
		byte[] ar = new byte[numberOfNs];
		for (int i = 0; i < ar.length; i++)
			ar[i] = 'N';
		GAP_SEQUENCE = new String(ar);
	}
}
