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
		nf.setMinimumIntegerDigits(4);
		nf.setMaximumIntegerDigits(4);
		int numSeqs = 0;
		while(it.hasNext()){
			cc = it.next();
			// Find the left-most segment in this connected component
			seg = cc.first();
			Set<ContigSegment> visited = new HashSet<ContigSegment>();
			fastaWriter.addNewSequence("scaffold"+nf.format(++numSeqs));
			while (!visited.contains(seg)) {
				System.out.println(seg.toString());
				if (visited.size() > 0 && GAP_SEQUENCE.length() > 0)
					fastaWriter.writeSequence(GAP_SEQUENCE);
				if (seg.inverted()){
					fastaWriter.writeSequenceInverted(sequence.get(seg.getContig().name).subSequence(seg.getStart()-1, seg.getEnd()));
					visited.add(seg);
					if (seg.getLeftSegments().size() > 0) {
						seg = seg.getLeftSegments().firstElement();
					}
				} else {
					fastaWriter.writeSequence(sequence.get(seg.getContig().name).subSequence(seg.getStart()-1, seg.getEnd()));
					visited.add(seg);
					if (seg.getRightSegments().size() > 0) {
						seg = seg.getRightSegments().firstElement();
					}
				}
			}
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
			blockL = tmpCtg.getLeftBlock();
			for (int i = 0; i < regions.size(); i++){
				tmpRegion = regions.get(i);
				if (tmpRegion.getMinScore() < MIN_PLP_SCORE){
					// get the end of this ContigSegment, and the block at that end
					end = tmpRegion.getMinPos();
					if (end == 585329)
						System.out.print("");
					blockR = tmpRegion.getLeftBlock();
					// create our ContigSegment object and this segments list of blocks
					if (end < start)
						System.out.print("");
					tmpSeg = new ContigSegment(tmpCtg, start, end);
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
					blockL = tmpRegion.getRightBlock();
					start = end+1;
					if (start == 585330)
						System.out.print("");
				}
			}
			end = tmpCtg.len;
			if (end < start)
				System.out.print("");
			tmpSeg = new ContigSegment(tmpCtg, start, end);
			tmpBlocks = new Vector<MisassemblyBlock>();
			segMap.put(tmpSeg, tmpBlocks);
			end = tmpCtg.len;
			blockR = tmpCtg.getRightBlock();
			// add the mapping between this segment and its right block
			if (blockR != null ){
				tmpBlocks.add(blockR);
				blockMap.put(blockR, tmpSeg);
			} 
			// add theB mapping between this segment and its left block
			if (blockL != null){
				tmpBlocks.add(blockL);
				blockMap.put(blockL, tmpSeg);
			} else {
				System.out.print("");
			}
		}
		
		Iterator<ContigSegment> segIt = segMap.keySet().iterator();
		while (segIt.hasNext())
			System.out.println(segIt.next().toString());
		
		segIt = segMap.keySet().iterator();
		
		
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
			blockIt = getLeftBlocks(tmpSeg, segMap.get(tmpSeg)).iterator();
			while(blockIt.hasNext()){
				tmpBlock = blockIt.next();
				if (!blockMap.containsKey(tmpBlock.getConnection()))
					continue;
				tmpCnct = blockMap.get(tmpBlock.getConnection());
				tmpSeg.addLeftConnection(tmpCnct, getOri(tmpBlock.getRev(),tmpBlock.getConnection().getRev()));
			}
			blockIt = getRightBlocks(tmpSeg, segMap.get(tmpSeg)).iterator();
			while(blockIt.hasNext()){
				tmpBlock = blockIt.next();
				if (!blockMap.containsKey(tmpBlock.getConnection()))
					continue;
				tmpCnct = blockMap.get(tmpBlock.getConnection());
				tmpSeg.addRightConnection(tmpCnct, getOri(tmpBlock.getRev(),tmpBlock.getConnection().getRev()));
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
		if (seg.getLeftSegments().size() > 0){
			Iterator<ContigSegment> it = seg.getLeftSegments().iterator();
			while (it.hasNext())
				buildConnectedComponent(it.next(),cc);
		}
		if (seg.getRightSegments().size() > 0){
			Iterator<ContigSegment> it = seg.getRightSegments().iterator();
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
