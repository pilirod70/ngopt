package org.halophiles.assembly.qc;

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
	
	public static void exportNewContigs(Map<String,Vector<MisassemblyRegion>> regions, Map<String,Contig> contigs){
		Set<SortedSet<ContigSegment>> conComps = buildGraph(regions, contigs);
		System.out.println("Found " + conComps.size()+" connected components");
	}
	
	/**
	 * Builds a graph of ContigSegments from the given MisassemblyRanges. This does so by  
	 * using the underlying MisassemblyBlocks that were used to call the MisassemblyRegions.
	 * 
	 * @param regions the MisassemblyRegions for each contig
	 * @param contigs the Contigs that this MisassemblyRanges exits
	 * @return a Set of connected components, where components are ContigSegments
	 */
	public static Set<SortedSet<ContigSegment>> buildGraph(Map<String,Vector<MisassemblyRegion>> regions, Map<String,Contig> contigs){
		Iterator<String> ctgIt = regions.keySet().iterator();
		Vector<MisassemblyRegion> reg = null;
		String ctgKey = null;
		Contig tmpCtg;
		Map<MisassemblyBlock, ContigSegment> blockMap = new HashMap<MisassemblyBlock, ContigSegment>();
		Map<ContigSegment, Vector<MisassemblyBlock>> segMap = new HashMap<ContigSegment, Vector<MisassemblyBlock>>();
		
		ContigSegment tmpSeg = null;
		MisassemblyRegion tmpReg = null;
		int start = 1;
		int end = -1;
		MisassemblyBlock blockL = null;
		MisassemblyBlock blockR = null;
		Vector<MisassemblyBlock> tmpBlocks = null;
		while (ctgIt.hasNext()) {
			ctgKey = ctgIt.next();
			reg = regions.get(ctgKey);
			tmpCtg = contigs.get(ctgKey);
			blockL = tmpCtg.getLeftBlock();
			for (int i = 0; i < reg.size(); i++){
				tmpReg = reg.get(i);
				if (tmpReg.getMinScore() < MIN_PLP_SCORE){
					// get the end of this ContigSegment, and the block at that end
					end = tmpReg.getMinPos();
					blockR = tmpReg.getLeftBlock();
					// create our ContigSegment object and this segments list of blocks
					tmpSeg = new ContigSegment(tmpCtg, start, end);
					tmpBlocks = new Vector<MisassemblyBlock>();
					segMap.put(tmpSeg, tmpBlocks);
					// add the mapping between this segment and its right block
					if (blockR == null)
						System.out.print("");
					tmpBlocks.add(blockR);
					blockMap.put(blockR, tmpSeg);
					// add the mapping between this segment and its left block
					if (blockL == null)
						System.out.print("");
					tmpBlocks.add(blockL);
					blockMap.put(blockL, tmpSeg);
					// update the left block on this segment
					blockL = tmpReg.getRightBlock();
					start = end+1;
				}
			}
			tmpSeg = new ContigSegment(tmpCtg, start, end);
			tmpBlocks = new Vector<MisassemblyBlock>();
			segMap.put(tmpSeg, tmpBlocks);
			end = tmpCtg.len;
			blockR = tmpCtg.getRightBlock();
			// add the mapping between this segment and its right block
			if (blockR == null)
				System.out.print("");
			tmpBlocks.add(blockR);
			blockMap.put(blockR, tmpSeg);
			// add the mapping between this segment and its left block
			if (blockL == null)
				System.out.print("");
			segMap.get(tmpSeg).add(blockL);
			blockMap.put(blockL, tmpSeg);
		}
		
		Iterator<ContigSegment> segIt = segMap.keySet().iterator();
		Iterator<MisassemblyBlock> blockIt = null;
		MisassemblyBlock tmpBlock = null;
		ContigSegment tmpCnct = null;
		/*
		 * for each ContigSegment, get its associated blocks, and add
		 * the ContigSegments that are connected via these blocks 
		 */
		while (segIt.hasNext()){
			tmpSeg = segIt.next();
			blockIt = getLeftBlocks(tmpSeg, segMap.get(tmpSeg)).iterator();
			while(blockIt.hasNext()){
				tmpBlock = blockIt.next();
				tmpCnct = blockMap.get(tmpBlock.getConnection());
				tmpSeg.addLeftConnection(tmpCnct, getOri(tmpBlock.getRev(),tmpBlock.getConnection().getRev()));
			}
			blockIt = getRightBlocks(tmpSeg, segMap.get(tmpSeg)).iterator();
			while(blockIt.hasNext()){
				tmpBlock = blockIt.next();
				tmpCnct = blockMap.get(tmpBlock.getConnection());
				tmpSeg.addRightConnection(tmpCnct, getOri(tmpBlock.getRev(),tmpBlock.getConnection().getRev()));
			}
		}
		
		/*
		 * BEGIN: Build connected components
		 */
		segIt = segMap.keySet().iterator();
		Set<SortedSet<ContigSegment>> connectedComponents = new HashSet<SortedSet<ContigSegment>>();
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
		if (seg == null)
			System.out.print("");
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
			if (tmp == null)
				System.out.print("");
			// if its closer to the left, call this a left block
			if (tmp.getLeft() - seg.getStart() < tmp.getRight() - seg.getEnd())
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
			if (tmp.getLeft() - seg.getStart() > tmp.getRight() - seg.getEnd())
				ret.add(tmp);
		}
		return ret;
	}
}
