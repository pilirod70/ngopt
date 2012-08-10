package org.halophiles.assembly.qc;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.Vector;

import org.halophiles.assembly.Contig;

public class ContigOrderer {

	public static final double MIN_PLP_SCORE = 0.9;
	
	public static void printOrder(ContigSegment start){
		
		LinkedList<ContigSegment> order = new LinkedList<ContigSegment>();
		ContigSegment curr = start;
		ContigSegment next = null;
		while(curr.getRightSegments() != null){
			if (curr.inverted())
				next = curr.getLeftSegments().get(0);
			else
				next = curr.getRightSegments().get(0);
			
			order.add(curr);
			curr = next;
		}
		ContigSegment prev = null;
		prev = curr;
		String tag = "";
		
		curr = order.pollFirst();
		while(!order.isEmpty()){
			curr = order.pollFirst();
			if (next == curr)
				tag = "circular";
			if (curr.inverted())
				System.out.println("inv["+curr+tag+"]");
			else 
				System.out.println(curr+tag);
			prev = curr;
			tag = "";
		}
	}
	
	/**
	 * Builds a graph of ContigSegments from the given MisassemblyRanges. This does so by using the underlying
	 * MisassemblyBlocks that were used to call the MisassemblyRanges
	 * 
	 * @param regions the MisassemblyRanges for each contig
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
		int left = 1;
		int right = -1;
		while (ctgIt.hasNext()) {
			ctgKey = ctgIt.next();
			reg = regions.get(ctgKey);
			tmpCtg = contigs.get(ctgKey);
			for (int i = 0; i < reg.size(); i++){
				tmpReg = reg.get(i);
				if (tmpReg.getMinScore() < MIN_PLP_SCORE){
					right = tmpReg.getMinPos();
					tmpSeg = new ContigSegment(tmpCtg, left, right);
					segMap.put(tmpSeg, new Vector<MisassemblyBlock>());
					segMap.get(tmpSeg).add(tmpReg.getLeftBlock());
					blockMap.put(tmpReg.getLeftBlock(), tmpSeg);
					if (i > 0){
						segMap.get(tmpSeg).add(reg.get(i-1).getRightBlock());
						blockMap.put(reg.get(i-1).getRightBlock(),tmpSeg);
					} else if (tmpCtg.hasEndBlock()) {
						segMap.get(tmpSeg).add(tmpCtg.getEndBlock());
						blockMap.put(tmpCtg.getEndBlock(), tmpSeg);
					}
				}
			}
		}
		
		Iterator<ContigSegment> segIt = segMap.keySet().iterator();
		Iterator<MisassemblyBlock> blockIt = null;
		MisassemblyBlock tmpBlock = null;
		ContigSegment tmpCnct = null;
		int tmpOri = -1;
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
				tmpSeg.addLeftConnection(tmpCnct, tmpOri);
			}
			blockIt = getRightBlocks(tmpSeg, segMap.get(tmpSeg)).iterator();
			while(blockIt.hasNext()){
				tmpBlock = blockIt.next();
				tmpCnct = blockMap.get(tmpBlock.getConnection());
				tmpSeg.addRightConnection(tmpCnct, tmpOri);
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
		
		// clear our visit marks
		segIt = segMap.keySet().iterator();
		while(segIt.hasNext())
			segIt.next().setVisited(false);
		
		return connectedComponents;
	}
	
	/**
	 * Recursive algorithm for 
	 * @param seg
	 * @param cc
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
