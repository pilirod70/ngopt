package org.halophiles.assembly.qc;

import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;
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
	
	public static void buildGraph(Map<String,Vector<MisassemblyRange>> regions, Map<String,Contig> contigs){
		Iterator<String> ctgIt = regions.keySet().iterator();
		Vector<MisassemblyRange> reg = null;
		String ctgKey = null;
		Contig tmpCtg;
		Map<MisassemblyBlock, ContigSegment> blockMap = new HashMap<MisassemblyBlock, ContigSegment>();
		Map<ContigSegment, Vector<MisassemblyBlock>> segMap = new HashMap<ContigSegment, Vector<MisassemblyBlock>>();
		
		ContigSegment tmpSeg = null;
		MisassemblyRange tmpReg = null;
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
		while (segIt.hasNext()){
			tmpSeg = segIt.next();
			blockIt = getLeftBlocks(tmpSeg, segMap.get(tmpSeg)).iterator();
			while(blockIt.hasNext()){
				tmpBlock = blockIt.next();
				tmpCnct = blockMap.get(tmpBlock.getConnection());
				
				tmpSeg.addLeftConnection(tmpCnct, tmpOri);
			}
		}
		
	}
	
	public static Vector<MisassemblyBlock> getLeftBlocks(ContigSegment seg, Vector<MisassemblyBlock> blocks) {
		Vector<MisassemblyBlock> ret = new Vector<MisassemblyBlock>();
		MisassemblyBlock tmp = null;
		for (int i = 0; i < blocks.size(); i++){
			tmp = blocks.get(i);
			if (tmp.getLeft() - seg.getStart() > tmp.getRight() - seg.getEnd())
				ret.add(tmp);
		}
		return ret;
	}
	
	public static Vector<MisassemblyBlock> getRightBlocks(ContigSegment seg, Vector<MisassemblyBlock> blocks) {
		Vector<MisassemblyBlock> ret = new Vector<MisassemblyBlock>();
		MisassemblyBlock tmp = null;
		for (int i = 0; i < blocks.size(); i++){
			tmp = blocks.get(i);
			if (tmp.getLeft() - seg.getStart() > tmp.getRight() - seg.getEnd())
				ret.add(tmp);
		}
		return ret;
	}
}
