package org.halophiles.assembly.dcj;

public class FastAccessTable {
	
	private Adjacency[] adjacencies;
	
	private int[][] locations;
	
	public FastAccessTable(Adjacency[] adjs, int[][] locs){
		adjacencies = adjs;
		locations = locs;
	}
	
	public Adjacency getHead(String block){
		return adjacencies[locations[Genome.blockIdMap.get(block)][1]];
	}
	
	public Adjacency getTail(String block){
		return adjacencies[locations[Genome.blockIdMap.get(block)][0]];
	}
	
	public Adjacency getAdjacency(String blockEnd){
		if (blockEnd.endsWith(BlockEnd.HEAD_TAG)){
			String block = blockEnd.substring(0,blockEnd.length()-2);
		 	return adjacencies[locations[Genome.blockIdMap.get(block)][1]];
		} else {
			String block = blockEnd.substring(0,blockEnd.length()-2);
			return adjacencies[locations[Genome.blockIdMap.get(block)][0]];
		}
	}
	

	/*
	private class Location {
		
		private int head;
		private int tail;
		
		public Location(int h, int t){
			head = h;
			tail = t;
		}
		
		public int getHead(){
			return head;
		}
		
		public int getTail(){
			return tail;
		}
	}
	*/
}

