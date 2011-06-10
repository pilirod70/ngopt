package org.halophiles.assembly;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Vector;

public class ContigTerminal {
	public static final String START = "_s";
	public static final String END = "_e";
	
	private String terminus;
	
	private Contig contig;
	
	private Map<ContigTerminal,Vector<Integer>> counts;
		
	private int totalConnections;
	
	public ContigTerminal(Contig c, String term){
		contig = c;
		counts = new HashMap<ContigTerminal, Vector<Integer>>();
		terminus = term;
		totalConnections = 0;
	}
	
	public void addLink(ContigTerminal ct, int dist){
		if (counts.containsKey(ct)) {
			counts.get(ct).add(dist);
		} else {
			Vector<Integer> tmp = new Vector<Integer>();
			tmp.add(dist);
			counts.put(ct, tmp);
		}
		totalConnections++;
	}
	
	public int hashCode(){
		return (contig.name+terminus).hashCode();
	}
	
	public int getDistance(ContigTerminal ct){
		if (!counts.containsKey(ct)) return -1;
		Vector<Integer> tmp = counts.get(ct);
		int ret = 0;
		Iterator<Integer> it = tmp.iterator();
		while (it.hasNext())
			ret += it.next();
		return ret/tmp.size();
	}
	
	public int nlinks(ContigTerminal ct){
		if (counts.containsKey(ct))
			return counts.get(ct).size();
		else 
			return 0;
	}
	
	public int numConnections(){
		return totalConnections;
	}
	
	public String getName(){
		return contig.name+terminus;
	}
	
	public String toString(){
		return getName();
	}
	
	public boolean isStart(){
		return terminus.equals(START);
	}
}
