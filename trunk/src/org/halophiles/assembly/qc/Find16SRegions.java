package org.halophiles.assembly.qc;

import java.io.FileReader;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Scanner;
import java.util.TreeSet;

public class Find16SRegions {

	private static int MINLEN = 100;
	
	private static double MINPID = 90.0;
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		if (args.length == 0){
			System.out.println("Usage: java Find16SRegions <m8_blast_output>");
			System.exit(-1);
		}
		Comparator<Location> comp = new Comparator<Location>(){
			@Override
			public int compare(Location bh1, Location bh2) {
				if (bh1.query.equals(bh2.query)){
					if (bh1.contig.equals(bh2.contig)){
						if ((bh1.left <= bh2.right && bh1.right >= bh2.left) || 
						   ((Math.abs(bh1.right-bh2.left) < 400 && Math.abs(bh1.left-bh2.right) < 3000) || 
							(Math.abs(bh1.left-bh2.right) < 400 && Math.abs(bh1.right-bh2.left) < 3000))){
							bh1.merge(bh2);
							bh2.merge(bh1);
							return 0;
						} else if (bh1.left < bh2.left) {
							return -1;
						} else if (bh1.left > bh2.left) {
							return 1;
						} else {
							if (bh1.right < bh2.right){
								return -1;
							} else if (bh1.right > bh2.right){
								return 1;
							} else {
								return 0;
							}
						}
				 	} else {
				 		return bh1.contig.compareTo(bh2.contig);
				 	}
				} else {
					return bh1.query.compareTo(bh2.contig);
				}
			}
			
		};
		TreeSet<Location> locs = new TreeSet<Location>(comp);
		Scanner in = null;
		try {
			in = new Scanner(new FileReader(args[0]));
		} catch (Exception e){
			e.printStackTrace();
			System.exit(-1);
		}
		String[] line = null;
		int len = 0;
		double pid = 0.0;
		while(in.hasNext()) {
			line = in.nextLine().split("\t");
			len = Integer.parseInt(line[3]);
			pid = Double.parseDouble(line[2]);
			if (len >= MINLEN && pid >= MINPID){
				locs.add(new Location(line));
			}
		}
		
		Iterator<Location> it = locs.iterator();
		while(it.hasNext()){
			System.out.println(it.next().toString());
		}

	}
	
	private static class Location {
		String query;
		String contig;
		int left;
		int right;
		boolean rev;
		
		public Location(String[] hit){
			this.query = hit[0].substring(0,hit[0].indexOf("|"));
			this.contig = hit[0].substring(hit[0].indexOf("|")+1);
			int qstart = Integer.parseInt(hit[6]);
			int qend = Integer.parseInt(hit[7]);
			this.left = Math.min(qstart, qend);
			this.right = Math.max(qstart, qend);
			rev = qstart > qend;
		}
		
		public void merge(Location loc) {
			this.left = Math.min(this.left, loc.left);
			this.right = Math.max(this.right, loc.right);
			if (this.rev != loc.rev) {
				System.err.println("Merging locations from different strands.");
			}
		}
	 
		public String toString(){
			String ret = query+"|"+contig+"\t"+(rev ? right+"\t"+left : left+"\t"+right);
			return ret;
		}
	}

}
