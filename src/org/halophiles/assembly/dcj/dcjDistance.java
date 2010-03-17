package org.halophiles.assembly.dcj;

import java.io.FileReader;
import java.io.IOException;
import java.util.Scanner;

public class dcjDistance {

	public static void main(String [] args)
	{
		
		if (args.length != 1){
			System.err.println("Usage: java -jar DCJ.jar <perm_file>\n"+
					           "  where perm_file contains one line per genome");
			System.exit(-1);
		}
		
		Scanner in = null;
		try {
			in = new Scanner(new FileReader(args[0]));
			
		} catch (IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
	
		String genomeX = in.nextLine();
		String genomeY = in.nextLine();
		
		// make some toys
		//String genomeX = "a,b,c,d,e,f,g,h*$";        
		//String genomeY = "a,c,-d,-g,-f,-e,h,b$";
		
		Genome x = new Genome(genomeX, "genomeX");
		Genome y = new Genome(genomeY, "genomeY");
		System.out.println("Computing DCJ distance between " + x.getName() + " and " + y.getName());
		x.printDesc(System.out);
		y.printDesc(System.out);
		
		int N = Genome.blockIdMap.size();
		System.out.println("Found N = " + N + " blocks");
		
		Adjacency[] adj = x.getAdjacencies();
		System.out.print("Adjacencies in X: {");
		for (int i = 0; i < adj.length; i++){
			System.out.print(adj[i].toString());
			if (i<adj.length-1)
				System.out.print(",");
			else 
				System.out.print("}\n");
		}
		adj = y.getAdjacencies();
		System.out.print("Adjacencies in Y: {");
		for (int i = 0; i < adj.length; i++){
			System.out.print(adj[i].toString());
			if (i<adj.length-1)
				System.out.print(",");
			else 
				System.out.print("}\n");
		}
		
		System.out.println("Creating Adjacency graph... ");
		AdjacencyGraph agXY = new AdjacencyGraph(x, y);
		int C = agXY.countCycles();
		int I = agXY.countOddPathes();
		System.out.println("Found C = " + C + " cycles");
		System.out.println("Found I = " + I + " odd pathes");
		int dcj = N - (C + I/2);
		System.out.println("dcj(X,Y) = " + dcj);
	//	AdjacencyGraph gag = new AdjacencyGraph(x,y);
	//	proposeFullDcjPath(gag);
	}
	
	public static int computeDCJ(Genome x, Genome y){
		AdjacencyGraph agXY = new AdjacencyGraph(x, y);
		int C = agXY.countCycles();
		int I = agXY.countOddPathes();
		return Genome.blockIdMap.size() - (C + I/2);
	}

	
}
