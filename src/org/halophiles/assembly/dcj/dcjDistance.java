package org.halophiles.assembly.dcj;

public class dcjDistance {

	public static void main(String [] args)
	{
		
		
		
		String genomeX = "a,b,c,d,e,f,g,h*$";
		String genomeY = "a,c,-d,-g,-f,-e,h,b$";
		Genome x = new Genome(genomeX, "genomeX");
		Genome y = new Genome(genomeY, "genomeY");
		
		System.out.println("Printing genomeX. Should be... "+genomeX);
		x.printDesc(System.out);
		System.out.println("Printing genomeY. Should be... "+genomeY);
		y.printDesc(System.out);
		
		System.out.println("Now let's try some fun stuff..... \n");
		
		System.out.println("We should find 8 blocks, 3 cycles, and 0 odd pathes");
		int N = Genome.blockIdMap.size();
		System.out.println("Found N = " + N + " blocks");
		
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

	
}
