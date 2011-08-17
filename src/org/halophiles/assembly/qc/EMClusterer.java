package org.halophiles.assembly.qc;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;

import org.halophiles.assembly.ReadPair;
import org.halophiles.assembly.ReadSet;

public class EMClusterer {
	private static long SEED = System.currentTimeMillis();
	
	private Random rand;
	
	private ReadPair[] reads;
	
	private double[][] P;
	
	private int[] C;
	
	private ReadSet[] gause;
	
	private double[] MU;
	
	private double[] SD;
	
	private double[] MU_last;
	
	private double[] SD_last;
	
	public EMClusterer(Collection<ReadPair> reads, int k){
		rand = new Random(SEED++);
		this.reads = new ReadPair[reads.size()];
		P = new double[reads.size()][k];
		C = new int[reads.size()];
		gause = new ReadSet[k];
		MU = new double[k];
		SD = new double[k];
		MU_last = new double[k];
		SD_last = new double[k];
		int i;
		for (i = 0; i < gause.length; i++)
			gause[i] = new ReadSet(i);
		
		Iterator<ReadPair> it = reads.iterator();
		i=0;
		int c = 0;
		while(it.hasNext()){
			this.reads[i] = it.next();
			c = rand.nextInt(k);
			C[i] = c;
			gause[c].add(this.reads[i]);
			i++;
		}
		maximize();
		/*for (i = 0; i < gause.length; i++){
			MU_last[i] = gause[i].mean();
			SD_last[i] = gause[i].sd();
		}*/
	}
	
	public int iterate(int I, double minDelta){
		/*
		System.out.print("initial: cluster0:"+gause[0].toString());
		for (int j = 1; j < gause.length; j++)
			System.out.print(" cluster"+j+gause[j].toString());
		System.out.println();
		*/
		double delta = 0;
		int i = 0;
		for (i = 0; i < I; i++){
			maximize();
			expect();
			/*
			System.out.print("iteration " + i + ": cluster0:"+gause[0].toString());
			for (int j = 1; j < gause.length; j++)
				System.out.print(" cluster"+j+gause[j].toString());
			NF.setMaximumFractionDigits(6);
			*/
			for (int j = 0; j < gause.length; j++){
				MU[j] = gause[j].mean();
				SD[j] = gause[j].sd();
			}
			delta = 0.0;
			for (int j = 0; j < gause.length; j++){
				delta += Math.abs(MU[j]-MU_last[j])/MU_last[j];
			}
			delta = delta/gause.length;
			if (delta < minDelta) {
				//System.out.print("*********");
				//System.out.println(" delta="+NF.format(delta));
				break;
			}
		}
		return i;
	}
	
	public Collection<ReadSet> getClusters(){
		Collection<ReadSet> clusters = new HashSet<ReadSet>();
		for (int i = 0; i < gause.length; i++){
			clusters.add(gause[i]);
		}
		return clusters;
	}
	
	private void expect(){
		double U = 0;
		int next = 0;
		for (int i = 0; i < reads.length; i++){
			U = rand.nextDouble();
			next = Arrays.binarySearch(P[i], U);
			if (next < 0)
				next = -(next + 1);
			if (C[i] != next){
				gause[C[i]].remove(reads[i]);
				gause[next].add(reads[i]);
				C[i] = next;
			}
		}
	}
	
	private void maximize(){
		for (int i = 0; i < gause.length; i++){
			MU_last[i] = gause[i].mean();
			SD_last[i] = gause[i].sd();
		}
		
		double total;
		for (int i = 0; i < reads.length; i++){
			total = 0.0;
			for (int j = 0; j < gause.length; j++){
				total += p(gause[j].mean(),gause[j].sd(),reads[i].getInsert());
				P[i][j] = total;	
			}
			for (int j = 0; j < gause.length; j++)
				P[i][j] = P[i][j]/total;
			
		}
	}
	
	private static double p(double mu, double sd, double x){
		double den = 2*sd*sd;
		return Math.exp(-Math.pow(mu-x,2)/den)/Math.sqrt(Math.PI*den);
	}
	
	
}
