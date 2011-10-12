package org.halophiles.assembly;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Map;
import java.util.Vector;

import org.halophiles.tools.SummaryStats;

public class ReadPair{
	public static boolean INWARD;
	public static boolean OUTWARD;
	
	public final String hdr;
	public int pos1 = 0;
	public boolean rev1 = false;
	public int len1 = 0;
	public int qual1 = -1;
	public String cig1;
	public Contig ctg1;
	public int pos2 = 0;
	public boolean rev2 = false;
	public int len2 = 0;
	public String cig2;
	public Contig ctg2;
	public int qual2 = -1;
	public boolean paired;
	public boolean outward;
	public boolean inward;
	private boolean endSpanning;
	private int ins = -1;
	
	public ReadPair(String hdr){
		this.hdr=hdr;
	}
	public String contigString(){
		return ctg1.name+"-"+ctg2.name;
	}
	public int addRead(int pos, boolean rev, int len, Contig ctg, int qual, String cig){
		if (pos1 == 0){
			pos1 = pos;
			rev1 = rev;
			len1 = len;
			ctg1 = ctg;
			cig1 = cig;
			qual1 = qual;
			paired = false;
			outward = false;
			inward = false;
			return 1;
		} else if (pos2 == 0){
			if (ctg1.getId() < ctg.getId()){ // alphabetize for consistency
				pos2 = pos;
				rev2 = rev;
				len2 = len;
				ctg2 = ctg;
				cig2 = cig;
				qual2 = qual;
			} else if (ctg1.getId() == ctg.getId()){ 
				if (pos1 > pos) { // order points for consistency
					pos2 = pos1;
					rev2 = rev1;
					len2 = len1;
					ctg2 = ctg1;
					cig2 = cig1;
					qual2 = qual1;
					pos1 = pos;
					rev1 = rev;
					len1 = len;
					ctg1 = ctg;
					cig1 = cig;
					qual1 = qual;
				} else {
					pos2 = pos;
					rev2 = rev;
					len2 = len;
					ctg2 = ctg;
					cig2 = cig;
					qual2 = qual;
				}
				outward = rev1 && !rev2;
				inward = !outward;
			} else {
				pos2 = pos1;
				rev2 = rev1;
				len2 = len1;
				ctg2 = ctg1;
				cig2 = cig1;
				qual2 = qual1;
				pos1 = pos;
				rev1 = rev;
				len1 = len;
				ctg1 = ctg;
				cig1 = cig;
				qual1 = qual;
			}
			paired = true;
			if (!ctg1.equals(ctg2)){
				System.out.print("");
			}
			return 2;
		} else {
			return -1;
		}
	}
	
	public void setEndSpanning(boolean endSpanning){
		this.endSpanning = endSpanning;
	}
	
	public int getInsert(){
		if (ins > -1){
			return ins;
		}
		if (this.ctg1 != null && this.ctg2 != null && this.ctg1.equals(this.ctg2)){
			if (endSpanning) {
				ins = this.ctg1.len-this.pos2+this.pos1+this.len1;
			} else 
				ins = this.pos2-this.pos1+this.len2;
			return ins;
		} else 
			return -1;
	}
	
	public int getQual(){
		if (pos1 == 0)
			return -1;
		else if (pos2 == 0)
			return qual1 % 255;
		else 
			return ((qual1 % 255)+ (qual2 % 255))/2;
	}
	
	public String toString(){
		/*
		String ret = null;
		if (pos1 != 0)
			ret = hdr+"/1\t"+(rev1?"-1":"1")+"\t"+ctg1.name+"\t"+pos1+"\t"+cig1+"\n";
		if (pos2 != 0)
			ret += hdr+"/2\t"+(rev2?"-1":"1")+"\t"+ctg2.name+"\t"+pos2+"\t"+cig2;
		return ret;
		*/
		return hdr;
	}
	
	public static void exportInsertSize(Collection<ReadPair> reads, PrintStream out) {
		Iterator<ReadPair> it = reads.iterator();
		ReadPair tmp = null;
		while(it.hasNext()){
			tmp = it.next();
			if (tmp.paired && tmp.ctg1.equals(tmp.ctg2))
				out.println(tmp.getInsert()+"\t"+(tmp.inward?1:-1));
		}
	}
	
	public static double[] estimateInsertSize(Collection<ReadPair> reads){
		Iterator<ReadPair> it = reads.iterator();
		ReadPair tmp = null;
		Vector<Double> vals = new Vector<Double>();
		while(it.hasNext()){
			tmp = it.next();
			if (tmp.paired && tmp.ctg1.equals(tmp.ctg2))
				vals.add(new Double(tmp.getInsert()));
		}
		Double[] arD = vals.toArray(new Double[vals.size()]);
		Arrays.sort(arD);
		// discard the upper and lower quartiles to prevent any outliers from distorting our estimate
		double[] dat = new double[arD.length];
		for (int i = 0; i < dat.length; i++)
			dat[i] = arD[i];
		double mean = SummaryStats.mean(dat);
		double stdev = Math.sqrt(SummaryStats.variance(dat,mean));
		double[] ret = {Math.round(mean),Math.round(stdev),dat.length};
		return ret;
	}
	/**
	 * 
	 * @param reads
	 * @param ins
	 * @param nSd
	 * @return a reference to the Map that was filtered
	 */
	public static Map<String,ReadPair> filterTailPairs(Map<String,ReadPair> reads, double alpha){
		ReadPair[] ar = new ReadPair[reads.size()];
		reads.values().toArray(ar);
		Arrays.sort(ar,new Comparator<ReadPair>() {
			public int compare(ReadPair o1, ReadPair o2) {
				return o1.getInsert() - o2.getInsert();
			}
		});
		
		Vector<String> rm = new Vector<String>();
		int l = (int) ((alpha/2) * ar.length);
		int r = (int) ((1-(alpha/2)) * ar.length);
		for (int i = 0; i < l; i++)
			rm.add(ar[i].hdr);
		for (int i = r; i < ar.length; i++)
			rm.add(ar[i].hdr);
		removeKeys(reads,rm);
		return reads;
	}
	

	public static double[] estimateInsertSizeIQR(Map<String,ReadPair> reads){
		Iterator<ReadPair> it = reads.values().iterator();
		ReadPair tmp = null;
		Vector<Double> vals = new Vector<Double>();
		while(it.hasNext()){
			tmp = it.next();
			if (tmp.paired && tmp.ctg1.equals(tmp.ctg2))
				vals.add(new Double(tmp.getInsert()));
		}
		Double[] arD = vals.toArray(new Double[vals.size()]);
		Arrays.sort(arD);
		// discard the upper and lower quartiles to prevent any outliers from distorting our estimate
		double[] dat = new double[(arD.length*8)/10];
		int j = arD.length/10;
		for (int i = 0; i < dat.length; i++)
			dat[i] = arD[j++];
		double mean = SummaryStats.mean(dat);
		double stdev = Math.sqrt(SummaryStats.variance(dat,mean));
		double[] ret = {Math.round(mean),Math.round(stdev),dat.length};
		return ret;
	}

	
	/**
	 * 
	 * @param reads
	 * @param ins
	 * @param nSd
	 * @return a reference to the Map that was filtered
	 */
	public static Map<String,ReadPair> filterDiagReads(Map<String,ReadPair> reads, double[] ins, int nSd){
		Iterator<String> it = reads.keySet().iterator();
		Vector<String> rm = new Vector<String>();
		String key = null;
		ReadPair r = null;
		while(it.hasNext()){
			key = it.next();
			r = reads.get(key);
			if (isDiag(r, ins, nSd)){
				rm.add(key);
			}
		}
		removeKeys(reads,rm);
		return reads;
	}
	/**
	 * 
	 * @param reads
	 * @param ins
	 * @param nSd
	 * @return a reference to the Map that was filtered
	 */
	public static Map<String,ReadPair> filterRange(Map<String,ReadPair> reads, double min, double max){
		Iterator<String> it = reads.keySet().iterator();
		Vector<String> rm = new Vector<String>();
		String key = null;
		ReadPair r = null;
		double ins = 0.0;
		while(it.hasNext()){
			key = it.next();
			r = reads.get(key);
			ins = r.getInsert();
			if (ins > min && ins < max){
				rm.add(key);
			}
		}
		removeKeys(reads,rm);
		return reads;
	}
	
	public static int getOrientation(Collection<ReadPair> reads){
		Iterator<ReadPair> it = reads.iterator();
		double NIN = 0;
		double NOUT = 0;
		ReadPair tmp = null;
		while(it.hasNext()){
			tmp = it.next();
			if (tmp.inward)
				NIN++;
			else if (tmp.outward)
				NOUT++;
		}
		if (NIN > NOUT)
			return 0;
		else
			return 1;
	}
	
	private static boolean isDiag(ReadPair r, double[] ins, int nSd){
		double dist = Math.abs(r.pos1 - r.pos2 );
		return (dist < ins[0]+nSd*ins[1] && dist > ins[0]-nSd*ins[1]) && r.ctg1.equals(r.ctg2);
	}
	
	private static <V> void removeKeys(Map<String,V> reads, Collection<String> torm){
		Iterator<String> it = torm.iterator();
		while(it.hasNext())
			reads.remove(it.next());
	}
}
