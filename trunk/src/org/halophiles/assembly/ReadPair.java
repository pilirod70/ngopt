package org.halophiles.assembly;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import org.halophiles.tools.SummaryStats;

public class ReadPair{
	public final String hdr;
	public int pos1 = 0;
	public boolean rev1 = false;
	public int len1 = 0;
	public int qual1 = -1;
	public String cig1;
	public Contig ctg1;
//	public String[] sam1;
	public int pos2 = 0;
	public boolean rev2 = false;
	public int len2 = 0;
	public String cig2;
	public Contig ctg2;
	public int qual2 = -1;
//	public String[] sam2;
	public boolean paired;
	public boolean outward;
	public boolean inward;
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
	
	public int getInsert(){
		if (this.ctg1 != null && this.ctg2 != null && this.ctg1.equals(this.ctg2))
			return this.pos2-this.pos1+this.len2;
		else 
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
		String ret = null;
		if (pos1 != 0)
			ret = hdr+"/1\t"+(rev1?"-1":"1")+"\t"+ctg1.name+"\t"+pos1+"\t"+cig1+"\n";
		if (pos2 != 0)
			ret += hdr+"/2\t"+(rev2?"-1":"1")+"\t"+ctg2.name+"\t"+pos2+"\t"+cig2;
		return ret;
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
	
}
