package org.halophiles.assembly;

public class ReadPair{
	public final String hdr;
	public int pos1 = 0;
	public boolean rev1 = false;
	public int len1 = 0;
	public Contig ctg1;
	public String[] sam1;
	public int pos2 = 0;
	public boolean rev2 = false;
	public int len2 = 0;
	public Contig ctg2;
	public String[] sam2;
	public boolean paired;
	public boolean outward;
	public boolean inward;
	public ReadPair(String hdr){
		this.hdr=hdr;
	}
	public String contigString(){
		return ctg1.name+"-"+ctg2.name;
	}
	public int addRead(int pos, boolean rev, int len, Contig ctg, String[] sam){
		if (pos1 == 0){
			pos1 = pos;
			rev1 = rev;
			len1 = len;
			ctg1 = ctg;
			sam1 = sam;
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
				sam2 = sam;
			} else if (ctg1.getId() == ctg.getId()){ 
				if (pos1 > pos) { // order points for consistency
					pos2 = pos1;
					rev2 = rev1;
					len2 = len1;
					ctg2 = ctg1;
					sam2 = sam1;
					pos1 = pos;
					rev1 = rev;
					len1 = len;
					ctg1 = ctg;
					sam1 = sam;
				} else {
					pos2 = pos;
					rev2 = rev;
					len2 = len;
					ctg2 = ctg;
					sam2 = sam;				 
				}
				outward = rev1 && !rev2;
				inward = !outward;
			} else {
				pos2 = pos1;
				rev2 = rev1;
				len2 = len1;
				ctg2 = ctg1;
				sam2 = sam1;
				pos1 = pos;
				rev1 = rev;
				len1 = len;
				ctg1 = ctg;
				sam1 = sam;
			}
			paired = true;
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
		if (sam1 == null)
			return -1;
		else if (sam2 == null)
			return Integer.parseInt(this.sam1[4]) % 255;
		else 
			return ((Integer.parseInt(this.sam1[4]) % 255)+ (Integer.parseInt(this.sam2[4]) % 255))/2;
	}
	
}
