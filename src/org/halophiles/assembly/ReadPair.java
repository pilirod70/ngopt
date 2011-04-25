package org.halophiles.assembly;

public class ReadPair{
	public final String hdr;
	public int pos1 = 0;
	public boolean rev1 = false;
	public Contig ctg1;
	public String[] sam1;
	public int pos2 = 0;
	public boolean rev2 = false;
	public Contig ctg2;
	public String[] sam2;
	public boolean paired;
	public ReadPair(String hdr){
		this.hdr=hdr;
	}
	public String contigString(){
		return ctg1.name+"-"+ctg2.name;
	}
	public int addRead(int pos, boolean rev, Contig ctg, String[] sam){
		if (pos1 == 0){
			pos1 = pos;
			rev1 = rev;
			ctg1 = ctg;
			sam1 = sam;
			paired = false;
			return 1;
		} else if (pos2 == 0){
			if (ctg1.getId() < ctg.getId()){ // alphabetize for consistency
				pos2 = pos;
				rev2 = rev;
				ctg2 = ctg;
				sam2 = sam;
			} else {
				pos2 = pos1;
				rev2 = rev1;
				ctg2 = ctg1;
				sam2 = sam1;
				pos1 = pos;
				rev1 = rev;
				ctg1 = ctg;
				sam1 = sam;
			}
			paired = true;
			return 2;
		} else {
			return -1;
		}
	}
}
