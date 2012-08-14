package org.halophiles.assembly.qc;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.Vector;

import net.sf.samtools.SAMFileReader;

import org.halophiles.assembly.Contig;
import org.halophiles.tools.HelperFunctions;

public class TestMBRefiner {
	public static void main(String[] args) {
		HelperFunctions.logInputs("MBRefiner", args);
		if (args.length == 4) {
			try {
				File bamFile = new File(args[0]);
				File bedFile = new File(args[1]);
				File connectionsFile = new File(args[2]);
				File ctgFile = new File(args[2]);
				File brokenFile = new File(args[3]);
				Map<String, Contig> contigs = MisassemblyBreaker.getContigs(new SAMFileReader(bamFile));
				Map<String, Vector<MisassemblyRegion>> ranges = MBRefiner.getRegions(bedFile, connectionsFile, contigs);
				MBRefiner.scoreAtBaseLevel(bamFile, bedFile, ctgFile, ranges, contigs);
				Map<String, int[]> junctions = MBRefiner.refine(ranges);
				MBRefiner.breakContigs(junctions, ctgFile.getAbsolutePath(), brokenFile.getAbsolutePath());

			} catch (IOException e) {
				e.printStackTrace();
				System.exit(-1);
			}	
		} else  
		/*if (args.length == 4) {
			try {
				String bedPath = args[0];
				String plpPath = args[1];
				String ctgPath = args[2];
				String brokenPath = args[3];
	
				Map<String, int[]> junctions = refine(bedPath, plpPath);
				breakContigs(junctions, ctgPath, brokenPath);
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(-1);
			}	
		} else*/ {
			System.out.println("Usage: MBRefiner <bam_file> <bed_file> <contig_file> <broken_ctgs_file>");
			System.exit(-1);
		}
	}
}
