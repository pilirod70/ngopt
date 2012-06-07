#!/usr/bin/perl -w
#$ -cwd
#$ -S /usr/bin/perl -w
use strict;

my $appath = "/home/koadman/software/ap/";
my $velvetpath = "/home/koadman/software/velvet_0.7.51/";
my $datapath = "/home/koadman/haloassembly/hfxIllumina/";
my $scratchpath = "/state/partition1/eisenlab/";
my $PRE = "ALLPATHS";
my $SPECIES = "Haloarcula";

my $covmin=30;
my $covmax=80;
my $covstep = 10;
my $pctlargemin=0;
my $pctlargemax=100;
my $pctlargestep = 20;
my $largeinsmin=2000;
my $largeinsmax=15000;
my $largeinsstep=2000;

my $readlen = 35;
my $shortins = 500;

my $velvet_hash_min = 21;
my $velvet_hash_max = 31;
my $velvet_hash_step = 2;
my $velvet_cov_cutoff_min = 0;
my $velvet_cov_cutoff_max = 15;
my $velvet_cov_cutoff_step = 2;
my $velvet_min_pair_count_min = 20;
my $velvet_min_pair_count_max = 400;
my $velvet_min_pair_count_step = 20;

open(COMLINES, ">commands.txt");
open(SCRIPTFILE, ">scriptblah.txt");
`mkdir -p jobscripts`;
my $jobid = 1; #SGE jobs start with ID 1
for(my $cov=$covmin; $cov <= $covmax; $cov += $covstep)
{
#	for(my $pctlarge=$pctlargemin; $pctlarge <= $pctlargemax; $pctlarge += $pctlargestep)
	for(my $pctlarge=$pctlargemax; $pctlarge >= $pctlargemin; $pctlarge -= $pctlargestep)
	{
		for(my $largeins=$largeinsmin; $largeins <= $largeinsmax; $largeins += $largeinsstep)
		{
			close SCRIPTFILE;
			open(SCRIPTFILE,">jobscripts/job.$jobid.sh");
			`chmod 755 jobscripts/job.$jobid.sh`;
			print SCRIPTFILE "#!/bin/bash\n";
			print SCRIPTFILE "export PATH=\$PATH:/home/koadman/bin/\n";
			simulateReadsMoreReal($cov,$pctlarge,$largeins, $shortins);
#			assembleReadsAllPaths($cov,$pctlarge,$largeins);
			for(my $vhm = $velvet_hash_min; $vhm < $velvet_hash_max; $vhm += $velvet_hash_step){
				# only creates the velvet hash, do this once per hash size
				assembleRealReadsVelvet(1,$cov,$pctlarge,$largeins,$vhm,3,100,20);
				for(my $vcc = $velvet_cov_cutoff_min; $vcc < $velvet_cov_cutoff_max; $vcc += $velvet_cov_cutoff_step){
					for(my $vmpc = $velvet_min_pair_count_min; $vmpc < $velvet_min_pair_count_max; $vmpc += $velvet_min_pair_count_step){
						assembleRealReadsVelvet(0,$cov,$pctlarge,$largeins,$vhm,$vcc,100,$vmpc);
					}
				}
			}
			$jobid++;
			last if $pctlarge==0;	# skip unnecessary loop iters
		}
	}
}

sub simulateReadsMoreReal
{
	my $cov = shift;
	my $pctlarge = shift;
	my $largeins = shift;
	my $shortins = shift;
	my $run = "run_cov$cov.pctl$pctlarge.li$largeins";
	my $simdir = $run;
	my $cleanup = "rm -rf $simdir && mkdir -p $simdir";
	executeCommand($cleanup, "$run.moveout", "$run.moveerr");
	my $shortcov = $cov * (100-$pctlarge)/100;
	my $largecov = $cov * $pctlarge/100;
	my $largesd = $largeins/10;
	my $shortsd = $shortins/10;
	my $simreads;
	$simreads = "makeMatePairList.pl $simdir $datapath/../s_6_sequence.fastq $shortins $shortsd $shortcov 35 $datapath/H_volcanii_DS2_chr3.fasta.mapping.filtered 6359 $datapath/H_volcanii_DS2_chr1.fasta.mapping.filtered 2847757 $datapath/H_volcanii_DS2_chr2.fasta.mapping.filtered 437906 $datapath/H_volcanii_DS2_chr4.fasta.mapping.filtered 85092 $datapath/H_volcanii_DS2_chr5.fasta.mapping.filtered 635786"; 
	executeCommand($simreads, "$run.simshortout", "$run.simshorterr");
	my $mover= "mv $simdir/simfastq $simdir/reads_short.fastq";
	executeCommand($mover, ">$run.moveout", ">$run.moveerr");

	$simreads = "makeMatePairList.pl $simdir $datapath/../s_6_sequence.fastq $largeins $largesd $largecov 35 $datapath/H_volcanii_DS2_chr3.fasta.mapping.filtered 6359 $datapath/H_volcanii_DS2_chr1.fasta.mapping.filtered 2847757 $datapath/H_volcanii_DS2_chr2.fasta.mapping.filtered 437906 $datapath/H_volcanii_DS2_chr4.fasta.mapping.filtered 85092 $datapath/H_volcanii_DS2_chr5.fasta.mapping.filtered 635786"; 
	executeCommand($simreads, "$run.simlargeout", "$run.simlargeerr");
	$mover = "mv $simdir/simfastq $simdir/reads_long.fastq";
	executeCommand($mover, ">$run.moveout", ">$run.moveerr");
	executeCommand("gzip -9 reads_*.fastq", ">$run.moveout", ">$run.moveerr");
	
	# copy to scratch storage
	my $mkdir = "mkdir -p $scratchpath/$simdir";
	my $cp = "cp reads_*.fastq.gz $scratchpath/$simdir/";
	executeCommand($mkdir, ">$run.moveout", ">$run.moveerr");
	executeCommand($cp, ">$run.moveout", ">$run.moveerr");
}

sub simulateReads
{
	my $cov = shift;
	my $pctlarge = shift;
	my $largeins = shift;
	my $run = "run_cov$cov.pctl$pctlarge.li$largeins";
	my $simdir = "$PRE/$SPECIES/$run"; 
	`rm -rf $simdir`;
	my $shortcov = $cov * (100-$pctlarge)/100;
	my $largecov = $cov * $pctlarge/100;
	my $shortparams = "N=$shortins,dev=2%,C=$shortcov";
	my $largeparams = "N=$largeins,dev=10%,C=$largecov";
	my $simreads = "$appath/SimulateReads PRE=$PRE DATA=$SPECIES RUN=$run K=20 LIBRARIES=";
	$simreads .= "\"n=$readlen,$shortparams:$largeparams\"";
	$simreads .= " CONSTRUCTION=D ERROR_GENERATOR_NAME=error_templates/template_2501.1_50";
	executeCommand($simreads, "$PRE/$SPECIES/$run.simout", "$PRE/$SPECIES/$run.simerr");
	my $fasta = "$appath/Fasta PRE=$PRE DATA=$SPECIES RUN=$run";
	executeCommand($fasta, "$PRE/$SPECIES/$run.fastaout", "$PRE/$SPECIES/$run.fastaerr");

	#now simulate the large and the small separately
	$run = "run_cov$cov.pctl$pctlarge.li$largeins.short";
	$simdir = "$PRE/$SPECIES/$run"; 
	`rm -rf $simdir`;
	$simreads = "$appath/SimulateReads PRE=$PRE DATA=$SPECIES RUN=$run K=20 LIBRARIES=";
	$simreads .= "\"n=$readlen,$shortparams\"";
	$simreads .= " CONSTRUCTION=D ERROR_GENERATOR_NAME=error_templates/template_2501.1_50";
	executeCommand($simreads, "$PRE/$SPECIES/$run.simout", "$PRE/$SPECIES/$run.simerr");
	$fasta = "$appath/Fasta PRE=$PRE DATA=$SPECIES RUN=$run";
	executeCommand($fasta, "$PRE/$SPECIES/$run.fastaout", "$PRE/$SPECIES/$run.fastaerr");

	$run = "run_cov$cov.pctl$pctlarge.li$largeins.large";
	$simdir = "$PRE/$SPECIES/$run"; 
	`rm -rf $simdir`;
	$simreads = "$appath/SimulateReads PRE=$PRE DATA=$SPECIES RUN=$run K=20 LIBRARIES=";
	$simreads .= "\"n=$readlen,$largeparams\"";
	$simreads .= " CONSTRUCTION=D ERROR_GENERATOR_NAME=error_templates/template_2501.1_50";
	executeCommand($simreads, "$PRE/$SPECIES/$run.simout", "$PRE/$SPECIES/$run.simerr");
	$fasta = "$appath/Fasta PRE=$PRE DATA=$SPECIES RUN=$run";
	executeCommand($fasta, "$PRE/$SPECIES/$run.fastaout", "$PRE/$SPECIES/$run.fastaerr");

}

sub assembleReadsAllPaths
{
	my $cov = shift;
	my $pctlarge = shift;
	my $largeins = shift;
	my $run = "run_cov$cov.pctl$pctlarge.li$largeins";
	
	my $apcom = "export PATH=\$PATH:$appath; RunAssemblyThroughUnipaths PRE=$PRE DATA=$SPECIES RUN=$run K=20 ERROR_TABLE_NAME=error_rates_2501.1_50 RECOVER_GAPS=True MAX_ERRORS=2 ec_ks=\"{24,20,16}\" ec_rounds=1 REMOVE_SUSPICIOUS_READS=True REMOVE_UNBUILDABLE_READS=True";
	executeCommand($apcom, "$PRE/$SPECIES/$run.apout", "$PRE/$SPECIES/$run.aperr");

	my $apcom2 = "export PATH=\$PATH:$appath; LocalizeReads PRE=$PRE DATA=$SPECIES RUN=$run SUBDIR=s";
	executeCommand($apcom2, "$PRE/$SPECIES/$run.ap2out", "$PRE/$SPECIES/$run.ap2err");

	my $fasta = "$appath/Fasta PRE=$PRE/$SPECIES/$run/ASSEMBLIES/ DATA=s RUN=run";
	executeCommand($fasta, "$PRE/$SPECIES/$run.fasta", "$PRE/$SPECIES/$run.fasta");
}

sub assembleReadsVelvet
{
	my $cov = shift;
	my $pctlarge = shift;
	my $largeins = shift;
	my $hsize = shift;
	my $cov_cutoff = shift;
	my $min_contig_length = shift;
	my $min_pair_count = shift;
	my $run = "run_cov$cov.pctl$pctlarge.li$largeins";
	my $runlarge = "run_cov$cov.pctl$pctlarge.li$largeins.large";
	my $runshort = "run_cov$cov.pctl$pctlarge.li$largeins.short";
	my $vh="$velvetpath/velveth $PRE/$SPECIES/$run $hsize -shortPaired -fasta $PRE/$SPECIES/$runshort/reads.fasta $PRE/$SPECIES/$runlarge/reads.fasta";
	executeCommand($vh, "$PRE/$SPECIES/$run.vhout", "$PRE/$SPECIES/$run.vherr");
	my $vg="$velvetpath/velvetg $PRE/$SPECIES/$run -exp_cov $cov -ins_length $shortins -ins_length2 $largeins -cov_cutoff $cov_cutoff -min_contig_lgth $min_contig_length -min_pair_count $min_pair_count";
	executeCommand($vg, "$PRE/$SPECIES/$run.vgout", "$PRE/$SPECIES/$run.vgerr");
}

sub assembleRealReadsVelvet
{
	my $hashing_only = shift;  # pass 1 to only create the hash and not assemble
	my $cov = shift;
	my $pctlarge = shift;
	my $largeins = shift;
	my $hsize = shift;
	my $cov_cutoff = shift;
	my $min_contig_length = shift;
	my $min_pair_count = shift;
	my $simrun = "run_cov$cov.pctl$pctlarge.li$largeins";
	my $run = "$simrun/velvet.h$hsize";
	my $rundir = "$scratchpath$run";
	my $vh = "";
	my $vg = "";
	if($pctlarge==100)
	{
		$vh="$velvetpath/velveth $rundir $hsize -shortPaired -fastqgz $simrun/reads_long.fastq.gz";
		$vg="$velvetpath/velvetg $rundir -exp_cov $cov -ins_length $largeins -cov_cutoff $cov_cutoff -min_contig_lgth $min_contig_length -min_pair_count $min_pair_count";
	}elsif($pctlarge==0){
		$vh="$velvetpath/velveth $rundir $hsize -shortPaired -fastqgz $simrun/reads_short.fastq.gz";
		$vg="$velvetpath/velvetg $rundir -exp_cov $cov -ins_length $shortins -cov_cutoff $cov_cutoff -min_contig_lgth $min_contig_length -min_pair_count $min_pair_count";
	}else{
		$vh="$velvetpath/velveth $rundir $hsize -shortPaired -fastqgz $simrun/reads_short.fastq.gz $simrun/reads_long.fastq.gz";
		$vg="$velvetpath/velvetg $rundir -exp_cov $cov -ins_length $shortins -ins_length2 $largeins -cov_cutoff $cov_cutoff -min_contig_lgth $min_contig_length -min_pair_count $min_pair_count";
	}
	executeCommand($vh, "$run.vhout", "$run.vherr") if( $hashing_only==1 );
	unless( $hashing_only==1 ){
		executeCommand($vg, "$run.vgout", "$run.vgerr");
		# copy results back to shared storage
		my $cp_cmd = "cp $rundir/contigs.fa $simrun/contigs.fa.h$hsize.cc$cov_cutoff.mcl$min_contig_length.mpc$min_pair_count";
		executeCommand($cp_cmd, ">$run.cpout", ">$run.cperr");
		# delete the corefile
		executeCommand("rm core.*", ">$run.cpout", ">$run.cperr");
	}
}

sub executeCommand {
  my $command = shift;
  my $stdout_file = shift;
  my $stderr_file = shift;
  $command .= " >$stdout_file 2>$stderr_file";
  print COMLINES "$command\n";
  if(defined(fileno(SCRIPTFILE))){
    print SCRIPTFILE "$command\n";
  }else{
    print "Executing $command\n";
    my $rval = system($command);
    `echo "Exited with code $rval" >> $stderr_file` if $rval != 0;
    return $rval;
  }
  return 0;
}


