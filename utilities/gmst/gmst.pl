#!/usr/bin/perl -w
#=============================================================
# Copyright Georgia Institute of Technology, Atlanta, Georgia, USA
# Distributed by GeneProbe Inc., Atlanta, Georgia, USA
#
# GeneMarkS-T
#
# Shiyuyun Tang and Mark Borodovsky 
# "GeneMarkS-T: identification of protein coding regions in RNA
# transcripts"
# 
# Besemer J., Lomsadze A. and Borodovsky M.
# Nucleic Acids Research, 2001, Vol. 29, No. 12, 2607-2618
# "GeneMarkS: a self-training method for prediction of
#  gene starts in microbial genomes. Implications for
#  finding sequence motifs in regulatory regions."
#
#
# Please report problems to Gena Tang at stang31@gatech.edu
#
#=============================================================

use strict;
use Getopt::Long;
use FindBin qw($RealBin);
use File::Spec;

my $VERSION = "5.1 March 2014";

#------------------------------------------------
# installation settings; required directories and programs

# GeneMarkS installation directory
my $gms_dir = $RealBin;

# directory with configuration and shared files
my $shared_dir = $gms_dir;

# GeneMark.hmm gene finding program <gmhmmp>; version 2.14
my $hmm = File::Spec->catfile( $gms_dir, "gmhmmp" );

# sequence parsing tool <probuild>; version 2.11c
my $build = File::Spec->catfile( $gms_dir, "probuild" );

# directopy with heuristic models for GeneMark.hmm; version 2.0
my $heu_dir_hmm = File::Spec->catdir( $shared_dir, "heuristic_mod" );

# directopy with heuristic models for GeneMark.hmm; version 2.5
my $heu_dir_hmm_2_5 = File::Spec->catdir( $shared_dir, "heuristic_mod_2_5" );

# gibbs sampler - NCBI software; for details see README in gibbs9_95
my $gibbs = File::Spec->catfile( $gms_dir, "gibbs" );

# gibbs sampler - from http://bayesweb.wadsworth.org/gibbs/gibbs.html ; version 3.10.001  Aug 12 2009
my $gibbs3 = File::Spec->catfile( $gms_dir, "Gibbs3" );

# MetaGeneMark model for initial prediction
my $meta_model = File::Spec->catfile( $gms_dir, "MetaGeneMark_v1.mod" );

#------------------------------------------------
# installation settings; optional directories and programs

# directory with heuristic models for GeneMark; version 2.0
my $heu_dir_gm = File::Spec->catdir( $shared_dir, "heuristic_mat" );

# codon table alteration for GeneMark
my $gm_4_tbl = File::Spec->catfile( $shared_dir, "gm_4.tbl" );
my $gm_1_tbl = File::Spec->catfile( $shared_dir, "gm_1.tbl" );

#------------------------------------------------

my $cmd_line = join(" ", @ARGV );

# command line parameters

my $output = '';
my $seqfile = '';
my $order = 4;
my $order_non = 2;
my $gcode = "1";
my $shape = "partial";
my $motif = 1;
my $width = 12;
my $prestart = 6;
my $identity = 0.99;
my $maxitr = 10;
my $fixmotif = 1;
my $offover = 1;
my $strand = "both";
my $bins = 0;
my $filter = 1;
my $par = '';
my $imod = '';
my $test = '';
my $verbose = '';
my $format = "LST";
my $faa = "";
my $fnn = "";
my $ext;
my $clean = 1;
my $prok = '';

my $gibbs_version = 3;   # Specifies which version of gibbs to use
my $heuristic_version = 2; 

my $version;

# "hard coded constants"
my $BINS = "1|2|3|0";
my $SHAPE_TYPE = "linear|circular|partial";
my $GENETIC_CODE = "11|4|1";
my $STRAND = "direct|reverse|both";
my $MIN_HEURISTIC_GC = 30;
my $MAX_HEURISTIC_GC = 70;
my $OUTPUT_FORMAT = "LST|GFF";
my $MIN_LENGTH = 10000; #minimal coding sequence requried for each bin

# "soft coded constants"
my $minGC = 30;
my $maxGC = 70;
my $metaout = "meta.lst";
my $logfile  = "gms.log";
my $seq = "sequence";
my $start_prefix = "startseq.";
my $gibbs_prefix = "gibbs_out.";
my $mod_prefix = "itr_";
my $mod_suffix = ".mod";
my $hmmout_prefix = "itr_";
my $hmmout_suffix = ".lst";
my $out_name = "GeneMark_hmm.mod";
my $out_name_heuristic = "GeneMark_hmm_heuristic.mod";
my $out_suffix = "_hmm.mod";
my $out_suffix_heu = "_hmm_heuristic.mod";
my $fnn_out = "";
my $faa_out = "";

my $meta_out = "initial.meta.lst";
my $gc_out = $meta_out.".feature";

#------------------------------------------------
# support for installation key in shared directory
# $hmm = $hmm . " -T $shared_dir  ";


#------------------------------------------------
# get program name
my $gms = $0; $gms =~ s/.*\///;

my $usage =
"GeneMarkS  version $VERSION
Usage: $gms [options] <sequence file name>

input sequence in FASTA format


Output options:
(output is in current working directory)

--output    <string> output file with predicted gene coordinates by GeneMarh.hmm
            and species parameters derived by GeneMarkS-T.
            (default: <sequence file name>.lst)
 
            GeneMark.hmm can be executed independently after finishing GeneMarkS training.
            This method may be the preferable option in some situations, as it provides accesses to GeneMarh.hmm options.
--format    <string> output coordinates of predicted genes in this format
            (default: $format; supported: LST and GFF)
--fnn       create file with nucleotide sequence of predicted genes
--faa       create file with protein sequence of predicted genes
--clean     <number> delete all temporary files
            (default: $clean; supported: 1 <true> and 0 <false>)
			
Run options:

--bins      <number> number of clusters for inhomogeneous genome
            (default: $bins; supported: 0(automatic clustering),1,2,3)
--filter    <number> keep at most one prediction per sequence
            (default: $filter; supported: 1 <true> and 0 <false>)
--strand    <string> sequence strand to predict genes in
            (default: '$strand'; supported: direct, reverse and both )
--order     <number> markov chain order
            (default: $order; supported in range: >= 0)
--order_non <number> order for non-coding parameters
            (default: $order_non)
--gcode     <number> genetic code
            (default: $gcode; supported: 11, 4 and 1)
--motif     <number> iterative search for a sequence motif associated with CDS start
            (default: $motif; supported: 1 <true> and 0 <false>)
--width     <number> motif width
            (default: $width; supported in range: >= 3)
--prestart  <number> length of sequence upstream of translation initiation site that presumably includes the motif
            (default: $prestart; supported in range: >= 0)
--fixmotif  <number> the motif is located at a fixed position with regard to the start; motif could overlap start codon
            (default: $fixmotif; supported: 1 <true> and 0 <false> if this option is on, it changes the meaning of --prestart 
            option which in this case will define the distance from start codon to motif start)
--offover   <number> prohibits gene overlap
            (default: $offover; supported: 1 <true> and 0 <false>)

Combined output and run options:

--prok      to run program on prokaryotic transcripts
            (this option is the same as:  --bins 1  --filter 0  --order 2  --order_non 2  --gcode 11 --width 6  --prestart 40 --fixmotif 0)


Test/developer options:

--par      <file name> custom parameters for GeneMarkS
           (default is selected based on gcode value: 'par_<gcode>.default' )
--gibbs    <number> version of Gibbs sampler software
           (default: $gibbs_version; supported versions: 1 and 3 ) 
--test     installation test
--identity  <number> identity level assigned for termination of iterations
            (default: $identity; supported in range: >=0 and <= 1)
--maxitr    <number> maximum number of iterations
            (default: $maxitr; supported in range: >= 1)
--verbose
--version

";


if ( $#ARGV == -1 ) { print $usage; exit 0; }

# parse command line
if ( !GetOptions
  (
  	'filter=i'    => \$filter,
	'bins=i'      => \$bins,
    'output=s'    => \$output,
    'order=i'     => \$order,
    'order_non=i' => \$order_non,
    'gcode=s'     => \$gcode,
    'motif=s'     => \$motif,
    'width=i'     => \$width,
    'prestart=i'  => \$prestart,
    'identity=f'  => \$identity,
    'par=s'       => \$par,
    'test'        => \$test,
    'verbose'     => \$verbose,
    'offover'     => \$offover,
    'strand=s'    => \$strand,
    'fixmotif=i'  => \$fixmotif,
    'maxitr=i'    => \$maxitr,
    'format=s'    => \$format,
    'faa'         => \$faa,
    'fnn'         => \$fnn,
    'ext=s'       => \$ext,
    'gibbs=i'     => \$gibbs_version,
    'version'     => \$version,
	'clean=i'     => \$clean,
	'prok'        => \$prok
  )
) { exit 1; }


#------------------------------------------------
# parse/check input/settings

if ( $test ) { &SelfTest(); exit 0; }
if ( defined $version ) { print "GeneMarkS version $VERSION\n"; exit 1; }

if ( $#ARGV == -1 ) { print "Error: sequence file name is missing\n"; exit 1; }
if ( $#ARGV > 0 ) { print "Error: more than one input sequence file was specified\n"; exit 1; }

$seqfile = shift @ARGV;
&CheckFile( $seqfile, "efr" )||exit 1;

if ( $prok )
{
  $bins = 1;
  $filter = 0;
  $order = 2;
  $order_non = 2;
  $offover  = '0';
  $gcode    = "11";
  $fixmotif = '0';
  $prestart = 40;
  $width    = 6;
  $fixmotif = 0;
}

if ( $OUTPUT_FORMAT !~ /\b$format\b/ )  { print "Error: format [$format] is not supported\n"; exit 1; }
if ( $order !~ /^\d\d?$/ )              { print "Error: Markov chain order [$order] is not supported\n"; exit 1; }
if ( $GENETIC_CODE !~ /\b$gcode\b/ )    { print "Error: genetic code [$gcode] is not supported\n"; exit 1 ; }
if (($motif ne '1')&&($motif ne '0'))   { print "Error: in value of motif option\n"; exit 1 ; }
if (($clean ne '1')&&($clean ne '0'))   { print "Error: in value of clean option\n"; exit 1 ; }
if (($filter ne '1')&&($filter ne '0')) { print "Error: in value of filter option\n"; exit 1 ; }
if (($fixmotif ne '1')&&($fixmotif ne '0')) { print "Error: in value of fixmotif option\n"; exit 1 ; }
if (($offover ne '1')&&($offover ne '0')) { print "Error: in value of offover option\n"; exit 1 ; }
if ( $BINS !~ /\b$bins\b/ )             { print "Error: in value of bin option\n"; exit 1 ; }
if ( $width !~ /^\d+$/ )                { print "Error: in value of motif width option\n"; exit 1; }
if ( $prestart !~ /^\d+$/ )             { print "Error: in value of prestart option\n"; exit 1; }
if (( $identity !~ /^\d?\.?\d+$/ )||( $identity > 1 ))  { print "Error: in value of identity option\n"; exit 1 ; }
if ( $maxitr !~ /^\d+$/ )               { print "Error: in value of maximum number of itteration\n"; exit 1; }
if ( $STRAND !~ /\b$strand\b/ )         { print "Error: strand [$strand] is not supported\n"; exit 1; }
if (($gibbs_version != 1)&&($gibbs_version != 3))   { print "Error: in specification of Gibbs version\n"; exit 1 ; }

if ( !$par )
  { $par = &GetParFileName( $gcode ); }
&CheckFile( $par, "efr" )||exit 1;


if ( ! $output )
{
  if ( $format eq "LST" )      { $output = File::Spec->splitpath( $seqfile ) .".lst";  }
  elsif ( $format eq "GFF" )   { $output = File::Spec->splitpath( $seqfile ) .".gff"; }

  if ($fnn) { $fnn_out = File::Spec->splitpath( $seqfile ) .".fnn"; }
  if ($faa) { $faa_out = File::Spec->splitpath( $seqfile ) .".faa"; }
}
else
{
  if ($fnn) { $fnn_out = $output . ".fnn"; }
  if ($faa) { $faa_out = $output . ".faa"; }
}

my $work_dir = `pwd`;
chomp $work_dir;
&CheckFile( $work_dir, "dwr" )||exit 1;


#------------------------------------------------
# start/append $logfile file
open( FILE, ">$logfile" )||die( "$!, $logfile," );
close FILE;
my $time = localtime();
my $logline =
"\n\n
-------------------------------------------------
start time            : $time
working directory     : $work_dir
command line          : $cmd_line
output file with predictions : $output
input sequence        : $seqfile
number of bins        : $bins
filter                : $filter
strand                : $strand
markov chain order    : cod = $order, non = $order_non
genetic code          : $gcode
search for motif      : $motif
motif width           : $width
prestart length       : $prestart
fixed motif position  : $fixmotif
off overlap           : $offover
identity threshold    : $identity
maximum iteration     : $maxitr
gene overlap off      : $offover
strand to predict on  : $strand
GeneMarkS parameters  : $par
initial hmm model     : $imod
output format         : $format
output nucleotides    : $fnn
output proteins       : $faa
gibbs_version         : $gibbs_version

         run starts here:
";

&Log($logline);

#------------------------------------------------
# more variables/settings

# use <probuild> with GeneMarkS parameter file <$par>
$build = $build . " --par " . $par;

# set options for <gmhmmp>

# switch gene overlap off in GeneMark.hmm; for eukaryotic intron-less genomes
if ( $offover )
  { $hmm = $hmm . " -p 0"; }

# set strand to predict
if ( $strand eq "direct" )
 { $hmm = $hmm . " -s d "; }
elsif ( $strand eq "reverse" )
 { $hmm = $hmm . " -s r "; }

# to run system calls
my $command;

# iteration counter
my $itr;

# model name in iteration cycle
my $mod;

# print prediction to file
my $next;

# file with prediction from previous step
my $prev;

# file with pre start sequence
my $start_seq;

# file with gibbs results
my $gibbs_out;

# difference in consecutive prediction to compare with identity threshold
my $diff;

# list of all temp files
my @list_of_temp;

# sequence G+C content
my $GC = 0;

my $do_iterations;

#------------------------------------------------
## tmp solution: get sequence size, get minimum sequence size from --par <file>
## compare, skip iterations if short

RunSystem( "$build --clean_join $seq --seq $seqfile --log $logfile", "prepare sequence\n" );
push @list_of_temp, $seq;

my $sequence_size = -s $seq;

$command = "grep MIN_SEQ_SIZE $par";
my $minimum_sequence_size = `$command`;
$minimum_sequence_size =~ s/\s*--MIN_SEQ_SIZE\s+//;

$do_iterations = 1;

if ( $sequence_size < $minimum_sequence_size )
{
        $do_iterations = 0;   
}

Log( "do_iterations = $do_iterations\n" );



#my $newseq = $seqfile.".GC";
#push @list_of_temp, $seqfile.".GC";
if($do_iterations){

#------------------------------------------------
# clustering

#initial prediction using MetaGeneMark model
#&RunSystem( "$hmm -m $meta_model -o $meta_out -d $seqfile -b" );
push @list_of_temp, $meta_out;
push @list_of_temp, $meta_out.".fna";

#&RunSystem( "$build --stat_fasta $gc_out --seq $meta_out.fna " );
#&RunSystem( "$build --stat_fasta --seq $meta_out.fna > $gc_out " );#get GC of coding region for each sequence
&RunSystem( "$build --stat_fasta --seq $seqfile > $gc_out " );#get GC of whole sequence for each sequence
push @list_of_temp, $gc_out;


#determine bin number and range
my ($bin_num, $cutoffs, $seq_GC) = cluster($gc_out, $bins);
&Log ("bin number = $bin_num\n");
&Log ( "GC range = ".join(",",@$cutoffs)."\n" );

#------------------------------------------------
# training

my $final_model;
my @seqs;
my @models; #n models
my %handles; # file handles for n bins.
if($bin_num == 1){
	$final_model = train($seqfile);

=cut	
	#-----------------------------
	#make a GC file for the input file
	#open NEWINPUT, ">", $newseq;
	#-----------------------------
	# read input sequences
	my $FA;
	open($FA, $seqfile) or die "can't open $seqfile: $!\n";
	my %read;
	while (read_fasta_seq($FA, \%read)) {
	
		if(!exists  $seq_GC->{$read{header}}){ #no coding region in the sequence
			$seq_GC->{$read{header}} = getGC($read{seq});
		}
		#print NEWINPUT ">$read{header}\t[gc=$seq_GC->{$read{header}}]\n$read{seq}\n";
	}
	#close NEWINPUT;
	#$seqfile = $newseq; #use seq with GC for final prediction
=cut
}
else{
	#open NEWINPUT, ">", $newseq;
	#-----------------------------
	#create sequence file for each bin
	#	
	for(my $i = 1; $i <= $bin_num; ++$i){
		my $fh;
		open ($fh, ">seq_bin_$i");
		push(@seqs, "seq_bin_$i");
		#push @list_of_temp, "seq_bin_$i";
		$handles{$i} = $fh;
	}
	#-----------------------------
	# read input sequences
	my $FA;
	open($FA, $seqfile) or die "can't open $seqfile: $!\n";
	my %read;
	while (read_fasta_seq($FA, \%read)) {
		if(!exists  $seq_GC->{$read{header}}){ #no coding region in the sequence
			$seq_GC->{$read{header}} = getGC($read{seq});
		}
		#--------------------------------------
		#decide which bin the sequence belongs to 
		#
		my $bin;
		if($bin_num == 2){
			if($seq_GC->{$read{header}} <= $cutoffs->[1]){
				$bin = 1;
			}
			else{
				$bin = 2;
			}
		}
		else{
			if( $seq_GC->{$read{header}} <= $cutoffs->[1] ){
				$bin = 1;
			}
			elsif( $seq_GC->{$read{header}} <= $cutoffs->[2]){
				$bin = 2;
			}
			else{
				$bin = 3;
			}
		}
		#output to corresponding output bin file
		print {$handles{$bin}} ">$read{header}\t[gc=$seq_GC->{$read{header}}]\n$read{seq}\n";
	}
	for(my $i = 1; $i <= $bin_num; ++$i){
		close ( $handles{$i} );
	}
	#train 
	for(my $i = 1; $i <= $bin_num; ++$i){
		$models[$i-1] = train( $seqs[$i-1] );
	}
	#combine individual models to make the final model file
	$final_model = combineModel( \@models, $cutoffs);
	
}#more than one bin 

&RunSystem( "cp $final_model $out_name", "output: $out_name\n" );
push @list_of_temp, $final_model if $out_name ne $final_model;

}#if $do_iterations, input sequence is long enough.

#------------------------------------------------
# final prediction

my $format_gmhmmp = "";
$format_gmhmmp = " -f G " if  $format eq "GFF" ;

$hmm .= " -b " if $filter eq 1;
#$hmm .= " -c " if $bin_num ne 1; #for multiple clusters, use GC-specified prediction
$hmm .= " -A $faa_out " if $faa;
$hmm .= " -D $fnn_out " if $fnn;
#$hmm .= " -a $faa_out " if $faa;
#$hmm .= " -d $fnn_out " if $fnn;

my $prediction_command;
my $output_with_seq = "with_seq.out";
push @list_of_temp, "with_seq.out";


if ( $do_iterations )
{
  if ( $motif )
  {
    $command = "$hmm -r -m $out_name -o $output $format_gmhmmp $seqfile";
    &RunSystem( $command, "predict genes using native model with motif\n" );
  }
  else
  {
    # no moitf option specified
    $command = "$hmm -m $out_name -o $output $format_gmhmmp $seqfile";
    &RunSystem( $command, "predict genes using native model and no motif\n" );
  }
}
else
{
   # no iterations - use heuristic only
   $command = "$hmm -m $meta_model -o $output $format_gmhmmp $seqfile";
   &RunSystem( $command, "predict genes using heuristic\n" );
}


#------------------------------------------------
# clean temp, finish

&RunSystem( "rm -f @list_of_temp" )  if $clean; 

$time = localtime();
&Log( "End: $time\n\n" );




#=============================================================
# sub section:
#=============================================================

# -----------------------------------------------
# gms training
# -----------------------------------------------
sub train{
	my $input_seq = shift;
	#------------------------------------------------
	# prepare sequence
	&RunSystem( "$build --clean_join $seq --seq $input_seq --log $logfile", "prepare sequence\n" );
	push @list_of_temp, $seq;

	#------------------------------------------------
	# tmp solution: get sequence size, get minimum sequence size from --par <file>
	# compare, skip iterations if short

	my $sequence_size = -s $seq;
	$command = "grep MIN_SEQ_SIZE $par";
	my $minimum_sequence_size = `$command`;
	$minimum_sequence_size =~ s/\s*--MIN_SEQ_SIZE\s+//;

	$do_iterations = 1;

	if ( $sequence_size < $minimum_sequence_size )
	{
	&RunSystem( "$build --clean_join $seq --seq $input_seq --log $logfile --MIN_CONTIG_SIZE 0 --GAP_FILLER ", "prepare sequence\n" );
	$do_iterations = 0;
	} 

	&Log( "do_iterations = $do_iterations\n" );

	#------------------------------------------------
	# run initial prediction
	$itr = 0;
	$next = &GetNameForNext( $itr );
	&RunSystem( "$hmm  $seq  -m $meta_model  -o $next", "run initial prediction\n" );
	push @list_of_temp, $next;

	#------------------------------------------------
	# enter iterations loop
    #print "(Liz) metamodel ".$meta_model." next ".$next."\n";
	&Log( "entering iteration loop\n" );

	while( $do_iterations )
	{
	$itr++;
	$mod  = GetNameForMod( $itr );

	if ( $motif && !($fixmotif) )
	{
		$start_seq = $start_prefix . $itr;
		$gibbs_out = $gibbs_prefix . $itr;
	}

	$command = "$build --mkmod $mod --seq $seq --geneset $next --ORDM $order --order_non $order_non --revcomp_non 1";

	if ( $motif && !$fixmotif )
	{ $command .= " --pre_start $start_seq --PRE_START_WIDTH $prestart"; }
	elsif ( $motif && $fixmotif )
	{ $command .= " --fixmotif --PRE_START_WIDTH $prestart --width $width --log $logfile"; }

	&RunSystem( $command, "build model: $mod for iteration: $itr\n" );
	push @list_of_temp, $mod;

	if ( $motif && !$fixmotif )
	{
		if ( $gibbs_version == 1 )
		{
			&RunSystem( "$gibbs $start_seq $width -n > $gibbs_out", "run gibbs sampler\n" );
		}
		elsif ( $gibbs_version == 3 )
		{
			#&RunSystem( "$gibbs3 $start_seq $width -o $gibbs_out -F -Z  -n -r -S 20  -y -x -m  -i 2000 -w 0.01", "run gibbs3 sampler\n" );
			&RunSystem( "$gibbs3 $start_seq $width -o $gibbs_out -F -Z  -n -r -y -x -m -s 1 -w 0.01", "run gibbs3 sampler\n" );
		}
    
		push @list_of_temp, $start_seq;

		&RunSystem( "$build --gibbs $gibbs_out --mod $mod --seq $start_seq --log $logfile", "make prestart model\n" );
		push @list_of_temp, $gibbs_out;
	}

	$prev = $next;
	$next = &GetNameForNext( $itr );

	$command = "$hmm  $seq  -m $mod  -o $next";
	if ( $motif )
		{ $command .= " -r"; }

	&RunSystem( $command, "prediction, iteration: $itr\n" );
	push @list_of_temp, $next;

	$command = "$build --compare --source $next --target $prev";
	&Log( "compare:\n" . $command . "\n" );

	$diff = `$command`;
	chomp( $diff );
	&Log( "compare $prev and $next: $diff\n" );

	if ( $diff >= $identity )
		{ &Log( "Stopped iterations on identity: $diff\n" ); last; }
	if ( $itr == $maxitr )
		{ &Log( "Stopped iterations on maximum number: $maxitr\n" ); last; }
	}
	#------------------------------------------------
	# create ouput
	
	
	if ( $do_iterations )
	{
		&RunSystem( "cp $mod $input_seq.mod", "create: $input_seq.mod\n" );
		return $input_seq.".mod";
		
	}
	else
	{
	   print "WARNING: unable to train model, using initialized system instead.\n";
		&RunSystem( "cp $next $input_seq.mod", "create: $input_seq.mod\n" );
		return $next;
	}

}

# -----------------------------------------------
# cluster sequences according to GC of coding regions
# -----------------------------------------------
sub cluster{
	my $feature_f = shift; #feature file from probuild
	my $clusters = shift; #user-defiend number of bins. Default=0
	#my $out = "/home/gena/_temp_out";
	#`~alexl/DISTR/src/probuild/probuild --stat_fasta  $out   --seq $seq_file`;
	
	my %gc_hash;
	my @cut_off_points;
	my ($min_GC, $max_GC, $one_third, $two_third, $one_half);
	my $num_of_seq = 0;
	my $total_length = 0;
	my %header_to_cod_GC;

	open (GC, "<", $feature_f) or die "can't open $feature_f: $!\n";
	while (<GC>){
		next if ($_ !~ /^>(.*?)\t(\d+)\s+(\d+)/);
		my $header = $1;
		my $length = $2;
		my $GC = $3;
		if($header =~ /^(.*?)\t/){
			$header = $1;
		}
		$header_to_cod_GC{$header} = $GC;
		$num_of_seq ++;
		$total_length += $length;
		$gc_hash{$GC} += $length;
		#$gc_hash{$GC} ++;
	
	}
	close GC;
	
	my @sorted_GC = sort {$a<=>$b} keys %gc_hash;
	$min_GC = $sorted_GC[0];
	$max_GC = $sorted_GC[-1];
	&Log ( "min_GC=$min_GC  max_GC=$max_GC total_seq_length=$total_length\n" );

	my $previous = 0;
	for my $key (@sorted_GC){
		$gc_hash{$key} += $previous;
	#	if($previous < $num_of_seq/3 && $gc_hash{$key} >= $num_of_seq/3){$one_third = $key};
	#	if($previous < $num_of_seq/3*2 && $gc_hash{$key} >= $num_of_seq/3*2){$two_third = $key};
	#	if($previous < $num_of_seq/2 && $gc_hash{$key} >= $num_of_seq/2){$one_half = $key};

		if($previous < $total_length/3 && $gc_hash{$key} >= $total_length/3){$one_third = $key};
		if($previous < $total_length/3*2 && $gc_hash{$key} >= $total_length/3*2){$two_third = $key};
		if($previous < $total_length/2 && $gc_hash{$key} >= $total_length/2){$one_half = $key};
		$previous = $gc_hash{$key};
	}
		&Log ("($one_third)->($gc_hash{$one_third})\n");
		&Log ("($one_half)->($gc_hash{$one_half})\n");
		&Log ("($two_third)->($gc_hash{$two_third})\n");

	if($clusters == 0){
		#cluster number is not specified by user
		#automatically choose cluster number.
		if( $two_third - $one_third > 3){
			$clusters = 3;
		}
		else {
			$clusters = 1;
		}
	}
	if($clusters == 3){
		if($two_third - $one_third < 1 || $max_GC - $two_third < 1 || $one_third - $min_GC < 1){
			&Log( "Total number of sequences is not enough for training in 3 clusters!\n" );
			$clusters = 1;
		}
		else{
			if($gc_hash{$one_third} > $MIN_LENGTH){
				push @cut_off_points, ($min_GC,$one_third,$two_third,$max_GC);}
			else{
				&Log( "Total length of sequences is not enough for training in 3 clusters!\n" );
				$clusters = 2;
			}
		}
	}
	if($clusters == 2){
		if($gc_hash{$one_half} > $MIN_LENGTH){
			push @cut_off_points, ($min_GC,$one_half,$max_GC);}
		else{
			&Log( "Total length of sequences is not enough for training in 2 clusters!\n" );
			$clusters = 1;
		}
	}

	if($clusters == 1){
		push @cut_off_points, ($min_GC, $max_GC);
	}
	#print $clusters," ",join(" ", @cut_off_points);
	return ($clusters,\@cut_off_points,\%header_to_cod_GC);
}

#--------------------------------------------------------------------
# read a fasta file on the fly. Sequence is not stored in the memory
# ------------------------------------------------------------------
sub read_fasta_seq {
   my ($fh, $seq_info) = @_;

   $seq_info->{seq} = undef; # clear out previous sequence

   # put the header into place
   $seq_info->{header} = $seq_info->{next_header} if $seq_info->{next_header};

   my $file_not_empty = 0; 
   while (<$fh>) {
      $file_not_empty = 1;
      next if /^\s*$/;  # skip blank lines
      chomp;    

      if (/^>/) { # fasta header line
         my $h = $_;    
         $h =~ s/^>//;  
         if ($seq_info->{header}) {
            $seq_info->{next_header} = $h;
            return $seq_info;   
         }              
         else { # first time through only
            $seq_info->{header} = $h;
         }              
      }         
      else {    
         s/\s+//;  # remove any white space
         $seq_info->{seq} .= $_;
      }         
   }    

   if ($file_not_empty) {
      return $seq_info;
   }    
   else {
      # clean everything up
      $seq_info->{header} = $seq_info->{seq} = $seq_info->{next_header} = undef;

      return;   
   }    
}

# -----------------------------------------------
# concatenate model files into one in MGM format:
# starts with "__GC" and ends with "end"
# -----------------------------------------------

sub combineModel{
	my $mod = $_[0];
	my @cut_offs = @{$_[1]};
#	my ($mod, $cut_offs) = @_;
	
	#change the min and max GC value of cut_offs to the minGC and maxGC used by gmhmmp
	$cut_offs[0] = $minGC;
	$cut_offs[scalar @cut_offs - 1] = $maxGC;
	
	my $b = 1;
	open MODEL, ">final_model";
	for(my $i = $minGC; $i <=$maxGC; $i ++){
		print MODEL "__GC".$i."\t\n";
		if($i == $cut_offs[$b] || $i == $maxGC){
			open my $fh, '<', $mod->[$b-1] or die;
			#my $data = do { local $/; <$fh> };
			#close $fh;
			my $data;
			while(my $line = <$fh>){
				if($line =~ /NAME/){
					chomp $line;
					$line .= "_GC<=$i\n";
				}
				$data .= $line;
			}
			close $fh;
			print MODEL $data;
			print MODEL "end \t\n\n";
			last if( $b > scalar(@$mod));
			$b ++;
		}
	}
	close MODEL;
	return "final_model";
}

#-----------------------------------------------
# $gm_1_tbl $gm_4_tbl
#-----------------------------------------------
sub CopyCodonTableForGeneMark
{
  my ( $code ) = @_;
  my $table_name = "";

  if ( $code eq "1" )
  {
    $table_name = $gm_1_tbl;
    $table_name =~ s/.*\///;
    &RunSystem( "cat $gm_1_tbl > $table_name", "copy codon table file $table_name for GeneMark $code\n" );
  }
  elsif ( $code eq "4" )
  {
    $table_name = $gm_4_tbl;
    $table_name =~ s/.*\///;
    &RunSystem( "cat $gm_4_tbl > $table_name", "copy codon table file $table_name for GeneMark $code\n" );
  }
}
#-----------------------------------------------
# $hmmout_prefix $hmmout_suffix
#-----------------------------------------------
sub GetNameForNext
{
  my ( $name ) = @_;
  return  $hmmout_prefix . $name . $hmmout_suffix;
}
#-----------------------------------------------
# mod_prefix $mod_suffix
#-----------------------------------------------
sub GetNameForMod
{
  my ( $name ) = @_;
  return  $mod_prefix . $name . $mod_suffix;
}
#-----------------------------------------------
sub RunSystem
{
  my( $com, $text ) = @_;
  if ( $text ) { &Log( $text ); }
  &Log( $com . "\n" );
  my $err_code = 0;
  if ( $err_code = system( $com ) )
  {
    &Log( "Error on last system call, error code $err_code\nAbort program!!!\n" );
    print "GeneMarkS: error on last system call, error code $err_code\nAbort program!!!\n";
    exit(1);
  }
  else
   { &Log( "system call done\n" ); }
}
#-----------------------------------------------
# $verbose $logfile
#-----------------------------------------------
sub Log
{
  my( $text ) = @_;
  $verbose && print $text;
  open( FILE, ">>$logfile" )||die( "error: $!, $logfile," );
  print FILE $text;
  close FILE;
}
#-----------------------------------------------
# $shared_dir
# "/par_" ".default"
#-----------------------------------------------
sub GetParFileName
{
  my( $code ) = @_;
  return $shared_dir . "/par_" . $code . ".default";
}
#------------------------------------------------
# $MIN_HEURISTIC_GC $MAX_HEURISTIC_GC
# "/heu_" "_"
#------------------------------------------------
sub GetHeuristicFileName
{
  my( $GC, $code, $dir, $ext, $ord ) = @_;
  $GC = int $GC;

  if( $GC < $MIN_HEURISTIC_GC ) { $GC = $MIN_HEURISTIC_GC; }
  if( $GC > $MAX_HEURISTIC_GC ) { $GC = $MAX_HEURISTIC_GC; }

  if ( defined $ord )
  {
  	return $dir ."/heu_" . $code . "_" . $GC . "_" . $ord . $ext;
  }
  else
  {
    return $dir ."/heu_" . $code . "_" . $GC . $ext;
  }
}
#------------------------------------------------
sub CheckFile
{
  my( $name, $option ) = @_;
  my @array = split( //, $option );
  my $result = 1;
  my $i;

  foreach $i ( @array )
  {
    if ( $i eq "e" )
      { ( -e $name )||( print("error: $!, $name\n"), $result = 0 ); }
    elsif ( $i eq "d" )
      { ( -d $name )||( print("error: not a dir $!, $name\n"), $result = 0 ); }
    elsif ( $i eq "f" )
      { ( -f $name )||( print("error: not a file: $name\n"), $result = 0 ); }
    elsif ( $i eq "x" )
      { ( -x $name )||( print("error: $!, $name\n"), $result = 0 ); }
    elsif ( $i eq "r" )
      { ( -r $name )||( print("error: $!, $name\n"), $result = 0 ); }
    elsif ( $i eq "w" )
      { ( -w $name )||( print("error: $!, $name not writable\n"), $result = 0 ); }
    else
      { die( "error, no support for file test option $i\n" ); }
    if( !$result ) { last; }
  }
  return $result;
}
#------------------------------------------------
sub SelfTest
{
  print "installation test ...\n";
  print "required ...\n";
  &SelfTestRequired();
  print "optional ...\n";
  &SelfTestOptional();
  print "done\n";
}
#------------------------------------------------
# $gms_dir $heu_dir_hmm $hmm $build $gibbs
# $GENETIC_CODE $MIN_HEURISTIC_GC $MAX_HEURISTIC_GC
# ".mod"
#------------------------------------------------
sub SelfTestRequired
{
  &CheckFile( $gms_dir, "edrx" );
  &CheckFile( $shared_dir, "edr" );
  &CheckFile( $heu_dir_hmm, "edr" );
  &CheckFile( $hmm , "efx");
  &CheckFile( $build, "efx" );
  &CheckFile( $gibbs, "efx" );

  my @array = split( /\|/, $GENETIC_CODE );
  my $code  = '';
  my $GC = 0;

  foreach $code ( @array )
    { &CheckFile( &GetParFileName($code), "efr" ); }

  # test files for GeneMark.hmm  ".mod"
  foreach $code ( @array )
  {
    $GC = $MIN_HEURISTIC_GC;
    while( $GC <= $MAX_HEURISTIC_GC )
    {
      &CheckFile( &GetHeuristicFileName( $GC, $code, $heu_dir_hmm, ".mod" ), "efr" );
      $GC++;
    }
  }
}

#------------------------------------------------
# get GC of sequence
#------------------------------------------------
sub getGC {
        my ($seq)  = @_;
        chomp($seq);
        my $gc     = 0;
        my $acgt = 0;
        my $length = length($seq);
        $gc = () =  $seq =~ /g|c/ig; 
        $acgt = () = $seq =~ /a|c|g|t/ig;
        if ( $length > 0 ) {
                if($acgt != $length){
                        print STDERR "letters that are not ACGT!\n";
                        return ( "0_noncanonical_letter") if $acgt eq 0;
                        return ( ($gc/$acgt*100)."_noncanonical_letter");
                }
                return ($gc/$length*100);
        }
        else {
			return -1;
                #return "It's an empty sequence";
        }
}

