#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use Getopt::Long;

### TITLE --- est2tss.pl
### AUTHOR --- KRISHNAKUMAR SRIDHARAN  Revised --> 12/21/2011 (Original --> 05/01/2010)
### PURPOSE --- To extract TSS-candidate Genome Positions using GeneSeqer spliced alignment output, by capturing 5'-EST ends aligning well to a very specific and short genomic local window (~5-10bp). This tool provides the positive training data input to machine-learning algorithms to predict TSSs. It also outputs a file with the Genomic coordinate and strand information of TSS Deserts (where no EST Matches to the genome are found. Currently hard-coded as atleast 800 bp in length)
### USAGE --- perl est2tss.pl -igsq <Input filename> -from <from-posn> -to <to-posn> -strand <forward/reverse/both> ...........

my ($infile_gsq,$fromposn,$toposn,$coveragethreshold,$simithreshold,$cdnaminimal,$cdnamaximal,$winbundling,$strand,$count_threshold,$outfile,$outfile_NO_est_match);

##################### Default values list ######################
$coveragethreshold = 0.9;
$simithreshold = 0.9;
$cdnaminimal = 1;
$cdnamaximal = 10;
$winbundling = 5; 
$strand = "both";
$count_threshold = 5; 
$outfile = "tss_occurence_and_strand_file.txt";
$outfile_NO_est_match = "genome_positions_with_NO_ESTs.txt";
################################################################

GetOptions (    "igsq=s" => \$infile_gsq,			### Input filename (.gsq.out from GeneSeqer)
                "from=i" => \$fromposn,        			### Position on genomic fragment to START analysis from
                "to=i" => \$toposn,              		### Position on genomic fragment to END analysis at
                "coverage=f" => \$coveragethreshold,            ### Coverage score threshold for 5'-trancript end extraction from GSQ-out
                "similarity=f" => \$simithreshold,             	### Similarity score threshold for 5'-transcript end extraction from GSQ-out
                "startmatchmin=i" => \$cdnaminimal,         	### MINimal position of transcript where alignment with genome should start
                "startmatchmax=i" => \$cdnamaximal,            	### MAXimal position of transcript where alignment with genome should start
		"winbundling=i" => \$winbundling,		### Length of genomic fragments within which aligned 5'-transcript ends will be combined/bundled as one
            	"strand=s" => \$strand,				### Strand in which analysis will be performed (Options= 'forward', 'reverse' or 'both')
		"estcountthresh=i" => \$count_threshold,	### Number of 5'-transcript ends that should evidence a "TSS"
		"outfile=s" => \$outfile,			### Output filename that will contain the TSS positon, number of 5'-transcript ends evidenced and the strand orientation information 
 		"outNOest=s" => \$outfile_NO_est_match );	### Output filename that will contain the genomic positions and orientations of atleast-800bp-length fragments with no 5'-transcript ends aligning

if (defined $infile_gsq && defined $fromposn && defined $toposn)
 {}
else
 {
  &optionlist;
  exit 1;
 }

open (INFILE_GSQ, $infile_gsq) or die ("Could not open $infile_gsq!!");
open (OUTFILE,">$outfile") or die ("Could not open $outfile!!");
open (OUTFILE_NO_EST,">$outfile_NO_est_match") or die ("Could not open $outfile_NO_est_match!!");

my ($exonposn1,$exonposn2,$cdnaposn1,$cdnaposn2,$simiscr,$temp_exon_count,$temp_cdna_count,$exoncount,$cdnacount,$covscore,$h,$tss_posn,$i,$j,$counter,$from_blast,$to_blast,$length_of_desert);
my (%est_orientation,%orientation,@exonposns,@est_posns,@temp_array,%countall)=();

while(<INFILE_GSQ>)
{
 if(/^\s*Exon\s+1\s+(\d+)\s+(\d+)\s+\(\s+(\d+)\s+n\);\s+cDNA\s+(\d+)\s+(\d+)\s+\(\s+(\d+)\s+n\);\s+score:\s+(\d*\.\d+)\s*\n$/)
 {
  #print "$1\t$2\t$3\t$4\t$5\t$6\t$7\n";  #630519  630703  185     1       186     186     0.941
  $exonposn1=$1;
  $exonposn2=$2;
  $cdnaposn1=$4;
  $cdnaposn2=$5;
  $simiscr=$7;
  $temp_exon_count=abs($exonposn2-$exonposn1)+1;
  $temp_cdna_count=abs($cdnaposn2-$cdnaposn1)+1;
  $exoncount= abs $temp_exon_count;
  $cdnacount=abs $temp_cdna_count;
  $covscore=$cdnacount/$exoncount;                                ## Coverage score for every exon
  
  ################ TSS Desert Finder Segment of EST2TSS ##########
  if ($exonposn1 > $exonposn2)   
  {
   $est_orientation{$exonposn1} = "reverse"; ## Hash to capture of orientation of each EST alignment ### First hash "%orientation" Key ->genomicDNA posn; Value ->strand orientation in alignment
  }
  elsif ($exonposn2 > $exonposn1)
  {
   $est_orientation{$exonposn1} = "forward";
  }
  ###############################################################
  if($covscore >= $coveragethreshold && $simiscr >= $simithreshold && $cdnaposn1>=$cdnaminimal && $cdnaposn1<=$cdnamaximal && $exonposn1>=$fromposn && $exonposn1<=$toposn)
## Threshold for Exon significance, Allowing for cdna posn of alignment to be specified
  {
    if ($exonposn1 > $exonposn2)
    {
     $orientation{$exonposn1} = "reverse"; ## First hash "%orientation" Key ->genomicDNA posn; Value ->strand orientation in alignment
    }
    elsif ($exonposn2 > $exonposn1)
    {
     $orientation{$exonposn1} = "forward";
    }
    
    if ($strand ne "both" && $orientation{$exonposn1} eq $strand)  ####Choose exon posns on the forward or reverse strand only
    { 
     push (@exonposns,$exonposn1);
    }
    elsif ($strand eq "both")
    {
     push (@exonposns,$exonposn1);
    }
  }    
 }
###################### Match lines and Search for the EST ID ##################
# if(/^MATCH\s+(\S+)\s+(\S+)\s+([0-9.]+)\s+(\d+)\s+([0-9.]+)[^PCG]+([PCG])/)              ## Match lines and Search for the EST ID
# {
#  $genid=$1;
#  $estid=$2;
#  $string=chop($estid);
#  print "ORIENTATION =  $string\n";
#  $alignmentscr=$3;
# }
###############################################################################
} #####End of while(<INFILE>) loop
close INFILE_GSQ;

#########################################################################################################
#--------------------------------------------------------------------------------------------------------
# Reading the file is done. Manipulation starts.
#--------------------------------------------------------------------------------------------------------
#########################################################################################################

my @sortedexonposns = sort {$a<=>$b} @exonposns;

################ TSS Desert Finder Segment of EST2TSS ##########
@est_posns = keys %est_orientation;
@est_posns = sort {$a<=>$b} @est_posns;
for(my $k=0;$k<$#est_posns;$k++)
{
 if($est_orientation{$est_posns[$k]} eq $est_orientation{$est_posns[$k+1]} && ($est_posns[$k+1]-$est_posns[$k])>=800) ####Find GP with no ESTs, look for segments 800bp(200+400+200) b/w GPs with atleast one EST match. Assumption = You want to extract a fragment [-200,+200] wrt to TSS
 {
   $from_blast = $est_posns[$k] + 200; ### Change this if you want overlapping windows with TSS +ve data. Assumption = You want to extract a fragment [-200,+200] wrt to TSS
   $to_blast = $est_posns[$k+1] - 200; 
   $length_of_desert = $to_blast-$from_blast+1;
   print OUTFILE_NO_EST "$from_blast,$to_blast,$est_orientation{$est_posns[$k]},$length_of_desert\n";
 }
}
################################################################

###########Information-file print commands##############
#print OUTFILE " @sortedexonposns\n\t ";
#print OUTFILE " \n\n ";
print OUTFILE " \nGENOME POSITIONS\t NO. OF 5' EST ends\t\tSTRAND ORIENTATION\n ";
print OUTFILE " \n----------------\t ------------------\t\t------------------\n ";
print OUTFILE " \n----------------\t ------------------\t\t------------------\n ";
########################################################
for ($i=0;$i<=$#sortedexonposns;$i++)
{
  $counter=1;
  push(@temp_array,$sortedexonposns[$i]);
  if($i==$#sortedexonposns)
  {
   if($orientation{$sortedexonposns[$#sortedexonposns]} eq $orientation{$sortedexonposns[$#sortedexonposns-1]} && 
 	$sortedexonposns[$#sortedexonposns]>($sortedexonposns[$#sortedexonposns-1]+$winbundling)) ####The last value is out of the window in which the second-last value is present
   {
    $counter = 1;
   }
   elsif($sortedexonposns[$#sortedexonposns]<=($sortedexonposns[$#sortedexonposns-1]+$winbundling))
   {
    last;   ###It has been counted as part of the previous loop itself
   }
  }

  elsif($i<$#sortedexonposns)
  {
   if($orientation{$sortedexonposns[$i+1]} eq $orientation{$sortedexonposns[$i]} && 
				$sortedexonposns[$i+1]>($sortedexonposns[$i]+$winbundling)) ####The i value is out of the window in which the i-1 value is present
   {
    $counter = 1;
   }
   elsif($sortedexonposns[$i+1]<=($sortedexonposns[$i]+$winbundling))
   {
    for($j=1;$j<=$#sortedexonposns;$j++)
    {
     if($i+$j > $#sortedexonposns)
     {
      last;
     }
     elsif($i+$j <= $#sortedexonposns)
     {
      if($orientation{$sortedexonposns[$i+$j]} eq $orientation{$sortedexonposns[$i]})
      {
       if($sortedexonposns[$i+$j]==$sortedexonposns[$i] || $sortedexonposns[$i+$j]<=($sortedexonposns[$i]+$winbundling))   ####The i value is in the same window in which the i-1 value is present
       {
        $counter++;
        push(@temp_array,$sortedexonposns[$i+$j]);
       }
       elsif($sortedexonposns[$i+$j]>($sortedexonposns[$i]+$winbundling))
       {
        last;
       }
      }
     }
    }
   }
  }
 if(scalar(@temp_array)>1)
 {
  $tss_posn = median(\@temp_array); ###Given more than one genomic position in a "winbundling" window, choose the middle value as TSS
 }
 else
 {
  $tss_posn = $sortedexonposns[$i];
 }

 $countall{$tss_posn}=$counter;  ## Second hash "%countall" Key ->genomicDNA posn; Value->Number of occurences of that posn in significant alignments

################## TEXT HISTOGRAM PRINTER segment ###############
# print HISTOGRAM "$tss_posn\t";
# for ( $his=0;$his<=$counter;$his++)
# {
#  print HISTOGRAM "\#";
# }
# print HISTOGRAM "\n";
#################################################################

 if($countall{$tss_posn} >= $count_threshold) ####Genomic positions/fragments only above a certain threshold of 5'-EST counts are counted
 { 
  print OUTFILE " \n$tss_posn\t\t\t\t $countall{$tss_posn}\t\t\t\t $orientation{$tss_posn}\n";
 }
 #if ($i+$j > $#sortedexonposns && $sortedexonposns[$i]==$sortedexonposns[$#sortedexonposns])
 #{
 # last;
 #}
@temp_array = ();
$i=$i+($counter-1); ###Go to the next genomic "winbundling" window/position to count
}### End of For loop

################# Defining a sub-routine for MEDIAN ####################################
	sub median 
	{
	 @_ == 1 or die ('Sub usage: $median = median(\@array);');
	 my ($array_ref) = @_;
	 my $count = scalar @$array_ref;
	 # Sort a COPY of the array, leaving the original untouched
	 my @array = sort { $a <=> $b } @$array_ref;
	 if ($count % 2) 
 	 {
	  return $array[($count/2)];
	 } 
         else 
	 {
	  return ($array[($count+1)/2]);
	 }
	}
######################################################################################

#########################################################################
###################### Options display ##################################
#########################################################################
sub optionlist
{
print"
	Program : est2tss.pl (perl interface for est2tss)
	Created by : Krishnakumar Sridharan (krish28\@iastate.edu)
	Usage : perl est2tss.pl -igsq <Input filename> -from <from-posn> -to <to-posn> -strand <forward/reverse/both> ...........

        ########################################################################################
        ################################### Option List ########################################
        ########################################################################################

        -igsq			Mandatory Option. String. Name of Input file to be parsed. Processes GeneSeqer output files.

        -from           	Mandatory Option. Integer. Position on genomic fragment to START analysis from. 

        -to             	Mandatory Option. Integer. Position on genomic fragment to END analysis at.

        -strand         	String. Specifies DNA strand to be analysed. Options are 'forward'- Original Strand, 'reverse' - Complementary Strand, 'both' - Both forward and reverse strands. Default= 'both'
        
        -coverage      		Floating point. Coverage score threshold for 5'-trancript end extraction from GeneSeqer output. Default= 0.9

        -similarity 		Floating point. Similarity score threshold for 5'-transcript end extraction from GSQ-out. Default= 0.9

        -startmatchmin 		Integer. MINimal position of transcript where alignment with genome should start. Default= 1

        -startmatchmax  	Integer. MAXimal position of transcript where alignment with genome should start. Default= 10 

        -winbundling 		Integer. Length of genomic fragments within which aligned 5'-transcript ends will be combined/bundled as one. Default= 5

        -estcountthresh		Integer. Number of 5'-transcript ends that should evidence a TSS. Default= 5

        -outfile		String. Output filename that will contain the TSS positon, number of 5'-transcript ends evidenced and the strand orientation information.
												Default=tss_occurence_and_strand_file.txt 

        -outNOest		String. Output filename that will contain the genomic positions and orientations of atleast-800bp-length fragments with no 5'-transcript ends aligning.
												Default=genome_positions_with_NO_ESTs.txt

        ###################     ####################    ##################      ###############
        ########################################################################################
\n
";
}

