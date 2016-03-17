#!/usr/bin/perl -w
use strict;

### snp2phylip.pl - converts SNP data input ###
### to blocks in PHYLIP format for analysis ###
### Author: Jeffrey A. Detras 		        ###
### Email: j.detras@irri.org		        ###
### Date: May 27, 2012 			            ###

### open input file and store in array ###
my $file = glob("*new.txt");
#my $file = $ARGV[0]; #imputed chromosome SNP data
my $row = 0; #
my @germplasm; #array to store SNP data input
open SNP, "< $file" or die "File missing! $!";
while (my $line = readline * SNP) {
	chomp $line;
	push(@germplasm, [split(/\t/,$line)]);
	$row++ #count number of rows 
}
close SNP;

### set gff variables ###
my $time = localtime; #timestamp the GFF3 output file
my $chr = substr($file, 0, 5); #get the chromosome number
my $source = "oryzaSNP_DB"; #source name for GFF3
my $feature = "SNP"; #feature name for GFF#
my $chrStart = 1; #set the start position for formatting
my $chrEnd = 16758485; #set the end position for formatting 

### statistics from the run ###
my $stats = "$chr".".stats";
open STATS, "> $stats" or die "$!";
print STATS "SIZE\tBLK_IN\tBLK_OUT\tTOTAL\n";

### set variables ###
my ($line, $germplasm, $sequence);
my (@snpGermplasm, @column);
my ($countGermplasm, $blockAllIn, $blockAllOut, $blockAll) = (0, 0, 0, 0);
my ($h, $j); #array reference call
my $block = "00001"; #initial block assignment
my $rejectedBlock = 1; #flag for block with > block filter size
my @windowSize = ("06", "12", "18", "24", "30"); #window block sizes
my @blockFilter = (100000, 200000, 300000, 400000, 500000); #set limit for a window block
my $last = $#windowSize; #number of array values
my $endPosition = 0;
my $i;

for ($i=0;$i<=$last;$i++){ #iteration based on window size
	### output file for GFF info ###
	my $gffOut = "$chr\_"."$windowSize[$i]".".gff3"; #output GFF3 file
	open GFF, "> $gffOut" or die "$!";

	### metadata information on GFF3 output ###
	print GFF "##gff-version 3\n";
	print GFF "##$time\n";
	print GFF "##sequence-region $chr $chrStart $chrEnd\n";

	### output file for info on discarded blocks ###
	my $discard = "$chr\_"."$windowSize[$i]"."_discard".".gff3";
	open DSCRD, "> $discard" or die "$!";

	### metada information on GFF3 output for discarded blocks ###
	print DSCRD "##gff-version 3\n";
	print DSCRD "##$time\n";
	print DSCRD "##sequence-region $chr $chrStart $chrEnd\n";
	
	my $firstSnp = 1; #first SNP in a window block
	my $lastSnp = $windowSize[$i]; #last SNP in a window block
	my $mid = $lastSnp/2; #midpoint SNP size
	
	until ($endPosition > $chrEnd) { #format onlyuntil specified position
	#until ($lastSnp >= $row) { #if doing whole chromosome
	
	my $startPosition = $germplasm[$firstSnp][2]; #start position of a block
	$endPosition = $germplasm[$lastSnp][2]; #end position of a block
	my $blockLength = $endPosition - $startPosition;#length of SNP window 
	
	if ($blockLength <= $blockFilter[$i]) { #accept block if less than filter
		### open file for output of formatted SNP block ###
		my $outfile = "seqboot_infile_"."$chr\_"."$windowSize[$i]\_"."$block".".phy";
		open PHY, "> $outfile" or die "$!";
		print PHY "27 $windowSize[$i]\n";
		
		### transposition and formatting of SNPs ###
		for $h (3 .. $#{$germplasm[0]}) {
			for $j (0, $firstSnp .. $lastSnp) {
				if ($j == $lastSnp) {
					print PHY "$germplasm[$j][$h]\n";
				} elsif ($germplasm[$j][$h] eq ""){
					print PHY "\n";
					last;
				} else {
					print PHY "$germplasm[$j][$h]";
				}
			}
		}
		
		### print out GFF file for blocks ###
		print GFF "$chr\t$source\t$feature\t$startPosition\t".
			"$endPosition\t.\t+\t.\tID=$chr\_$windowSize[$i]\_$block\n";
	
		### create folder for output file by window and block ###
		system("mkdir -p $chr\/$chr\_$windowSize[$i]\/$chr\_$windowSize[$i]\_$block");
		system("mv $outfile $chr\/$chr\_$windowSize[$i]\/$chr\_$windowSize[$i]\_$block");
		### set next block window parameters ###
		$block++;
		$firstSnp = $firstSnp + $mid;
		$lastSnp = $lastSnp + $mid;	
		#$startPosition = $germplasm[$firstSnp][2]; #start position of a block
		#$endPosition = $germplasm[$lastSnp][2]; #end position of a block

	} else { ###for rejected blocks###
		### print out GFF file for discarded blocks ###
		print DSCRD "$chr\t$source\t$feature\t$startPosition\t".
			"$endPosition\t.\t+\t.\tID=$chr\_$windowSize[$i]\_$block\n";
		
		$rejectedBlock++; #counts rejected block
		$firstSnp = $firstSnp + $mid; #assign next first SNP 
		$lastSnp = $lastSnp + $mid; #assign next last SNPi
		#$startPosition = $germplasm[$firstSnp][2]; #start position of a block
		#$endPosition = $germplasm[$lastSnp][2]; #end position of a block

	}
	close PHY;
	}
	
	### remove extra count ###
	$block = $block - 1;
	$rejectedBlock = $rejectedBlock - 1;
	
	### statistics on blocks ###
	my $blockTotal = $block + $rejectedBlock;
	$blockAllIn = $blockAllIn + $block;
	$blockAllOut = $blockAllOut + $rejectedBlock;
	$blockAll = $blockAllIn + $blockAllOut;
	print STATS "$windowSize[$i]\t$block\t$rejectedBlock\t$blockTotal\n";
	
	### reset to initial for next window size iteration
	$block = "00001"; #reassign block number 
	$endPosition = 0; #reset end position
}
system("mkdir $chr\_gff3");
system("mv *.gff3 $chr\_gff3");

print STATS "00\t$blockAllIn\t$blockAllOut\t$blockAll\n";

close GFF;
close DSCRD;
close STATS;

exit;
