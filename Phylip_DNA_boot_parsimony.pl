#!/usr/bin/perl
use warnings;
use strict;
use Cwd;

# Local variable definitions
##############################################################################
my (@phylip_file, @seqboot_out, $replicate_number, $phylip_file_nmbr, $random_nmbr, $outgroup_nmbr, $jumble);

# Initialize variable values
##############################################################################
$replicate_number = 100;   # Number of bootstraps in seqboot
$outgroup_nmbr = 1;        # Outgroup number for neighbor program
$jumble = 5;               # Number of times to jumble data
##############################################################################

# See the randon number generator.
# Combines the current time with the current process id in a
# weak attempt to come up with a random seed
srand(time|$$);
#system("rm *.txt");

&seqboot;	#####added subroutine to execute seqboot
&dnapars;	#####added subroutine to execute dnapars
&consense;	#####added subroutine to execute consense

sub seqboot {
# Commence Phylip program: seqboot
##############################################################################
# Read .phy file
@phylip_file = glob("*.phy") || print "\nNo phylip file found\n";
$phylip_file_nmbr = ($#phylip_file + 1);
if ($phylip_file_nmbr > 1) {
	die "Too many phylip files" }

# Copy .phy file - rename as infile
#link("$phylip_file[0]", "seqboot_infile.txt"); 	###***This is the original line to input***### 
system("cp $phylip_file[0] seqboot_infile.txt");	###***used this instead. changed others as well***###
system("sync"); # Finish disk write            		###***I think this has no use since copying does it***### 

#Select random number for seqboot
##############################################################################
$random_nmbr = int(rand 100000);
# If random number is even, subtract 1 to ensure odd random number
if (($random_nmbr % 2) == 0) {
	$random_nmbr = $random_nmbr - 1;
	}
#print "seqboot random number is $random_nmbr\n";

# Make input file for seqboot
open(OUT_FILE, ">", "seqboot_input.txt")
	|| die "Cannot open filehandle OUT_FILE to make input file for seqboot";
print OUT_FILE ("seqboot_infile.txt\n"); 
print OUT_FILE ("R\n");
print OUT_FILE ("$replicate_number\n");
print OUT_FILE ("Y\n");
print OUT_FILE ("$random_nmbr");
close(OUT_FILE);

# Run seqboot
system("seqboot < seqboot_input.txt > seqboot_screenout.txt");	###***change original directory***###
system("sync"); # Finish disk write

# Change names of seqboot output files
#@seqboot_out = glob("outfile") || print "No seqboot outfile present";
#link("outfile", "dnapars_infile.txt");
system("cp outfile dnapars_infile.txt");
system("sync"); # Finish disk write
#system("cp outfile dnapars_infile.txt");
rename("outfile", "seqboot_outfile.txt");
system("sync"); # Finish disk write

}

sub dnapars {

# Make input file for dnapars
##############################################################################
#Select random number for dnapars
$random_nmbr = int(rand 100000);
# If random number is even, subtract 1 to ensure odd random number
if (($random_nmbr % 2) == 0) {
	$random_nmbr = $random_nmbr - 1;
	}
#print "dnapars random number is $random_nmbr\n";

open(OUT_FILE, ">", "dnapars_input.txt")
	|| die "Cannot open filehandle OUT_FILE to make input file for dnapars";
print OUT_FILE ("dnapars_infile.txt\n");

print OUT_FILE ("S\n");
print OUT_FILE ("Y\n");

print OUT_FILE ("J\n"); # Randomize order of sequences
print OUT_FILE ("$random_nmbr\n"); # Random number for dnapars
print OUT_FILE ("$jumble\n"); # Number of times to jumble data

print OUT_FILE ("O\n"); # Use outgroup, need to fix this to choose specific
print OUT_FILE ("$outgroup_nmbr\n");

print OUT_FILE ("M\n"); # Analyze Multiple datasets
print OUT_FILE ("D\n"); # Selects Multiple datasets option
print OUT_FILE ("$replicate_number\n"); # Number of datasets

print OUT_FILE ("Y\n"); # Settings correct - ready to run program - dnapars
close(OUT_FILE);
system("sync"); # Finish disk write

# Run dnapars (use output from Seqboot as input)
system("dnapars < dnapars_input.txt > dnapars_screenout.txt");
system("sync"); # Finish disk write

# Change names of protpars output files
#link("outtree", "consense_intree.txt");
system("cp outtree consense_intree.txt");
system("sync"); # Finish disk write
#system("cp outtree consense_intree.txt");
rename("outtree", "dnapars_outtree.txt");
system("sync"); # Finish disk write
rename("outfile", "dnapars_outfile.txt");
system("sync"); # Finish disk write

}

sub consense {

# Make input file for consense
##############################################################################
open(OUT_FILE, ">", "consense_input.txt")
	|| die "Cannot open filehandle OUT_FILE to make input file for dnadist";
print OUT_FILE ("consense_intree.txt\n");
print OUT_FILE ("O\n");
print OUT_FILE ("27\n");
print OUT_FILE ("Y\n");
close(OUT_FILE);
system("sync"); # Finish disk write

# Run Consense (use tree file from dnapars as its input)
system("consense < consense_input.txt > consense_screenout.txt");
system("sync"); # Finish disk write

}

print "Finished Phylip_DNA_boot_parsimony.pl script\n";

