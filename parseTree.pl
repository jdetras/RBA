#!/usr/bin/perl -w
use strict;
use Bio::TreeIO;

#Usage: ./parseTree.pl "outtree file"
#Requires "germplasm_lines" 

##########################################
#Calculate the distance between two nodes#
#and make a matrix of distance values    #
#Author: Jeffrey A. Detras		         #
#Date: March 5, 2010-March 9, 2010	     #
##########################################

#Input file: Germplasm lines and Newick-formatted tree

my $file = 'germplasm_lines'; #text file listing germplasm lines 
my @germplasm = ();
my %group = ();
open GERM, "< $file" or die "file missing! $!";
while (my $line = readline*GERM) {
	#chomp $line;
	push @germplasm, $line; #store germplasm in array
}

#my $treefile = $ARGV[0];
#my $treefile = glob ("*bs-parsconpars-outtree.txt"); #newick-formatted tree
my $treefile = glob ("*usertree_outtree.txt"); #newick-formatted tree
#my $out = $treefile;
my $tri = $treefile;
#$out =~ s/.txt/_distance_list.txt/;
$tri =~ s/.txt/_distance_matrix.txt/;
#open DIST, "> $out" or die "Cannot create file $!"; #print out tab values of distances 
open TRI, "> $tri" or die "Cannot create file $!"; #print out distances in matrix form
print TRI "Variety\t";

my $lastelement = 'ZS97B_impu'; #last array element of germplasm lines
foreach my $line (@germplasm) {
	$line =~ s/\r//g;
	chomp($line);
	if ($line =~ /$lastelement/) {
		print TRI "$line\n";
	} else {
		print TRI "$line\t";
	}
}

my $treeio = Bio::TreeIO->new(-format => 'newick',
				-file => $treefile);
my @node_list = ();
while (my $tree = $treeio->next_tree) {	
	my @taxa = $tree->get_leaf_nodes; #get the varieties or leaf nodes
	for my $vars (@germplasm) { #outer loop to compare germplasm
		for my $nodes ($tree->get_leaf_nodes) { #inner loop to get node value
			my $id = $nodes->id; 
			if ($vars =~ /$id/) { #compare germplasm with the node id
				push (@node_list, $nodes); #store germplasm as node 
			}
		}
	}
	my @variety = @node_list; #duplicate the node list
	my $k = 0; #zero for upper and lower matrix; 1 if upper matrix only
	my $row = 26; #assign second to the last row for matrix; 25 if upper matrix only
	my $col = 26; #assign last column for the matrix
	my $tab = "\t";
	my $node = '';
	my $other_node = '';
	my %var_dist = ();
	for my $i (0..$row) {
		for my $j ($k..$col) {
			my $distance = $tree->distance(-nodes => [$node_list[$i], $variety[$j]]);
			$node = $node_list[$i]->id; 
			$other_node = $variety[$j]->id;
			push @{$var_dist{$node}}, $distance;
			#print DIST "$node", "\t", "$other_node", "\t", "$distance\n";
		}
		#$k++; #matrix delimiter to get upper triangle only; comment to get upper and lower
	}
	foreach my $germs (@germplasm) {
		foreach my $variety (keys %var_dist) {
			if ($germs =~ /$variety/) {
			print TRI "$variety\t";
			my @dist = @{$var_dist{$variety}};
			print TRI join "\t", @dist;
			print TRI "\n";
			}
		}
	}
}

exit;
