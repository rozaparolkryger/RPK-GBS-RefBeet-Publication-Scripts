##########################################
### Roza Parol-Kryger                  ###
### parolkrr@cebitec.uni-bielefeld.de  ###
### v0.1                               ###
##########################################


#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;

# global variables
our($opt_h,$opt_f);
getopts('h:f:');

if (!$opt_h) {
    print "\nERROR: No ABH file given!\n";
    print "\nUsage: perl ABH_covered_position.pl -h <ABH file> -f <pileup file> ???\n\n";
    exit 1;
};
if (!$opt_f) {
    print "\nERROR: No pileup file given!\n";
    print "\nUsage: perl ABH_covered_position.pl -h -f ???\n\n";
    exit 1;
};

# Open file
open(FILE1, $opt_h) || die "can't open ABH file: $opt_h\n";
open(FILE2, $opt_f) || die "can't open pileup file: $opt_f\n";

my %abh_values;
my $header;

while(<FILE1>){
	my $line = $_;
	my (@ABH) = split("\t",$line);

	if($line=~/Scaffold\s.+Position/){
		$header = $line;
	}
	else{
		my $ABH_scaffold = $ABH[0];
		my $ABH_position = $ABH[1];
	
		my @alleles;
		foreach(my $i=2;$i<$#ABH;$i=$i+1){
			push(@alleles,$ABH[$i]);
		}
		$abh_values{$ABH_scaffold."\t".$ABH_position} = \@alleles;
	}
		
}
close(FILE1);

#"Gen"	"Pos"	1	2	3	4	5	6	7	8	9
# "Bv_00110_aycw"	3710	"A"	"A"	"A"	"A"	"A"	"A"	"A"	"B"	"A"
print STDERR "finished file 1\n";

#print Dumper %abh_values;
my $count_update = 0;
while(<FILE2>){
	my $line = $_;
	my (@pileup) = split("\t",$line);
	my $pileup_scaffold = $pileup[0];
	my $pileup_position = $pileup[1];
	my $pileup_base = $pileup[2];

 	my @genotypes;
	foreach(my $i=3;$i<=$#pileup;$i=$i+3){
		push(@genotypes,$pileup[$i]);
	}

	if($abh_values{$pileup_scaffold."\t".$pileup_position}){
		my @alleles = @{$abh_values{$pileup_scaffold."\t".$pileup_position}};
		for (my $j=0;$j<=$#alleles;$j++){
			if($alleles[$j] eq "-"){
				if($genotypes[$j] > 0){
					$alleles[$j] = "A";
					$count_update++;
				}
			}
		}
		# Update abh entries
		$abh_values{$pileup_scaffold."\t".$pileup_position} = \@alleles;
	}
}
close(FILE2);

print STDERR "finished file 2\n";



#print HEADER for ABH-Table, same as in preABH-table
#print "Scaffold\tPosition\t1K_2945\t1K_2948\t6FB2002\t6FB2003\t6FB2004\t6FB2005\t6FB2006\t6FB2008\t6FB2009\t9FB2134\t9FB2135\t9FB2138";
print "$header\n";

# print updated abh table
foreach my $position (sort{$a cmp $b} keys (%abh_values)){
	# get updated entries
	my @abh = @{$abh_values{$position}};

	# print values for Bvchr1
	#if($position=~m/Bvchr1/){
	#	print "$position\t".join("\t",@abh)."\n";
	#}
	# print all values
	print "$position\t".join("\t",@abh)."\n";
}

#print Dumper @genotype;

# # Bvchr1.sca001   13236   G       4       ....    IIGI    6       ......  IEEIHI  6       ......  IIIIII  1       .       H       2       ..      GG      0       *       *       0       *       *       0       *       *       2       ..      ?I      5       .....   ?IBI:   2
# #        ..      FI      2       ..      II
# 
# if ($ABH_scaffold == $pileup_scaffold){
# 	if ($ABH_position == $pileup_position){
# 		if ($ABH_g1 =="A"){
# 			if ($pileup_g1>0)
# 			print"g1 covered by" $pileup_g1 "reads";
# 		}
# 		else
# 	}
# }
# 
# 
# close(FILE1,FILE2);
print STDERR "\nFinished, changed entries $count_update\n\n";
