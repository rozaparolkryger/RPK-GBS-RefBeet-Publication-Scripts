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
our($opt_d);
getopts('d:');

if (!$opt_d) {
    print "\nERROR: No directory given!\n";
    print "\nUsage: perl Create_ABH_on_freq.pl -d ???\n\n";
    exit 1;
};

# Open directory
opendir(DIR,$opt_d) || die "can't open directory: $opt_d\n";
my $count_file = 0;

my $count_falseSNP = 0;
my $count_B_SNP = 0;
my $count_b_SNP = 0;
my $count_H_SNP = 0;
my $count_a_SNP = 0;
my $count_A_SNP = 0;
my $count_A_allele = 0;


my %save;
my %files;

my %test;

my @files_sorted = sort{$a cmp $b} readdir(DIR);

# Analyze each file that is stored in the directory of intrest
#while (my $file = (sort {$a cmp $b }(readdir(DIR)))) {
foreach my $file (@files_sorted){ 
	# Only use txt files
    if ($file=~m/.*.txt/){

        print STDERR "\nAnalyzing $file \n";
        $files{$count_file} = $file;

        # Open file and fetch data
        open(FILE,"< $opt_d/$file") || die "can't open file: $opt_d/$file\n";
        my $count_snp = 0;
    
        while(<FILE>){
            my $line = $_;
           
	    my ($name,$ref_pos,$Variant_type,$reference,$allele_variants,$z,$x,$y,$zygosity,$a,$b,$frequency) = split("\t",$line);

#"Chromosome","Region","Type","Reference","Allele","Reference allele","Length", "Linkage","Zygosity","Count","Coverage","Frequency","Probability","Forward read count","Reverse read count","Forward/reverse balance","Average quality"
#"Bvchr1.sca001","220322","SNV","G","C","No","","Heterozygous","31","66","46.97","1","31","0","0","38.742"

            $name =~s/"//g;
            #$name =~s/ mapping//g
            $ref_pos =~s/"//g;
			#$Variant_type =~s/"//g;
	    	$reference =~s/"//g;
            $allele_variants =~s/"//g;
	    	$zygosity =~s/"//g;
			$frequency=~s/"//g;

#print STDERR "$name\t$ref_pos\t$reference\t$allele_variants\t$zygosity\n";

            if ($Variant_type=~m/SNV/){
                $count_snp++;

				# is a position on a certain scaffold defined in hash %save?
                if($save{$name."\t".$ref_pos}){
					# store array with already existing entries
                    my @variants = @{$save{$name."\t".$ref_pos}};

		    		# extend the array by integrating further information of another file
                    $variants[$count_file] = $frequency;
                    $save{$name."\t".$ref_pos} = \@variants;
                }
				# create a new hash entry on a certain position, if nothing is known before
                else{
                    my @variants;
                    $variants[$count_file] = $frequency;
                    $save{$name."\t".$ref_pos} = \@variants;
                }
            }
        }
        $count_file++;
        print STDERR "Number of SNPs $count_snp\n";
        close(FILE);
    }
}

#print Dumper %save;

print "Scaffold\tPosition\t";
foreach my $f (sort{$a <=> $b} keys(%files)){
	# find file name prefix
	my $old_filename = $files{$f};
	# name in brackets will be stored
	#my ($short_filename) = ($old_filename =~/(.+)_filtered/); 
	my ($short_filename) = ($old_filename =~/(.+)_SNPs/); 
	print $short_filename."\t";
}
print "\n";

# # print with nucleotides
# foreach my $position (sort {$a cmp $b} keys(%save)){
# 
#     my @values = @{$save{$position}};
#     print $position."\t";
#     foreach (my $i=0;$i<$count_file;$i++){
#         if ($values[$i]){
#             print $values[$i]."\t";
#         }
#         else{
#             print "-\t";
#         }
#     }
#     print "\n";
# }

#print as ABH
foreach my $position (sort {$a cmp $b} keys(%save)){
    my @values = @{$save{$position}};
    print $position."\t";
    foreach (my $i=0;$i<$count_file;$i++){

		# is the entry defined for a certain file on a certian scaffold and position?
        if ($values[$i]){
            my $allele = $values[$i];
            if ($allele > 79){
                $count_B_SNP++;
				print "B\t";
            }
			#elsif ($allele >=67 && $allele <80){
			#	$count_b_SNP++;
                #print "b\t";
            #}
	    	elsif ($allele >33 && $allele <67){
				$count_H_SNP++;
                print "H\t";
            }
			#elsif ($allele >20 && $allele <=33){
				#$count_a_SNP++;
                #print "a\t";
            #}
			#elsif ($allele <21){
				#$count_A_SNP++;
                #print "A\t";
            #}
			else{
				$count_falseSNP++;
				print "?\t";
				#print STDERR "Wrong allele type found: $allele\n";
			}
    	}
		# if not, print "-"
		else {	
				$count_A_allele++;
            	print "-\t";
        	}
    }
    print "\n";
}
print STDERR "\nFinished\n";
print STDERR "\nNumber of B SNPs: ".$count_B_SNP;
#print STDERR "\nNumber of b SNPs: ".$count_b_SNP;
print STDERR "\nNumber of H SNPs: ".$count_H_SNP;
#print STDERR "\nNumber of a SNPs: ".$count_a_SNP;
#print STDERR "\nNumber of A SNPs: ".$count_A_SNP;
print STDERR "\nNumber of unassign SNPs:".$count_falseSNP."\n";
print STDERR "\nNumber of undiscovered alleles: ".$count_A_allele."\n";
print STDERR "\n";

close(DIR);
