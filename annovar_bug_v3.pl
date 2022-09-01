#!/usr/bin/perl -w
use strict;
die "perl $0 <annovar output vcf> <annovar output txt> <repair bug output>" unless @ARGV==3;
my ($vcfin,$txtin,$out)=@ARGV;
my %chr;
open TXT, $txtin or die $!;
while(<TXT>){
	chomp;
	my @c=split/\t/;
	next if ($c[0]=~/Chr/);
	my $id="$c[1]\t$c[3]\t$c[4]";#pos ref alt
	$chr{$id}=$c[0];
}
close TXT;

open OUT, ">$out";
open VCF, $vcfin or die $!;
my $column;
while(<VCF>){
	chomp;
	#my $column;
	if($_=~/^##/){print OUT "$_\n";}
	elsif($_=~/^#CHROM/){
		print OUT "$_\n";
		my @column=split/\t/;
		$column=scalar@column;
	}
	else{
		my @c=split/\t/;#43609969
		if($c[0]=~/^chr/){
			print OUT "$_\n";
		}
		elsif($c[0]=~/^[0-9]/){
			my $info=$c[6];
			my $format;
			my $id="$c[0]\t$c[2]\t$c[3]";
			if(exists $chr{$id}){
				if($c[7]=~/(GT:.+)(;ANNOVAR_DATE=.+;ALLELE_END$)/){
					$format=$1;
					$info.=$2;
					if($column==10){print OUT "$chr{$id}\t$c[0]\t$c[1]\t$c[2]\t$c[3]\t$c[4]\t$c[5]\t$info\t$format\t$c[8]\n";}
					elsif($column==11){print OUT "$chr{$id}\t$c[0]\t$c[1]\t$c[2]\t$c[3]\t$c[4]\t$c[5]\t$info\t$format\t$c[8]\t$c[9]\n";}
				}
				else{
					print "FORMAT ERROR in column 8!\n";
				}
			}
			else{
				print "ERROR!!!VCF id $id does not exist in TXT!!!\n$vcfin\n$txtin\n";
				exit;	
			}
		}
		else{
			print "ERROR: Unknown line $_\n";
			exit;
		}
	}
}
close VCF;
close OUT;
				
			
		

