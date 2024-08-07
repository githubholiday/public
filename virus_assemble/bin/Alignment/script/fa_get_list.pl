#!/usr/bin/perl
use strict;
use warnings;


my $fa_file=shift;
my $list=shift;
die "perl $0 fa_file list \n" unless $fa_file;

my %hash1=();
open IN,$list or die $!;
while (<IN>){
	chomp;
        my @cut=split;
        $hash1{$cut[0]} = 1;
}
close(IN);

open IN,$fa_file or die $!;
$/=">";$/=<IN>;$/="\n";
while (<IN>){
	chomp;
#	my ($id,$cvg)=split;
	(my $id=$_)=~s/\s+.*$//;
	
	$/=">";
	my $seq=<IN>;
	chomp $seq;
#	$seq=~s/\s+//g;
	$/="\n"; #取第二个id  

	if(exists $hash1{$id})
	{
#		print ">$id\t$cvg\n$seq";
		print ">$id\n$seq";
	}
	
#	my $length = length($seq);
#	$seq =~ s/N//g;
#	my $seq_len = length($seq);
#	print "$id\t$length\t$seq_len\n";
#	print "$id\t1\t$length\n";
}
close IN;
