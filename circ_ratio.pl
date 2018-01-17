################################################
#File Name: circ_ratio.pl
#Author: Wanqiu Ding 
################################################

#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use vars qw($circ $junc $out);

my ($circ,$junc,$out);
GetOptions(
				"circ|c:s"		=>		\$circ,
				"junc|j:s"		=>		\$junc,
				"out|o:s"		=>		\$out,
);

open CIRC,"$circ" or die "can't open the file!";
open JUNC,"$junc" or die "can't open the file!";
open OUT,">$out" or die "can't open the new file!";

my (%juncleft,%juncright);
while(<JUNC>){
	chomp;
	my @field = split "\t",$_;
	$juncleft{"$field[0]:$field[1]"}->{"$field[0]:$field[2]"}=$field[4];
	$juncright{"$field[0]:$field[2]"}->{"$field[0]:$field[1]"}=$field[4];
	}
=cut
foreach my $keys1 (keys %junc){
	my $hash2 = $junc{$key1};
	foreach my $key2 (sort{$hash2->{$b}<=>$hash2->{$a}}keys %$hash2){
		print $key1."\t".$key2."\t".$hash2->{$key2}."\n";
		}
		}
=cut
while(<CIRC>){
	chomp;
	my @line = split "\t",$_;
	my ($right,$left)=(0,0);
	if(exists $juncleft{"$line[0]:$line[2]"}){
		my $hash2 = $juncleft{"$line[0]:$line[2]"};
		foreach my $key2 (sort{$hash2->{$b}<=>$hash2->{$a}}keys %$hash2)
		{
		$right+=$hash2->{$key2};
		}}else{
			$right=0;}
#		print OUT join "\t",(@line[0..2],$right,$line[5],$line[6]);
#		print OUT "\n";
	if(exists $juncright{"$line[0]:$line[1]"}){
		my $hash2 = $juncright{"$line[0]:$line[1]"};
		foreach my $key2 (sort{$hash2->{$b}<=>$hash2->{$a}}keys %$hash2){
			$left+=$hash2->{$key2};
			}}else{
				$left=0;}
	print OUT join "\t",(@line,$left,$right);
	print OUT "\n";
	}
