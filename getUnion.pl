#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my $file1 = '';
my $file2 = '';
my $key = '';
my $outfile = '';
my ($key1,$key2,$point,@key1,@key2,@point);

GetOptions ('file1|f1=s'=>\$file1, 'key1|k1=s'=>\$key1, 'key2|k2=s'=>\$key2,
			'point|p=s'=>\$point, 'file2|f2=s'=>\$file2, 'key|k=s'=>\$key, 
			'outfile|o=s'=>\$outfile, 'h|help'=>sub{usage()})||usage();

open (FILE1,"$file1") or die "can't open the file!";
open (FILE2,"$file2") or die "can't open the file!";
open (KEY,"$key") or die "can't open the key file!";
open (OUT,">$outfile") or die "can't open the new file!";

@key1 = split ',',$key1;
@key2 = split ',',$key2;
@point = split ',',$point;

my (%hash1);
while(<FILE1>){
	chomp;
	my @field = split "\t",$_;
	my $file1_key = join (':',@field[@key1]);
	$hash1{$file1_key} = join "\t",@field[@point];
	}

my (%hash2);
while(<FILE2>){
	chomp;
	my @field = split "\t",$_;
	my $file2_key = join (':',@field[@key1]);
	$hash2{$file2_key} = join "\t",@field[@point];
	}	

while(<KEY>){
	chomp;
	my @field = split "\t",$_;
	my $line_key = join (':',@field[@key2]);
	my ($value1,$value2);
	if(exists $hash1{$line_key}){
		$value1 = $hash1{$line_key};
			}else{
					$value1 = join "\t", map {0*$_} @point;
				}
	if(exists $hash2{$line_key}){
		$value2 = $hash2{$line_key};
			}else{
				$value2 = join "\t", map {0*$_} @point;
				}
	print OUT join "\t",(@field,$value1,$value2);
	print OUT "\n";
}
sub usage{
	my $scriptName = basename $0;
print <<HELP;
Usage: to get the union value of two sets based on a common key to file1
	perl $scriptName -f1 file1 -f2 file2 -k1 1,2,4 -k2 2,4,5 -k file_with_key -p keypoint -o outfile}
perl $scriptName --file1=/mnt/share/dingwq/ADAR_KD/293-total-NC_fusion/293-NC-long_range_circ --file2=/mnt/share/dingwq/ADAR_KD/293-total-KD_fusion/293-KD-long_range_circ --key=uniq_addjunc_293_NC_KD_longsplice_fil_psi_NCedit_only --key1=0,1,2,3 --point=4 --key2=10,11,12,20,13 --outfile=pltest
Options:
	'file1|fi=s'=>\$file1			FILE	file1 which is base on k1 to subtract interested p(keypoint)
	'key1|k1=s'=>\$key1						comm-separated key of file1
	'key2|k2=s'=>\$key2						comm-separated key of file2
	'point|p=i'=>\$point					field which is interested
	'file2|f2=s'=>\$file2			FILE	file2 which is base on k2 to subtract interested p(keypoint)
	'key|k=s'=>\$key				FILE	file with key and is the common 		
	'outfile|o=s'=>\$outfile		FILE	out file
	'h|help'=>sub{usage()}					print this help information
HELP
	exit(-1);
}
