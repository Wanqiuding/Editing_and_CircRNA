#!/usr/bin/perl 

use strict;
use Getopt::Long;
use File::Basename;
use lib dirname $0;
use Bio::DB::Fasta;

#circular RNA annalysis toolkit.

my ($fusion, $genome, $anno, $mapped, $fasta, $output, $pair, $version, $bin);
sub usage{
	my $scriptName = basename $0;
	print <<HELP;
Usage: to identify circRNA using NGS data with the tophat-fusion result. The sequencing data can be single-end or paired-end(with the option -p <mapped reads bam file>). 
For single-end reads, we detect backsplice reads indicating circRNAs. For paired-end reads, we choose the reads: two paired reads are both backsplice reads and indicate the same circRNA
or one of a pair is backspliced and the other is on the same transcript.
	perl $scriptName OPTION
Options:
			-h|--help	=>	sub{usage()},
			-v|--version	=>	\$version,
			-f|--fusion	=>	\$fusion,
			-g|--genome	=>	\$fasta,
			-r|--ref	=>	\$anno,
			-b|--bin	=>\$bin,
			-p|--paired 	=>	\$pair,
			-o|--output	=>	\$output,
HELP
	exit(-1);
}

GetOptions(
			'h|help'	=>	sub{usage()},
			'v|version'	=>	\$version,
			'f|fusion=s'	=>	\$fusion,
			'g|genome=s'	=>	\$fasta,
			'r|ref=s'	=>	\$anno,
			'b|bin'	=>\$bin,
			'p|paired=s' 	=>	\$pair,
			'm|mapped'	=>	\$mapped,
			'o|output=s'	=>	\$output,
);

open ANN,"$anno" or die "can't open the gpe file!";

my %annotation;
while (<ANN>){
	chomp;
	my @fields = split "\t",$_;
	shift @fields if defined $bin;
	$annotation{$fields[0]}={
											start => $fields[3],
											end	=>	$fields[4],
											chr => $fields[1],
											strand => $fields[2],
											cdsSt	=>	$fields[5],
											cdsEnd	=>	$fields[6],
											exonStarts	=>	$fields[8],
											exonEnds	=>	$fields[9],
											gene	=>	$fields[11],
											};
}

if (-B $fusion){
	open FUSION, "samtools view -h $fusion|" or die "can't open $fusion: $!";
	}else{
		open FUSION, "$fusion" or die "can't open $fusion: $!";
		}

my (%fusion,%fusionName);		
while(<FUSION>){
	chomp;
	my $line = $_;
	my @field = split "\t",$line;
	next if($field[1] == 256);   #filter secondary alignment
	next if($line !~ /.*XF.*/);  #filter no-fusion reads
	my ($chr1,$chr2);
	$line=~/chr(\d+)-chr(\d+)/;
	$chr1 = $1;
	$chr2 = $2;
	next if(defined($chr1) && ($chr1 ne $chr2));
	my $strand = $field[1] & 0b10000 ? '-' : '+';
	if($line =~ /.*XF:Z:1.*/){
		if($line =~ /(\d+)m(\d+)F(\d+)m/){
			if($1>=8 && $3>=8){
				my $frag1_s = $field[3]-1;
				my $frag2_s =  $frag1_s + $1;
				my $fus_start = $2 - $3;
				my $fus_site = $2;
				$fusion{$field[2]}->{$line}={
					start1 => $frag1_s,
					end1 => $frag2_s,
					start2 => $fus_start,
					end2 => $fus_site,
					strand =>  $strand,
					chr => $field[2],
					name => $field[0],
				};
				$fusionName{$field[0]}=0;
			}
		}
	}
}		
my %mapPair;
if(defined $pair){
	open MAP, "samtools view -h $pair|" or die "can't open $pair: $!";
my $number = 0;
while(<MAP>){
	next if /^@/; 
	chomp;
	my @fields = split "\t", $_;
	my ($readPair, $flag, $chr, $start, $mapQ, $cigar) = @fields[0..5];
	my $strand = $flag & 0b10000 ? '-' : '+';
	$start--;
	$number++;
if (exists $fusionName{$readPair}){
	my (@blockStarts, @blockEnds);
	my ($blockStarts, $blockEnds) = cigar2block($start, $cigar);
	foreach my $factor(@$blockStarts){
		push @blockStarts,$factor;
		}
	foreach my $factor(@$blockEnds){
		push @blockEnds,$factor;
		}
	$mapPair{$readPair}->{$number}={
				chrPair	=>	$chr,
				startPair	=>	$blockStarts[0],
				endPair	=>	$blockEnds[-1],
				};
}
}
}
open OUT, ">$output" or die "can't open the new file!";
my $db = Bio::DB::Fasta->new($fasta);
for my $transName (keys %annotation){
	my $transV = $annotation{$transName};
	my @exonS = split ',',$transV->{exonStarts};
	my @exonE = split ',',$transV->{exonEnds};
	my $readH = $fusion{$transV->{chr}};
	my (@circS,@circE,@leftIndex,@rightIndex,@index,@readName,@circStrand,@readSeq);
	for my $read ( keys %$readH ){
		my $mapSeq;
		if($read =~/[Mm]\s+(\w+)/){
			$mapSeq = $1;
			}
        my $readV=$readH->{$read};
		if ((grep {($readV->{start1}) == $_} @exonS) && (grep {($readV->{end2}) == $_} @exonE)){    #exactly mapped to the annotation
			push @circS,$readV->{start1};
			push @circE,$readV->{end2};
			push @readName,$readV->{name};
			push @circStrand,$readV->{strand};
			push @readSeq,$mapSeq;
			@index = grep {($exonS[$_]>=($readV->{start1})) && ($exonS[$_]<($readV->{end2}))} (0..@exonS);
			push @leftIndex,$index[0];
			push @rightIndex,$index[-1];
		}else{
			my @errorIndex;
			my ($errorS,$errorE)=(0,0);
			foreach my $i (0..$#exonS){  #set 10 bp fixation realignment
				if((($exonS[$i]-10)<=($readV->{start1})) && (($exonS[$i]+10)>=($readV->{start1}))){
					$errorS = $exonS[$i];
					push @errorIndex,$i;
				}
				if((($exonE[$i]-10)<=($readV->{end2})) && (($exonE[$i]+10)>=($readV->{end2}))){
					$errorE = $exonE[$i];
					push @errorIndex,$i;
				}
			}
			if($errorS !=0 && $errorE !=0 && $errorS !=($readV->{start1}) && $errorE !=($readV->{end2})){
				my $anchor1 = $db->seq($readV->{chr},$readV->{start1}+1,$errorS);
				my $anchor2 = $db->seq($readV->{chr},$readV->{end2}+1,$errorE);
				if(&diffStr($anchor1,$anchor2)){
					push @circS,$errorS;
					push @circE,$errorE;
					push @readName,$readV->{name};
					push @circStrand,$readV->{strand};
					push @readSeq,$mapSeq;
					push @leftIndex,$errorIndex[0];
					push @rightIndex,$errorIndex[1];#输出intron信息
				 }
			}	
		}
	}
my @location;
foreach my $j (0..$#circS){
	push @location,join('-',$circS[$j],$circE[$j]);
	}
my @out;
	foreach my $i (0..$#location){
		push @out,join('|',$location[$i],$leftIndex[$i],$rightIndex[$i]);
		}
		
if(defined $pair){
	my (@pairRead);
	my %count;
	@pairRead = grep{++$count{$_}>=2} @readName;   
	my (@spliceIndex,@dupIndex);
	foreach my $dup (@pairRead){
		my (@dupName,@dupSeq,@dupLindex,@dupLocus);
		my @indexR;
		foreach my $n (0..@readName){
			if($readName[$n] eq $dup){
				push @dupIndex,$n;
				push @dupSeq,$readSeq[$n];
				push @indexR, $n;
			}
		}
		my @uniqDup = &uniq(@dupSeq);
		if($#uniqDup==0){   #one read is mapped to different position(multiple primary alignments)
			foreach my $index (@indexR){
				my $tag = 0;
				if(exists $mapPair{$readName[$index]}){
					my $mapH = $mapPair{$readName[$index]};
					my @locus = split '-',$location[$index];
					foreach my $mapRead (keys %$mapH){
					my $mapV = $mapH->{$mapRead};
						if(($mapV->{chrPair} ==$transV->{chr}) && ($mapV->{startPair} >=$locus[0]) && ($mapV->{endPair} <=$locus[1])){
						$tag = 1;
						push @dupLindex,$index;
						push @dupLocus,$location[$index];
						}
					}
				}
				if(!$tag){
					push @spliceIndex,$index;
				}
			}
			my @uniqdupLocus = uniq(@dupLocus); 
			if($#uniqdupLocus ==0){   
				pop @dupLindex;
			}
			foreach my $i (@dupLindex){
				push @spliceIndex,$i;
			}
		}else{   #read pairs are both fusion reads
			my (@one,@two,@sameOne,@sameTwo,@diff,@pairLoca);
			foreach my $i (0..$#dupSeq){
				if($dupSeq[$i] eq $uniqDup[0]){
					push @one,$indexR[$i];
					}else{
						push @two,$indexR[$i];
					}
			}
			my $samePair=0;
			foreach my $j (@one){
				foreach my $k (@two){
					if($location[$k]==$location[$j]){
						$samePair++;
						push @pairLoca,$location[$k];
						push @sameOne,$j;
						push @sameTwo,$k;
						}else{
							push @diff,$k;
							push @diff,$j;
							}
				}
			}
			my @uniqPair = uniq(@pairLoca);
			if($#uniqPair==0){
				pop @sameOne;
				foreach my $m(@sameOne){
					push @spliceIndex,$m;
					}
				foreach my $m(@sameTwo){
					push @spliceIndex,$m;
					}
				foreach my $m(@diff){
					push @spliceIndex,$m;
					}
			}
			if(($samePair == 0) || ($#uniqPair>=1)){
				foreach my $m (@indexR){
					push @spliceIndex,$m;
				}
			}
		}
	}
	foreach my $m (0..@readName){   #a read in a pair is fusion reads and the other is not
		if((!grep{$m==$_}@dupIndex)){
		my $tag =0;
		if(exists $mapPair{$readName[$m]}){
			my $mapH = $mapPair{$readName[$m]};
			my @locus = split '-',$location[$m];
			foreach my $mapRead (keys %$mapH){
					my $mapV = $mapH->{$mapRead};
			if(($mapV->{chrPair} ==$transV->{chr}) && ($mapV->{startPair} >=$locus[0]) && ($mapV->{endPair} <=$locus[1])){
				$tag = 1;
			}
			}
		}
		if(!$tag){
			push @spliceIndex,$m;
		}
		}
	}
	foreach my $index (@spliceIndex){
		undef($out[$index]);
		}
	@out = grep { defined($_) } @out;
	my %pos;
	my $num = 0;
	map {$pos{$_}++} @out;
	foreach my $key(keys %pos){
		my @code = split '\|',$key;
		my @site = split '-',$code[0];
		if(($code[1] != 0) && ($code[2] != $#exonS)){
			print OUT join "\t",($transV->{chr},$site[0],$site[1],$transName,$pos{$key},$transV->{strand},$code[1],$code[2],$exonE[$code[1]-1],$exonS[$code[2]+1],$transV->{gene});
			print OUT "\n";
		}elsif(($code[1] == 0) && ($code[2] != $#exonS)){
			print OUT join "\t",($transV->{chr},$site[0],$site[1],$transName,$pos{$key},$transV->{strand},$code[1],$code[2],"none",$exonS[$code[2]+1],$transV->{gene});
			print OUT "\n";
			}elsif(($code[1] != 0) && ($code[2] == $#exonS)){
				print OUT join "\t",($transV->{chr},$site[0],$site[1],$transName,$pos{$key},$transV->{strand},$code[1],$code[2],$exonE[$code[1]-1],"none",$transV->{gene});
				print OUT "\n";
				}else{
					print OUT join "\t",($transV->{chr},$site[0],$site[1],$transName,$pos{$key},$transV->{strand},$code[1],$code[2],"none","none",$transV->{gene});
					print OUT "\n";}
	}
 }else{
	my %pos;
	my %count;
	my $num = 0;
	my @multiMap = grep{++$count{$_}>=2} @readName;
	foreach my $h (0..$#multiMap){
	my (@multiPos,@multiIndex);
		foreach my $r (0..$#readName){
			if($readName[$r] eq $multiMap[$h]){
				push @multiPos,$location[$r];
				push @multiIndex,$r;
			}
		}
		my @uniq = uniq(@multiPos);
		if($#uniq == 0){
			pop @multiIndex;
			}
		foreach my $index (@multiIndex){
			undef($out[$index]);
			}
		@out = grep {defined($_)} @out;
		}
	 map {$pos{$_}++} @out;
	foreach my $key(keys %pos){
		my @code = split '\|',$key;
		my @site = split '-',$code[0];
		if(($code[1] != 0) && ($code[2] != $#exonS)){
			print OUT join "\t",($transV->{chr},$site[0],$site[1],$transName,$pos{$key},$transV->{strand},$code[1],$code[2],$exonE[$code[1]-1],$exonS[$code[2]+1],$transV->{gene});
			print OUT "\n";
		}elsif(($code[1] == 0) && ($code[2] != $#exonS)){
			print OUT join "\t",($transV->{chr},$site[0],$site[1],$transName,$pos{$key},$transV->{strand},$code[1],$code[2],"none",$exonS[$code[2]+1],$transV->{gene});
			print OUT "\n";
			}elsif(($code[1] != 0) && ($code[2] == $#exonS)){
				print OUT join "\t",($transV->{chr},$site[0],$site[1],$transName,$pos{$key},$transV->{strand},$code[1],$code[2],$exonE[$code[1]-1],"none",$transV->{gene});
				print OUT "\n";
				}else{
					print OUT join "\t",($transV->{chr},$site[0],$site[1],$transName,$pos{$key},$transV->{strand},$code[1],$code[2],"none","none",$transV->{gene});
					print OUT "\n";}
	}
	}
}
	
sub diffStr{
	my ($str1,$str2) = @_;
if(length($str1)==length($str2)){
	my @str1 = split '',$str1;
	my @str2 = split '',$str2;
	my $diffNum = 0;
	foreach my $j (0..$#str1){
		$diffNum++ if $str1[$j] ne $str2[$j];
		}
	if($diffNum ==0){    #set new fault-tolerant parameter
	return 1;
	}else{
		return 0;
		}
}else{
	return 0;
	}
}
sub get_seq{
         my ($ref,$chr,$start,$end,$strand) = @_;
		 my $db = Bio::DB::Fasta->new($ref);
         my $seq;
		 $seq = $db->seq("$chr:$start-$end");
		 return $seq;
     }
sub uniq {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}
sub cigar2block {
	my ($blockStart, $cigar);
	($blockStart, $cigar) = @_;
	my (@blockStarts, @blockEnds);
	while($cigar =~ s/([^N]+?)(\d+)N//){
		my ($match1, $match2) = ($1, $2);
        my $blockEnd = $blockStart;
        push @blockStarts, $blockStart;
        $blockEnd += $_ for($match1 =~ /(\d+)[MD=X]/g);
        $blockStart = $blockEnd + $match2;
        push @blockEnds, $blockEnd;
		}
	my $blockEnd = $blockStart;
    push @blockStarts, $blockStart;
    $blockEnd += $_ for($cigar =~ /(\d+)[MD=X]/g);
    push @blockEnds, $blockEnd;
    return (\@blockStarts, \@blockEnds);
}