#!/bin/sh
#######################
###Name: Ding Wanqiu###
#######################

usage(){
    cat << EOF
Description: Calculate the circular ratio [circular num/(circular num + maximum linear num)] and PSI of flanking junctions in each circRNA.
Dependencies: circ.uniq.num file through command 
"cut -f1-3,5-6,9-10 $query|awk '$6!= "none" && $7!= "none"'|sort|uniq|awk -F "\t" 'OFS="\t"{print $1,$2,$3,"circ_"NR,$4,$5,$6,$7}'" 
and uniqjunc.bed12 from sam2junc.pl script.
Options:
	c	FILE	circ.bed6+ file through circRNA identification pipeline
	j	FILE	junction file through sam2junc.pl script
	s	FILE	scan.es.bed12 file in a sample
	n	CHR	sample name such as R_brain
Usage: sh $0 -c circ.bed6+ -j uniqJunc.bed12+ -s scan.es.bed12 -n R_brain
EOF
    exit 0
}
[ $1 ] || usage

while getopts "hc:j:n:s:" OPTION
do
    case $OPTION in
        h) usage;;
	c) circ=$OPTARG;;
	j) junc=$OPTARG;;
	n) name=$OPTARG;;
	s) scanes=$OPTARG;;
	?) usage
	   exit 0
	   ;;
    esac
done

cut -f1-3,5-6,9-10 $circ|awk '$6!= "none" && $7!= "none"'|sort|uniq|awk -F "\t" 'OFS="\t"{print $1,$2,$3,"circ_"NR,$4,$5,$6,$7}' >${name}.uniq.num.bed6+
grep -v '#' $scanes|cut -f1-3,5- |sort|uniq|awk -F "\t" 'OFS="\t"{print $1,$2,$3,"es_"NR,$0}'|cut -f1-4,8- >${name}.scan.es.uniq.bed12+

perl /mnt/share/dingwq/bin/circ_ratio.pl -c ${name}.uniq.num.bed6+ -j $junc -o ${name}.uniq.num.max.bed6+
awk -F "\t" '{if($9>$10)print $0"\t"$9"\t"$5/($5+$9);else print $0"\t"$10"\t"$5/($5+$10)}' ${name}.uniq.num.max.bed6+ |cut -f1-8,11,12 >${name}.uniq.ratio.bed6+

perl /mnt/share/dingwq/bin/calculate_psi.pl --file=${name}.scan.es.uniq.bed12+ --key=${name}.uniq.ratio.bed6+ --junc=$junc --outfile=${name}.uniq.ratio.psi.bed6+

sed -i '1i#chr\tstart\tend\tcirc_id\tbacksplice_supp\tstrand\tleft_junc\tright_junc\tmax_juncNum\tcirc_ratio\tpsi\tflank_juncNum\ttotal_juncNum\tlabel' ${name}.uniq.ratio.psi.bed6+