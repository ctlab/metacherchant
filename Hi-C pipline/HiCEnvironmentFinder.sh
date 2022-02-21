#!/bin/bash

metacherchant=""
wgs_reads=""
seq_path=""
work_dir=""
hi_c_r1=""
hi_c_r2=""
k=31
coverage=5
maxradius=100000

echo "Program parameters:"
while [ -n "$1" ]
do
	case "$1" in
	--metacherchant) metacherchant="$2"
	echo "--metacherchant = $metacherchant"
	shift ;;
	--reads ) wgs_reads="$2"
	echo "--reads = $wgs_reads"
	shift ;;
	--seq) seq_path="$2"
	echo "--seq = $seq_path"
	shift ;;
	--work-dir) work_dir="$2"
	echo "--work-dir = $work_dir"
	shift;;
	--k ) k="$2"
	echo "--k = $k"
	shift ;;
	--coverage) coverage="$2"
	echo "--coverage = $coverage"
	shift ;;
	--maxradius) maxradius="$2"
	echo "--maxradius = $maxradius"
	shift;;
	--hi-c-r1) hi_c_r1="$2"
	echo "--hi-c-r1= $hi_c_r1"
	shift ;;
	--hi-c-r2) hi_c_r2="$2"
	echo "--hi-c-r2= $hi_c_r2"
	shift;;
	*) echo "$1 is not an option";;
	esac
	shift
done

if [[ ${wgs_reads} = "" || ${seq_path} = "" || ${work_dir} = "" || ${hi_c_r1} = "" || ${hi_c_r2} = "" || ${metacherchant} = "" ]]
then
echo "Please add all mandatory parameters: --reads --seq --work-dir --metacherchant --hi-c-r1 --hi-c-r2"
exit -1
fi

work_dir="${work_dir}/"

java -jar $metacherchant --k $k --coverage $coverage --reads $wgs_reads --seq $seq_path --output "${work_dir}output/1" --work-dir "${work_dir}workDir/1" --maxradius $maxradius --bothdirs False --chunklength 10

mkdir "${work_dir}1"
mkdir "${work_dir}2"

bwa index "${work_dir}output/1/seqs.fasta"
bwa mem -t 12 "${work_dir}output/1/seqs.fasta" "${hi_c_r1}" "${hi_c_r2}" > "${work_dir}1/all_hic_reads.sam"
samtools view -f 0x5 -F 0x908 -o "${work_dir}1/filteredHiC_1.bam" "${work_dir}1/all_hic_reads.sam"
samtools view "${work_dir}1/filteredHiC_1.bam" | awk ' {print ">"1"\n"$10} ' - > "${work_dir}1/selected_reads.fasta"

java -jar $metacherchant --k $k --coverage $coverage --reads $wgs_reads --seq $seq_path --output "${work_dir}output/2" --work-dir "${work_dir}workDir/2" --maxradius $maxradius --bothdirs False --chunklength 10 --merge true --hicseq "${work_dir}1/selected_reads.fasta"


bwa index "${work_dir}output/2/merged/seqs.fasta"
bwa mem -t 12 "${work_dir}output/2/merged/seqs.fasta" "${hi_c_r1}" "${hi_c_r2}" > "${work_dir}2/filteredHiC_2.sam"

samtools view -@ 12 -f 1 -F 2060 -b -o "${work_dir}2/filteredHiC_pair_map.bam" "${work_dir}2/filteredHiC_2.sam"
samtools view "${work_dir}2/filteredHiC_pair_map.bam" | awk '($3!=$7 && $7!="=")' - > "${work_dir}2/filtered_HiC_diff_chr.sam"

python hic_map.py "${work_dir}2/" "filtered_HiC_diff_chr.sam"
