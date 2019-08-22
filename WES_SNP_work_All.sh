#!/bin/bash
dir=$1
sample=$2
control=$3
cat $sample|while read sp 
do
echo "
source /datg/wangrui/Immune_Variant/Exon/WES_SNP_pipeline.sh $dir $sp $control

do_QC
do_bwa_mapping
do_MarkDuplicates
do_RealignerTargetCreator
do_check_GATK RealignerTargetCreator
do_IndelRealigner
do_check_GATK IndelRealigner
do_BaseRecalibrator
do_check_GATK BaseRecalibrator
do_PrintReads
do_check_GATK PrintReads
do_HaplotypeCaller
do_HaplotypeCaller_chrM
do_check_GATK HaplotypeCaller
do_Haplo_VariantFilter
do_snpEff_Anno_Hap
do_SIFT_Polyphen_Anno_Hap
#do_ExtractVariant_Hap

if [[ $sp -ne $control ]];
then
	do_Mutect2
	do_Mutect2_chrM
	do_snpEff_Anno_Mutect
	do_SIFT_Polyphen_Anno_Mutect
#	do_ExtractVarint_Mutect
elif
	exit 0
fi

" > $sp.work.tmp.sh
qsub -cwd -l vf=30g,p=7 -V $sp.work.tmp.sh
done
