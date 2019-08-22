#!/bin/bash

# setting 
dir=$1
indir=$dir/00.raw_data
outpath=$dir/01.clean_data
outdir=$dir/02.bam
mkdir -p $indir
mkdir -p $outdir 
mkdir -p $outpath
sp=$2
control=$3

# script
QC_script=/datg/wangrui/bin/QC_plus_rm_primer_polyA_T.pl

#software
bwa=/datc/wangrui/software/bwa-0.7.12/bwa
SAMTOOLS=/datc/wangrui/software/samtools-1.2/samtools
PICARD="java -jar /datc/wangrui/software/picard-tools-1.130/picard.jar"
GATK="java -jar /datb1/wangrui/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar"
#GATK="java -jar /datc/wangrui/software/GenomeAnalysisTK/GenomeAnalysisTK.jar" # older version
snpEff=/datg/wangrui/Software_NewVersion/snpEff_v4.3/snpEff/snpEff.jar
SnpSift=/datg/wangrui/Software_NewVersion/snpEff_v4.3/snpEff/SnpSift.jar
snpEff_config=/datg/wangrui/Software_NewVersion/snpEff_v4.3/snpEff/snpEff.config
SnpSift_db=/datb1/wangrui/Database/dbNSFP/hg19/dbNSFP.txt.gz
#Database
#H_fai=/datb1/wangrui/Database/human/hg19/hg19.fa.fai
#Human_ref=/datb1/wangrui/Database/human/hg19/hg19.fa
H_fai=/datb1/wangrui/Database/IBM_Database/GATK/genome.fa.fai
Human_ref=/datb1/wangrui/Database/IBM_Database/GATK/genome.fa # include chrM
vcf=/datc/wangrui/Database/indel_annotation/Mills_and_1000G_gold_standard.indels.hg19.vcf
Exon_target=/datg/wangrui/Immune_Variant/Exon/Human_All_Exon_V6/S07604514_Regions.bed
#knownSites_1=/datc/wangrui/Database/dbSNP/dbsnp_135.hg19.vcf
knownSites_1=/datg/wangrui/Colon_Liver_addSample/Exon/bin/dbsnp_135.hg19.vcf
knownSites_2=/datc/wangrui/Database/indel_annotation/Mills_and_1000G_gold_standard.indels.hg19.vcf
knownSites_3=/datc/wangrui/Database/dbSNP/1000G_phase1.indels.hg19.vcf

function do_QC {
mkdir -p $outpath/$sp
perl $QC_script -indir $indir -outdir $outpath -sample $sp -end 2 -scRNA 0

}



function do_bwa_mapping {

mkdir -p $outdir/$sp

$bwa mem -M -t 2 -R "@RG\tID:$sp\tLB:$sp\tPL:ILLUMINA\tPM:HISEQ\tSM:$sp" $Human_ref  $outpath/$sp/$sp.R1.clean.fq.gz $outpath/$sp/$sp.R2.clean.fq.gz| $SAMTOOLS view -u -b -S -t $H_fai -| $SAMTOOLS sort -m 200000000 - $outdir/$sp/$sp.sort

}



function do_MarkDuplicates {


$PICARD MarkDuplicates I=$outdir/$sp/$sp.sort.bam O=$outdir/$sp/$sp.Mdu.sort.bam M=$indir/$sp/marked_dup_metrics.txt CREATE_INDEX=TRUE

}



function do_RealignerTargetCreator {

ref=$Human_ref

$GATK -T RealignerTargetCreator -R $ref -I $outdir/$sp/$sp.Mdu.sort.bam -o $outdir/$sp/$sp.forindelRealigner.intervals --known $vcf -L $Exon_target
}


function do_IndelRealigner {

ref=$Human_ref

$GATK -T IndelRealigner -R $ref -I $outdir/$sp/$sp.Mdu.sort.bam -targetIntervals $outdir/$sp/$sp.forindelRealigner.intervals -known $vcf -o $outdir/$sp/$sp.realignedBam.bam

}


function do_BaseRecalibrator {

ref=$Human_ref

$GATK -T BaseRecalibrator -R $ref -I $outdir/$sp/$sp.realignedBam.bam -L $Exon_target -knownSites $knownSites_1  -knownSites $knownSites_2 -knownSites $knownSites_3  -o $outdir/$sp/$sp.recal.grp

}


function do_PrintReads {

ref=$Human_ref

$GATK -T PrintReads -R $ref -I $outdir/$sp/$sp.realignedBam.bam -BQSR $outdir/$sp/$sp.recal.grp -o $outdir/$sp/$sp.realn_Recal.bam

}


function do_HaplotypeCaller {

ref=$Human_ref

$GATK -T HaplotypeCaller -R $ref -I  $outdir/$sp/$sp.realn_Recal.bam -L $Exon_target -stand_call_conf 30.0 -o $outdir/$sp/$sp.HaplotypeCaller.variants_result.vcf

}


function do_HaplotypeCaller_chrM {

ref=$Human_ref

$GATK -T HaplotypeCaller -R $ref -I  $outdir/$sp/$sp.realn_Recal.bam -L chrM -stand_call_conf 30.0 -ERC BP_RESOLUTION -o $outdir/$sp/$sp.HaplotypeCaller.variants_result.chrM.gvcf

}


function do_Mutect2 {

ref=$Human_ref

$GATK -T MuTect2 -L $Exon_target -I:tumor $outdir/$sp/$sp.realn_Recal.bam -I:normal $outdir/$control/$control.realn_Recal.bam --dbsnp $knownSites_1 --output_mode EMIT_VARIANTS_ONLY -o $outdir/$sp/$sp.Directed_mutect2.vcf.gz -R $ref

}


function do_Mutect2_chrM {

ref=$Human_ref

$GATK -T MuTect2 -L chrM -I:tumor $outdir/$sp/$sp.realn_Recal.bam -I:normal $outdir/$control/$control.realn_Recal.bam --dbsnp $knownSites_1 --output_mode EMIT_VARIANTS_ONLY -o $outdir/$sp/$sp.Directed_mutect2.chrM.vcf.gz -R $ref

}

function do_Haplo_VariantFilter {

ref=$Human_ref

$GATK -T VariantFiltration -R $ref -V $outdir/$sp/$sp.HaplotypeCaller.variants_result.vcf  -filterName FS -filter "FS > 30.0" -filterName DP -filter "DP>30" -o  $outdir/$sp/$sp.HaplotypeCaller.filted.vcf
}


function do_snpEff_Anno_Mutect {

java -Xmx4g -jar $snpEff -c $snpEff_config  hg19  $outdir/$sp/$sp.Directed_mutect2.vcf.gz >$outdir/$sp/$sp.mutect2.ann.vcf

}


function do_snpEff_Anno_Hap {


java -Xmx4g -jar $snpEff -c $snpEff_config  hg19 $outdir/$sp/$sp.HaplotypeCaller.filted.vcf >$outdir/$sp/$sp.HaplotypeCaller.ann.vcf

}


function do_SIFT_Polyphen_Anno_Hap {

java -jar $SnpSift dbnsfp -v -db $SnpSift_db $outdir/${sp}/${sp}.HaplotypeCaller.ann.vcf >$outdir/${sp}/${sp}.HaplotypeCaller.ann.dbNSFP.vcf
}



function do_SIFT_Polyphen_Anno_Mutect {

indir=$1
sp=$2

java -jar $SnpSift dbnsfp -v -db $SnpSift_db $outdir/${sp}/${sp}.mutect2.ann.vcf > $outdir/${sp}/${sp}.mutect2.ann.dbNSFP.vcf

}


function do_ExtractVariant_Hap {


paste <(java -jar $SnpSift extractFields -s "," -e "." $outdir/${sp}/${sp}.HaplotypeCaller.ann.dbNSFP.vcf  CHROM POS REF ALT FILTER "ANN[*].EFFECT" "ANN[*].GENE" "ANN[*].IMPACT" "ANN[*].BIOTYPE" "ANN[*].HGVS" dbNSFP_SIFT_pred dbNSFP_Polyphen2_HDIV_pred dbNSFP_Polyphen2_HVAR_pred) <(cut -f 10 $outdir/${sp}/${sp}.HaplotypeCaller.ann.dbNSFP.vcf|grep -v "^#"|cut -d ":" -f 2 ) >$outdir/${sp}/${sp}.HaplotypeCaller.extract.txt

}



function do_ExtractVarint_Mutect {

paste <(java -jar $SnpSift extractFields -s "," -e "." $outdir/${sp}/${sp}.mutect2.ann.dbNSFP.vcf CHROM POS REF ALT FILTER "ANN[*].EFFECT" "ANN[*].GENE" "ANN[*].IMPACT" "ANN[*].BIOTYPE" "ANN[*].HGVS" dbNSFP_SIFT_pred dbNSFP_Polyphen2_HDIV_pred dbNSFP_Polyphen2_HVAR_pred ) <(cut -f 10 $outdir/${sp}/${sp}.mutect2.ann.dbNSFP.vcf|grep -v "^#" |cut -d ':' -f 1,2,3) <(cut -f 7 $outdir/${sp}/${sp}.mutect2.ann.dbNSFP.vcf|grep -v "^#" )> $outdir/${sp}/${sp}.mutect2.extract.txt
}


# get jobID
jobID=`grep $sp.work.tmp.sh  JobID_file.e | cut -d ' ' -f 3`
echo $jobID

function do_check_GATK {
step=$1
error_num=`grep 'ERROR' $sp.work.tmp.sh.e${jobID}|wc -l`
if [[ $error_num -gt 0 ]]
then
	echo "$step Error" >> $sp.e
	exit 0 
else 
	echo "$step OK"  >> $sp.e
fi
}
