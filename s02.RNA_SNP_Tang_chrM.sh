#!/bin/bash
indir=$1  ## Path to TopHat Result
sample=$2  ## one sample name per line
outdir=$3  ## StatInfo/02.SNP_Calling

ref=/datb1/wangrui/Database/IBM_Database/GATK/genome.fa
PICARD="java -jar /datc/wangrui/software/picard-tools-1.130/picard.jar"
GATK="java -jar /datb1/wangrui/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar"

knownSites_1=/datg/wangrui/Colon_Liver_addSample/Exon/bin/dbsnp_135.hg19.vcf
knownSites_2=/datc/wangrui/Database/indel_annotation/Mills_and_1000G_gold_standard.indels.hg19.vcf
knownSites_3=/datc/wangrui/Database/dbSNP/1000G_phase1.indels.hg19.vcf
vcf=/datc/wangrui/Database/indel_annotation/Mills_and_1000G_gold_standard.indels.hg19.vcf

function do_Remove_ERCC {
inpath=$1
#spike_in=/Share/BP/wangrui/Database/GATK/ERCC_RGC.list
#spike_in=/WPSnew/wangrui/Database/Share_Database/Database/GATK/ERCC_RGC.list
spike_in=/datg/wangrui/Lung_Organoid/RNA/FurtherAnalysis/01.chrM_Trance/bin/ERCC_RGC.list
sp=$2
outpath=$3

mkdir -p $outpath/$sp
cat $spike_in| grep -f - -v  <(samtools view -h $inpath/$sp/accepted_hits.bam) >  $outpath/$sp/$sp.NPI.sam 

}

function do_ReorderSam {
inpath=$1
PICARD=$PICARD
sp=$2
ref=$ref
outpath=$inpath

$PICARD ReorderSam I=$inpath/$sp/$sp.NPI.sam O=$outpath/$sp/$sp.reorder.sam REFERENCE=$ref 

}


function do_samTobam {
inpath=$1
outpath=$inpath
sp=$2
samtools view -bS $inpath/$sp/$sp.reorder.sam -o $outpath/$sp/$sp.reorder.bam
}


function do_RG_added_sorted {
inpath=$1
PICARD=$PICARD
sp=$2
outpath=$inpath

$PICARD AddOrReplaceReadGroups I=$inpath/$sp/$sp.reorder.bam O=$outpath/$sp/$sp.rg_added_sorted.bam SO=coordinate RGID=$sp RGLB=$sp RGPL=illumina RGPU=$sp RGSM=$sp

}

function do_mark_dup {
path=$1
PICARD=$PICARD
sp=$2
		$PICARD MarkDuplicates I=$path/$sp/$sp.rg_added_sorted.bam O=$path/$sp/$sp.dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$path/$sp/$sp.metrics
}

function do_splitNTrim {
path=$1
GATK=$GATK
sp=$2
ref=$ref
		$GATK -T SplitNCigarReads -R $ref -I $path/$sp/$sp.dedupped.bam -o $path/$sp/$sp.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS -L chrM

}

function do_RealignerTargetCreator {
inpath=$1
outpath=$inpath
GATK=$GATK
sp=$2

		$GATK -T RealignerTargetCreator -R $ref -I $inpath/$sp/$sp.split.bam -o $outpath/$sp/$sp.forIndelRealigner.intervals --known $vcf -L chrM
}


function do_IndelRealigner {
inpath=$1
outpath=$inpath
GATK=$GATK
sp=$2


		$GATK -T IndelRealigner -R $ref -I $inpath/$sp/$sp.split.bam -targetIntervals $outpath/$sp/$sp.forIndelRealigner.intervals -o $outpath/$sp/$sp.realignedBam.bam -known $vcf -L chrM

}

function do_BaseRecalibrator {
inpath=$1
outpath=$inpath
GATK=$GATK
sp=$2
ref=$ref
knownSites_1=$knownSites_1
knownSites_2=$knownSites_2
knownSites_3=$knownSites_3
		$GATK -T BaseRecalibrator -R $ref -I $inpath/$sp/$sp.realignedBam.bam -knownSites $knownSites_1 -knownSites $knownSites_2 -knownSites $knownSites_3 -o $outpath/$sp/$sp.recal.grp -L chrM
}

function do_BaseRecalibrator_BQSR {
inpath=$1
outpath=$inpath
GATK=$GATK
sp=$2
ref=$ref
knownSites_1=$knownSites_1
knownSites_2=$knownSites_2
knownSites_3=$knownSites_3
		$GATK -T BaseRecalibrator -R $ref -I $inpath/$sp/$sp.realignedBam.bam -knownSites $knownSites_1 -knownSites $knownSites_2 -knownSites $knownSites_3 -BQSR $inpath/$sp/$sp.recal.grp -o $outpath/$sp/$sp.pos_recal.grp -L chrM
}

function do_PrintReads {
inpath=$1
outpath=$inpath
GATK=$GATK
sp=$2
ref=$ref
knownSites_1=$knownSites_1
knownSites_2=$knownSites_2
knownSites_3=$knownSites_3
		$GATK -T PrintReads -R $ref -I $inpath/$sp/$sp.realignedBam.bam -BQSR $inpath/$sp/$sp.recal.grp -o $outpath/$sp/$sp.realn_Recal.bam -L chrM

}


function do_HaplotypeCaller {
inpath=$1
outpath=$inpath
GATK=$GATK
sp=$2
ref=$ref
		$GATK -T HaplotypeCaller -R $ref -I $inpath/$sp/$sp.realn_Recal.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -oi $outpath/$sp/$sp.variants_result.vcf -L chrM
}

function do_HaplotypeCaller_GVCF {
inpath=$1
outpath=$inpath
GATK=$GATK
sp=$2
ref=$ref
        $GATK -T HaplotypeCaller -R $ref -I $inpath/$sp/$sp.realn_Recal.bam -G StandardAnnotation --emitRefConfidence GVCF -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -variant_index_type LINEAR -variant_index_parameter 128000 -o $outpath/$sp/$sp.variants_result.AS.g.vcf -L chrM

		#$GATK -T UnifiedGenotyper -R $ref -I $inpath/$sp/$sp.realn_Recal.bam -stand_call_conf 20.0 -stand_emit_conf 20.0 -o $outpath/$sp/$sp.variants_result.UnifiedGenotyper.vcf
}


function do_HaplotypeCaller_GVCF_V2 {
inpath=$1
outpath=$inpath
GATK=$GATK
sp=$2
ref=$ref
        $GATK -T HaplotypeCaller -R $ref -I $inpath/$sp/$sp.realn_Recal.bam -L chrM -G StandardAnnotation --emitRefConfidence BP_RESOLUTION -dontUseSoftClippedBases -variant_index_type LINEAR -variant_index_parameter 128000 -o $outpath/$sp/$sp.variants_result_V2.AS.g.vcf

        #$GATK -T UnifiedGenotyper -R $ref -I $inpath/$sp/$sp.realn_Recal.bam -stand_call_conf 20.0 -stand_emit_conf 20.0 -o $outpath/$sp/$sp.variants_result.UnifiedGenotyper.vcf
}


function do_VarintFilter {
inpath=$1
outpath=$inpath
GATK=$GATK
sp=$2
ref=$ref
		$GATK -T VariantFiltration -R $ref -V $inpath/$sp/$sp.variants_result.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o $outpath/$sp/$sp.filtered_variants_result.vcf
}


function do_SnpEFF_annatation {
inpath=$1
outpath=$inpath
GATK=$GATK
sp=$2
#snpEff=/Share/BP/wangrui/software/snpEff/snpEff.jar
#snpEff=/WPSnew/wangrui/Software/snpEff_v4.3/snpEff/snpEff.jar
snpEff=/datg/wangrui/Software_NewVersion/snpEff_v4.3/snpEff/snpEff.jar
#		java -Xmx4g -jar $snpEff -c /Share/BP/wangrui/software/snpEff/snpEff.config hg19 $inpath/${sp}/${sp}.filtered_variants_result.vcf >$outpath/${sp}/${sp}.ann.vcf
#		java -Xmx10g -jar $snpEff -c /WPSnew/wangrui/Software/snpEff/snpEff.config hg19 $inpath/${sp}/${sp}.filtered_variants_result.vcf >$outpath/${sp}/${sp}.ann.vcf
		java -Xmx10g -jar $snpEff -c /datg/wangrui/Software_NewVersion/snpEff_v4.3/snpEff/snpEff.config hg19 $inpath/${sp}/${sp}.filtered_variants_result.vcf >$outpath/${sp}/${sp}.ann.vcf

}

function do_SnpEff_annotation_gVCF {
inpath=$1
outpath=$inpath
sp=$2
#snpEff=/Share/BP/wangrui/software/snpEff/snpEff.jar
#snpEff=/WPSnew/wangrui/Software/snpEff/snpEff.jar
#snpEff=/WPSnew/wangrui/Software/snpEff_v4.3/snpEff/snpEff.jar
snpEff=/datg/wangrui/Software_NewVersion/snpEff_v4.3/snpEff/snpEff.jar
#java -Xmx4g -jar $snpEff -c /Share/BP/wangrui/software/snpEff/snpEff.config hg19 $inpath/${sp}/${sp}.variants_result.AS.g.vcf > $outpath/${sp}/${sp}.ann.g.vcf
#java -Xmx10g -jar $snpEff -c /WPSnew/wangrui/Software/snpEff/snpEff.config hg19 $inpath/${sp}/${sp}.variants_result.AS.g.vcf > $outpath/${sp}/${sp}.ann.g.vcf
java -Xmx10g -jar $snpEff -c /datg/wangrui/Software_NewVersion/snpEff_v4.3/snpEff/snpEff.config hg19 $inpath/${sp}/${sp}.variants_result.AS.g.vcf > $outpath/${sp}/${sp}.ann.g.vcf

}

function do_UnifiedGenotyper {
inpath=$1
outpath=$inpath
GATK=$GATK
sp=$2
ref=$ref
        $GATK -T UnifiedGenotyper  -R $ref -I $inpath/$sp/$sp.realn_Recal.bam -o $outpath/$sp/$sp.variants_result.Unified.raw.vcf --output_mode EMIT_ALL_SITES


}



do_Remove_ERCC $indir  $sample $outdir #step 00
do_ReorderSam $outdir $sample          #step 01
do_samTobam $outdir $sample            #step 02
do_RG_added_sorted $outdir $sample     #step 03
do_mark_dup $outdir $sample            #step 04
do_splitNTrim $outdir $sample          #step 05
do_RealignerTargetCreator $outdir $sample #step 06
do_IndelRealigner $outdir $sample      #step 07
do_BaseRecalibrator $outdir $sample    #step 08
do_BaseRecalibrator_BQSR $outdir $sample #step 09
do_PrintReads $outdir $sample          #step 10
#do_HaplotypeCaller $outdir $sample     #step 11
#do_VarintFilter $outdir $sample        #step 12
#do_SnpEFF_annatation $outdir $sample   #step 13


#Option
#do_HaplotypeCaller_GVCF $outdir $sample
do_HaplotypeCaller_GVCF_V2 $outdir $sample
#do_SnpEff_annotation_gVCF $outdir $sample
#do_UnifiedGenotyper $outdir $sample

