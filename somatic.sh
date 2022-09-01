#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'


################################################################################
# Input Parameters
################################################################################
export PATH="/storeData/project/user/xuzhenyu/software/micromamba/envs/snpsift/bin:${PATH}"
export PATH="/storeData/project/user/xuzhenyu/software/micromamba/envs/cdx/bin:${PATH}"
REF="/storeData/project/user/xuzhenyu/database/genome/hg19/hg19.fa"
dbsnp="/storeData/project/user/xuzhenyu/cdx/database/gatk_bundle/common_all_20180423.vcf.gz"
Mills="/storeData/project/user/xuzhenyu/cdx/database/gatk_bundle/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
g1000="/storeData/project/user/xuzhenyu/cdx/database/gatk_bundle/1000G_phase1.snps.high_confidence.hg19.sites.vcf"
cosmic_vcf="/storeData/project/user/xuzhenyu/cdx/database/cosmic/CosmicCodingMuts.chr.vcf.gz"
annovar_db="/storeData/project/user/xuzhenyu/cdx/database/annovar/humandb"

# data.list
DATA_LIST="${1}"
SM=$(tail -1 "${DATA_LIST}" | awk '{print $1}')
R1=$(tail -n +2 "${DATA_LIST}" | awk '{print $2}' | sed ':t;N;s/\n/ /;b t' -)
R2=$(tail -n +2 "${DATA_LIST}" | awk '{print $3}' | sed ':t;N;s/\n/ /;b t' -)
IFS=' ' read -ra R2_args <<< "${R2}"
IFS=' ' read -ra R1_args <<< "${R1}"


OUT_DIR=${2}
NCPU=${3:-40}
cd "${OUT_DIR}" || exit

# check input
printf 'DATA_LIST:\t%s\n' "${DATA_LIST}"
printf 'REF:\t%s\n' "${REF}"
printf 'OUT_DIR:\t%s\n' "${OUT_DIR}"
printf 'NCPU:\t%s\n' "${NCPU}"

# directories
DATA_DIR="${OUT_DIR}/1.data"
ALIGN_DIR="${OUT_DIR}/2.alignment"
QC_DIR="${OUT_DIR}/3.qc"
SOMATIC_DIR="${OUT_DIR}/4.somatic"
ANNOTATION_DIR="${OUT_DIR}/5.annotation"

[ -d "${DATA_DIR}" ] || mkdir -p "${DATA_DIR}"
[ -d "${ALIGN_DIR}" ] || mkdir -p "${ALIGN_DIR}"
[ -d "${QC_DIR}" ] || mkdir -p "${QC_DIR}"
[ -d "${SOMATIC_DIR}" ] || mkdir -p "${SOMATIC_DIR}"
[ -d "${ANNOTATION_DIR}" ] || mkdir -p "${ANNOTATION_DIR}"


################################################################################
# QC with fastp
################################################################################
R1_CLEAN="${DATA_DIR}/${SM}.1.clean.fq.gz"
R2_CLEAN="${DATA_DIR}/${SM}.2.clean.fq.gz"
fastp \
    --stdin \
    -i <(zcat "${R1_args[@]}") \
    -I <(zcat "${R2_args[@]}") \
    --thread 16 \
    -f 0 \
    -F 0 \
    -o "${R1_CLEAN}" \
    -O "${R2_CLEAN}" \
    -h "${DATA_DIR}/${SM}.qc.fastp.html" \
    -j "${DATA_DIR}/${SM}.qc.fastp.json"


################################################################################
# mapping with bwa
################################################################################
rawbam="${ALIGN_DIR}/${SM}.raw.bam"
addGRPbam="${ALIGN_DIR}/${SM}.rg.bam"
sortbam="${ALIGN_DIR}/${SM}.sorted.bam"
bwa mem \
    -M \
    -t ${NCPU} \
    -K 10000000 \
    -R '@RG\tID:'"${SM}"'\tPL:Illumina\tLB:Target\tSM:'"${SM}" \
    "${REF}" \
    "${R1_CLEAN}" "${R2_CLEAN}" | \
    samtools view -@ ${NCPU} -b -o "$rawbam" -
# Add groups
gatk AddOrReplaceReadGroups \
    -I "${rawbam}" \
    -O "$addGRPbam" \
    -RGLB hg19 \
    -RGPL Illumina \
    -RGSM "${SM}" \
    -RGPU new \
    --VALIDATION_STRINGENCY SILENT
# sort
gatk SortSam \
    -I "$addGRPbam" \
    -O "$sortbam" \
    -SO coordinate \
    --VALIDATION_STRINGENCY SILENT


################################################################################
# mark duplicate with picard
################################################################################
markdupbam="${ALIGN_DIR}/${SM}.md.bam"
markdupmetrics="${ALIGN_DIR}/${SM}.markdup.metrics"
gatk MarkDuplicates \
    -I "$sortbam" \
    -O "$markdupbam" \
    -M "$markdupmetrics" \
    --VALIDATION_STRINGENCY SILENT
samtools index -@ "${NCPU}" "$markdupbam"


################################################################################
# recalibrate base quality score
################################################################################
bqsrrecal="${ALIGN_DIR}/${SM}.bqsr.metrics"
bqsrrecalbam="${ALIGN_DIR}/${SM}.bqsr.bam"
gatk BaseRecalibrator \
    -R "$REF" \
    -I "$markdupbam" \
    -O "$bqsrrecal" \
    --known-sites $dbsnp \
    --known-sites $Mills \
    --known-sites $g1000 \
    --known-sites $cosmic_vcf \
    --interval-padding 100
gatk ApplyBQSR \
    -R "$REF" \
    -I "$markdupbam" \
    -O "$bqsrrecalbam" \
    --bqsr-recal-file "$bqsrrecal" \
    --static-quantized-quals 10 \
    --static-quantized-quals 20 \
    --static-quantized-quals 30 \
    --add-output-sam-program-record \
    --create-output-bam-md5 \
    --use-original-qualities \
    --interval-padding 100


################################################################################
# gatk Mutect2
################################################################################
m2_raw_vcf="${SOMATIC_DIR}/${SM}.m2.raw.vcf"
m2_fix_vcf="${SOMATIC_DIR}/${SM}.m2.fix.vcf"
m2_filter_vcf="${SOMATIC_DIR}/${SM}.m2.filter.vcf"
gatk Mutect2 \
    --max-mnp-distance 2 \
    --native-pair-hmm-threads "${NCPU}" \
    -R "$REF" \
    -I "${bqsrrecalbam}" \
    --max-reads-per-alignment-start 0 \
    --kmer-size 10 \
    --kmer-size 25 \
    --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
    -A AlleleFraction \
    -A BaseQuality \
    -A ChromosomeCounts \
    -A ClippingRankSumTest \
    -A Coverage \
    -A DepthPerAlleleBySample \
    -A DepthPerSampleHC \
    -A ExcessHet \
    -A FisherStrand \
    -A InbreedingCoeff \
    -A LikelihoodRankSumTest \
    -A MappingQuality \
    -A MappingQualityRankSumTest \
    -A QualByDepth \
    -A RMSMappingQuality \
    -A ReadPosRankSumTest \
    -A ReadPosition \
    -A OrientationBiasReadCounts \
    -A StrandBiasBySample \
    -A StrandOddsRatio \
    -A TandemRepeat \
    -A UniqueAltReadCount \
    --interval-padding 30 \
    -O "${m2_raw_vcf}"
# remove INFO tag AS_UNIQ_ALT_READ_COUNT
bcftools annotate \
    -x 'INFO/AS_UNIQ_ALT_READ_COUNT' \
    "${m2_raw_vcf}" \
    > "${m2_fix_vcf}"
gatk FilterMutectCalls \
    -R "${REF}" \
    -V "${m2_fix_vcf}" \
    -O "${m2_filter_vcf}" \
    --stats "${m2_raw_vcf}.stats"


################################################################################
# indel realign
################################################################################
m2_norm_vcf="${SOMATIC_DIR}/${SM}.m2.norm.vcf"
m2_aligned_vcf="${SOMATIC_DIR}/${SM}.m2.aligned.vcf"
bcftools norm \
    -m -both \
    "${m2_filter_vcf}" \
    > "${m2_norm_vcf}"
gatk LeftAlignAndTrimVariants \
    -R "${REF}" \
    -V "${m2_norm_vcf}" \
    -O "${m2_aligned_vcf}"


################################################################################
# snpEff Annotation
################################################################################
m2_snpeff_vcf="${ANNOTATION_DIR}/${SM}.m2.snpeff.vcf"
snpEff -Xmx16G ann -noStats -no PROTEIN_STRUCTURAL_INTERACTION_LOCUS \
    -no-downstream -no-intergenic -no-upstream -no NEXT_PROT \
    -hgvs1LetterAa \
    -hgvsTrId hg19 "${m2_aligned_vcf}" | \
    SnpSift -Xmx16G annotate -id "$cosmic_vcf" - | \
    SnpSift -Xmx16G annotate -id "$dbsnp" - \
    > "${m2_snpeff_vcf}"


################################################################################
# ANNOVAR Annotation
################################################################################
m2_annovar_prefix="${ANNOTATION_DIR}/${SM}.m2.snpeff.annovar"
m2_annovar_vcf="${ANNOTATION_DIR}/${SM}.m2.snpeff.annovar.hg19_multianno.vcf"
table_annovar.pl \
    "${m2_snpeff_vcf}" \
    "${annovar_db}" \
    --buildver hg19 \
    --remove --polish \
    --protocol refGene,cytoBand,esp6500siv2_all,1000g2015aug_all,exac03,gnomad_genome,gnomad_exome,avsnp150,cosmic70,clinvar_20220320,dbnsfp42a,icgc28 \
    --nastring . \
    --operation g,r,f,f,f,f,f,f,f,f,f,f \
    --vcfinput \
    --outfile "${m2_annovar_prefix}" \
    --argument "--splicing_threshold 3 --hgvs,--hgvs,--hgvs,--hgvs,--hgvs,--hgvs,--hgvs,--hgvs,--hgvs,--hgvs,--hgvs,--hgvs"


################################################################################
# vcf conversion
################################################################################
m2_annovar_txt="${ANNOTATION_DIR}/${SM}.m2.snpeff.annovar.txt"
m2_annovar_xlsx="${ANNOTATION_DIR}/${SM}.m2.snpeff.annovar.xlsx"
bcftools query -H \
    -f'%CHROM\t%TYPE\t%POS\t%REF\t%ALT[\t%AF][\t%AD][\t%DP]\t%FILTER\t%AAChange.refGene\t%Func.refGene\t%ExonicFunc.refGene\t%cosmic70\t%Gene.refGene\t%ANN\n' \
    "${m2_annovar_vcf}" | \
    sed '1c CHROM\tTYPE\tPOS\tREF\tALT\tAF\tAD\tDP\tFILTER\tAAChange\tFunc\tExonicFunc\tCOSMIC\tGENE\tANN' \
    > "${m2_annovar_txt}"
# vcf to xlsx
python /storeData/project/user/xuzhenyu/cdx/pipeline/txt2excel.py \
    "${m2_annovar_txt}" "${m2_annovar_xlsx}"


################################################################################
# split vcf
################################################################################
split_xlsx="${ANNOTATION_DIR}/${SM}.report.VCF.VUS.xlsx"
python /storeData/project/user/xuzhenyu/cdx/st_reporter/utilities/split_vcf.py \
    "${m2_annovar_xlsx}" \
    "${split_xlsx}" \
    --THERAPIES /storeData/project/user/xuzhenyu/cdx/st_reporter/blob/MCG_DB/raw_st/MCG_ST_V4/therapies.json \
    --TRIALS /storeData/project/user/xuzhenyu/cdx/st_reporter/blob/MCG_DB/raw_st/MCG_ST_V4/trials.json \
    --COSM_FILE /storeData/project/user/xuzhenyu/cdx/st_reporter/blob/CosmicMutantExportNM.txt \
    --EXON_INDEL_HGVS /storeData/project/user/xuzhenyu/cdx/st_reporter/blob/Exon.INDELs.HGVS.xlsx