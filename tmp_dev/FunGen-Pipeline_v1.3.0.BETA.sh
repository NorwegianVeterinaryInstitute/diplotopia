#!/bin/bash
#SBATCH --account=nn9305k

##memory specs
#SBATCH --mem=38G     ## How much memory does your job need
#SBATCH --cpus-per-task=12    ## total number of cpus, max on SAGA is 40
#SBATCH --time=12:00:00
#SBATCH --job-name=FunGen

set -o errexit ##Recommended for easier debugging, it kills the job when something goes wrong.
set -o nounset ## Treat any unset variables as an error

function Usage ()
{
cat <<-ENDOFMESSAGE
#####################################################
##                 FunGen-Pipeline                 ##
##                Version 1.3.0_BETA               ##
#####################################################

Usage

$0 [OPTIONS] -R <reference.fasta> -1 <file_R1.fastq> -2 <file1_R2.fastq> -gff <reference.gff> -p <your_prefix> -P <ploidy>

Flag description:

	-R --reference	your fasta reference (required)
	-1 --fq1		your fastq_R1 file (required)
	-2 --fq2		your fastq_R2 file (required)
	-gff --gff		your reference gff file (required)
	-p --prefix		prefix for the output files (default = myprefix)

Options:
	-P --ploidy		Organism ploidy (default = 1)
	-h --help		display this message

ENDOFMESSAGE
	exit 1
}

function Die ()
{
	echo "$*"
	exit 1
}

function GetOpts() {
	reference=""
	fq1=""
	fq2=""
	gff=""
	prefix="myprefix"
	ploidy=1
	# argv=()
	while [ $# -gt 0 ]
	do
		opt=$1
		shift
		case ${opt} in
			-R|--reference)
				arg_has_reference=1
				if [ $# -eq 0 -o "${1:0:1}" = "-" ]; then
					Die "The ${opt} option requires an argument."
				fi
				reference="$1"
				shift
				;;
			-1|--fq1)
				arg_has_fq1=1
				if [ $# -eq 0 -o "${1:0:1}" = "-" ]; then
					Die "The ${opt} option requires an argument."
				fi
				fq1="$1"
				shift
				;;
			-2|--fq2)
				arg_has_fq2=1
				if [ $# -eq 0 -o "${1:0:1}" = "-" ]; then
					Die "The ${opt} option requires an argument."
				fi
				fq2="$1"
				shift
				;;
			-p|--prefix)
				arg_has_prefix=1
				if [ $# -eq 0 -o "${1:0:1}" = "-" ]; then
					Die "The ${opt} option requires an argument."
				fi
				prefix="$1"
				shift
				;;
#			-gff|--gff)
#				arg_has_gff=1
#				if [ $# -eq 0 -o "${1:0:1}" = "-" ]; then
#					Die "The ${opt} option requires an argument."
#				fi
#				gff="$1"
#				shift
#				;;
			-P|--ploidy)
				arg_has_ploidy=1
				if [ $# -eq 0 -o "${1:0:1}" = "-" ]; then
					Die "The ${opt} option requires an argument."
				fi
				ploidy="$1"
				shift
				;;
			-h|--help)
				Usage;;
			*)
				if [ "${opt:0:1}" = "-" ]; then
					Die "${opt}: unknown option."
				fi
#					argv+=(${opt});;
		esac
	done
	
		
		if [ $arg_has_reference -ne 0 ]; then
				GENOME=$reference
		else
			echo "Error: -R (reference.fasta) absent"
				Usage
		fi
		
		if [ $arg_has_fq1 -ne 0 ]; then
				R1=$fq1
		else
			echo "Error: fq1 file absent"
				Usage
		fi
		
		if [ $arg_has_fq2 -ne 0 ]; then
				R2=$fq2
		else
			echo "Error: fq2 file absent"
				Usage
		fi
		
#		if [ $arg_has_gff -ne 0 ]; then
#				GFF=$gff
#		else
#			echo "Error: gff file absent"
#				Usage
#		fi
		
		if [ $arg_has_ploidy -ne 0 ]; then
				PLOIDY=$ploidy
		else
			echo "Error: ploidy file absent"
				Usage
		fi
		
		if [ $arg_has_prefix -ne 0 ]; then
				SAMPLE=$prefix
		else
			echo "Error: prefix file absent"
				Usage
#		if [ -z "$prefix" ]; then
#				SAMPLE="default"

#echo $SAMPLE
		exit 1
	fi
}

GetOpts $*

#GENOME=${reference}
#R1=$fq1
#R2=$fq2
#SAMPLE=$prefix
#SAMPLE=BG2
echo "Reference genome file used: $GENOME"
echo "Fastq1 file used: $R1"
echo "Fastq2 file used: $R2"
echo "Prefix used: $SAMPLE"
#echo "GFF# file used: $GFF"
echo "Ploidy used: $PLOIDY"
dict="${GENOME%.*}"


########################
### STEP-1 - Reconstruction of short reads using a ref_genome

echo "######### STEP-1 - Reconstruction of short reads using a ref_genome #########"

module load BWA/0.7.18-GCCcore-12.3.0

bwa mem -M $GENOME $R1 $R2 > $SAMPLE.sam

# cleaning up standard environment

echo "### Clean up std_evn ###"

module purge

### STEP-2 - conversion, coordinate sorting and indexing

echo "######### STEP-2 - conversion, coordinate sorting and indexing #########"

module load SAMtools/1.19.2-GCC-13.2.0
module load picard/3.0.0-Java-17
module list
samtools view -bS $SAMPLE.sam > $SAMPLE.bam
samtools sort $SAMPLE.bam -o $SAMPLE.sorted.bam
samtools index $SAMPLE.sorted.bam

### STEP-3 - Fixing the read groups

echo "######### STEP-3 - Fixing the read groups #########"

java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups -I $SAMPLE.sorted.bam -O $SAMPLE.fixed.sorted.bam --SORT_ORDER coordinate --RGID 241011_LH00534 --RGLB dnaseq --RGPL illumina --RGSM 'WGS' --CREATE_INDEX TRUE --RGPU Unknown --VALIDATION_STRINGENCY SILENT

echo "### Clean up std_env ###"

module purge

### STEP-4 Raw variant calling

echo "######### STEP-4 Raw variant calling #########"

module load picard/3.0.0-Java-17
module load GATK/4.5.0.0-GCCcore-12.3.0-Java-17
module list

echo "### Picard - MarkDuplicates ###"

java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I $SAMPLE.fixed.sorted.bam -O $SAMPLE.sorted.marked.bam --CREATE_INDEX TRUE --METRICS_FILE picard_info_$SAMPLE.txt --REMOVE_DUPLICATES false --ASSUME_SORTED true --VALIDATION_STRINGENCY SILENT

echo "### GATK - HaplotypeCaller ###"

gatk HaplotypeCaller -R $GENOME -I $SAMPLE.sorted.marked.bam -ploidy 1 -O $SAMPLE.raw_variants.vcf

### STEP-5 - Base Quality Score Recalibration - First pass

echo "######### STEP-5 - Base Quality Score Recalibration - First pass #########"
echo "### GATK - BaseRecalibrator ###"

gatk BaseRecalibrator -R $GENOME -I $SAMPLE.sorted.marked.bam --known-sites $SAMPLE.raw_variants.vcf -O $SAMPLE.recal_data.table

echo "### GATK - ApplyBQSR ###"

gatk ApplyBQSR -R $GENOME -I $SAMPLE.sorted.marked.bam -bqsr-recal-file $SAMPLE.recal_data.table -O $SAMPLE.recal_reads.bam

### STEP-6 - Base Quality Score Recalibration - Second pass

echo "######### STEP-6 - Base Quality Score Recalibration - Second pass #########"
echo "### GATK - BaseRecalibrator ###"

gatk BaseRecalibrator -R $GENOME -I $SAMPLE.recal_reads.bam --known-sites $SAMPLE.raw_variants.vcf -O $SAMPLE.post_recal_data.table

echo "### GATK - ApplyBQSR ###"
gatk ApplyBQSR -R $GENOME -I $SAMPLE.recal_reads.bam -bqsr-recal-file $SAMPLE.post_recal_data.table -O $SAMPLE.post_recal_reads.bam

### STEP-7 - Quality Control

echo "######### STEP-7 - Quality Control #########"

java -jar $EBROOTPICARD/picard.jar CollectWgsMetrics -R $GENOME -I $SAMPLE.post_recal_reads.bam -O $SAMPLE.metrics.txt

echo "### SAMtools - clean std_env and start module ###"

module purge
module load SAMtools/1.19.2-GCC-13.2.0

echo "### Get coverage data ###"

samtools depth $SAMPLE.post_recal_reads.bam > $SAMPLE.coverage.txt

### STEP-8 - Obtaining mapping stats

echo "######### STEP-8 - Obtaining mapping stats #########"

samtools flagstat $SAMPLE.post_recal_reads.bam > $SAMPLE.flagstat.txt

### STEP-9 - Calling high quality variants

echo "######### STEP-9 - Calling high quality variants #########"
echo "### Clean up std_env ###"

module purge
module load GATK/4.5.0.0-GCCcore-12.3.0-Java-17

echo "### GATK HaplotypeCaller ###"
gatk HaplotypeCaller -R $GENOME -I $SAMPLE.post_recal_reads.bam -O $SAMPLE.raw_variants_recal.vcf -ERC GVCF --pcr-indel-model NONE -ploidy $PLOIDY -stand-call-conf 30 -mbq 20 -A QualByDepth

echo "### GATK - GenotypeGVCFs ###"
gatk GenotypeGVCFs -R $GENOME -V $SAMPLE.raw_variants_recal.vcf -O $SAMPLE.genotyped_variants_recal.vcf

###
### TEST OF AN EXTENDED FILTRATION, MERGING AND LowConf EXCLUDING FOLLOWED BY ALTERNATING REFERENCE GENOME
###

### STEP-10 - Filtering SNPs based on quality

echo "######### STEP-10 - Filtering SNPs based on quality #########"
echo "### GATK SelectVariants - SNPs ###"
gatk SelectVariants -R $GENOME -V $SAMPLE.genotyped_variants_recal.vcf -O $SAMPLE.raw_snps_recal.vcf --select-type-to-include SNP -select "vc.getGenotype('WGS').getAD().get(1) > 0 && vc.getGenotype('WGS').getDP() > 0 && vc.getGenotype('WGS').getAD().get(1) * 1.0 / vc.getGenotype('WGS').getDP() > 0.90"

echo "### GATK - SelectVariants -indels ###"
gatk SelectVariants -R $GENOME -V $SAMPLE.raw_variants_recal.vcf -O $SAMPLE.raw_indels_recal.vcf --select-type-to-include INDEL

echo "### GATK - VariantFiltration on SNP file###"
gatk VariantFiltration -R $GENOME -V $SAMPLE.raw_snps_recal.vcf -O $SAMPLE.filtered_snps_final.vcf -filter "QD < 2.0" --filter-name "LowConf" -filter "FS > 60.0" --filter-name "LowConf" -filter "MQ < 40.0" --filter-name "LowConf" -filter "MQRankSum < -12.5" --filter-name "LowConf" -filter "ReadPosRankSum < -8.0" --filter-name "LowConf" -filter "SOR > 4.0" --filter-name "LowConf" -filter "DP < 5" --filter-name "LowConf" -G-filter "GQ < 50" -G-filter-name "FILTER_GQ-50"

echo "### GATK - VariantFiltration on INDELs file###"
gatk VariantFiltration -R $GENOME -V $SAMPLE.raw_indels_recal.vcf -O $SAMPLE.filtered_indels_final.vcf -filter "QD < 2.0" --filter-name "LowConf" -filter "FS > 60.0" --filter-name "LowConf" -filter "MQ < 40.0" --filter-name "LowConf" -filter "MQRankSum < -12.5" --filter-name "LowConf" -filter "ReadPosRankSum < -8.0" --filter-name "LowConf" -filter "SOR > 4.0" --filter-name "LowConf" -filter "DP < 5" --filter-name "LowConf" -G-filter "GQ < 50" -G-filter-name "FILTER_GQ-50"


### Merging SNP and INDEL files
echo "Merging SNP and INDEL files"
gatk MergeVcfs -I "$SAMPLE.filtered_snps_final.vcf" -I "$SAMPLE.filtered_indels_final.vcf" -O "$SAMPLE.merged_variants.vcf"

### Line count in filtered SNP file
echo "Merging SNP and INDEL files"
snp_lines=$(grep -v "#" "$SAMPLE.filtered_snps_final.vcf" | wc -l)

### Line count in filtered INDEL file
echo "Line count in filtered INDEL file"
indel_lines=$(grep -v "#" "$SAMPLE.filtered_indels_final.vcf" | wc -l)

### Line count in merged files
echo "Line count in merged files"
merged_lines=$(grep -v "#" "$SAMPLE.merged_variants.vcf" | wc -l)

### Merge effect control
echo "########################################################"
echo "###############                          ###############"
echo "###############   Merge effect control   ###############"
echo "###############                          ###############"
echo "###############   $SAMPLE                ###############"
echo "########################################################"

    if [ "$((snp_lines + indel_lines))" -eq "$merged_lines" ]; then
      echo "SNP + INDEL were paired correctly."
    else
      echo "ERROR: SNP + INDEL not merged correctly."
      echo "SNP linie: $snp_lines"
      echo "INDEL linie: $indel_lines"
      echo "Merged file: $merged_lines"
    fi
    
echo "######### STEP-11 - Filtration number control #########"
echo "Genotyped variant positions (#) prior to quality filtration:"
grep -v "#" "$file" | wc -l
    
echo "Final variants - post-filtration:"
grep -v "#" "$SAMPLE.merged_variants.vcf" | grep -v "LowConf" | grep -v "FILTER_GQ-50" | wc -l

echo "Preparation of a no-LowConf file for genome reconstruction and phylogenomics"
grep -v "LowConf" "$SAMPLE.merged_variants.vcf" > "$SAMPLE.noLowQual.vcf"

echo "Indexing refined noLowConf file"
gatk IndexFeatureFile -I "$SAMPLE.noLowQual.vcf"

echo "NOW - GET A FASTA FORMATTED ALTERNATIVE GENOME"
gatk FastaAlternateReferenceMaker -R "$GENOME" -V "$SAMPLE.noLowQual.vcf" -O "$SAMPLE.consensus_genome.fasta"


echo "########################################################"
echo "###############                          ###############"
echo "###############   $SAMPLE                ###############"
echo "###############   process FINISHED       ###############"
echo "########################################################"


###
### THE BELOW FRAGMENT IS TO BE REMOVED AFTER TESTING
###

### STEP-10 - Filtering SNPs based on quality

echo "######### STEP-10 - Filtering SNPs based on quality #########"
echo "### GATK SelectVariants - SNPs ###"
gatk SelectVariants -R $GENOME -V $SAMPLE.genotyped_variants_recal.vcf -O $SAMPLE.raw_snps_recal.vcf --select-type-to-include SNP -select "vc.getGenotype('WGS').getAD().get(1) > 0 && vc.getGenotype('WGS').getDP() > 0 && vc.getGenotype('WGS').getAD().get(1) * 1.0 / vc.getGenotype('WGS').getDP() > 0.90"

echo "### GATK - SelectVariants -indels ###"
gatk SelectVariants -R $GENOME -V $SAMPLE.raw_variants_recal.vcf -O $SAMPLE.raw_indels_recal.vcf --select-type-to-include INDEL

echo "### GATK - VariantFiltration ###"
gatk VariantFiltration -R $GENOME -V $SAMPLE.raw_snps_recal.vcf -O $SAMPLE.filtered_snps_final.vcf -filter "QD < 2.0" --filter-name "LowConf" -filter "FS > 60.0" --filter-name "LowConf" -filter "MQ < 40.0" --filter-name "LowConf" -filter "MQRankSum < -12.5" --filter-name "LowConf" -filter "ReadPosRankSum < -8.0" --filter-name "LowConf" -filter "SOR > 4.0" --filter-name "LowConf" -filter "DP < 5" --filter-name "LowConf" -G-filter "GQ < 50" -G-filter-name "FILTER_GQ-50"

### STEP-11 - Filtration number control
echo "######### STEP-11 - Filtration number control #########"
echo "Genotyped variant positions (#) prior to quality filtration:"
grep -v "#" $SAMPLE.genotyped_variants_recal.vcf | wc -l
echo "Final variants - post-filtration:"
grep -v "#" $SAMPLE.filtered_snps_final.vcf | grep -v "LowConf" | grep -v "FILTER_GQ-50" | wc -l
