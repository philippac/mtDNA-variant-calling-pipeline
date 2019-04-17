

#####################

#!/bin/bash
samtoolsVIEW="/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/samtools-1.8/samtools view"
SamToFastq="java -Xmx2g -jar /homes/athosnew/Genetics_Centre_Bioinformatics/resourses/picard-2.815/picard.jar SamToFastq"	
rCRS="/homes/athosnew/Genetics_Centre_Bioinformatics/mtDNA_ampliconSeq/mtDNA_reference/rCRS.fasta"
bwa_mem="/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/bwa/bwa mem"
samtoolsSORT="/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/samtools-1.8/samtools sort"
samtoolsINDEX="/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/samtools-1.8/samtools index"
MarkDuplicates="/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/picard-2.815/picard.jar MarkDuplicates"
Mutect2="/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/GATK_v4_BETA_5/GenomeAnalysisTK.jar Mutect2"
FilterMutectCalls="/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/GATK_v4_BETA_5/GenomeAnalysisTK.jar FilterMutectCalls"
SelectVariants="/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/GATK_v4_BETA_5/GenomeAnalysisTK.jar SelectVariants"
CollectRawWgsMetrics="/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/picard-2.815/picard.jar CollectRawWgsMetrics"

myPROJECT=$1
masterDIR="/homes/athosnew/Genetics_Centre_Bioinformatics/mtDNA_exomeSeq/philippa"
inputDIR="/homes/athosnew/Genetics_Centre_Bioinformatics/Exomes/Aligned"

mkdir $masterDIR/$myPROJECT
projectDIR=$masterDIR/$myPROJECT

echo $projectDIR

samples=`ls $inputDIR/$myPROJECT`
echo $samples

mkdir $projectDIR/Bam_Files
mkdir $projectDIR/VCF_Files
mkdir $projectDIR/Coverage_data
mkdir $projectDIR/ReadCount_data
##extract mtdna

for sample in $samples; do

	mkdir $projectDIR/Bam_Files/$sample
	mkdir $projectDIR/VCF_Files/$sample
	mkdir $projectDIR/Coverage_data/$sample
	mkdir $projectDIR/ReadCount_data/$sample
	
		$samtoolsVIEW \
		-b $inputDIR/$myPROJECT/$sample/${sample}_sorted_unique_recalibrated.bam MT \
		> $projectDIR/Bam_Files/$sample/${sample}_raw_mtdna.bam
		
		echo "extraction done"
		
##bam to fastq		
		$SamToFastq \
		I=$projectDIR/Bam_Files/$sample/${sample}_raw_mtdna.bam \
		F=$projectDIR/Bam_Files/$sample/${sample}_mtdna_r1.fastq \
		F2=$projectDIR/Bam_Files/$sample/${sample}_mtdna_r2.fastq \
		INCLUDE_NON_PRIMARY_ALIGNMENTS=false \
		INCLUDE_NON_PF_READS=true \
		COMPRESSION_LEVEL=0 \
		VALIDATION_STRINGENCY=SILENT \
		QUIET=true
		echo "fastq done"

##align to reference

		$bwa_mem \
		-t 4 \
		-R "@RG\\tID:${sample}\\tSM:${sample}\\tPL:ILLUMINA" \
		$rCRS \
		$projectDIR/Bam_Files/$sample/${sample}_mtdna_r1.fastq \
		$projectDIR/Bam_Files/$sample/${sample}_mtdna_r2.fastq \
		>$projectDIR/Bam_Files/$sample/${sample}_mtdna.sam
		echo "aligment done"
	
##post aligment proccessing
##convert sam to bam

		$samtoolsVIEW \
		-bh $projectDIR/Bam_Files/$sample/${sample}_mtdna.sam \
		>$projectDIR/Bam_Files/$sample/${sample}_mtdna.bam
		echo "sam to bam done"
		
##sort bam

		$samtoolsSORT \
		$projectDIR/Bam_Files/$sample/${sample}_mtdna.bam \
		>$projectDIR/Bam_Files/$sample/${sample}_mtdna.sorted.bam
				
##index bam

		$samtoolsINDEX \
		$projectDIR/Bam_Files/$sample/${sample}_mtdna.sorted.bam
		echo "bam sort and index done"
	
##markduplicates

		java -jar $MarkDuplicates \
		I=$projectDIR/Bam_Files/$sample/${sample}_mtdna.sorted.bam \
		O=$projectDIR/Bam_Files/$sample/MD_${sample}.bam \
		M=$projectDIR/Bam_Files/$sample/${sample}duplication_metrics.txt

##index markduplicates

		$samtoolsINDEX \
		$projectDIR/Bam_Files/$sample/MD_${sample}.bam
		echo "markduplicates and index done"
		
##variant calling

		java -jar $Mutect2 \
		-R $rCRS \
		-I $projectDIR/Bam_Files/$sample/MD_${sample}.bam \
		-O $projectDIR/VCF_Files/$sample/VC_${sample}.vcf \
		-tumor $sample
		echo "variant calling done"
		
		
##filter called variants

		java -jar $FilterMutectCalls \
		-V $projectDIR/VCF_Files/$sample/VC_${sample}.vcf \
		-O $projectDIR/VCF_Files/$sample/VC_${sample}_filtered.vcf \
		--max_germline_posterior 1000 \
		-maxEventsInHaplotype 1000
	
	
		java -Xmx20g -jar $SelectVariants \
		-R $rCRS \
		-V $projectDIR/VCF_Files/$sample/VC_${sample}_filtered.vcf \
		-O $projectDIR/VCF_Files/$sample/${sample}_passed.vcf \
		--select "vc.isNotFiltered()"
		echo "Filtering of called variants done"

##add readcount
	
		$samtoolsVIEW -c \
		$projectDIR/Bam_Files/$sample/MD_${sample}.bam \
		> $projectDIR/ReadCount_data/$sample/MD_${sample}_readcount.bam
		echo "Read count done"
##Summarymetrics		
		
		java -jar $CollectRawWgsMetrics \
		R=$rCRS \
		INPUT=$projectDIR/Bam_Files/$sample/MD_${sample}.bam \
		OUTPUT=$projectDIR/Coverage_data/$sample/MD_${sample}_WgsMetrics.txt \
		INCLUDE_BQ_HISTOGRAM=true
		echo "wgsmetrics done"	
done

exit