#
# Purpose


# best practise protocol version 4 of the Broad Institute 
(http://www.broadinstitute.org/gsa/wiki/index.php/ Best_Practice_Variant_Detection_with_the_GATK_v4);

We begin by mapping the sequence reads to the reference genome to produce a file in SAM/BAM format sorted by coordinate. Next, we mark duplicates to mitigate biases introduced by data generation steps such as PCR amplification.

Finally, we recalibrate the base quality scores, because the variant calling algorithms rely heavily on the quality scores assigned to the individual base calls in each sequence read.


```BASH
EXPORT=/LUSTRE/apps/bioinformatica/bwa
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/STAR-2.7.11b/bin/Linux_x86_64_static/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/gatk-4.6.0.0
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/RSEM/bin/samtools-1.3/
export PATH=$PATH:$EXPORT

```

Prior to mapping step, you may considere first construct SuperTranscripts, which attempt to collapse overlapping contigs within de novo assembly (See freedman et al., 2019 and Trinity protocol), 

# 1) Map to reference
Tools involved: BWA, MergeBamAlignments

Mapping the sequence reads to the reference Genome/Transcriptome.

Samples were individually mapped against the transcriptome assembly using the Burrows–Wheeler Aligner (Li & Durbin 2009), and default parameters (except n = 0.005 and k = 5).

*Note*
For RNAseq short variant discovery (SNPs + Indels), Gatk recommed use STAR aligner because it increased sensitivity compared to TopHat (especially for INDELS).

```BASH
reference=evigene_transcripts.cdna

mkdir -p bwa_index

if [ ! -f "bwa_index/${reference}.bwt" ]; then
    bwa index $reference
    
    mv ${reference}.* bwa_index

else
    echo "index for '$reference' already exists."
fi

# Alignment for DNAseq

bwa mem -t 20 INDEX/evigene_transcripts.cdna BLA_C1_1P.fq BLA_C1_2P.fq

# Alignment for RNAseq using STAR

mkdir -p STAR_index

if [ ! -f "STAR_index/Genome" ]; then
    
    STAR --runMode genomeGenerate --genomeDir STAR_index --genomeFastaFiles $reference --limitGenomeGenerateRAM 45284502453 --runThreadN 20 --genomeSAindexNbases 12

else
    echo "index for '$reference' already exists."
fi




# Align

# use comma separated list for read1, followed by space, followed by comma separaeted list for read2, e.g.: s1read1.fq,s2read1.fq s1read2.fq,s2read2.fq

# The default MAPQ=255 for the unique mappers maybe changed with --outSAMmapqUnique parameter (integer 0 to 255) to ensure compatibility with downstream tools such as GATK.

# uT SAM tag indicates reason for not mapping:



gatk MergeBamAlignment \
            --REFERENCE_SEQUENCE ${reference} \
            --UNMAPPED_BAM unmapped.bam \ # <-- requiered !!!
            --ALIGNED_BAM Aligned.sortedByCoord.out.bam \
            --OUTPUT MergeBamAlignment.bam \
            --INCLUDE_SECONDARY_ALIGNMENTS false \
            --VALIDATION_STRINGENCY SILENT

#-limitBAMsortRAM ${star_mem+"000000000"}


# --outSAMmapqUnique 255
# --outSAMunmapped Within \
# --outSAMattributes Standard 

# STAR --runMode alignReads --readFilesIn BLA_C1_1P.fq,LOL_C3_1P.fq BLA_C1_2P.fq,LOL_C3_2P.fq

# slurm


star_mem=100
thread_count=20 # $SLURM_NPROCS

for i in $(ls *1P.fq);
do
base_name="${i%_1P.fq}"

left_file=${base_name}_1P.fq
right_file=${base_name}_2P.fq

STAR --runMode alignReads \
--genomeDir STAR_index \
--runThreadN $thread_count \
--readFilesIn $left_file $right_file \
--outFileNamePrefix ${base_name}. \
--outSAMunmapped Within \
--outSAMtype BAM SortedByCoordinate \
--twopassMode Basic \
--limitBAMsortRAM ${star_mem}"000000000"


SORTEDBAM=${base_name}.Aligned.sortedByCoord.out.bam

samtools view -b -f 4 -@ $thread_count $SORTEDBAM > ${SORTEDBAM%.bam}.al.bam
samtools view -b -F 4 -@ $thread_count $SORTEDBAM > ${SORTEDBAM%.bam}.uT.bam

done

gatk MergeBamAlignment \
            --REFERENCE_SEQUENCE ${reference} \
            --UNMAPPED_BAM ${unaligned_bam} \
            --ALIGNED_BAM ${star_bam} \
            --OUTPUT ${base_name}.bam \
            --INCLUDE_SECONDARY_ALIGNMENTS false \
            --VALIDATION_STRINGENCY SILENT
 

gatk MarkDuplicates \
 	        --INPUT ${input_bam} \
 	        --OUTPUT ${base_name}.bam  \
 	        --CREATE_INDEX true \
 	        --VALIDATION_STRINGENCY SILENT \
 	        --METRICS_FILE ${base_name}.metrics
```
# 2) Mark Duplicates
Tools involved: MarkDuplicatesSpark / MarkDuplicates + SortSam
```BASH

```

# Base (Quality Score) Recalibration
Tools involved: BaseRecalibrator, Apply Recalibration, AnalyzeCovariates (optional)

De Wit & Palumbi (2012) omit the Base Quality Score Recalibration step for red abalone transcriptome, as this step requires known variant sites as input. 


# 3)
```BASH
gatk AnyTool toolArgs
```