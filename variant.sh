#!/bin/sh
## Directivas
#SBATCH --job-name=gatk
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20
#SBATCH --propagate=STACK

# In addition to --propagate, 
# increase the max number of open files with Linux command ulimit before running STAR
ulimit -s unlimited

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/STAR-2.7.11b/bin/Linux_x86_64_static/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/RSEM/bin/samtools-1.3/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/gatk-4.6.0.0
export PATH=$PATH:$EXPORT


star_mem=100 # gb
thread_count=$SLURM_NPROCS
reference=$1

# Alignment for RNAseq using STAR

mkdir -p STAR_index

if [ ! -f "STAR_index/Genome" ]; then
    STAR --runMode genomeGenerate --genomeDir STAR_index --genomeFastaFiles $reference --limitGenomeGenerateRAM 45284502453 --runThreadN $thread_count --genomeSAindexNbases 12
else
    echo "index for '$reference' already exists."
fi


for i in $(ls *1P.fq);
do
base_name="${i%_1P.fq}"

left_file=${base_name}_1P.fq
right_file=${base_name}_2P.fq

# echo ${left_file} ${right_file} $base_name | column -t > ${base_name}_manifest.tsv

STAR --runMode alignReads \
--genomeDir STAR_index \
--runThreadN $thread_count \
--readFilesIn $left_file $right_file \
--outFileNamePrefix ${base_name}. \
--outSAMunmapped Within \
--outSAMtype BAM SortedByCoordinate \
--twopassMode Basic \
--limitBAMsortRAM ${star_mem}"000000000" \
--outSAMattrRGline ID:${base_name}
# --outSAMmapqUnique 1

# use outSAMattrRGline or readFilesManifes to avoid read group issue in MergeBamAlignment step
# --readFilesManifest ${base_name}_manifest.tsv

# after alignment w/ STAR, error w/ processing bam for MergeBamAlignment <----

SORTEDBAM=${base_name}.Aligned.sortedByCoord.out.bam

samtools view -b -f 4 -@ $thread_count $SORTEDBAM | samtools sort -@ $thread_count -o ${SORTEDBAM%.bam}.al.bam
samtools view -b -F 4 -@ $thread_count $SORTEDBAM | samtools sort -@ $thread_count -o ${SORTEDBAM%.bam}.uT.bam

# samtools view -b -f 4 -@ $thread_count $SORTEDBAM > ${SORTEDBAM%.bam}.al.bam
# samtools view -b -F 4 -@ $thread_count $SORTEDBAM > ${SORTEDBAM%.bam}.uT.bam


gatk SortSam --I ${SORTEDBAM%.bam}.al.bam  --O ${SORTEDBAM%.bam}.al.sort.bam --SO queryname
gatk SortSam --I ${SORTEDBAM%.bam}.uT.bam  --O ${SORTEDBAM%.bam}.uT.sort.bam --SO queryname


# Step 2, pre-process bam for gatk


if [ ! -f "${reference%.*}.dict" ]; then
    gatk CreateSequenceDictionary -R ${reference}
else
    echo "gatk dir for '$reference' already exists."
fi

gatk MergeBamAlignment \
            --REFERENCE_SEQUENCE ${reference} \
            --ALIGNED_BAM ${SORTEDBAM%.bam}.al.bam \
            --UNMAPPED_BAM ${SORTEDBAM%.bam}.uT.bam \
            --SORT_ORDER coordinate \
            --OUTPUT ${base_name}.merged.bam \
            --INCLUDE_SECONDARY_ALIGNMENTS false \
            --VALIDATION_STRINGENCY SILENT
 
            
done

exit

gatk MarkDuplicates \
 	        --INPUT ${base_name}.merged.bam \
 	        --OUTPUT ${base_name}.merged.dup.bam  \
 	        --CREATE_INDEX true \
 	        --VALIDATION_STRINGENCY SILENT \
 	        --METRICS_FILE ${base_name}.metrics