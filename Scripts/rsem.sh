#!/bin/sh
## Directivas
#SBATCH --job-name=quant
#SBATCH -N 2
#SBATCH --mem=70GB
#SBATCH --ntasks-per-node=20

EXPORT=/LUSTRE/apps/bioinformatica/bowtie2/bin/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/RSEM/bin/samtools-1.3/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/RSEM/bin/
export PATH=$PATH:$EXPORT

# 1) INDEX THE REFERENCE

REFERENCE=$1 # reference.fasta
ref_index_name=idx # ${reference%.cdna}

thread_count=$SLURM_NPROCS

mkdir -p INDEX

# Build Bowtie2 index if not already present
if [ ! -f "INDEX/$ref_index_name.1.bt2" ]; then
    bowtie2-build --threads $thread_count $REFERENCE INDEX/$ref_index_name
else
    echo "index '$ref_index_name' already exists."
fi

echo "Continue with step two"

# Parameters for alignment

aligner_params="--no-mixed --no-discordant --gbar 1000 --end-to-end -k 200"
read_type="-q"
max_ins_size=800

mkdir -p S1_BOWTIE2_BAM_FILES


for i in $(ls *_1.P.qtrim.fq.gz)
do
bs="${i%*_1.P.qtrim.fq.gz}"

left_file=${bs}_1.P.qtrim.fq.gz
right_file=${bs}_2.P.qtrim.fq.gz

bam_file=${bs}.sorted.bam
met_file=${bs}.met.txt

# # Test if the alignment was previously done!

if [ ! -f "S1_BOWTIE2_BAM_FILES/$bam_file" ]; then
    
    echo "Aligning reads back to reference"

    bowtie2 --met-file $met_file $aligner_params $read_type -X $max_ins_size \
    -x INDEX/$ref_index_name -1 $left_file -2 $right_file -p $thread_count 2> ${bs}.stderr | \
    samtools view -@ $thread_count -F 4 -S -b | samtools sort -@ $thread_count -n -o S1_BOWTIE2_BAM_FILES/$bam_file
else
    echo "File $bam_file already exists in S1_BOWTIE2_BAM_FILES directory"

fi

unlink $left_file
unlink $right_file

done

mkdir -p stats
mv *.met.txt stats
mv *.stderr stats
mv stats S1_BOWTIE2_BAM_FILES

# grep 'overall alignment rate' S1_BOWTIE2_BAM_FILES/stats/*.stderr

echo "Continue with step three"

rsem_prefix=idx # conveniently as $ref_index_name

if [ ! -f "$rsem_prefix.seq" ]; then
    rsem-prepare-reference $REFERENCE $rsem_prefix
else
    echo "index '$rsem_prefix' already exists."
fi

# 2.1) Convert bam to rsem 

mkdir -p S2_RSEM_CALCULATION_FILES

for i in $(ls S1_BOWTIE2_BAM_FILES/*.sorted.bam);
do

bs=`basename ${i%.bam}`
bam_for_rsem=${bs}.rsem
ls
if [ ! -f "S2_RSEM_CALCULATION_FILES/$bs.rsem.bam" ]; then
    
    echo "convert-sam-for-rsem fo file $i"

   convert-sam-for-rsem -p $thread_count $i S2_RSEM_CALCULATION_FILES/$bam_for_rsem
else
    echo "File $bam_for_rsem already exists in S2_RSEM_CALCULATION_FILES directory"
fi
done

fragment_length=200
fragment_std=80
fraglength_info_txt="--fragment-length-mean $fragment_length --fragment-length-sd $fragment_std"

paired_flag_text="--paired-end"  
no_qualities_string=""
keep_intermediate_files_opt="--keep-intermediate-files"

SS_opt="--forward-prob 1.0"
rsem_bam_flag="--no-bam-output"

# Estimate abundance:
echo "Continue with final step"

for i in $(ls S2_RSEM_CALCULATION_FILES/*.rsem.bam);
do
output_prefix=${i%.sorted.rsem.bam}
rsem-calculate-expression $no_qualities_string $paired_flag_text -p $thread_count $fraglength_info_txt $keep_intermediate_files_opt $SS_opt $rsem_bam_flag --bam $i $rsem_prefix $output_prefix

done

# OUTPUT COUNT MATRIX TO DIFFERENTIAL EXPRESSION ANALYSIS\

rsem-generate-data-matrix S2_RSEM_CALCULATION_FILES/*.isoforms.results > ${ref_index_name}_isoforms.matrix
rsem-generate-data-matrix S2_RSEM_CALCULATION_FILES/*.genes.results > ${ref_index_name}_genes.matrix
#


for i in $(ls *.sorted.bam)
do
samtools flagstat $i > ${i%.sorted.bam}.flagstats.txt
done

exit

