# Quantification step
## Using rsem

```bash
# 0) LOAD TOOLS
# git clone https://github.com/deweylab/RSEM.git
# make
# make install

EXPORT=/LUSTRE/apps/bioinformatica/RSEM/bin/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/bowtie2/bin/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/RSEM/bin/samtools-1.3/
export PATH=$PATH:$EXPORT
```

## Build index
```bash
#!/bin/sh
## Directivas
#SBATCH --job-name=bbuild
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20

# 1) INDEX THE REFERENCE

reference=$1 # evigene_transcripts.cdna
db_index_name=idx # ${reference%.cdna}

# Build Bowtie2 index if not already present
if [ ! -f "$index_name.1.bt2" ]; then
    bowtie2-build --threads $SLURM_NPROCS $reference $db_index_name
else
    echo "index '$index_name' already exists."
fi

rsem-prepare-reference $reference $db_index_name

exit
```

## Bowtie alignment

sbatch rsem.sh transcripts.fasta

```bash
#!/bin/sh
## Directivas
#SBATCH --job-name=rsem
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20


# 1.1) BOWTIE PARAMETERS

aligner_params="--no-mixed --no-discordant --gbar 1000 --end-to-end -k 200"
read_type="-q"
thread_count=$SLURM_NPROCS
max_ins_size=800

reference=$1
db_index_name=idx


# 2) ALIGN READS BACK TO THE REFERENCE
#if [ ! -f "${base}.sorted.bam" ]; then
#    left_file=${base}_1.P.qtrim.fq.gz
#    right_file=${base}_2.P.qtrim.fq.gz
#else
#    echo "Alignment for '${base}.' samples already exists."
#fi

for i in $(ls *1.P.qtrim.fq.gz);
do
base="${i%_1.P.qtrim.fq.gz}"
left_file=${base}_1.P.qtrim.fq.gz
right_file=${base}_2.P.qtrim.fq.gz

bam_file=${base}.sorted.bam
met_file=${base}.met.txt

bowtie2 --met-file $met_file $aligner_params $read_type -X $max_ins_size -x $db_index_name -1 $left_file -2 $right_file -p $thread_count | samtools view -@ $thread_count -F 4 -S -b | samtools sort -@ $thread_count -n -o $bam_file

done

# 2.1) GENERATE RSEM FILES:

for i in $(ls *.sorted.bam);
do
bam_for_rsem=${i%.bam}.rsem
convert-sam-for-rsem -p $thread_count $i $bam_for_rsem
done

# 3) RUN ESTIMATION ABUNDANCE
# 3.1) RSEM PARAM

fragment_length=200
fragment_std=80
fraglength_info_txt="--fragment-length-mean $fragment_length --fragment-length-sd $fragment_std"

paired_flag_text="--paired-end"  
no_qualities_string=""
keep_intermediate_files_opt="--keep-intermediate-files"

SS_opt="--forward-prob 1.0"
rsem_bam_flag="--no-bam-output"


rsem_prefix=idx


# ESTIMATE ABUNDANCE:

for i in $(ls *.rsem.bam);
do
output_prefix=${i%.sorted.rsem.bam}

rsem-calculate-expression $no_qualities_string $paired_flag_text -p $thread_count $fraglength_info_txt $keep_intermediate_files_opt $SS_opt $rsem_bam_flag --bam $i $rsem_prefix $output_prefix

done

exit


# 4) OUTPUT COUNT MATRIX TO DIFFERENTIAL EXPRESSION ANALYSIS\

ls -ltrh *isoforms.results
rsem-generate-data-matrix *isoforms.results > isoforms.counts.matrix

```

## hisat alignment 
```bash
for i in $(ls *1.P.qtrim.fq.gz); do sbatch hisat_align.sh transcripts.fasta $i; done
```