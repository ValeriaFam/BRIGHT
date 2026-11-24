#!/bin/bash

#PBS -S /bin/sh
#PBS -N Modkit
#PBS -l select=1:ncpus=20:mem=10G
#PBS -l walltime=120:00:00
#PBS -q workq

cd /projects/CGS_shared/vfama/MODIFICATIONS/MODKIT/NICASSIO

source /home/vfama/miniconda3/bin/activate /home/vfama/miniconda3/envs/MINIMAP2
# Percorsi
MODKIT="/projects/CGS_shared/vfama/MODIFICATIONS/MODKIT/dist_modkit_v0.6.1_481e3c9/modkit"
TRANSCRIPTOME_FA="/projects/CGS_shared/vfama/MODIFICATIONS/MODKIT/gencode.v48.transcripts.fa"
OUTDIR="/projects/CGS_shared/vfama/MODIFICATIONS/MODKIT/NICASSIO/"
THREADS=24

mkdir -p ${OUTDIR}/{logs,pileup,histograms}

# Loop su tutti i BAM nella cartella corrente
for bam in IVT.bam; do
    filename=$(basename "$bam")
    prefix=${filename%.bam}  # rimuove solo .bam

    #/projects/CGS_shared/vfama/MODIFICATIONS/MODKIT/dorado-1.3.0-linux-x64/bin/dorado aligner --mm2-opts "-x map-ont" --threads $THREADS $TRANSCRIPTOME_FA "$bam" | samtools view -hb -F 2324| samtools sort -@ $THREADS -o "${prefix}_aligned.bam"

    #echo "Indexing ${prefix}_aligned.bam..."
    #samtools index "${prefix}_aligned.bam"

    #BAM="${prefix}_aligned.bam"  # BAM da passare a Modkit
    BAM="$bam"

    # --- MODKIT ---
    echo "Step 1: summary - $prefix"
    "$MODKIT" summary "$BAM" --threads "$THREADS" \
        > "$OUTDIR/logs/${prefix}_summary.txt"

   echo "Step 3: pileup m6A all contexts - $prefix"
    "$MODKIT" pileup "$BAM" --mod-thresholds a:0.98 \
        "$OUTDIR/pileup/${prefix}_m6A_all.bed" \
        --ref "$TRANSCRIPTOME_FA" \
        --threads "$THREADS" \
        --log-filepath "$OUTDIR/logs/${prefix}_m6A_all.log" \
        --high-depth

done
