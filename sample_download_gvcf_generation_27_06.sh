#!/usr/bin/bash

#SBATCH --account=project0039        # project ID
#SBATCH --job-name=schisto_gvcf      # descriptive job name
#SBATCH --partition=nodes            # CPU-only partition (node[02-10] idle)
#SBATCH --nodes=1                    # one physical node
#SBATCH --ntasks=1                   # one Slurm task
#SBATCH --cpus-per-task=16           # 16 logical cores for Spark / bwa / etc.
#SBATCH --mem=64G                    # ample RAM for 16 threads
#SBATCH --time=3-00:00:00            # 3 days (well under 7-day limit)
#SBATCH --output=%x.out              # stdout file:  <jobname>.out
#SBATCH --error=%x.err               # stderr file:  <jobname>.err

set -euo pipefail

###############################################################################
# TIMING HELPERS  (fixed)
###############################################################################
hms() {                                  # convert seconds → HH:MM:SS
    printf '%02d:%02d:%02d' $(( $1/3600 )) $(( ($1/60)%60 )) $(( $1%60 ))
}

step_begin() {                           # $1 = step-id   $2 = label
    printf '\n[%s] BEGIN %s\n' "$1" "$2" >&2
    date +%s
}

step_end() {                             # $1 = step-id   $2 = label   $3 = start-ts
    local elapsed=$(( $(date +%s) - $3 ))
    printf '[%s] END %s  (%s)\n' "$1" "$2" "$(hms $elapsed)" >&2
}

run_start=$(date +%s)

###############################################################################
# CONFIG
###############################################################################
REF_GZ="ref_files/schistosoma_mansoni.PRJEA36577.WBPS19.genomic.fa.gz"
WINDOW_SIZE=10000
CONTAINERS="../containers"
SAMPLES="samples.txt"
THREADS=${SLURM_CPUS_PER_TASK:-8}
export OMP_NUM_THREADS=$THREADS

###############################################################################
# ENV & DIRS
###############################################################################
module load apps/apptainer/1.3.4
mkdir -p reads work bam stats depth depth_stats windows window_cov gvcf tmp
export TMPDIR=$PWD/tmp

###############################################################################
# 0) Decompress reference
###############################################################################
step_t0=$(step_begin "0/5" "decompress-reference")
REF="${REF_GZ%.gz}"
if [[ -s $REF ]]; then
    echo "[0/5] ⏩ $REF already present"
else
    echo "[0/5] gunzip → $REF"
    gunzip -c "$REF_GZ" > "$REF"
fi
step_end "0/5" "decompress-reference" "$step_t0"

###############################################################################
# 1) Download reads
###############################################################################
step_t1=$(step_begin "1/5" "download-reads")
while read -r SAMP; do
    [[ -z $SAMP ]] && continue
    R1="reads/${SAMP}_1.fastq.gz"
    R2="reads/${SAMP}_2.fastq.gz"

    if [[ -s $R1 && -s $R2 ]]; then
        echo "SKIP: $SAMP FASTQs already present"
        continue
    fi

    echo "FETCH: $SAMP"
    curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SAMP}&result=read_run&fields=fastq_ftp&format=tsv" \
      | awk 'NR>1 {print $NF}' \
      | tr ';' '\n' \
      | while read -r URL; do
            [[ -z $URL ]] && continue
            # prefer HTTPS if given FTP
            if [[ $URL == ftp://* ]]; then
                URL="https://${URL#ftp://}"
            fi
            OUT="reads/$(basename "$URL")"
            if [[ -s $OUT ]]; then
                echo "  SKIP: $OUT"
            else
                # use wget with resume and retries
                wget --continue --tries=5 --timeout=30 --waitretry=10 -O "$OUT" "$URL"
                if [[ $? -ne 0 ]]; then
                    echo "  ERROR: failed to download $URL" >&2
                    exit 1
                fi
            fi
        done
done < "$SAMPLES"
step_end "1/5" "download-reads" "$step_t1"

###############################################################################
# 2) Index reference + bwa
###############################################################################
step_t2=$(step_begin "2/5" "index-reference")
DICT="${REF%.fa}.dict"
need_index=false
for f in "$REF"{.fai,.amb,.ann,.bwt,.pac,.sa}; do
    [[ ! -s $f ]] && need_index=true
done
if $need_index; then
    echo "[2/5] Building samtools .fai and bwa index"
    apptainer exec "${CONTAINERS}/samtools_1.9.simg" samtools faidx "$REF"
    apptainer exec "${CONTAINERS}/bwa.0.7.17.simg" bwa index "$REF"
else
    echo "[2/5] SKIP: samtools/bwa indices already present"
fi
if [[ ! -s $DICT ]]; then
    echo "[2/5] Creating sequence dictionary ($DICT)"
    apptainer exec "${CONTAINERS}/samtools_1.9.simg" samtools dict "$REF" -o "$DICT"
else
    echo "[2/5] SKIP: $DICT already exists"
fi
step_end "2/5" "index-reference" "$step_t2"

###############################################################################
# 3) Build fixed windows
###############################################################################
step_t3=$(step_begin "3/5" "windows")
WIN_BED="windows/${WINDOW_SIZE}bp_windows.bed"
if [[ -s $WIN_BED ]]; then
    echo "[3/5] ⏩ $WIN_BED exists"
else
    echo "[3/5] bedtools makewindows"
    apptainer exec "${CONTAINERS}/bedtools.2.31.simg" \
        bedtools makewindows -g "${REF}.fai" -w "$WINDOW_SIZE" > "$WIN_BED"
fi
step_end "3/5" "windows" "$step_t3"

###############################################################################
# 4) Per-sample pipeline
###############################################################################
step_t4=$(step_begin "4/5" "per-sample-processing")
: > gvcf_list.txt
while read -r SAMP; do
    [[ -z $SAMP ]] && continue
    BAM=bam/${SAMP}.bam
    GVCF=gvcf/${SAMP}.g.vcf

    echo -e "\n===== SAMPLE: $SAMP ====="

    # a) align & sort
    if [[ -s work/${SAMP}.sorted.bam ]]; then
        echo "  SKIP a) align+sort"
    else
        echo "  RUN  a) align+sort"
        apptainer exec "${CONTAINERS}/bwa.0.7.17.simg" \
            bwa mem -t "$THREADS" "$REF" reads/${SAMP}_1.fastq.gz reads/${SAMP}_2.fastq.gz \
          | apptainer exec "${CONTAINERS}/samtools_1.9.simg" \
            samtools sort -@ "$THREADS" -o work/${SAMP}.sorted.bam -
    fi

    # b) mark duplicates
    if [[ -s work/${SAMP}.dedup.bam ]]; then
        echo "  SKIP b) markdup"
    else
        echo "  RUN  b) markdup"
        apptainer exec "${CONTAINERS}/gatk_4.1.3.0.simg" gatk MarkDuplicates \
            --INPUT work/${SAMP}.sorted.bam \
            --OUTPUT work/${SAMP}.dedup.bam \
            --METRICS_FILE work/${SAMP}.markdup.metrics.txt
    fi

    # c) add read groups
    if [[ -s $BAM ]]; then
        echo "  SKIP c) addRG"
    else
        echo "  RUN  c) addRG"
        apptainer exec "${CONTAINERS}/picard_2.27.5.simg" \
          java -jar /usr/picard/picard.jar AddOrReplaceReadGroups \
          I=work/${SAMP}.dedup.bam O=$BAM RGID=$SAMP RGLB=l1 RGPL=ILLUMINA RGPU=unit1 RGSM=$SAMP
    fi

    # d) index & QC
    if [[ -s ${BAM}.bai ]]; then
        echo "  SKIP d) index"
    else
        echo "  RUN  d) index"
        apptainer exec "${CONTAINERS}/samtools_1.9.simg" samtools index $BAM
    fi
    if [[ -s stats/${SAMP}.flagstat && -s stats/${SAMP}.idxstats ]]; then
        echo "  SKIP d) QC stats"
    else
        echo "  RUN  d) QC stats"
        apptainer exec "${CONTAINERS}/samtools_1.9.simg" samtools flagstat $BAM > stats/${SAMP}.flagstat
        apptainer exec "${CONTAINERS}/samtools_1.9.simg" samtools idxstats $BAM > stats/${SAMP}.idxstats
    fi

    # e) depth & coverage
    if [[ -s depth/${SAMP}.depth ]]; then
        echo "  SKIP e) depth"
    else
        echo "  RUN  e) depth"
        apptainer exec "${CONTAINERS}/samtools_1.9.simg" \
            samtools depth -aa "$BAM" > depth/${SAMP}.depth
    fi

    # f) window coverage
    if [[ -s window_cov/${SAMP}.window.cov ]]; then
        echo "  SKIP f) win_cov"
    else
        echo "  RUN  f) win_cov"
        apptainer exec "${CONTAINERS}/bedtools.2.31.simg" \
            bedtools coverage -sorted -a "$WIN_BED" -b $BAM > window_cov/${SAMP}.window.cov
    fi

    # g) haplotype-caller → GVCF (full‐genome blocks)
    if [[ -s $GVCF ]]; then
        echo "  SKIP g) hapcaller"
    else
        echo "  RUN  g) hapcaller"
        apptainer exec "${CONTAINERS}/gatk_4.1.3.0.simg" gatk HaplotypeCaller \
            --native-pair-hmm-threads "$THREADS" \
            --emit-ref-confidence GVCF \
            --output-mode EMIT_ALL_SITES \
            -I "$BAM" \
            -R "$REF" \
            -O "$GVCF"
    fi

    echo "$GVCF" >> gvcf_list.txt
    echo "===== SAMPLE: $SAMP done ====="
done < "$SAMPLES"
step_end "4/5" "per-sample-processing" "$step_t4"

###############################################################################
# 5) Merge per-sample GVCFs
###############################################################################
step_t5=$(step_begin "5/5" "merge-gvcfs")
COMBINED=combined.g.vcf.gz
if [[ -s $COMBINED ]]; then
    echo "[5/5] ⏩ $COMBINED exists"
else
    apptainer exec "${CONTAINERS}/gatk_4.1.3.0.simg" gatk CombineGVCFs \
        -R "$REF" \
        $(awk '{ printf "--variant %s ", $1 }' gvcf_list.txt) \
        -O "$COMBINED"
fi
step_end "5/5" "merge-gvcfs" "$step_t5"

###############################################################################
# 6) Joint‐genotype gVCF → final cohort VCF (incl. non-variant sites)
###############################################################################
echo -e "\n[6/6] BEGIN genotype-gvcfs"
step_t6=$(step_begin "6/6" "genotype-gvcfs")
COHORT_VCF=cohort.fullmatrix.vcf.gz
if [[ ! -s $COHORT_VCF ]]; then
    apptainer exec "${CONTAINERS}/gatk_4.1.3.0.simg" \
        gatk --java-options "-Xmx60G" GenotypeGVCFs \
        -R "$REF" \
        -V "$COMBINED" \
        --include-non-variant-sites \
        -O "$COHORT_VCF"
else
    echo "[6/6] ⏩ $COHORT_VCF exists"
fi
step_end "6/6" "genotype-gvcfs" "$step_t6"

###############################################################################
# WRAP-UP
###############################################################################
run_elapsed=$(( $(date +%s) - run_start ))
printf '\n[✔] Pipeline complete ➟ %s (total runtime %s)\n' \
       "$COHORT_VCF" "$(hms $run_elapsed)"
