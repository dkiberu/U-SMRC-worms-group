#!/bin/bash

# Set up variables
INPUT_VCF="/home/dkiberu/combined_fibroschot_raw.all.vcf.gz"
RESULTS_DIR="results_o4"
META_FILE="combined_meta"
CHROM_MAP="chrom_mapping.txt"
SUBSET_SIZE=50000
SMALL_SUBSET_SIZE=10000
START_TIME=$(date +%s)

# Create results directory
mkdir -p $RESULTS_DIR

# Step 1: Filter SNPs and by MAF using bcftools (faster than vcftools)
if [ ! -f "$RESULTS_DIR/maf_filtered.recode.vcf.gz" ]; then
  echo "Step 1: Filtering SNPs by MAF using bcftools..."
  bcftools view -v snps -i 'MAF[0] > 0.05' -Oz -o $RESULTS_DIR/maf_filtered.recode.vcf.gz $INPUT_VCF || { echo "Error in SNP/MAF filtering"; exit 1; }
  bcftools index $RESULTS_DIR/maf_filtered.recode.vcf.gz
  echo "SNP + MAF filtering completed."
else
  echo "Step 1 skipped: MAF filtered file already exists."
fi

# Step 2: Filter by Missingness
if [ ! -f "$RESULTS_DIR/filtered.recode.vcf.gz" ]; then
  echo "Step 2: Filtering by missingness..."
  vcftools --gzvcf $RESULTS_DIR/maf_filtered.recode.vcf.gz --max-missing 0.9 --recode --recode-INFO-all --out $RESULTS_DIR/filtered || { echo "Error in missingness filtering"; exit 1; }
  bgzip $RESULTS_DIR/filtered.recode.vcf
  bcftools index $RESULTS_DIR/filtered.recode.vcf.gz
  echo "Missingness filtering completed."
else
  echo "Step 2 skipped: Missingness filtered file already exists."
fi

# Step 3: LD Pruning
if [ ! -f "$RESULTS_DIR/ld_pruned.vcf.gz" ]; then
  echo "Step 3: LD pruning..."
  bcftools annotate --set-id '%CHROM:%POS' -Oz -o $RESULTS_DIR/filtered_with_ids.vcf.gz $RESULTS_DIR/filtered.recode.vcf.gz
  tabix -p vcf $RESULTS_DIR/filtered_with_ids.vcf.gz

  plink --vcf $RESULTS_DIR/filtered_with_ids.vcf.gz \
        --indep-pairwise 50 5 0.2 \
        --double-id --allow-extra-chr \
        --out $RESULTS_DIR/pruned || { echo "Error in pruning"; exit 1; }

  bcftools view -i "ID=@$RESULTS_DIR/pruned.prune.in" \
    $RESULTS_DIR/filtered_with_ids.vcf.gz \
    -Oz -o $RESULTS_DIR/ld_pruned.vcf.gz
  tabix -p vcf $RESULTS_DIR/ld_pruned.vcf.gz
  echo "LD pruning completed."
else
  echo "Step 3 skipped: LD-pruned file exists."
fi

# Step 4: Subset 50k SNPs
if [ ! -f "$RESULTS_DIR/subset_main_50k.vcf.gz" ]; then
  echo "Step 4: Subsetting 50k SNPs..."
  bcftools query -f '%CHROM\n' $RESULTS_DIR/ld_pruned.vcf.gz | sort | uniq -c > $RESULTS_DIR/variant_counts_per_chrom.txt
  awk -v subset_size=$SUBSET_SIZE '{proportion = $1 / total; snps = int(proportion * subset_size); print $2, snps}' \
    total=$(awk '{sum += $1} END {print sum}' $RESULTS_DIR/variant_counts_per_chrom.txt) \
    $RESULTS_DIR/variant_counts_per_chrom.txt > $RESULTS_DIR/snp_proportions_50k.txt
  while read chrom snps; do
    bcftools view -r $chrom $RESULTS_DIR/ld_pruned.vcf.gz | \
    bcftools query -f '%CHROM\t%POS\n' | \
    shuf -n $snps >> $RESULTS_DIR/subset_snps_main_50k.txt
  done < $RESULTS_DIR/snp_proportions_50k.txt
  bcftools view -T $RESULTS_DIR/subset_snps_main_50k.txt $RESULTS_DIR/ld_pruned.vcf.gz -Oz -o $RESULTS_DIR/subset_main_50k.vcf.gz
  bcftools index $RESULTS_DIR/subset_main_50k.vcf.gz
  echo "50k SNP subsetting completed."
else
  echo "Step 4 skipped: 50k SNP subset already exists."
fi

# Step 5: Subset 10k SNPs
if [ ! -f "$RESULTS_DIR/subset_small_10k.vcf.gz" ]; then
  echo "Step 5: Subsetting 10k SNPs..."
  awk -v subset_size=$SMALL_SUBSET_SIZE '{proportion = $1 / total; snps = int(proportion * subset_size); print $2, snps}' \
    total=$(awk '{sum += $1} END {print sum}' $RESULTS_DIR/variant_counts_per_chrom.txt) \
    $RESULTS_DIR/variant_counts_per_chrom.txt > $RESULTS_DIR/snp_proportions_10k.txt
  while read chrom snps; do
    bcftools view -r $chrom $RESULTS_DIR/ld_pruned.vcf.gz | \
    bcftools query -f '%CHROM\t%POS\n' | \
    shuf -n $snps >> $RESULTS_DIR/subset_snps_small_10k.txt
  done < $RESULTS_DIR/snp_proportions_10k.txt
  bcftools view -T $RESULTS_DIR/subset_snps_small_10k.txt $RESULTS_DIR/ld_pruned.vcf.gz -Oz -o $RESULTS_DIR/subset_small_10k.vcf.gz
  bcftools index $RESULTS_DIR/subset_small_10k.vcf.gz
  echo "10k SNP subsetting completed."
else
  echo "Step 5 skipped: 10k SNP subset already exists."
fi

# Step 6: Calculate Pi
if [ ! -f "$RESULTS_DIR/pi_non_overlapping.windowed.pi" ] || [ ! -f "$RESULTS_DIR/pi_overlapping.windowed.pi" ]; then
  echo "Step 6: Calculating Pi..."
  vcftools --gzvcf $RESULTS_DIR/subset_main_50k.vcf.gz --window-pi 10000 --out $RESULTS_DIR/pi_non_overlapping || { echo "Error in non-overlapping Pi calculation"; exit 1; }
  vcftools --gzvcf $RESULTS_DIR/subset_main_50k.vcf.gz --window-pi 10000 --window-pi-step 5000 --out $RESULTS_DIR/pi_overlapping || { echo "Error in overlapping Pi calculation"; exit 1; }
  echo "Pi calculation completed."
else
  echo "Step 6 skipped: Pi results already exist."
fi

# Step 7: FST Calculation
if [ ! -f "$RESULTS_DIR/fst_villages.weir.fst" ]; then
  echo "Step 7: FST Calculation..."
  awk -F"," '$3 == "Buhirigi" {if (NR > 1) print $1}' $META_FILE > $RESULTS_DIR/buhirigi_samples.txt
  awk -F"," '$3 == "Kaiso" {if (NR > 1) print $1}' $META_FILE > $RESULTS_DIR/kaiso_samples.txt
  vcftools --gzvcf $RESULTS_DIR/subset_main_50k.vcf.gz \
           --weir-fst-pop $RESULTS_DIR/buhirigi_samples.txt \
           --weir-fst-pop $RESULTS_DIR/kaiso_samples.txt \
           --out $RESULTS_DIR/fst_villages || { echo "Error in FST calculation for villages"; exit 1; }
  echo "FST calculation completed."
else
  echo "Step 7 skipped: FST results already exist."
fi

# Step 8: PCA
if [ ! -f "$RESULTS_DIR/pca.eigenvec" ]; then
  echo "Step 8: PCA analysis..."
  plink --vcf $RESULTS_DIR/subset_main_50k.vcf.gz --pca 20 --out $RESULTS_DIR/pca --double-id --allow-extra-chr || { echo "Error in PCA analysis"; exit 1; }
  echo "PCA analysis completed."
else
  echo "Step 8 skipped: PCA results already exist."
fi

# Step 9: KING kinship analysis
if [ ! -f "$RESULTS_DIR/king_kinship.kin0" ]; then
  echo "Step 9: KING kinship analysis..."
  bcftools annotate --rename-chrs $CHROM_MAP -Oz -o $RESULTS_DIR/subset_main_50k_chr_renamed.vcf.gz $RESULTS_DIR/subset_main_50k.vcf.gz
  bcftools index $RESULTS_DIR/subset_main_50k_chr_renamed.vcf.gz
  plink --vcf $RESULTS_DIR/subset_main_50k_chr_renamed.vcf.gz --double-id --allow-extra-chr --make-bed --out $RESULTS_DIR/king_data || { echo "Error in KING preparation"; exit 1; }
  king -b $RESULTS_DIR/king_data.bed --kinship --prefix $RESULTS_DIR/king_kinship || { echo "Error in KING kinship analysis"; exit 1; }
  echo "KING analysis completed."
else
  echo "Step 9 skipped: KING results already exist."
fi

# Step 10: ngsRelate
if [ ! -f "$RESULTS_DIR/ngsrelate_results.res" ]; then
  echo "Step 10: Running ngsRelate..."
  ngsRelate -h $RESULTS_DIR/subset_main_50k.vcf.gz -O $RESULTS_DIR/ngsrelate_results.res || { echo "Error in ngsRelate"; exit 1; }
  echo "ngsRelate completed."
else
  echo "Step 10 skipped: ngsRelate results already exist."
fi

# Step 11: Phylogenetic Tree Construction
if [ ! -f "$RESULTS_DIR/phylo.treefile" ]; then
  echo "Step 11: Constructing phylogenetic tree..."
  vcf2phylip.py -i $RESULTS_DIR/subset_small_10k.vcf.gz -f || { echo "Error in VCF to phylip conversion"; exit 1; }
  iqtree2 -s subset_small_10k.min4.fasta -st DNA -m GTR+G -nt AUTO --alrt 1000 -B 1000 --prefix $RESULTS_DIR/phylo || { echo "Error in tree construction"; exit 1; }
  echo "Phylogenetic tree construction completed."
else
  echo "Step 11 skipped: Tree already exists."
fi

# Step 12c: ADMIXTURE for K=2,3,4
for K in 2 3 4; do
  if [ ! -f "$RESULTS_DIR/admixture_K${K}.Q" ]; then
    echo "Step 12c: Performing ADMIXTURE analysis for K=${K}..."
    admixture $RESULTS_DIR/king_data.bed $K > $RESULTS_DIR/admixture_K${K}.log || { echo "Error in ADMIXTURE analysis for K=${K}"; exit 1; }
    echo "ADMIXTURE analysis for K=${K} completed."
  else
    echo "Step 12c skipped: ADMIXTURE results for K=${K} already exist."
  fi
done

# Timing
END_TIME=$(date +%s)
ELAPSED_TIME=$((END_TIME - START_TIME))
ELAPSED_MINUTES=$((ELAPSED_TIME / 60))
echo "Pipeline completed in $ELAPSED_MINUTES minutes."
