#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  01 - Quality control and filtering of caribou SNP data
#
#  Companion code for:
#    Bourbon et al. (2026) Inbreeding and Demographic History of Caribou
#    (Rangifer tarandus) in Western Canada Inferred From Genome-Wide SNP Data.
#    Evolutionary Applications. doi:10.1111/eva.70311
#
#  This script leads to the merged, chromosome-aligned SNP dataset (available on
#  Dryad: doi:10.5061/dryad.8pk0p2p3x) and applies the quality-control pipeline
#  described in the Methods, producing the final dataset of 33,346 autosomal
#  SNPs for 759 individuals used in all downstream analyses.
#
#  No MAF or LD pruning is applied. Following Meyermans et al. (2020), such
#  filtering can bias ROH detection and LD-based Ne estimation by removing
#  informative rare variants and reducing the linkage signal.
#
#  REQUIREMENTS
#    Software : PLINK v1.9 and PLINK2 on the system PATH (called via system()).
#    R        : dplyr
#
#  INPUTS (place in the directory set by `data_dir`)
#    caribou_snps.bed/.bim/.fam  Merged, chromosome-aligned genotypes (Dryad).
#    no_call_SNPs.txt            Poorly-genotyped markers (Carrier et al. 2022).
#    monomorphic_SNPs.txt        Monomorphic markers (Carrier et al. 2022).
#    sample_herd_ids.csv         Columns: old_FID, old_IID, new_FID, new_IID.
#                                Maps each individual to its subpopulation (herd).
#                                Herd label "JNP_MISC" marks individuals with no
#                                herd assignment, which are removed.
#
#  OUTPUT
#    caribou_759_autosomes.*     Final analysis dataset (BED/PED/VCF),
#                                33,346 SNPs x 759 individuals, autosomes only.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
data_dir <- "data"      # input files (see INPUTS above)
out_dir  <- "output"    # results are written here

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(data_dir)

library(dplyr)

# Write outputs to out_dir while keeping the input directory clean
op <- function(f) file.path("..", out_dir, f)

# The caribou assembly has 34 autosomes + 1 sex chromosome (chr 35).
# --chr-set 35 and --allow-extra-chr are required by PLINK for this genome.


# ---------------------------------------------------------------------------
# 1. Remove no-call and monomorphic markers
# ---------------------------------------------------------------------------
# These marker lists were reported by Carrier et al. (2022) during array design.

system("plink --bfile caribou_snps --allow-extra-chr --chr-set 35 --exclude no_call_SNPs.txt --make-bed --out caribou_step1")
system("plink --bfile caribou_step1 --allow-extra-chr --chr-set 35 --exclude monomorphic_SNPs.txt --make-bed --out caribou_step2")


# ---------------------------------------------------------------------------
# 2. Missingness filters
# ---------------------------------------------------------------------------
# Remove individuals (--mind) then markers (--geno) with >5% missing genotypes.

system("plink --bfile caribou_step2 --mind 0.05 --allow-extra-chr --chr-set 35 --make-bed --out caribou_mind")
system("plink --bfile caribou_mind --geno 0.05 --allow-extra-chr --chr-set 35 --make-bed --out caribou_geno")


# ---------------------------------------------------------------------------
# 3. Identify and remove duplicate individuals (IBS/PI_HAT > 0.95)
# ---------------------------------------------------------------------------
system("plink --bfile caribou_geno --genome --allow-extra-chr --chr-set 35 --out caribou_ibs")

genome <- read.table("caribou_ibs.genome", header = TRUE, sep = "")

duplicate_pairs <- genome %>%
  filter(PI_HAT > 0.95) %>%
  select(FID1, IID1, FID2, IID2, PI_HAT) %>%
  arrange(desc(PI_HAT))

cat("Duplicate pairs found (PI_HAT > 0.95):", nrow(duplicate_pairs), "\n")

# Full pair info (kept for the genotyping-error estimate in section 8)
write.table(duplicate_pairs, op("duplicate_pairs.txt"),
            sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)

# One individual per duplicate pair, to be dropped
duplicates_to_remove <- duplicate_pairs %>% select(FID1, IID1)
write.table(duplicates_to_remove, "duplicates_to_remove.txt",
            sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)

system("plink --bfile caribou_geno --remove duplicates_to_remove.txt --allow-extra-chr --chr-set 35 --make-bed --out caribou_nodup")


# ---------------------------------------------------------------------------
# 4. Assign herd (subpopulation) IDs and drop unassigned individuals
# ---------------------------------------------------------------------------
# sample_herd_ids.csv maps individuals to subpopulations. Individuals labelled
# "JNP_MISC" have no herd assignment and are removed.

herd_ids <- read.csv("sample_herd_ids.csv", header = TRUE)

# Individuals with a valid herd, in PLINK --update-ids format:
# old_FID  old_IID  new_FID  new_IID
assigned <- herd_ids[herd_ids$new_FID != "JNP_MISC", ]
write.table(assigned, "herd_id_update.txt",
            sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)

# Individuals to drop (no herd assignment)
unassigned <- herd_ids[herd_ids$new_FID == "JNP_MISC", 1:2]
write.table(unassigned, "unassigned_to_remove.txt",
            sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)

system("plink --bfile caribou_nodup --remove unassigned_to_remove.txt --allow-extra-chr --chr-set 35 --make-bed --out caribou_assigned")
system("plink --bfile caribou_assigned --update-ids herd_id_update.txt --allow-extra-chr --chr-set 35 --make-bed --out caribou_herds")


# ---------------------------------------------------------------------------
# 5. Hardy-Weinberg equilibrium filter (p < 1e-6)
# ---------------------------------------------------------------------------
system("plink --bfile caribou_herds --hwe 1e-6 --allow-extra-chr --chr-set 35 --make-bed --out caribou_hwe")
# Result at this stage: 33,531 SNPs x 759 individuals


# ---------------------------------------------------------------------------
# 6. Final dataset: autosomes only (exclude sex chromosome, chr 35)
# ---------------------------------------------------------------------------
# Sex chromosomes are excluded because they differ from autosomes in
# inheritance, Ne, and recombination rate. This yields the final analysis set.

system("plink --bfile caribou_hwe --chr 1-34 --allow-extra-chr --chr-set 35 --make-bed --out caribou_759_autosomes")
system("plink --bfile caribou_759_autosomes --allow-extra-chr --chr-set 34 --recode --out caribou_759_autosomes")
system("plink --bfile caribou_759_autosomes --allow-extra-chr --chr-set 34 --recode vcf --out caribou_759_autosomes")
# Final analysis dataset: 33,346 SNPs x 759 individuals

# A version with close relatives removed (KING kinship cutoff 0.177), used where
# relatedness could bias an estimate.
system("plink2 --bfile caribou_759_autosomes --allow-extra-chr --chr-set 34 --king-cutoff 0.177 --make-bed --out caribou_759_autosomes_unrelated")


# ---------------------------------------------------------------------------
# 7. Summary statistics on the final dataset
# ---------------------------------------------------------------------------
# Per-individual heterozygosity
system("plink --bfile caribou_759_autosomes --allow-extra-chr --chr-set 34 --het --out caribou_759_het")

# Missingness
system("plink --bfile caribou_759_autosomes --allow-extra-chr --chr-set 34 --missing --out caribou_759_missing")

lmiss <- read.table("caribou_759_missing.lmiss", header = TRUE)
imiss <- read.table("caribou_759_missing.imiss", header = TRUE)
cat("Mean locus missingness:      ", mean(lmiss$F_MISS), "\n")
cat("Mean individual missingness: ", mean(imiss$F_MISS), "\n")


# ---------------------------------------------------------------------------
# 8. Empirical genotyping error rate from duplicate pairs
# ---------------------------------------------------------------------------
# Duplicate (re-genotyped) individuals are compared marker-by-marker to estimate
# an empirical error rate, used to assess ROH robustness to genotyping error.

dups <- read.table(op("duplicate_pairs.txt"), header = FALSE)
colnames(dups) <- c("FID1", "IID1", "FID2", "IID2", "PI_HAT")

total_mismatches  <- 0
total_comparisons <- 0

for (i in seq_len(nrow(dups))) {

  write.table(data.frame(dups$FID1[i], dups$IID1[i]), "tmp_id1.txt",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(data.frame(dups$FID2[i], dups$IID2[i]), "tmp_id2.txt",
              row.names = FALSE, col.names = FALSE, quote = FALSE)

  system("plink --bfile caribou_step2 --keep tmp_id1.txt --allow-extra-chr --chr-set 35 --recode --out tmp_geno1")
  system("plink --bfile caribou_step2 --keep tmp_id2.txt --allow-extra-chr --chr-set 35 --recode --out tmp_geno2")

  ped1 <- read.table("tmp_geno1.ped", header = FALSE)
  ped2 <- read.table("tmp_geno2.ped", header = FALSE)

  geno1 <- as.character(unlist(ped1[1, 7:ncol(ped1)]))
  geno2 <- as.character(unlist(ped2[1, 7:ncol(ped2)]))

  valid       <- geno1 != "0" & geno2 != "0"   # exclude missing calls
  total_mismatches  <- total_mismatches  + sum(geno1[valid] != geno2[valid])
  total_comparisons <- total_comparisons + sum(valid)
}

file.remove(Filter(file.exists,
  c("tmp_id1.txt", "tmp_id2.txt",
    paste0("tmp_geno1", c(".ped", ".map", ".log")),
    paste0("tmp_geno2", c(".ped", ".map", ".log")))))

error_rate <- total_mismatches / total_comparisons
cat(sprintf("Empirical genotyping error rate: %.4f (%.2f%%)\n",
            error_rate, error_rate * 100))
