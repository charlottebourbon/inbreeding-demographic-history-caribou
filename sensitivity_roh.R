#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  ROH robustness to genotyping error (parameter sensitivity)
#
#  Companion code for:
#    Bourbon et al. (2026) Inbreeding and Demographic History of Caribou
#    (Rangifer tarandus) in Western Canada Inferred From Genome-Wide SNP Data.
#    Evolutionary Applications. doi:10.1111/eva.70311
#
#  Monte Carlo assessment of how ROH detection under the consecutive method
#  tolerates genotyping error (Methods 2.3; Supporting Information S1). For each
#  combination of detection parameters (minSNP, maxOppRun, maxMissRun), a run of
#  n_snp SNPs is simulated many times with per-SNP error at the observed array
#  rate (~1.3%; array reproducibility 98.7%, Carrier et al. 2022). A run is
#  "recovered" if its number of erroneous calls does not exceed the combined
#  heterozygous + missing tolerance. Recovery is the fraction of runs retained.
#
#  Reports the recovery rate for the parameter set used in the paper:
#  minSNP = 20, maxOppRun = 2, maxMissRun = 2.
#
#  Base R only.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

set.seed(1)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
out_dir <- "output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

error_rate <- 0.013     # per-SNP genotyping error (1 - 0.987 reproducibility)
n_snp      <- 200       # SNPs in the simulated run
n_sim      <- 1000      # simulations per parameter combination

minSNP_values  <- c(20, 30, 40, 50)
maxOpp_values  <- 0:4   # heterozygous-SNP tolerance
maxMiss_values <- 0:2   # missing-genotype tolerance

# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------
# One simulated run: n_snp Bernoulli draws, each an erroneous (non-homozygous)
# call with probability error_rate.
simulate_roh <- function(n_snp, error_rate) {
  rbinom(n_snp, 1, error_rate)
}

# A run is retained if it meets the minimum length and its error count does not
# exceed the combined opposite + missing tolerance.
is_recovered <- function(roh, min_snp, max_opp, max_miss) {
  if (length(roh) < min_snp) return(FALSE)
  sum(roh) <= (max_opp + max_miss)
}

results <- expand.grid(
  minSNP  = minSNP_values,
  maxOpp  = maxOpp_values,
  maxMiss = maxMiss_values
)
results$recovery <- NA_real_

for (i in seq_len(nrow(results))) {
  keep <- replicate(n_sim, {
    roh <- simulate_roh(n_snp, error_rate)
    is_recovered(roh, results$minSNP[i], results$maxOpp[i], results$maxMiss[i])
  })
  results$recovery[i] <- mean(keep)
}

# ---------------------------------------------------------------------------
# Report and save
# ---------------------------------------------------------------------------
chosen <- subset(results, minSNP == 20 & maxOpp == 2 & maxMiss == 2)
cat(sprintf("Recovery at chosen parameters (minSNP=20, maxOppRun=2, maxMissRun=2): %.3f\n",
            chosen$recovery))

write.table(results, file.path(out_dir, "roh_error_sensitivity.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
