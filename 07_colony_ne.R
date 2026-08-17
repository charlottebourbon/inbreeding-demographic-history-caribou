#!/bin/bash
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=1-00:00:00
#SBATCH --job-name=colony_Ne
#SBATCH --output=logs/colony_%j.log
#
#  07 - Metapopulation Ne cross-check with COLONY (sibship reconstruction)
#
#  Companion code for:
#    Bourbon et al. (2026) Inbreeding and Demographic History of Caribou
#    (Rangifer tarandus) in Western Canada Inferred From Genome-Wide SNP Data.
#    Evolutionary Applications. doi:10.1111/eva.70311
#
#  Structure-insensitive, LD-independent cross-check on the currentNe2 estimates
#  (Methods 2.5; Supporting Information S1). For each metapopulation:
#    1. Select the 500 most informative SNPs: MAF > 0.4, low missingness, and
#       LD-pruned to r^2 < 0.1 (PLINK) for marker independence.
#    2. Format a COLONY input (.dat) implementing the Full Likelihood method,
#       polygamous mating for both sexes, and non-random mating (inbreeding
#       accounted for), with conservative error (0.01) and missing (0.02) rates.
#    3. Run COLONY and record the sibship-based Ne estimate.
#
#  Requires: PLINK v1.9, Python 3, COLONY v2.0.7.2. Set paths in CONFIGURATION.
#
#  INPUT   : data/metapopulations/<name>.ped/.map  (one per metapopulation)
#  OUTPUT  : results/colony/<name>/<name>_out.Ne   (COLONY Ne output)
#            results/colony/colony_Ne_summary.txt  (collated Ne estimates)
# ===========================================================================

set -uo pipefail

# ---------------------------------------------------------------------------
# CONFIGURATION  -- edit for your environment
# ---------------------------------------------------------------------------
BASE="."
PLINK="plink"
COLONY="colony2s.gnu.out"          # COLONY executable

INDIR="${BASE}/data/metapopulations"
OUTBASE="${BASE}/results/colony"

# Metapopulations to analyse (all six)
METAPOPS=("North-western" "North-eastern" "Itcha-Ilgachuz" \
          "Central-eastern" "Jasper-Banff" "South-eastern")

N_SNP=500                          # informative SNPs retained per metapopulation
SEED=98765                         # COLONY random seed

mkdir -p "$OUTBASE" "${BASE}/logs"

SUMMARY="${OUTBASE}/colony_Ne_summary.txt"
{
    echo "COLONY sibship-based Ne summary"
    echo "Run date: $(date)"
    echo "Settings: ${N_SNP} SNPs (MAF>0.4, r^2<0.1), Full Likelihood,"
    echo "          run length 3 (very long), precision 2 (high),"
    echo "          polygamy both sexes, non-random mating (inbreeding on)."
    echo "------------------------------------------------------------"
} > "$SUMMARY"

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
for POP in "${METAPOPS[@]}"; do
    echo "==== ${POP} ===="
    OUTDIR="${OUTBASE}/${POP}"
    mkdir -p "$OUTDIR"

    IN_STUB="${INDIR}/${POP}"
    if [ ! -f "${IN_STUB}.ped" ]; then
        echo "  WARNING: ${IN_STUB}.ped not found, skipping"; continue
    fi

    # --- Step 1: select 500 informative, LD-pruned SNPs --------------------
    echo "  Selecting ${N_SNP} informative SNPs..."

    # LD prune to independent markers with MAF > 0.4
    "$PLINK" --file "$IN_STUB" --allow-extra-chr --chr-set 34 \
        --maf 0.4 --geno 0.05 --indep-pairwise 50 10 0.1 \
        --out "${OUTDIR}/${POP}_candidates" > /dev/null 2>&1

    # Frequency and missingness for ranking the pruned markers
    "$PLINK" --file "$IN_STUB" --allow-extra-chr --chr-set 34 \
        --extract "${OUTDIR}/${POP}_candidates.prune.in" \
        --freq --missing --out "${OUTDIR}/${POP}_stats" > /dev/null 2>&1

    # Rank by lowest missingness, then highest MAF; keep the top N_SNP
    OUTDIR="$OUTDIR" POP="$POP" N_SNP="$N_SNP" python3 - << 'EOF'
import os
outdir, pop, n_snp = os.environ["OUTDIR"], os.environ["POP"], int(os.environ["N_SNP"])

lmiss_path = f"{outdir}/{pop}_stats.lmiss"
frq_path   = f"{outdir}/{pop}_stats.frq"
if not os.path.exists(lmiss_path):
    raise SystemExit(f"missing {lmiss_path}")

miss = {}
with open(lmiss_path) as f:
    next(f)
    for line in f:
        p = line.split()
        if p:
            miss[p[1]] = float(p[4])

data = []
with open(frq_path) as f:
    next(f)
    for line in f:
        p = line.split()
        if p:
            snp = p[1]
            data.append({"snp": snp, "miss": miss[snp], "maf": float(p[4])})

data.sort(key=lambda x: (x["miss"], -x["maf"]))
with open(f"{outdir}/{pop}_colony_snps.txt", "w") as out:
    for item in data[:n_snp]:
        out.write(item["snp"] + "\n")
EOF

    # Extract the final SNP set to PED/MAP
    "$PLINK" --file "$IN_STUB" --allow-extra-chr --chr-set 34 \
        --extract "${OUTDIR}/${POP}_colony_snps.txt" \
        --recode --out "${OUTDIR}/${POP}_final" > /dev/null 2>&1

    # --- Step 2: build the COLONY .dat input -------------------------------
    echo "  Formatting COLONY input..."
    OUTDIR="$OUTDIR" POP="$POP" SEED="$SEED" python3 - << 'EOF'
import os
outdir, pop, seed = os.environ["OUTDIR"], os.environ["POP"], os.environ["SEED"]

# PLINK A/C/G/T alleles -> COLONY integer codes; missing (0) stays 0.
code = {"A": "1", "C": "2", "G": "3", "T": "4", "0": "0"}
def convert(line):
    p = line.split()
    return f"{p[1]} " + " ".join(code.get(g, "0") for g in p[6:])

with open(f"{outdir}/{pop}_final.map") as f:
    markers = [l.split()[1] for l in f if l.strip()]
with open(f"{outdir}/{pop}_final.ped") as f:
    offspring = [convert(l) for l in f if l.strip()]

with open(f"{outdir}/{pop}_colony.dat", "w") as out:
    out.write(f"'{pop}'\n'{pop}_out'\n")
    out.write(f"{len(offspring)}\n{len(markers)}\n{seed}\n")
    out.write("0\n")     # 0 = do not update allele frequencies
    out.write("1\n")     # 1 = dioecious
    out.write("1\n")     # 1 = non-random mating (inbreeding accounted for)
    out.write("0\n")     # 0 = diploid
    out.write("1 1\n")   # polygamy for males and females
    out.write("0\n")     # 0 = no clones
    out.write("1\n")     # 1 = scale sibship (for Ne estimation)
    out.write("0\n")     # 0 = no sibship size prior
    out.write("0\n")     # 0 = allele frequencies unknown
    out.write("1\n")     # number of runs
    out.write("3\n")     # run length: 3 = very long
    out.write("0\n")     # 0 = no monitoring
    out.write("1000\n")  # monitor interval
    out.write("0\n")     # 0 = non-Windows version
    out.write("2\n")     # analysis method: 2 = Full Likelihood
    out.write("2\n")     # precision: 2 = high

    out.write(" ".join(markers) + "\n")
    out.write(" ".join(["0"]    * len(markers)) + "\n")  # 0 = codominant markers
    out.write(" ".join(["0.02"] * len(markers)) + "\n")  # dropout (allelic) rate
    out.write(" ".join(["0.01"] * len(markers)) + "\n")  # genotyping error rate
    for row in offspring:
        out.write(row + "\n")
    # No candidate parents / known relationships blocks
    out.write("0.0 0.0\n0 0\n0 0\n0 0\n0\n0\n0\n0\n0\n0\n")
EOF

    # --- Step 3: run COLONY ------------------------------------------------
    echo "  Running COLONY (Full Likelihood; this is slow)..."
    ( cd "$OUTDIR" && "$COLONY" IFN:"${POP}_colony.dat" > /dev/null 2>&1 )

    # --- Step 4: record result --------------------------------------------
    if [ -f "${OUTDIR}/${POP}_out.Ne" ]; then
        echo "  Ne output found."
        {
            echo ">> ${POP}"
            cat "${OUTDIR}/${POP}_out.Ne"
            echo "------------------------------------------------------------"
        } >> "$SUMMARY"
    else
        echo "  ERROR: COLONY produced no Ne output for ${POP}"
    fi
done

echo "COLONY analysis complete: $(date)"
echo "Summary: $SUMMARY"
