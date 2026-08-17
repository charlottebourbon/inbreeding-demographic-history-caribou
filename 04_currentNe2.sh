#!/bin/bash
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=1-00:00:00
#SBATCH --job-name=currentNe2
#SBATCH --output=logs/currentNe2_%j.log
#
#  04 - Contemporary effective population size (Ne) with currentNe2
#
#  Companion code for:
#    Bourbon et al. (2026) Inbreeding and Demographic History of Caribou
#    (Rangifer tarandus) in Western Canada Inferred From Genome-Wide SNP Data.
#    Evolutionary Applications. doi:10.1111/eva.70311
#
#  Workflow (see Methods 2.5):
#    Step 1. Build pooled genetic units by merging subpopulations (PLINK).
#    Step 2. Estimate Ne with currentNe2 at three hierarchical levels:
#              - subpopulations   (panmictic model)
#              - genetic units    (panmictic model)
#              - metapopulations  (panmictic AND structure-aware "-x" models,
#                                  compared per Santiago et al. 2025)
#            The full dataset is also run under both models as a reference.
#
#  Recombination rate: r = 1.0 cM/Mb (primary), 0.9 and 1.1 (sensitivity).
#
#  Requires: PLINK v1.9, currentNe2 v2.0. Set the paths in CONFIGURATION.
#
#  INPUT LAYOUT (relative to $BASE)
#    data/subpopulations/*.ped/.map    one PED/MAP per subpopulation (n > 5)
#    data/metapopulations/*.ped/.map   one PED/MAP per metapopulation
#    data/caribou_759_autosomes.ped    full dataset (34 autosomes)
#
#  OUTPUT
#    results/ne_subpopulation/   Ne per subpopulation, per rate
#    results/ne_genetic_unit/    Ne per pooled genetic unit, per rate
#    results/ne_metapopulation/  Ne per metapopulation (panmixia + structure)
#    logs/                       one log per run
# ===========================================================================

set -uo pipefail

# ---------------------------------------------------------------------------
# CONFIGURATION  -- edit these for your environment
# ---------------------------------------------------------------------------
BASE="."                                   # project root
PLINK="plink"                              # PLINK v1.9 executable
CURRENTNE2="currentne2"                    # currentNe2 v2.0 executable
THREADS=8
REC_RATES=(0.9 1.0 1.1)                    # primary 1.0; sensitivity 0.9, 1.1

SUBPOP_DIR="${BASE}/data/subpopulations"
METAPOP_DIR="${BASE}/data/metapopulations"
FULL_DATASET="${BASE}/data/caribou_759_autosomes.ped"

UNIT_DIR="${BASE}/data/genetic_units"      # merged units are written here
RES_SUBPOP="${BASE}/results/ne_subpopulation"
RES_UNIT="${BASE}/results/ne_genetic_unit"
RES_METAPOP="${BASE}/results/ne_metapopulation"

mkdir -p "$UNIT_DIR" "$RES_SUBPOP" "$RES_UNIT" "$RES_METAPOP" "${BASE}/logs"

echo "currentNe2 analysis started: $(date)"

# ---------------------------------------------------------------------------
# STEP 1. Build pooled genetic units by merging subpopulations
# ---------------------------------------------------------------------------
# Genetic units aggregate subpopulations that share the finest-scale cluster in
# Deakin et al. (2026). The first argument is the output unit name; the rest are
# the subpopulation file stems (basename without extension) in SUBPOP_DIR.
merge_unit() {
    local output_name=$1; shift
    local members=("$@")

    echo "Building genetic unit: ${output_name}  (${members[*]})"

    if [ "${#members[@]}" -eq 1 ]; then
        # Single-member "unit": just copy through PLINK for consistent formatting
        "$PLINK" --ped "${SUBPOP_DIR}/${members[0]}.ped" \
                 --map "${SUBPOP_DIR}/${members[0]}.map" \
                 --chr-set 34 --recode \
                 --out "${UNIT_DIR}/${output_name}"
        return
    fi

    # Merge list: every member after the first
    local merge_list="${UNIT_DIR}/tmp_${output_name}_merge.txt"
    rm -f "$merge_list"
    for i in "${!members[@]}"; do
        [ "$i" -gt 0 ] || continue
        echo "${SUBPOP_DIR}/${members[$i]}.ped ${SUBPOP_DIR}/${members[$i]}.map" >> "$merge_list"
    done

    "$PLINK" --ped "${SUBPOP_DIR}/${members[0]}.ped" \
             --map "${SUBPOP_DIR}/${members[0]}.map" \
             --merge-list "$merge_list" \
             --chr-set 34 --recode \
             --out "${UNIT_DIR}/${output_name}"

    rm -f "$merge_list"
}

# Pooled genetic units (Table 1 / Table S2)
merge_unit "Carcross_Atlin"                         "Carcross" "Atlin"
merge_unit "Little_Rancheria_Level_Kawdy_Horseranch_Tsenaglode" \
                                                    "Little_Rancheria" "Level_Kawdy" "Horseranch" "Tsenaglode"
merge_unit "Telkwa_Tweedsmuir"                      "Telkwa" "Tweedsmuir"
merge_unit "Calendar_Maxhamish"                     "Calendar" "Maxhamish"
merge_unit "Snake_Sahtaneh_Hay_River"               "Snake_Sahtaneh" "Hay_River"
merge_unit "Westside_Fort_Nelson_Chinchaga"         "Westside_Fort_Nelson" "Chinchaga"
merge_unit "Frog_Gataga"                            "Frog" "Gataga"
merge_unit "Finlay_East_Williston_Graham"           "Finlay" "East_Williston" "Graham"
merge_unit "Narraway_Quintette"                     "Narraway" "Quintette"
merge_unit "Narrow_Lake_North_Cariboo_Barkerville"  "Narrow_Lake" "North_Cariboo" "Barkerville"
merge_unit "Maligne_Brazeau_Banff"                  "Maligne" "Brazeau" "Banff"

# ---------------------------------------------------------------------------
# Helper: run currentNe2 over all recombination rates
#   run_ne <input.ped> <output_dir> <label> [extra currentNe2 flags]
# ---------------------------------------------------------------------------
run_ne() {
    local ped=$1 outdir=$2 label=$3; shift 3
    local extra=("$@")
    for r in "${REC_RATES[@]}"; do
        "$CURRENTNE2" "${extra[@]}" -r "$r" -t "$THREADS" \
            -o "${outdir}/${label}_r${r}" "$ped" \
            > "${BASE}/logs/${label}_r${r}.log" 2>&1 \
            && echo "  ok  ${label} r=${r}" \
            || echo "  FAIL ${label} r=${r} (see logs)"
    done
}

# ---------------------------------------------------------------------------
# STEP 2a. Subpopulations (panmictic model)
# ---------------------------------------------------------------------------
echo "== Subpopulations =="
for ped in "$SUBPOP_DIR"/*.ped; do
    [ -f "$ped" ] || continue
    run_ne "$ped" "$RES_SUBPOP" "$(basename "$ped" .ped)"
done

# ---------------------------------------------------------------------------
# STEP 2b. Pooled genetic units (panmictic model)
# ---------------------------------------------------------------------------
echo "== Genetic units =="
for ped in "$UNIT_DIR"/*.ped; do
    [ -f "$ped" ] || continue
    run_ne "$ped" "$RES_UNIT" "$(basename "$ped" .ped)"
done

# ---------------------------------------------------------------------------
# STEP 2c. Metapopulations (panmictic AND structure-aware "-x")
# ---------------------------------------------------------------------------
# Compare the two models per metapopulation: overlapping 90% CIs => effectively
# panmictic (keep panmixia); separated CIs => structured (keep "-x").
echo "== Metapopulations =="
for ped in "$METAPOP_DIR"/*.ped; do
    [ -f "$ped" ] || continue
    b=$(basename "$ped" .ped)
    run_ne "$ped" "$RES_METAPOP" "${b}_panmixia"
    run_ne "$ped" "$RES_METAPOP" "${b}_structure" -x
done

# Full dataset, both models (reference)
if [ -f "$FULL_DATASET" ]; then
    echo "== Full dataset (n=759) =="
    run_ne "$FULL_DATASET" "$RES_METAPOP" "ALL_759_panmixia"
    run_ne "$FULL_DATASET" "$RES_METAPOP" "ALL_759_structure" -x
fi

echo "currentNe2 analysis complete: $(date)"
