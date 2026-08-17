#!/bin/bash
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00:00
#SBATCH --job-name=GONE2_historical
#SBATCH --output=logs/GONE2_%j.log
#
#  06 - Historical effective population size (Ne) with GONE2
#
#  Companion code for:
#    Bourbon et al. (2026) Inbreeding and Demographic History of Caribou
#    (Rangifer tarandus) in Western Canada Inferred From Genome-Wide SNP Data.
#    Evolutionary Applications. doi:10.1111/eva.70311
#
#  Reconstructs Ne over the last ~100 generations with GONE2 v2.0 under the
#  panmictic model (Methods 2.7). For each dataset:
#    1. Remove chromosome 33 (genetic length < 20 cM; below the threshold for
#       reliable Ne estimation), leaving 31,859 SNPs.
#    2. Run 100 independent seeds per recombination rate (0.9, 1.0, 1.1 cM/Mb;
#       primary rate 1.0). Ne is later summarized as the mean +/- SD across seeds.
#    3. Consolidate the per-seed Ne trajectories into one table (generation x seed).
#
#  Datasets are analysed at three hierarchical levels: individual subpopulations
#  (n > 5), pooled genetic units, and metapopulations (plus the full dataset).
#  Pooled genetic units are built here by merging subpopulations with PLINK.
#
#  Requires: PLINK v1.9, GONE2 v2.0. Set paths in CONFIGURATION.
#
#  INPUT LAYOUT (relative to $BASE)
#    data/subpopulations/*.ped/.map     one per subpopulation
#    data/metapopulations/*.ped/.map    one per metapopulation
#    data/caribou_759_autosomes.ped/.map  full dataset
#
#  OUTPUT
#    results/gone2/<dataset>_r<rate>_Ne.txt   generation x 100-seed Ne table
# ===========================================================================

set -uo pipefail

# ---------------------------------------------------------------------------
# CONFIGURATION  -- edit for your environment
# ---------------------------------------------------------------------------
BASE="."
PLINK="plink"
GONE2="gone2"

N_SEEDS=100
REC_RATES=(0.9 1.0 1.1)     # primary 1.0; sensitivity 0.9, 1.1
MAF=0.0                     # no MAF filtering (see Methods)
LOWER_C=0.001               # GONE2 lower recombination bin
UPPER_C=0.05                # GONE2 upper recombination bin

SUBPOP_DIR="${BASE}/data/subpopulations"
METAPOP_DIR="${BASE}/data/metapopulations"
UNIT_DIR="${BASE}/data/genetic_units"          # merged units written here
FULL_PED="${BASE}/data/caribou_759_autosomes.ped"
FULL_MAP="${BASE}/data/caribou_759_autosomes.map"

RES_DIR="${BASE}/results/gone2"
TMP_ROOT="${BASE}/tmp_gone2"

mkdir -p "$UNIT_DIR" "$RES_DIR" "$TMP_ROOT" "${BASE}/logs"

echo "GONE2 historical Ne started: $(date)"


# ---------------------------------------------------------------------------
# STEP 1. Build pooled genetic units (PLINK merge)
# ---------------------------------------------------------------------------
# Arg 1 = output unit name; remaining args = subpopulation file stems in SUBPOP_DIR.
merge_unit() {
    local output_name=$1; shift
    local members=("$@")

    if [ "${#members[@]}" -eq 1 ]; then
        "$PLINK" --ped "${SUBPOP_DIR}/${members[0]}.ped" \
                 --map "${SUBPOP_DIR}/${members[0]}.map" \
                 --chr-set 34 --recode \
                 --out "${UNIT_DIR}/${output_name}" > /dev/null 2>&1
        return
    fi

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
             --out "${UNIT_DIR}/${output_name}" > /dev/null 2>&1
    rm -f "$merge_list"
    echo "  built genetic unit: ${output_name}"
}

echo "== Building genetic units =="
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
# Helper: remove chromosome 33 from a PED/MAP pair
# ---------------------------------------------------------------------------
# Chr 33 falls below the 20 cM length threshold for reliable Ne estimation.
filter_chr33() {
    local ped=$1 map=$2 out_prefix=$3

    awk '$1 != 33' "$map" > "${out_prefix}.map"
    awk '$1 == 33 {print NR}' "$map" > "${out_prefix}_chr33pos.txt"

    if [ -s "${out_prefix}_chr33pos.txt" ]; then
        # Drop the genotype columns (two per SNP) belonging to chr-33 markers.
        awk -v posfile="${out_prefix}_chr33pos.txt" '
            BEGIN { while ((getline < posfile) > 0) pos[$1] = 1; close(posfile) }
            {
                printf "%s %s %s %s %s %s", $1, $2, $3, $4, $5, $6
                for (i = 1; i <= (NF-6)/2; i++)
                    if (!(i in pos)) printf " %s %s", $(6 + i*2 - 1), $(6 + i*2)
                printf "\n"
            }' "$ped" > "${out_prefix}.ped"
    else
        cp "$ped" "${out_prefix}.ped"
    fi
    rm -f "${out_prefix}_chr33pos.txt"
}


# ---------------------------------------------------------------------------
# Helper: consolidate per-seed GONE2 Ne outputs into one table
# ---------------------------------------------------------------------------
# Produces: Generation  Seed0001  Seed0002 ... (one Ne column per seed).
consolidate() {
    local prefix=$1 out_file=$2 n_seeds=$3
    local first="" seeds_found=0

    for s in $(seq 1 "$n_seeds"); do
        [ -f "${prefix}_seed${s}_GONE2_Ne" ] || continue
        [ -z "$first" ] && first="${prefix}_seed${s}_GONE2_Ne"
        ((seeds_found++))
    done

    if [ "$seeds_found" -eq 0 ]; then
        echo "  ERROR: no seed outputs for $(basename "$prefix")"; return 1
    fi

    # Generation column from the first available seed
    tail -n +2 "$first" | awk '{print $1}' > "${out_file}.gen"

    local cols=("${out_file}.gen")
    for s in $(seq 1 "$n_seeds"); do
        local f="${prefix}_seed${s}_GONE2_Ne"
        if [ -f "$f" ]; then
            tail -n +2 "$f" | awk '{print $2}' > "${out_file}.s${s}"
        else
            tail -n +2 "$first" | awk '{print "NA"}' > "${out_file}.s${s}"
        fi
        cols+=("${out_file}.s${s}")
    done

    { printf "Generation"; for s in $(seq 1 "$n_seeds"); do printf "\tSeed%04d" "$s"; done; printf "\n"; } > "$out_file"
    paste "${cols[@]}" >> "$out_file"
    rm -f "${out_file}.gen" "${out_file}".s*

    echo "  consolidated $(basename "$out_file")  (${seeds_found}/${n_seeds} seeds)"
}


# ---------------------------------------------------------------------------
# Helper: full GONE2 run for one dataset (chr33 filter + all rates + consolidate)
# ---------------------------------------------------------------------------
run_gone2() {
    local ped=$1 map=$2 label=$3
    local tmp="${TMP_ROOT}/${label}"
    mkdir -p "$tmp"

    filter_chr33 "$ped" "$map" "${tmp}/${label}_noChr33"
    local fped="${tmp}/${label}_noChr33.ped"

    for r in "${REC_RATES[@]}"; do
        echo "  ${label}  r=${r}: running ${N_SEEDS} seeds"
        for s in $(seq 1 "$N_SEEDS"); do
            "$GONE2" -S "$s" -r "$r" -t 1 -M "$MAF" -l "$LOWER_C" -u "$UPPER_C" \
                -o "${tmp}/${label}_r${r}_seed${s}" "$fped" > /dev/null 2>&1
        done
        consolidate "${tmp}/${label}_r${r}" \
                    "${RES_DIR}/${label}_r${r}_Ne.txt" "$N_SEEDS"
    done

    rm -rf "$tmp"
}


# ---------------------------------------------------------------------------
# STEP 2. Run GONE2 at each hierarchical level
# ---------------------------------------------------------------------------
echo "== Subpopulations =="
for ped in "$SUBPOP_DIR"/*.ped; do
    [ -f "$ped" ] || continue
    b=$(basename "$ped" .ped)
    run_gone2 "$ped" "${ped%.ped}.map" "$b"
done

echo "== Genetic units =="
for ped in "$UNIT_DIR"/*.ped; do
    [ -f "$ped" ] || continue
    b=$(basename "$ped" .ped)
    run_gone2 "$ped" "${ped%.ped}.map" "$b"
done

echo "== Metapopulations =="
for ped in "$METAPOP_DIR"/*.ped; do
    [ -f "$ped" ] || continue
    b=$(basename "$ped" .ped)
    run_gone2 "$ped" "${ped%.ped}.map" "$b"
done

if [ -f "$FULL_PED" ] && [ -f "$FULL_MAP" ]; then
    echo "== Full dataset (n=759) =="
    run_gone2 "$FULL_PED" "$FULL_MAP" "ALL_759"
fi

echo "GONE2 historical Ne complete: $(date)"
