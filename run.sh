#!/usr/bin/env bash
# Reproducible run recipes for pantax-dbg-analysis.
#
# Every run is expressed as explicit command-line flags; no config file has to
# be edited and no parameter block has to be copied around.  The effective
# config of each pipeline run is saved next to its results
# (<outdir>/configs/<scenario>.R) and the full command is written to the log.
#
# Usage:
#   ./run.sh <scenario> [extra run_pipeline.R flags]
#
# Scenarios:
#   spiec_rc       read counts + SPIEC-EASI + 10% prevalence      -> res_rc
#   eras           config/eras.yaml + config/eras.tsv + time ordering -> res_rc
#   relative       relative abundance + SPIEC-EASI                -> res_raclr_no_filter
#   flashweave     read counts + FlashWeave (prevalence skipped)  -> res_flashweave_rc
#   sparcc         read counts + SparCC                           -> res_sparcc_rc
#   filter_weight  read counts + SPIEC-EASI + |weight| >= 0.1     -> res_rc_filter
#   filter_conf    read counts + SPIEC-EASI + confidence >= 0.8   -> res_raclr_econ_filter
#   merged         merged-cohort table (set MERGED_TABLE)         -> res_rc
#   network_only   build Control networks only                    -> all_country/
#   core_network   core-species Control networks per country      -> all_country_core/
#   hadza          Hadza five- and three-timepoint season networks
#   analysis       downstream tables and figures on a result tree
#   core_species   core-species role analysis, figures, summary
#   full           both cohort tables into one tree, then analysis
#
# Every scenario writes under $OUT_ROOT (default: ./runs), never into the
# original high_rare_net4 directory.  Cohort networks go to
# <scenario>/network/ and analysis outputs to <scenario>/analysis/.
# Useful environment variables:
#   OUT_ROOT=/path            root for all scenario outputs
#   COHORTS=/path             non-ERAS cohort table
#   CORE_COHORTS=/path        core-species country table
#   HADZA_COHORTS=/path       Hadza cohort table
#   ERAS_COHORTS=/path        ERAS cohort table
#   MERGED_TABLE=/path        merged-cohort table
#   ANALYSIS_RESULTS=/path    analysis output root (default <OUT_ROOT>/full/analysis)
#   ANALYSIS_NETWORK=/path    cohort network root (default <OUT_ROOT>/full/network)
#   ANALYSIS_SUFFIX=res_rc    suffix to analyse
#   ANALYSIS_STEP=all|net_compare,rho,...  steps for the analysis scenario
#   CORE_NETWORK=/path        network root for the core-species analysis
#   CORE_STEP=all             analysis, figures, manuscript (comma-separated)
#   DRY_RUN=1                 print the commands without running them
#
# Examples:
#   ./run.sh spiec_rc
#   ./run.sh spiec_rc --only PRJNA398089_598_CD
#   ./run.sh eras
#   ./run.sh full
#   DRY_RUN=1 ./run.sh full
#   ANALYSIS_RESULTS=/data/net4 ./run.sh analysis
#
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SRC=${SRC:-/home/work/wenhai/meta/high_rare_net4}
R=${R:-/home/wenhai/miniconda3/envs/r_4.3/bin/Rscript}

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-/home/wenhai/miniconda3/envs/r_4.3/lib}
export R_LIBS=${R_LIBS:-/home/wenhai/miniconda3/envs/r_4.3/lib/R/library:/home/wenhai/application/Qmatey/tools/R}

OUT_ROOT=${OUT_ROOT:-$ROOT/runs}
COHORTS=${COHORTS:-$ROOT/config/cohorts.tsv}
CORE_COHORTS=${CORE_COHORTS:-$ROOT/config/cohorts_core.tsv}
HADZA_COHORTS=${HADZA_COHORTS:-$ROOT/config/hadza.tsv}
ERAS_COHORTS=${ERAS_COHORTS:-$ROOT/config/eras.tsv}
MERGED_TABLE=${MERGED_TABLE:-$SRC/scripts/merged_network_metadata.tsv}
ANALYSIS_RESULTS=${ANALYSIS_RESULTS:-$OUT_ROOT/full/analysis}
ANALYSIS_NETWORK=${ANALYSIS_NETWORK:-$OUT_ROOT/full/network}
ANALYSIS_SUFFIX=${ANALYSIS_SUFFIX:-res_rc}
ANALYSIS_STEP=${ANALYSIS_STEP:-all}
CORE_NETWORK=${CORE_NETWORK:-$OUT_ROOT/full/network}
CORE_STEP=${CORE_STEP:-all}
DRY_RUN=${DRY_RUN:-0}

usage() {
  awk 'NR > 1 && /^#/ { sub(/^# ?/, ""); print; next }
       NR > 1 { exit }' "${BASH_SOURCE[0]}"
}

scenario=${1:-}
if [ -z "$scenario" ] || [ "$scenario" = "-h" ] || [ "$scenario" = "--help" ]; then
  usage
  exit 0
fi
shift

if [ ! -f "$ROOT/run_pipeline.R" ]; then
  echo "run_pipeline.R not found under $ROOT" >&2
  exit 1
fi

ensure_dirs() {
  if [ "$DRY_RUN" != "1" ]; then
    mkdir -p "$OUT_ROOT/logs"
  fi
}

# Run one pipeline invocation into $2.  Remaining arguments are passed
# verbatim to run_pipeline.R; the output directory and saved config are
# appended here so every recipe is logged the same way.
run_pipeline_to() {
  local name=$1
  local outdir=$2
  shift 2
  local config_file="$outdir/configs/$name.yaml"
  local log_file="$OUT_ROOT/logs/$name.log"
  ensure_dirs
  if [ "$DRY_RUN" != "1" ]; then
    mkdir -p "$outdir/configs"
  fi

  local -a cmd=(
    "$R" "$ROOT/run_pipeline.R" "$@"
    --outdir "$outdir"
    --save-config "$config_file"
  )
  printf '+ '; printf '%q ' "${cmd[@]}"; printf '\n'
  if [ "$DRY_RUN" = "1" ]; then
    return 0
  fi
  "${cmd[@]}" 2>&1 | tee "$log_file"
}

run_analysis_to() {
  local outdir=$1
  shift
  ensure_dirs
  local log_file="$OUT_ROOT/logs/analysis_$ANALYSIS_SUFFIX.log"
  local -a cmd=(
    "$R" "$ROOT/analysis/run_analysis.R"
    --results "$outdir" --suffix "$ANALYSIS_SUFFIX" "$@"
  )
  printf '+ '; printf '%q ' "${cmd[@]}"; printf '\n'
  if [ "$DRY_RUN" = "1" ]; then
    return 0
  fi
  "${cmd[@]}" 2>&1 | tee "$log_file"
}

# The analysis steps expect scripts/disease_name_mapping.tsv and a
# scripts/final_figures directory inside the result tree.
prepare_analysis_root() {
  local outdir=$1
  if [ "$DRY_RUN" = "1" ]; then
    echo "+ mkdir -p $outdir/scripts/final_figures/main_split"
    echo "+ cp $SRC/scripts/disease_name_mapping.tsv $outdir/scripts/"
    return 0
  fi
  mkdir -p "$outdir/scripts/final_figures/main_split"
  if [ ! -f "$outdir/scripts/disease_name_mapping.tsv" ]; then
    cp "$SRC/scripts/disease_name_mapping.tsv" "$outdir/scripts/"
  fi
}

# Core-species analysis reads networks and writes tables/figures/summary under
# $OUT_ROOT/core_species/.  The effective config is saved next to the results.
run_core_species_to() {
  local outdir=$1
  shift
  ensure_dirs
  local log_file="$OUT_ROOT/logs/core_species.log"
  local -a cmd=(
    "$R" "$ROOT/core_species/run_core_species.R"
    --base-dir "$CORE_NETWORK"
    --outdir "$outdir/analysis"
    --figure-dir "$outdir/analysis/figures"
    --manuscript-dir "$outdir/analysis/manuscript_stats"
    --step "$CORE_STEP"
    --save-config "$outdir/analysis/configs/core_species.yaml"
    "$@"
  )
  if [ "$DRY_RUN" != "1" ]; then
    mkdir -p "$outdir/analysis/configs"
  fi
  printf '+ '; printf '%q ' "${cmd[@]}"; printf '\n'
  if [ "$DRY_RUN" = "1" ]; then
    return 0
  fi
  "${cmd[@]}" 2>&1 | tee "$log_file"
}

case "$scenario" in
  # ---- read counts, SPIEC-EASI, 10% prevalence ----------------------------
  spiec_rc)
    run_pipeline_to spiec_rc "$OUT_ROOT/spiec_rc/network" \
      --cohorts "$COHORTS" \
      --data-type read_count \
      --tool spiec-easi \
      --prevalence 0.1 \
      "$@"
    ;;

  # ---- ERAS follow-up cohort ----------------------------------------------
  eras)
    if [ "$DRY_RUN" != "1" ] && [ ! -f "$ERAS_COHORTS" ]; then
      echo "ERAS cohort table not found: $ERAS_COHORTS" >&2
      exit 1
    fi
    run_pipeline_to eras "$OUT_ROOT/eras/network" \
      --config "$ROOT/config/eras.yaml" \
      --cohorts "$ERAS_COHORTS" \
      "$@"
    ;;

  # ---- relative abundance, SPIEC-EASI -------------------------------------
  relative)
    run_pipeline_to relative "$OUT_ROOT/relative/network" \
      --cohorts "$COHORTS" \
      --data-type relative_abundance \
      --tool spiec-easi \
      --spiec-easi-ra TRUE \
      --prevalence 0.1 \
      "$@"
    ;;

  # ---- FlashWeave, read counts -------------------------------------------
  # The legacy read-count FlashWeave config set Notprefilter = TRUE.
  flashweave)
    run_pipeline_to flashweave "$OUT_ROOT/flashweave/network" \
      --cohorts "$COHORTS" \
      --data-type read_count \
      --tool flashweave \
      --prevalence 0.1 \
      --skip-prevalence-filter TRUE \
      "$@"
    ;;

  # ---- SparCC, read counts -----------------------------------------------
  sparcc)
    run_pipeline_to sparcc "$OUT_ROOT/sparcc/network" \
      --cohorts "$COHORTS" \
      --data-type read_count \
      --tool sparcc \
      --prevalence 0.1 \
      "$@"
    ;;

  # ---- edge filters -------------------------------------------------------
  filter_weight)
    run_pipeline_to filter_weight "$OUT_ROOT/filter_weight/network" \
      --cohorts "$COHORTS" \
      --data-type read_count \
      --tool spiec-easi \
      --prevalence 0.1 \
      --filter-weight TRUE \
      "$@"
    ;;

  filter_conf)
    run_pipeline_to filter_conf "$OUT_ROOT/filter_conf/network" \
      --cohorts "$COHORTS" \
      --data-type read_count \
      --tool spiec-easi \
      --prevalence 0.1 \
      --filter-confidence TRUE \
      "$@"
    ;;

  # ---- merged-cohort table ------------------------------------------------
  merged)
    if [ "$DRY_RUN" != "1" ] && [ ! -f "$MERGED_TABLE" ]; then
      echo "Merged cohort table not found: $MERGED_TABLE" >&2
      echo "Set MERGED_TABLE=/path/to/merged_network_metadata.tsv" >&2
      exit 1
    fi
    run_pipeline_to merged "$OUT_ROOT/merged/network" \
      --cohorts "$MERGED_TABLE" \
      --merged TRUE \
      --data-type read_count \
      --tool spiec-easi \
      --prevalence 0.1 \
      "$@"
    ;;

  # ---- Control networks only ---------------------------------------------
  network_only)
    run_pipeline_to network_only "$OUT_ROOT/network_only/network" \
      --cohorts "$COHORTS" \
      --data-type read_count \
      --tool spiec-easi \
      --prevalence 0.1 \
      --network-only TRUE \
      "$@"
    ;;

  # ---- core-species Control networks only --------------------------------
  core_network)
    run_pipeline_to core_network "$OUT_ROOT/core_network/network" \
      --config "$ROOT/config/core_network.yaml" \
      --cohorts "$CORE_COHORTS" \
      "$@"
    ;;

  # ---- Hadza season networks ---------------------------------------------
  hadza)
    run_pipeline_to hadza "$OUT_ROOT/hadza" \
      --config "$ROOT/config/hadza.yaml" \
      --cohorts "$HADZA_COHORTS" \
      "$@"
    ;;

  # ---- core-species role analysis on an existing network tree ------------
  core_species)
    if [ "$DRY_RUN" != "1" ] && [ ! -d "$CORE_NETWORK" ]; then
      echo "Network root not found: $CORE_NETWORK" >&2
      echo "Run ./run.sh full first, or set CORE_NETWORK=/path/to/network" >&2
      exit 1
    fi
    run_core_species_to "$OUT_ROOT/core_species" "$@"
    ;;

  # ---- downstream analysis on an existing result tree --------------------
  analysis)
    if [ "$DRY_RUN" != "1" ] && [ ! -d "$ANALYSIS_RESULTS" ]; then
      echo "Analysis root not found: $ANALYSIS_RESULTS" >&2
      echo "Run ./run.sh full first, or set ANALYSIS_RESULTS=/path/to/analysis" >&2
      exit 1
    fi
    if [ "$DRY_RUN" != "1" ] && [ ! -d "$ANALYSIS_NETWORK" ]; then
      echo "Network root not found: $ANALYSIS_NETWORK" >&2
      echo "Set ANALYSIS_NETWORK=/path/to/network" >&2
      exit 1
    fi
    prepare_analysis_root "$ANALYSIS_RESULTS"
    run_analysis_to "$ANALYSIS_RESULTS" \
      --network-dir "$ANALYSIS_NETWORK" \
      --step "$ANALYSIS_STEP" "$@"
    ;;

  # ---- full reproduction: 73 cohorts + ERAS, then all analysis steps -----
  full)
    full_out="$OUT_ROOT/full"
    full_network="$full_out/network"
    full_analysis="$full_out/analysis"
    run_pipeline_to full_cohorts "$full_network" \
      --cohorts "$COHORTS" \
      --data-type read_count \
      --tool spiec-easi \
      --prevalence 0.1
    run_pipeline_to full_eras "$full_network" \
      --config "$ROOT/config/eras.yaml" \
      --cohorts "$ERAS_COHORTS"
    prepare_analysis_root "$full_analysis"
    run_analysis_to "$full_analysis" \
      --network-dir "$full_network" \
      --step "$ANALYSIS_STEP"
    ;;

  *)
    echo "Unknown scenario: $scenario" >&2
    usage >&2
    exit 1
    ;;
esac
