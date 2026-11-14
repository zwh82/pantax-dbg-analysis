#!/bin/bash
# Reproduce the 2026-11-14 full run (73 non-ERAS cohorts + ERAS, analysis,
# core-species Control networks and core-species role analysis).
#
# One command:
#   bash reproduce.sh
#
# It launches the 12 recorded cohort shards on node001/node003, waits for all
# 13 pipeline workers, checks their exit codes, then runs analysis, builds the
# 20 core-species Control networks and runs the core-species analysis on
# node001. Outputs are written to:
#   network/       cohort networks
#   analysis/      downstream tables and figures
#   core_network/  20 Control-only core-species networks
#   core_species/  core-species tables, figures and manuscript summary
#
# Re-running overwrites network/, analysis/, core_network/ and core_species/
# (the original /home/work/wenhai/meta/high_rare_net4 tree is never touched).
#
# Useful overrides:
#   NODE_A=node001 NODE_B=node003
#   SKIP_ANALYSIS=1        only rebuild the networks
#   SKIP_CORE=1            skip core_network/ and core_species/
#   DRY_RUN=1              print the ssh commands without running anything
#   SNAPSHOT_ONLY=1        only save the executed code/config snapshot
#
# Before launching, the exact code, configs, cohort tables and launchers used
# for the run are copied to control/executed_files_<UTC timestamp>/ together
# with the git commit and a COMMANDS.txt transcript.
#
# The per-worker commands are the same as the recorded launchers in
# control/launchers/; each worker saves its effective config under
# control/configs/.

set -euo pipefail

RUN=${RUN:-/home/work/wenhai/meta/pantax-dbg-analysis_results_20261114}
REPO=${REPO:-/home/work/wenhai/meta/pantax-dbg-analysis}
RSCRIPT=${RSCRIPT:-/home/wenhai/miniconda3/envs/r_4.3/bin/Rscript}
NODE_A=${NODE_A:-node001}
NODE_B=${NODE_B:-node003}
SKIP_ANALYSIS=${SKIP_ANALYSIS:-0}
SKIP_CORE=${SKIP_CORE:-0}
SNAPSHOT_ONLY=${SNAPSHOT_ONLY:-0}
MAX_WAIT=${MAX_WAIT:-21600}   # seconds; 6 hours
DRY_RUN=${DRY_RUN:-0}

WORKERS_A=(w01_node001 w02_node001 w03_node001 w04_node001 w05_node001 w06_node001 e13_eras)
WORKERS_B=(w07_node003 w08_node003 w09_node003 w10_node003 w11_node003 w12_node003)
EXPECTED=$(( ${#WORKERS_A[@]} + ${#WORKERS_B[@]} ))

launch_worker() {
  local node=$1 worker=$2
  local launcher="$RUN/control/launchers/$worker.sh"
  if [ ! -f "$launcher" ]; then
    echo "Missing launcher: $launcher" >&2
    exit 1
  fi
  if [ "$DRY_RUN" = "1" ]; then
    printf 'ssh %s "setsid bash %q" &\n' "$node" "$launcher"
    return 0
  fi
  echo "launch $worker on $node"
  ssh -o BatchMode=yes "$node" \
    "setsid bash '$launcher' </dev/null >/dev/null 2>&1 & echo \$!" \
    > "$RUN/control/status/$worker.pid"
}

wait_pipeline() {
  local expected=$1
  local deadline=$((SECONDS + MAX_WAIT))
  while :; do
    local done_n
    done_n=$(find "$RUN/control/status" -maxdepth 1 \
      \( -name 'w*.exit' -o -name 'e13_eras.exit' \) | wc -l)
    echo "$(date -Is) pipeline workers finished: $done_n/$expected"
    [ "$done_n" -ge "$expected" ] && break
    if [ "$SECONDS" -ge "$deadline" ]; then
      echo "Timed out after ${MAX_WAIT}s; see control/logs/." >&2
      exit 1
    fi
    sleep 30
  done
}

check_pipeline() {
  local fail=0 f value
  for f in "$RUN"/control/status/w*.exit "$RUN"/control/status/e13_eras.exit; do
    value=$(cat "$f")
    if [ "$value" != "0" ]; then
      echo "FAILED: $(basename "$f" .exit) exit=$value" >&2
      fail=1
    fi
  done
  return "$fail"
}

# Copy every file that the run executes (code, configs, cohort tables,
# launchers) into the result tree, with the git commit and a command
# transcript, so the run can be reconstructed after the code changes.
snapshot_executed_files() {
  local stamp dir
  stamp=$(date -u +%Y%m%dT%H%M%SZ)
  dir="$RUN/control/executed_files_$stamp"
  mkdir -p "$dir/analysis" "$dir/control"
  cp -a "$REPO/run_pipeline.R" "$REPO/run.sh" "$REPO/R" "$REPO/config" "$dir/"
  cp -a "$REPO/analysis/run_analysis.R" "$REPO/analysis/R" \
        "$REPO/analysis/config" "$REPO/analysis/verification" "$dir/analysis/"
  cp -a "$REPO/core_species" "$dir/"
  cp -a "$REPO/verification" "$dir/"
  cp -a "$REPO/README.md" "$REPO/README.zh.md" "$REPO/LICENSE" "$dir/" 2>/dev/null || true
  cp -a "$RUN/control/launchers" "$RUN/control/cohort_sets" \
        "$RUN/control/configs" "$RUN/control/disease_name_mapping.tsv" "$dir/control/"
  cp -a "$0" "$dir/reproduce.sh"
  if [ -d "$REPO/.git" ]; then
    (cd "$REPO" && git rev-parse HEAD > "$dir/GIT_COMMIT" &&
                   git status --short > "$dir/GIT_STATUS")
  fi
  {
    echo "snapshot_utc=$stamp"
    echo "run_root=$RUN"
    echo "repo=$REPO"
    echo "host=$(hostname)"
    echo "r_version=$("$RSCRIPT" -e 'cat(R.version.string)' 2>/dev/null || true)"
    echo
    for f in "$RUN"/control/launchers/*.sh; do
      echo "### $f"
      cat "$f"
      echo
    done
  } > "$dir/COMMANDS.txt"
  echo "$dir" > "$RUN/control/last_snapshot.txt"
  echo "Saved executed files: $dir"
}

run_analysis() {
  local launcher="$RUN/control/launchers/analysis.sh"
  if [ "$DRY_RUN" = "1" ]; then
    printf 'ssh %s "setsid bash %q" &\n' "$NODE_A" "$launcher"
    return 0
  fi
  rm -f "$RUN/control/status/analysis.exit" "$RUN/control/status/analysis.state"
  echo "launch analysis on $NODE_A"
  ssh -o BatchMode=yes "$NODE_A" \
    "setsid bash '$launcher' </dev/null >/dev/null 2>&1 & echo \$!" \
    > "$RUN/control/status/analysis.pid"
  local deadline=$((SECONDS + 3600))
  while [ ! -f "$RUN/control/status/analysis.exit" ]; do
    if [ "$SECONDS" -ge "$deadline" ]; then
      echo "Analysis timed out; see control/logs/analysis.log" >&2
      exit 1
    fi
    sleep 20
  done
  local status
  status=$(cat "$RUN/control/status/analysis.exit")
  echo "analysis exit=$status"
  [ "$status" = "0" ]
}

# Launch a recorded launcher on NODE_A and wait for its .exit file.
run_recorded() {
  local label=$1
  local launcher="$RUN/control/launchers/$label.sh"
  local timeout=$2
  if [ "$DRY_RUN" = "1" ]; then
    printf 'ssh %s "setsid bash %q" &\n' "$NODE_A" "$launcher"
    return 0
  fi
  rm -f "$RUN/control/status/$label.exit" "$RUN/control/status/$label.state"
  echo "launch $label on $NODE_A"
  ssh -o BatchMode=yes "$NODE_A" \
    "setsid bash '$launcher' </dev/null >/dev/null 2>&1 & echo \$!" \
    > "$RUN/control/status/$label.pid"
  local deadline=$((SECONDS + timeout))
  while [ ! -f "$RUN/control/status/$label.exit" ]; do
    if [ "$SECONDS" -ge "$deadline" ]; then
      echo "$label timed out; see control/logs/$label.log" >&2
      exit 1
    fi
    sleep 20
  done
  local status
  status=$(cat "$RUN/control/status/$label.exit")
  echo "$label exit=$status"
  [ "$status" = "0" ]
}

if [ "$DRY_RUN" != "1" ]; then
  mkdir -p "$RUN/control/status"
  if [ "$SNAPSHOT_ONLY" != "1" ]; then
    rm -f "$RUN"/control/status/w*.exit "$RUN/control/status/e13_eras.exit" \
          "$RUN/control/status/w*.state" "$RUN/control/status/e13_eras.state" \
          "$RUN/control/status/analysis.exit" "$RUN/control/status/analysis.state" \
          "$RUN/control/status/core_network.exit" "$RUN/control/status/core_network.state" \
          "$RUN/control/status/core_species.exit" "$RUN/control/status/core_species.state"
  fi
  snapshot_executed_files
fi

if [ "$SNAPSHOT_ONLY" = "1" ]; then
  echo "Snapshot only; nothing was launched."
  exit 0
fi

for worker in "${WORKERS_A[@]}"; do
  launch_worker "$NODE_A" "$worker"
done
for worker in "${WORKERS_B[@]}"; do
  launch_worker "$NODE_B" "$worker"
done

if [ "$DRY_RUN" = "1" ]; then
  echo "+ bash $RUN/control/launchers/analysis.sh"
  if [ "$SKIP_CORE" != "1" ]; then
    echo "+ bash $RUN/control/launchers/core_network.sh"
    echo "+ bash $RUN/control/launchers/core_species.sh"
  fi
  exit 0
fi

wait_pipeline "$EXPECTED"
if ! check_pipeline; then
  echo "Some pipeline workers failed; see control/logs/." >&2
  exit 1
fi

if [ "$SKIP_ANALYSIS" != "1" ]; then
  run_analysis
fi

if [ "$SKIP_CORE" != "1" ]; then
  run_recorded core_network 7200
  run_recorded core_species 21600
fi

echo "Done."
echo "  networks:      $RUN/network"
echo "  analysis:      $RUN/analysis"
echo "  core networks: $RUN/core_network/all_country_core"
echo "  core species:  $RUN/core_species"
echo "  logs:          $RUN/control/logs"
