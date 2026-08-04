#!/usr/bin/env bash
# Launch the multi-start neutron optimization, then the follow-up analysis.
#
#   nohup scripts/run_neutron_optimization.sh > results/logs/fit.log 2>&1 &
#
# Renamed from run_neutron_optimization_weekend.sh: with the measured budget it is a
# ~16 h overnight job at 4 chains x 100 evaluations, not a weekend one.
#
# One process per Nelder-Mead start at 4 threads each. Intra-process q-threading is only
# ~73% efficient at 4 threads on this box and saturates near 3.4x at 32, so N processes x 4
# threads beats one process x 32 threads severalfold. Total threads are kept at or below the
# physical core count deliberately: on the DGX, oversubscribing past the physical cores
# reduced throughput.
#
# Every start writes its own log and its own CSV and is independent of the others, so one
# start dying costs one start. The follow-up reads whatever best_*.csv exist, so it still
# produces results if some starts failed.
#
# Env overrides: YZGO_N_STARTS, YZGO_THREADS_PER_START, YZGO_ANALYSIS_THREADS,
# YZGO_NEUTRON_OPT_CONTROLS.

set -u
cd "$(dirname "$0")/.." || exit 1
mkdir -p results/logs

# Defaults are the LATENCY optimum, which is what a sequential optimizer needs: Nelder-Mead
# cannot be parallelised internally, so per-chain wall time is set by latency, and the DGX
# measured 32 threads as the minimum (818 s at 225 q, against 1773 s at 8 threads and 875 s at
# 64). Four concurrent chains cost only 24% per-chain latency, so 4 starts come for 1.24x the
# price of one.
#
# These deliberately REPLACE an earlier 8 x 4 default, which was the THROUGHPUT optimum. The two
# optimise different quantities and 8 x 4 was simply the wrong one for this optimizer.
#
# On a 128-core box use the defaults. On a 32-core box, 1 x 32 gives one start at minimum
# latency; more chains there trade latency for starts and do worse on both.
N=${YZGO_N_STARTS:-4}
T=${YZGO_THREADS_PER_START:-32}

echo "=== $(date) launching $N starts at $T threads each ($((N * T)) total) ==="
julia --version

pids=()
for k in $(seq 0 $((N - 1))); do
    YZGO_START_INDEX=$k nohup julia -t "$T" --project=. \
        scripts/optimize_neutron_neldermead.jl \
        > "results/logs/opt_start_${k}.log" 2>&1 &
    pids+=($!)
    echo "  start $k -> pid $! -> results/logs/opt_start_${k}.log"
done

echo "=== waiting for ${#pids[@]} starts ==="
fail=0
for p in "${pids[@]}"; do
    if wait "$p"; then :; else
        echo "  pid $p exited non-zero"
        fail=$((fail + 1))
    fi
done
echo "=== $(date) all starts done ($fail non-zero exits) ==="

# The follow-up runs even if some starts failed: it reads whatever best_*.csv exist, and a
# partial multi-start still yields a usable optimum plus profiles.
echo "=== $(date) follow-up analysis ==="
# The follow-up does REAL KPM work (profile likelihood, seed validation), so -t auto is the
# worst available choice: the measured latency optimum is 32 threads and 64 is already worse.
# On a 256-thread box -t auto is badly wrong, and it contradicted this file's own advice.
AT=${YZGO_ANALYSIS_THREADS:-32}
julia -t "$AT" --project=. scripts/analyze_neutron_optimum.jl \
    > results/logs/opt_followup.log 2>&1
echo "=== $(date) finished (follow-up exit $?) ==="
