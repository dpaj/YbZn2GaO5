#!/usr/bin/env bash
# Launch the multi-start neutron optimization, then the follow-up analysis.
#
#   nohup scripts/run_neutron_optimization_weekend.sh > results/logs/weekend.log 2>&1 &
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
# Env overrides: YZGO_N_STARTS, YZGO_THREADS_PER_START, YZGO_NEUTRON_OPT_CONTROLS.

set -u
cd "$(dirname "$0")/.." || exit 1
mkdir -p results/logs

N=${YZGO_N_STARTS:-8}
T=${YZGO_THREADS_PER_START:-4}

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
julia -t auto --project=. scripts/analyze_neutron_optimum.jl \
    > results/logs/opt_followup.log 2>&1
echo "=== $(date) finished (follow-up exit $?) ==="
