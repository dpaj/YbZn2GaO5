#!/usr/bin/env bash
# Drive the multi-process scaling test. One process per NUMA node, pinned with
# taskset to that node's 16 PHYSICAL cores (SMT siblings excluded to avoid
# hyperthread contention masquerading as bandwidth saturation).
#
# numactl is not installed, so memory policy relies on Linux first-touch: with all
# of a process's threads on one node's CPUs, its allocations land on that node's
# memory. That is effectively --membind for this workload, but it is an assumption
# rather than an enforcement, which is why an UNPINNED control run is included.
set -u
SP="${SP:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd /home/vdp/repos/YbZn2GaO5

NSPEC=${NSPEC:-6}
NQSIDE=${NQSIDE:-3}
THREADS=${THREADS:-16}

# Physical cores per NUMA node on this box (2x EPYC 7742, 8 nodes x 16 cores).
NODE_CPUS=(0-15 16-31 32-47 48-63 64-79 80-95 96-111 112-127)

run_config() {
  local N=$1 pinned=$2 label=$3
  local B=$SP/mpbar_${label}; rm -rf "$B"; mkdir -p "$B"
  local pids=()
  for ((k=0; k<N; k++)); do
    if [ "$pinned" = "yes" ]; then
      taskset -c "${NODE_CPUS[$k]}" env MP_TAG="n$k" MP_NSPEC=$NSPEC MP_NQ_SIDE=$NQSIDE \
        MP_BARRIER="$B" MP_NPROC=$N \
        julia -t $THREADS --project=. $SP/benchmark_process_scaling.jl > "${OUT:-/tmp}/mp_${label}_n$k.log" 2>&1 &
    else
      env MP_TAG="n$k" MP_NSPEC=$NSPEC MP_NQ_SIDE=$NQSIDE \
        MP_BARRIER="$B" MP_NPROC=$N \
        julia -t $THREADS --project=. $SP/benchmark_process_scaling.jl > "${OUT:-/tmp}/mp_${label}_n$k.log" 2>&1 &
    fi
    pids+=($!)
  done
  for p in "${pids[@]}"; do wait "$p"; done

  # Aggregate over the barrier-synchronised timed phase.
  awk -v N=$N -v label="$label" -v nspec=$NSPEC '
    /^RESULT/ {
      for (i=1;i<=NF;i++) {
        split($i,a,"="); if (a[1]=="compute_s") c+=a[2];
        if (a[1]=="t0") { if (t0=="" || a[2]<t0) t0=a[2] }
        if (a[1]=="t1") { if (a[2]>t1) t1=a[2] }
        if (a[1]=="per_spec_s") { ps+=a[2]; n++ }
        if (a[1]=="rss_MB") { rss+=a[2] }
        if (a[1]=="chunks") { ch=a[2] }
      }
    }
    END {
      span = t1 - t0
      printf "%-10s N=%-3d chunks=%-4s mean_per_spec=%6.2fs  span=%7.2fs  spectra=%3d  throughput=%6.3f spec/s  rel=%s  totalRSS=%.1f GB\n",
             label, N, ch, ps/n, span, N*nspec, (N*nspec)/span, "", rss/1024
    }' ${OUT:-/tmp}/mp_${label}_n*.log
}

echo "=== per-process: $THREADS threads, $NSPEC spectra, grid side $NQSIDE ==="
echo "free RAM before: $(free -g | awk '/^Mem:/{print $7" GB"}')"
for N in 1 2 4 8; do
  run_config $N yes "pin$N"
  echo "  free RAM now: $(free -g | awk '/^Mem:/{print $7" GB"}')"
done
echo "=== unpinned control at N=8 (separates NUMA locality from mere multi-process) ==="
run_config 8 no "nopin8"
echo "DRIVER_DONE"
