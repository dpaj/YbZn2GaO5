# Working from a third machine (laptop, home, a different Claude session)

Short answer: **yes, this works, and it is not problematic for the things a laptop is good at.**
It becomes problematic in three specific ways, all avoidable, and all of them are about *state*
rather than about compute.

The existing setup is two machines that share state only through GitHub — a Windows desktop and
`neutrons-dgx01.ornl.gov` — with separate Claude sessions that cannot see each other. **You are the
only channel between them.** A third machine does not change that structure; it adds one more node
that can be stale, and one more place line endings and environments can diverge.

## The one rule that matters most

**`git fetch && git pull --ff-only` before doing or discussing anything.**

This has already bitten us once, and expensively: the DGX spent a session reporting that the new
magnetization and AC data "are not in the repo" and treating that as a coordination gap. The data
had been on `origin/main` for nine commits. `git status` said "up to date with 'origin/main'"
because that compares HEAD to the *locally cached* remote ref, not to GitHub.

With three nodes the chance of one being stale roughly doubles, and a stale node produces confident
wrong statements rather than errors. So: fetch first, and if a session tells you something is
missing, have it show you `git log --oneline -1` and `git ls-tree -r origin/main -- <path>` before
believing it.

## What a laptop is genuinely good for

Everything tagged **free** or **cheap** in `scripts/README.md` — which is more than it sounds,
because most of the recent scientific results came from data-only scripts:

| script | what it gives you |
|---|---|
| `plot_mpms3_centering_correction.jl` | the 1.49x centring correction, from the raw file |
| `plot_ac_susceptibility_plateau_test.jl` | AC to 18 T, both orientations, T1 as drift reference |
| `plot_neutron_1d_two_incident_energies.jl` | Ei = 3.32 vs 4.65 overlay |
| `check_thermal_slope_disorder.jl` | the thermal-slope calculation, no Sunny at all |
| `plot_background_variance_effect.jl` | reads a CSV, no recompute |
| `plot_gamma_first_scan.jl` | reads a CSV, no recompute |
| `crystal_field_van_vleck.jl` | level scheme, g-tensor, chi_VV |

Plus all of `docs/`, `CLAUDE.md`, reading and writing code, and reviewing the argument. A laptop is
a fine place to think, and a bad place to run KPM.

**Do not run the heavy fits from a laptop.** Not because it breaks anything, but because one
six-cut neutron evaluation is ~460-520 s on 32 threads and the noise floor needs n >= 34
realizations to resolve a chi2_red improvement of 1.0. That is DGX work.

## The three real traps

### 1. `results/` is gitignored, so figures do NOT travel

Every figure lives under `results/`, which is deliberately ignored — all of it is regenerable. So a
laptop clone has **no figures**, including any made on the desktop earlier that day. Consequences:

- To review a **free**-tagged figure from home, just regenerate it; it costs seconds to a minute.
- To review a **medium** or **heavy** figure (neutron cuts, 2D maps), you cannot. Either copy the
  PNG out of band, or ask for it before leaving.
- Do not ask a laptop session "what does the current fit look like" and expect it to look. It will
  either regenerate (slow, or impossible) or, worse, describe a figure it has not seen. Ask it to
  read the **feature tables** instead — those are also under `results/` and also ignored, so the
  same caveat applies. The durable record is `CLAUDE.md` and the commit messages.

### 2. Julia environment and `Manifest.toml` churn

`Manifest.toml`, `envs/sunny-main/Manifest.toml` and `envs/sunny-kpm-gpu/Manifest.toml` are
tracked, carry 167-185 `_jll` dependencies, and pin `julia_version = "1.12.3"`. A laptop with a
different Julia patch version, or a different OS, will rewrite them on `Pkg.instantiate()`.

    julia --version                 # compare against 1.12.3
    julia --project=. -e 'using Pkg; Pkg.instantiate()'
    git checkout -- Manifest.toml   # DO NOT commit the churn

If you commit it, it ping-pongs between machines forever. The CUDA artifacts in the GPU env resolve
differently on Linux, so that one is worse.

### 3. Line endings — the trap most likely to make a mess

There is **no `.gitattributes`**, and 111 tracked files carry CRLF in the index. A laptop whose git
has a different `core.autocrlf` will produce a diff touching hundreds of files that has nothing to
do with your work, and it is easy to commit by accident.

Before touching anything on a new machine:

    git config --global user.name  "Daniel Pajerowski"
    git config --global user.email "daniel@pajerowski.com"
    git config --global core.autocrlf false      # match the index; do NOT let git rewrite

Setting the identity matters independently: without it git invents one from `username@hostname`,
which is how `vdp@neutrons-dgx01.ornl.gov` got into the history. `.mailmap` canonicalizes it for git
tooling but GitHub does not honour `.mailmap` for account linking.

Then `git status` should be clean. **If it shows a large diff you did not make, stop** — that is the
line-ending renormalization, not your work. Do not commit it. Renormalizing properly is a one-time
noisy commit touching all 111 files and has to be done while every machine is clean and pulled
immediately afterwards.

## What is NOT in the repo, so a laptop clone will not have it

- **`../references/`** — 13 MB of published PDFs, a sibling of the repo, deliberately untracked
  because they are copyrighted and the DGX is a shared box.
- **`YZGO/CNCS_data/`** — outside the repo, and what the three `scripts/legacy/*.jl` want. Those
  three also have hardcoded `C:\Users\vdp\...` paths and will fail anywhere else; the fix pattern is
  `get(ENV, "YZGO_DATA_DIR", raw"C:\...")`, already used at
  `plot_yzgo_2d_data_vs_model_legacy.jl:5910`.
- **`.claude/settings.local.json`** — machine-specific and ignored, by design.
  `.claude/settings.json` and `CLAUDE.md` ARE shared config; push and pull them like code.

## Ordering, when two machines both have work

**Pull first, then commit.** Committing first makes the box simultaneously ahead and behind for no
reason. `--ff-only` is preferred because it *refuses* rather than silently creating a merge commit
when two sides have diverged — which is what you want when you cannot see the other machine.

If a laptop session and the desktop session both edited `CLAUDE.md`, expect a conflict there before
anywhere else; it is the file every session wants to write to.

## Briefing a fresh session

A new Claude on a new machine reads `CLAUDE.md` and nothing else automatically. It will not know
what happened today. The cheapest useful briefing is:

    git log --oneline -15

The commit messages in this repo are written to carry the reasoning, including retractions, so the
log is a genuine handoff and not just a changelog. Point a fresh session at it rather than
re-explaining.
