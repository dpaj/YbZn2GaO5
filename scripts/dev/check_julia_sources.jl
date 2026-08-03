#!/usr/bin/env julia
# Pre-flight check for Julia sources: parse them, and catch the mistakes that survive a parse
# check and only fail at RUN time -- which in this repo has meant failing after the compute.
#
#   julia scripts/dev/check_julia_sources.jl [paths...]     # default: scripts/ and src/
#
# Run this before executing any script that costs real time. It needs no packages and no
# project, so it starts in well under a second.
#
# WHY THIS EXISTS
#
# `@printf`/`@sprintf` require a LITERAL format string. Writing
#
#     @printf("a long format " * "continued", x)
#
# throws "First argument to @printf must be a format string" at MACRO EXPANSION during load,
# so `Meta.parseall` accepts it happily and the failure only appears when the script runs.
# That has now happened FOUR times in this repo, once costing 950 s of KPM and once a full
# overnight window. Writing the gotcha into CLAUDE.md did not prevent recurrence; a check that
# is actually executed does.
#
# Fix by putting the ARGUMENTS on continuation lines and keeping the format on one line, or by
# using string interpolation instead of a format string for long prose.

const PATTERNS = [
    (r"@s?printf\s*\(\s*(?:[A-Za-z_][A-Za-z0-9_]*\s*,\s*)?\".*\"\s*\*\s*$",
     "@printf/@sprintf format string is CONCATENATED; it must be a single literal"),
]

function check_file(path)
    problems = Tuple{Int,String}[]
    src = read(path, String)

    # 1. Does it parse at all?
    ex = try
        Meta.parseall(src)
    catch err
        push!(problems, (0, "does not parse: " * first(split(sprint(showerror, err), '\n'))))
        return problems
    end
    for a in ex.args
        if a isa Expr && (a.head === :error || a.head === :incomplete)
            push!(problems, (0, "parse error / incomplete expression"))
        end
    end

    # 2. Line patterns that parse fine but die on load or at run time.
    for (i, line) in enumerate(split(src, '\n'))
        stripped = strip(line)
        startswith(stripped, "#") && continue
        for (re, msg) in PATTERNS
            occursin(re, line) && push!(problems, (i, msg))
        end
    end
    return problems
end

function collect_files(paths)
    files = String[]
    if isempty(paths)
        root = normpath(joinpath(@__DIR__, "..", ".."))
        for d in ("scripts", "src", "test")
            isdir(joinpath(root, d)) || continue
            for (dir, _, fs) in walkdir(joinpath(root, d))
                for f in fs
                    endswith(f, ".jl") && push!(files, joinpath(dir, f))
                end
            end
        end
    else
        for p in paths
            if isdir(p)
                for (dir, _, fs) in walkdir(p)
                    for f in fs
                        endswith(f, ".jl") && push!(files, joinpath(dir, f))
                    end
                end
            else
                push!(files, p)
            end
        end
    end
    return sort(unique(files))
end

function main(args)
    files = collect_files(args)
    root = normpath(joinpath(@__DIR__, "..", ".."))
    nbad = 0
    for f in files
        ps = check_file(f)
        isempty(ps) && continue
        nbad += 1
        rel = relpath(f, root)
        for (ln, msg) in ps
            println(ln == 0 ? "  $rel: $msg" : "  $rel:$ln: $msg")
        end
    end
    println(nbad == 0 ? "OK -- $(length(files)) files, no problems" :
                        "FAILED -- $nbad of $(length(files)) files have problems")
    return nbad
end

exit(main(ARGS) == 0 ? 0 : 1)
