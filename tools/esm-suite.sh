#!/usr/bin/env bash
# Validate, units-check, and test a set of .esm files with the repo's esm CLI.
#
#   tools/esm-suite.sh                 # every stage2/*.esm
#   tools/esm-suite.sh stage2/cam_*.esm
#   tools/esm-suite.sh --no-test ...   # validate + units only (cheap)
#
# Sets ESD_ROOT (needed by the CICE files) if it is not already set, and
# runs `esm validate` and `esm units` file by file (they take one file each;
# see EarthSciML/EarthSciAST#643). Then runs `esm test` once on the explicit
# file list (rather than a directory; see EarthSciAST#442). The last line is
#   SUITE: <P>P/<F>F/<E>E  validate-failures=<n>  units-issues=<n>
# which is the P/F/E to quote in `stage2: … — <P>P/<F>F/<E>E` commit messages.
# The exit status is nonzero if any validation, units check, or test fails.
#
# Run it inside a slurm job, not on a login node.
set -uo pipefail

repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
main=/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM
export ESD_ROOT=${ESD_ROOT:-/projects/illinois/eng/cee/ctessum/ctessum/code/EarthSciDiscretizations}

esm=./esm
[ -x "$esm" ] || esm="$repo/esm"
[ -x "$esm" ] || esm="$main/esm"
[ -x "$esm" ] || { echo "esm-suite: no esm binary (copy $main/esm here)" >&2; exit 2; }

run_tests=1
if [ "${1:-}" = "--no-test" ]; then run_tests=0; shift; fi
if [ $# -eq 0 ]; then set -- stage2/*.esm; fi

vfail=0; ufail=0
printf '%-44s %-10s %s\n' FILE VALIDATE UNITS
for f in "$@"; do
    if "$esm" validate "$f" >/dev/null 2>&1; then v=ok; else v=FAIL; vfail=$((vfail+1)); fi
    rep=$("$esm" units "$f" 2>&1); chk=$("$esm" units --check "$f" 2>&1); urc=$?
    nlit=$(grep -c 'units_dimensionless_literal_mismatch' <<<"$rep")
    nmis=$(grep -oE '[0-9]+ mismatched' <<<"$chk" | grep -oE '^[0-9]+' || echo 0)
    if [ "$urc" -ne 0 ] || [ "${nmis:-0}" -ne 0 ] || [ "$nlit" -ne 0 ]; then
        u="${nmis:-?} mismatched, $nlit literal"; ufail=$((ufail+1))
    else
        u=ok
    fi
    printf '%-44s %-10s %s\n' "$f" "$v" "$u"
done

P=0; F=0; E=0; trc=0
if [ "$run_tests" -eq 1 ]; then
    out=$("$esm" test "$@" 2>&1); trc=$?
    echo "$out" | sed -n '/^ *-\{5,\}/,$p' | head -n 200
    read -r P F E < <(echo "$out" | awk '$1=="TOTAL"{print $2, $3, $4}')
    P=${P:-0}; F=${F:-0}; E=${E:-0}
    if [ -z "$(echo "$out" | awk '$1=="TOTAL"')" ]; then
        echo "$out" | tail -n 20
        echo "esm-suite: no TOTAL row from esm test (exit $trc)" >&2
        E=${E:-0}; trc=${trc:-1}; [ "$trc" -eq 0 ] && trc=1
    fi
fi

echo "SUITE: ${P}P/${F}F/${E}E  validate-failures=$vfail  units-issues=$ufail"
[ "$vfail" -eq 0 ] && [ "$ufail" -eq 0 ] && [ "$F" -eq 0 ] && [ "$E" -eq 0 ] && [ "$trc" -eq 0 ]
