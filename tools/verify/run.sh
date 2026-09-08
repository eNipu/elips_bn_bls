#!/usr/bin/env bash
# Run every independent verification. Works locally and in CI.
#
#   tools/verify/run.sh              everything
#   tools/verify/run.sh params       PARI curve parameters only
#   tools/verify/run.sh proofs       symbolic proofs only
#   tools/verify/run.sh crosscheck   py_ecc pairing only
#
# See tools/verify/README.md for what each one proves, and what it does not.
set -u

cd "$(dirname "$0")/../.." || exit 1
WANT="${1:-all}"
VENV="${ELIPS_VERIFY_VENV:-build/verify-venv}"
rc=0

hdr() { printf '\n\033[1m== %s ==\033[0m\n' "$1"; }
miss() { printf '  SKIPPED: %s\n' "$1"; }

run_gp() {
    local name=$1 file=$2
    hdr "$name"
    if ! command -v gp >/dev/null 2>&1; then
        miss "PARI/GP not found. apt-get install --no-install-recommends pari-gp,
           or use the container: tools/verify/run-docker.sh"
        return 0   # absence of the tool is not a failed proof
    fi
    gp -q < "$file" || rc=1
}

# if-blocks rather than [ ] || [ ] && cmd: that chain returns the status of the
# last test when it does not run the command, which is a good way to lose an
# exit code.
if [ "$WANT" = all ] || [ "$WANT" = params ]; then
    run_gp "Curve parameters, counted independently by PARI" tools/verify/curve_params.gp
fi
if [ "$WANT" = all ] || [ "$WANT" = proofs ]; then
    run_gp "Symbolic proofs in Z[X], for every seed" tools/verify/identities.gp
fi

if [ "$WANT" = all ] || [ "$WANT" = crosscheck ]; then
    hdr "BLS12-381 pairing against py_ecc"
    BIN=build/dump_pairing
    if [ ! -x "$BIN" ]; then
        printf '  building %s\n' "$BIN"
        cmake -S . -B build -DCMAKE_BUILD_TYPE=Release >/dev/null 2>&1
        cmake --build build --target dump_pairing -j"$(nproc 2>/dev/null || echo 4)" >/dev/null 2>&1
    fi
    if [ ! -x "$BIN" ]; then
        miss "could not build $BIN (is GMP installed?)"
    else
        PY=python3
        if ! $PY -c 'import py_ecc' >/dev/null 2>&1; then
            if [ ! -d "$VENV" ]; then
                printf '  creating venv at %s\n' "$VENV"
                python3 -m venv "$VENV" >/dev/null 2>&1
            fi
            if [ -x "$VENV/bin/python" ]; then
                "$VENV/bin/pip" install -q -r tools/verify/requirements.txt >/dev/null 2>&1
                PY="$VENV/bin/python"
            fi
        fi
        if $PY -c 'import py_ecc' >/dev/null 2>&1; then
            $PY tools/verify/crosscheck_pyecc.py "$BIN" || rc=1
        else
            miss "py_ecc unavailable and the venv could not be created.
           pip install -r tools/verify/requirements.txt"
        fi
    fi
fi

printf '\n'
if [ "$rc" -eq 0 ]; then printf 'verification: everything that ran, passed\n'
else printf 'verification: FAILURES above\n'; fi
exit "$rc"
