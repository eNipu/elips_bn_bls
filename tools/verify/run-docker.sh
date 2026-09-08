#!/usr/bin/env bash
# Run the verification in a clean container, so the result does not depend on
# what happens to be installed on this machine.
set -eu
cd "$(dirname "$0")/../.."
docker build -q -t elips-verify -f tools/verify/Dockerfile . >/dev/null
exec docker run --rm elips-verify
