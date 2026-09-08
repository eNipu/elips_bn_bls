# Reproducible environment for the independent verification.
#
# The point of tools/verify/ is that it does not share anything with the
# library. That argument is weaker if the tool versions drift, so this pins the
# base image and the py_ecc version, and installs PARI/GP from the distribution
# rather than building it here.
#
#   docker build -t elips-verify -f tools/verify/Dockerfile .
#   docker run --rm elips-verify
#
# NOTE ON HOW FAR THIS WAS TESTED. The image build itself has not been run:
# the environment it was written in has the docker client but no daemon. What
# WAS verified there is every ingredient -- the apt set resolves, the pinned
# py_ecc is the version the cross-check actually passes against, and
# tools/verify/run.sh passes on an ubuntu:24.04 userland with exactly these
# packages. The CI job runs run.sh directly for the same reason, so a broken
# image cannot silently stop the verification from happening.
#
# or just:  tools/verify/run-docker.sh
FROM ubuntu:24.04

# --no-install-recommends matters: pari-gp otherwise pulls in an X11 plotting
# stack, and one of those packages 404s on the current archive.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential cmake pkg-config libgmp-dev \
        pari-gp python3 python3-pip python3-venv ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY tools/verify/requirements.txt /tmp/requirements.txt
RUN pip3 install --no-cache-dir --break-system-packages -r /tmp/requirements.txt

COPY . /src
RUN cmake -S . -B build -DCMAKE_BUILD_TYPE=Release >/dev/null \
    && cmake --build build --target dump_pairing -j"$(nproc)" >/dev/null

CMD ["tools/verify/run.sh"]
