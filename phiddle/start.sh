#!/bin/bash
# Start the CrystalShiftAPI Julia backend, wait for it to come up, then launch the GUI.
# Works from any directory (including double-clicking). Arguments are passed to phiddle.py.

cd "$(dirname "$0")" || exit 1

PORT=8080
TIMEOUT=300  # seconds; first launch after a Julia update can spend minutes precompiling

server_up() { curl -s -o /dev/null "http://127.0.0.1:$PORT/phase_names"; }

if server_up; then
    # Leave an already-running server alone; we didn't start it, so we don't stop it.
    echo "Using CrystalShiftAPI already running on port $PORT."
else
    julia --project=../CrystalShiftAPI --threads=6 ../CrystalShiftAPI/src/CrystalShiftAPI.jl &
    pid=$!
    trap 'kill $pid 2>/dev/null' EXIT INT TERM

    echo "Waiting for CrystalShiftAPI on port $PORT..."
    for ((i = 0; i < TIMEOUT; i++)); do
        server_up && break
        if ! kill -0 $pid 2>/dev/null; then
            echo "CrystalShiftAPI exited during startup." >&2
            exit 1
        fi
        sleep 1
    done
    if ((i == TIMEOUT)); then
        echo "CrystalShiftAPI did not respond within ${TIMEOUT}s." >&2
        exit 1
    fi
fi

# Prefer uv (uses the locked env in ../.venv); fall back to whatever python is active.
if command -v uv >/dev/null 2>&1; then
    uv run python phiddle.py "$@"
else
    python phiddle.py "$@"
fi
