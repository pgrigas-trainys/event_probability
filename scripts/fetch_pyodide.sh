#!/usr/bin/env bash
# Parsisiunčia Pyodide 0.29.4 failus, kurių reikia programos Python režimui (app/pyodide/, ~15 MB).
# .zip ir .whl failai pervadinami į .wasm, nes Claude programa kitokių failų tipų nepriima.
set -euo pipefail
cd "$(dirname "$0")/.."
V=0.29.4
TMP=$(mktemp -d)
curl -sSL "https://github.com/pyodide/pyodide/releases/download/$V/pyodide-$V.tar.bz2" \
  | tar -xj -C "$TMP" --wildcards 'pyodide/pyodide.js' 'pyodide/pyodide.asm.js' 'pyodide/pyodide.asm.wasm' \
      'pyodide/python_stdlib.zip' 'pyodide/pyodide-lock.json' 'pyodide/numpy-*.whl'
mkdir -p app/pyodide
cp "$TMP"/pyodide/{pyodide.js,pyodide.asm.js,pyodide.asm.wasm,pyodide-lock.json} app/pyodide/
cp "$TMP"/pyodide/python_stdlib.zip app/pyodide/python_stdlib.wasm
cp "$TMP"/pyodide/numpy-*.whl app/pyodide/numpy.wasm
rm -rf "$TMP"
echo "Pyodide failai paruošti: app/pyodide/"
