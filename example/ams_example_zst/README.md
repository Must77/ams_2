# ams_example compressed data

This directory stores the compressed copy of the external AMS example data that
was provided under:

```text
/home/akama/ams/ams_example
```

The files are kept on the `dev` branch as debugging and handoff data. They are
not on `main` because the solver behavior is still under investigation.

## Contents

Each source file is compressed independently with `zstd`:

```bash
zstd -19 -T0 --long=27 <source-file> -o <source-file>.zst
```

The original dataset contains:

| File | Meaning |
| --- | --- |
| `CSR_I.txt.zst` | Full CSR row pointer, 1-based. |
| `CSR_J.txt.zst` | Full CSR column indices, 1-based. |
| `csr_K1.txt.zst` | Full CSR real values. |
| `csr_K2.txt.zst` | Full CSR imaginary values. |
| `jcol.txt.zst` | Upper-triangle CSR row pointer, 1-based. |
| `irw.txt.zst` | Upper-triangle CSR column indices, 1-based. |
| `ssor_k1.txt.zst` | Upper-triangle real values. |
| `ssor_k2.txt.zst` | Upper-triangle imaginary values. |
| `edges.txt.zst` | Edge id and endpoint node ids, 1-based. |
| `node.txt.zst` | Node id and 3D coordinates, 1-based ids. |
| `right-hand.txt.zst` | Right-hand side data. |
| `slove.txt.zst` | Direct-solver reference solution, real and imaginary columns. |
| `main.f90.zst` | Caller example from the collaborator. |

## Verify

```bash
cd example/ams_example_zst
sha256sum -c SHA256SUMS
for f in *.zst; do zstd -t "$f"; done
```

## Decompress

To restore the dataset into a separate working directory:

```bash
mkdir -p /tmp/ams_example
for f in example/ams_example_zst/*.zst; do
  zstd -d -f "$f" -o "/tmp/ams_example/$(basename "$f" .zst)"
done
```

## Size note

The largest compressed files are `csr_K1.txt.zst` and `csr_K2.txt.zst`.
They are below GitHub's 100 MiB hard limit for ordinary Git files, but above the
50 MiB warning threshold. Keep future replacements compressed and verify file
sizes before pushing.
