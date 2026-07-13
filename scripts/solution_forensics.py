#!/usr/bin/env python3
"""解代回残差法证：验收迭代解的唯一可靠标准。

用法: python3 scripts/solution_forensics.py <data_dir> <solution.txt> [more_solutions...]

对每个解文件计算 ||K z - S|| / ||S||（K = Kr + i*Ki 由上三角 CSR 对称展开，
复数方程 K z = S 与 upper_example_check.cpp 的 reference 残差定义等价）。
同时给出各解之间的差异分解：低 σ 边（|diag Ki| 相对极小，空气层）占比——
与 slove.txt 的全局 L2 差集中在低 σ 边属规范自由度假象，不是错误
（见 docs/project_handoff.md 2026-06-13 节）。

数据文件格式（与 UpperExampleCheck 相同）：
  jcol.txt   上三角 row_ptr（1-based，带表头）
  irw.txt    上三角 col_idx（1-based，带表头）
  ssor_k1/ssor_k2.txt  实部/虚部值（带表头，"idx value" 两列）
  right-hand.txt       右端项（带表头，"idx re im x x" 五列）
解文件：每行 "re im" 两列，无表头（ams_solution.txt / slove.txt 格式）。
"""
import sys

import numpy as np
import scipy.sparse as sp


def load_indexed(path, usecols):
    return np.loadtxt(path, skiprows=1, usecols=usecols)


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    data_dir, sol_paths = sys.argv[1].rstrip("/"), sys.argv[2:]

    rp = load_indexed(f"{data_dir}/jcol.txt", 1).astype(np.int64) - 1
    ci = load_indexed(f"{data_dir}/irw.txt", 1).astype(np.int64) - 1
    k1 = load_indexed(f"{data_dir}/ssor_k1.txt", 1)
    k2 = load_indexed(f"{data_dir}/ssor_k2.txt", 1)
    rhs = load_indexed(f"{data_dir}/right-hand.txt", (1, 2))
    n = len(rp) - 1
    S = rhs[:, 0] + 1j * rhs[:, 1]

    rows = np.repeat(np.arange(n), np.diff(rp))
    upper = sp.coo_matrix((k1 + 1j * k2, (rows, ci)), shape=(n, n)).tocsr()
    strict = sp.triu(upper, k=1)
    K = upper + strict.T  # 对称展开

    ki_diag = np.abs(upper.diagonal().imag)
    low_sigma = ki_diag < 1e-6 * ki_diag.max()
    print(f"n_edges {n}  nnz_upper {len(k1)}  "
          f"low_sigma_edges {low_sigma.sum()} ({100 * low_sigma.mean():.1f}%)")

    s_norm = np.linalg.norm(S)
    sols = []
    for path in sol_paths:
        z2 = np.loadtxt(path)
        z = z2[:, 0] + 1j * z2[:, 1]
        relres = np.linalg.norm(K @ z - S) / s_norm
        print(f"{path}: true relres {relres:.4e}")
        sols.append((path, z))

    for i in range(len(sols)):
        for j in range(i + 1, len(sols)):
            (pa, za), (pb, zb) = sols[i], sols[j]
            d = za - zb
            rel = np.linalg.norm(d) / np.linalg.norm(zb)
            frac_low = (np.linalg.norm(d[low_sigma]) / np.linalg.norm(d)) ** 2
            hi = ~low_sigma
            rel_hi = np.linalg.norm(d[hi]) / np.linalg.norm(zb[hi])
            print(f"diff {pa} vs {pb}: rel_l2 {rel:.4e}  "
                  f"diff_energy_in_low_sigma {frac_low:.4f}  "
                  f"rel_l2_high_sigma_only {rel_hi:.4e}")


if __name__ == "__main__":
    main()
