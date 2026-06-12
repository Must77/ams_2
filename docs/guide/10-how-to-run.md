# 10 操作手册：构建、回归、运行、工具

> 前置：04（不读也能照抄命令，但读了才看得懂输出）。
> 环境：MaginaLaptop，Ubuntu 26.04，PETSc 3.24（apt 安装，带 Hypre），项目在 `~/ams_2`。

## 1. 依赖与构建

```bash
# 依赖（已装好；新机器重来一遍即可）
sudo apt install petsc-dev libmetis-dev libeigen3-dev libnanoflann-dev gfortran cmake g++ zstd

# 配置 + 编译（产物在 build_ubuntu/）
cd ~/ams_2
cmake -S src_ams -B build_ubuntu -DCMAKE_BUILD_TYPE=Release \
  -DEIGEN_DIR=/usr -DMETIS_DIR=/usr -DNANOFLANN_DIR=/usr \
  -DPETSC_DIR=/usr/lib/petsc -DPETSC_ARCH=
cmake --build build_ubuntu -j$(nproc)
```

产物：`FemAms`（旧文件式入口）、`AdapterSelfCheck`（单元测试）、
`UpperExampleCheck`（大算例检查）、`UpperDumpSystem` + `AmsLab`（诊断工具）。

## 2. 改完代码必跑的两个回归（30 秒）

```bash
# ① 单元测试：1基转换/上三角展开/符号逻辑（静默+退出码0=通过）
cd ~/ams_2 && ./build_ubuntu/AdapterSelfCheck && echo PASS

# ② eg_1 黄金回归：结果必须逐字节复现这个哈希
cd ~/ams_2/exmaple/eg_1 && ../../build_ubuntu/FemAms > /dev/null 2>&1
sha256sum result.txt
# 必须是: 342042dbcff12d00a13339367327a851158f0dc137f6af5caa1ee6779e9fd266
```

哈希变了 = 改动影响了默认路径 = 停下来想清楚。
（注意 ① 会顺带修改 eg_1 目录下的日志文件，提交前 `git checkout` 还原之。）

## 3. 大算例数据准备（已做过；换机器时参考）

```bash
# 压缩数据在仓库 example/ams_example_zst/（dev 分支），校验+解压：
cd ~/ams_2/example/ams_example_zst
sha256sum -c SHA256SUMS
mkdir -p ~/ams_example
for f in *.zst; do zstd -d -f --long=27 "$f" -o ~/ams_example/"$(basename "$f" .zst)"; done
# 解压后 1.7 GB
```

## 4. 标准运行命令与输出解读

```bash
cd ~/ams_2
./build_ubuntu/UpperExampleCheck ~/ams_example \
  -fortran_upper_unscale_edge_len \
  -ams_beta_mass_poisson \
  -ams_beta_mass_shift 1e-6 \
  -B_ksp_type gmres -B_ksp_max_it 15 \
  -em_outer_max_it 60 -em_outer_rtol 1e-7
```

单核约 35 分钟（前 ~5 分钟是 1.7G 文本解析，静默，别以为卡死了）。输出解读：

```text
reference_relative_residual 6.81e-07   ← slove.txt 代回方程的残差。这是数据自洽性
                                          检查，和我们的求解器无关；它变了说明数据坏了
solve_relative_l2 …                    ← 与 slove.txt 的全局 L2 差。因规范自由度
                                          (08 篇)，预期 ~0.69，**这不是失败指标**
ams_solution.txt                       ← 解向量落盘（每行: 实部 虚部），真正的验收
                                          用它做代回残差/分区域比对（见第 6 节）
```

## 5. 选项速查表

| 选项 | 默认 | 作用 / 何时动它 |
| --- | --- | --- |
| `-fortran_upper_unscale_edge_len` | 关 | 修复 1：对方数据必开；Whitney 约定数据必关 |
| `-ams_beta_mass_poisson` | 关 | 修复 2：8 位精度数据必开；16 位数据可试关 |
| `-ams_beta_mass_shift <ε>` | 1e-6 | 修复 3：1e-8~1e-2 等效，基本不用动 |
| `-B_ksp_type gmres` | fcg | 内层算法：标准命令固定 gmres（FCG 有浮点误报坑） |
| `-B_ksp_max_it <n>` | 100 | 内层步数上限：15 是性价比点（外层多走几步换内层不磨洋工） |
| `-B_ksp_rtol <t>` | 1e-2 | 内层容差：一般不动 |
| `-em_outer_max_it <n>` | 100 | 外层上限（本轮新增选项） |
| `-em_outer_rtol <t>` | 1e-9 | 外层容差：8 位数据下 1e-9 不可达，用 1e-7 |
| `-A_ksp_converged_reason` | - | 打印外层收敛原因，排查必备 |
| `-B_ksp_converged_reason` | - | 打印每次内层收敛原因（输出多，短跑用） |

## 6. 验收数字怎么复算（不依赖程序自报）

```python
# python3, 需 numpy/scipy (apt: python3-numpy python3-scipy)
# 加载矩阵/rhs/两个解 → 代回残差 + 分区域比对
# (完整脚本逻辑见 06 篇 9 节, 当时版本在远程 /tmp/solution_forensics.py)
import numpy as np, scipy.sparse as sp
D = "/home/magina/ams_example/"
jcol = np.loadtxt(D+"jcol.txt", skiprows=1, dtype=np.int64)[:,1]
irw  = np.loadtxt(D+"irw.txt",  skiprows=1, dtype=np.int64)[:,1]
k1   = np.loadtxt(D+"ssor_k1.txt", skiprows=1)[:,1]
k2   = np.loadtxt(D+"ssor_k2.txt", skiprows=1)[:,1]
n = len(jcol)-1
rr = np.repeat(np.arange(n), np.diff(jcol)); cc = irw-1
def full(v):
    U = sp.coo_matrix((v,(rr,cc)),shape=(n,n)).tocsr()
    return U + U.T - sp.diags(U.diagonal())
A = full(k1) + 1j*full(k2)
rh = np.loadtxt(D+"right-hand.txt", skiprows=1); b = rh[:,1] + 1j*rh[:,2]
x = np.loadtxt("/home/magina/ams_2/ams_solution.txt"); x = x[:,0] + 1j*x[:,1]
print("代回相对残差:", np.linalg.norm(A@x-b)/np.linalg.norm(b))   # 预期 ~7e-8
```

## 7. 诊断工具：UpperDumpSystem + AmsLab

预条件器/参数实验**不要**用 UpperExampleCheck 反复跑（每轮 5 分钟读文本）。正确姿势：

```bash
# 一次性: 把内层系统 (B, M, G, 坐标, rhs) 存成 PETSc 二进制 (~2.5 分钟)
./build_ubuntu/UpperDumpSystem ~/ams_example ~/ams_dump_whitney -fortran_upper_unscale_edge_len

# 之后: 每个配置数秒出结果
./build_ubuntu/AmsLab ~/ams_dump_whitney -lab_beta_mass -lab_mass_shift 1e-6 \
    -lab_ksp_type gmres                      # 复现"3 步收敛"
./build_ubuntu/AmsLab ~/ams_dump_whitney -lab_pc jacobi      # 对照: 烂预条件器
./build_ubuntu/AmsLab ~/ams_dump_whitney -lab_pc none        # 对照: 无预条件
# 输出一行: reason CONVERGED_RTOL its 3 rnorm ... xnorm ...
```

AmsLab 专属选项：`-lab_pc ams|jacobi|none`、`-lab_beta_mass`（修复 2）、
`-lab_mass_shift <ε>`（修复 3）、`-lab_random_rhs`（随机右端项；注意随机 rhs 会激发
近奇异模态，FCG 必炸，属预期）、其余 `-lab_ksp_*` 透传 PETSc。

## 8. 已知小坑

- 迭代日志文件名带乱码前缀（`*group-000.log`）：历史遗留的未初始化字符串，无害，
  已加入 .gitignore；
- 日志内外层行混杂、缩进动态变化，**勿盲目脚本解析**（04 篇 5 节）；
- `Solving for dual mode` 会把同一系统解两遍（CSEM 多源结构遗留），耗时翻倍，
  是已知的待优化项；
- 远程写 C++ 文件时若走 ssh 单引号包裹，代码里的 `'/'` 等字符会破坏引号配对，
  用 base64 中转。
