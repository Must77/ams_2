# AMS solver project handoff

本文档面向 0 上下文接手的开发者或 agent，记录项目背景、当前代码状态、数据格式、运行方式、已验证证据、当前问题和建议下一步。

## 项目背景

对方原有地球物理程序中存在一个直接求解器。直接求解器结果可信，但内存消耗过大，因此希望改用 AMS 迭代求解器。

当前仓库中的 C/C++ 版本可以理解为一个 PETSc/Hypre AMS 求解器接口原型。对方后续希望它能被 C/Fortran 像普通函数一样调用，并希望输入数据从内存数组传入，而不是从文件系统读取。

当前迁移目标不是立刻重写成 Fortran，而是先把 C++ adapter 跑通，确认输入格式、符号约定、AMS 所需网格数据、求解结果都正确，再考虑绑定或重写。

## 代码位置

- 代码仓库：`/home/akama/ams/ams_2`
- 对方新算例数据：`/home/akama/ams/ams_example`
- 旧自带回归算例：`/home/akama/ams/ams_2/exmaple/eg_1`
- 主要源码目录：`/home/akama/ams/ams_2/src_ams`

重要文件：

| 文件 | 作用 |
| --- | --- |
| `src_ams/main.cpp` | 旧文件入口。读取当前目录下的 `edgesN.txt/nodes.txt/rowPtr.txt/colIdx.txt/rhs.txt/val_real_imag.txt`，调用 `solve_eg1()`，写 `result.txt`。 |
| `src_ams/PetscAMGInterface.h` | 暴露 solver、adapter、结果写出等函数声明。 |
| `src_ams/PetscAMGInterface.cc` | 核心实现：输入数组加载、AMS 初始化、矩阵/右端项构造、MatShell 乘法、KSP 求解、结果拷贝。 |
| `src_ams/adapter_self_check.cpp` | 小型 adapter 回归检查。 |
| `src_ams/upper_example_check.cpp` | 读取 `ams_example` 上三角 1 基数据并调用上三角 adapter 的检查程序。 |
| `src_ams/CMakeLists.txt` | 构建 `FemAms`、`AdapterSelfCheck`、`UpperExampleCheck`。 |

## 当前代码状态

当前代码已经从“`solve_eg1()` 自己读文件并写结果”改为“调用者准备数组并接收输出数组”：

```text
main.cpp
  read_eg1_input_files()
  solve_eg1(... arrays ..., out_real, out_imag)
  write_result_arrays(...)
```

已经加入的关键接口：

| 函数 | 作用 |
| --- | --- |
| `solve_eg1()` | 当前主求解函数。假设 PETSc 已由调用侧初始化；输入为数组；输出写到 `out_real/out_imag`。 |
| `load_context_from_arrays()` | 把调用者数组加载到 `EMContext`。 |
| `copy_result_to_arrays()` | 从 PETSc 解向量拷贝到调用者输出数组。 |
| `write_result_arrays()` | 把输出数组写成旧格式 `result.txt`。 |
| `expand_upper_triangle_1based_to_full_csr()` | 把 Fortran 1 基上三角 CSR 展开为当前 solver 使用的 0 基全量 CSR。 |
| `prepare_fortran_upper_1based_inputs()` | 处理 Fortran/对方输入：1 基边节点、1 基上三角 CSR、符号转换、rhs 转换。 |
| `solve_eg1_fortran_upper_1based()` | 当前上三角 adapter 入口。内部调用 `prepare_fortran_upper_1based_inputs()` 和 `solve_eg1()`。 |

注意：还没有做 Fortran 绑定；当前只是 C++ adapter。

## 依赖和构建

当前 `build_ubuntu` 是在 Ubuntu 包路线下构建出来的。`CMakeCache.txt` 中的依赖路径为：

```text
EIGEN_DIR=/usr
METIS_DIR=/usr
NANOFLANN_DIR=/usr
PETSC_DIR=/usr/lib/petsc
PETSC_ARCH=
```

可执行文件：

```text
build_ubuntu/FemAms
build_ubuntu/AdapterSelfCheck
build_ubuntu/UpperExampleCheck
```

典型构建命令：

```bash
cd /home/akama/ams/ams_2
cmake -S src_ams -B build_ubuntu \
  -DEIGEN_DIR=/usr \
  -DMETIS_DIR=/usr \
  -DNANOFLANN_DIR=/usr \
  -DPETSC_DIR=/usr/lib/petsc \
  -DPETSC_ARCH=
cmake --build build_ubuntu -j
```

依赖含义：

| 依赖 | 作用 |
| --- | --- |
| PETSc | Krylov 迭代器、向量/矩阵、PC 接口；当前代码通过 PETSc 调 Hypre AMS。 |
| Hypre AMS | 真正的 AMS 预条件器，通过 PETSc `PCHYPRESetType(..., "ams")` 使用。 |
| MPI | PETSc 构建和运行需要。当前测试主要用单进程。 |
| Eigen | 代码依赖的头文件库。 |
| METIS | 当前 CMake 链接依赖。 |
| nanoflann | 头文件库，当前 CMake 要求存在 `nanoflann.hpp`。 |

## 运行方式

### 旧自带算例回归

```bash
cd /home/akama/ams/ams_2/exmaple/eg_1
../../build_ubuntu/FemAms
```

此前已验证：当前 `FemAms` 输出与原始 `result.txt` 字节级一致，sha256 都是：

```text
342042dbcff12d00a13339367327a851158f0dc137f6af5caa1ee6779e9fd266
```

这说明“把 `solve_eg1()` 改成数组输入输出”没有破坏旧自带算例。

### Adapter 小回归

```bash
cd /home/akama/ams/ams_2
./build_ubuntu/AdapterSelfCheck
```

用途：检查上三角展开和符号转换的基本逻辑。

### 对方上三角算例

```bash
cd /home/akama/ams/ams_2
./build_ubuntu/UpperExampleCheck /home/akama/ams/ams_example
```

这个程序读取：

```text
jcol.txt
irw.txt
ssor_k1.txt
ssor_k2.txt
edges.txt
node.txt
right-hand.txt
slove.txt
```

它会先计算 `slove.txt` 代回 A/b 的残差，再调用 `solve_eg1_fortran_upper_1based()` 求解并比较结果。

注意：对方上三角算例目前没有收敛到 `slove.txt`；不要把它当作已通过测试。

## 数据文件语义

### 旧 `FemAms` 文件入口需要的文件

| 文件 | 含义 | 基准 |
| --- | --- | --- |
| `edgesN.txt` | 每条边连接的两个 node id，用于构造 AMS 离散梯度矩阵 G。 | 当前旧入口期望 0 基。 |
| `nodes.txt` | 节点坐标，每行 x/y/z，用于 `PCSetCoordinates()`。 | 无 id 列。 |
| `rowPtr.txt` | CSR 行指针。 | 0 基。 |
| `colIdx.txt` | CSR 列编号。 | 0 基。 |
| `val_real_imag.txt` | CSR 每个非零元的实部和虚部。 | 与 `colIdx` 一一对应。 |
| `rhs.txt` | 右端项实部和虚部。 | 与边/未知量顺序一致。 |
| `result.txt` | 输出解向量，旧格式为 `i real imag`。 | `i` 为输出序号。 |

### 对方 `ams_example` 文件

| 文件 | 含义 | 当前确认 |
| --- | --- | --- |
| `node.txt` | 节点坐标，头部为 `97901 3 0 0`，后续为 `node_id x y z`。 | node id 是 1 基。 |
| `edges.txt` | 边和两端节点，头部为边数，后续为 `edge_id node1 node2`。 | edge id 和 node id 都是 1 基。 |
| `right-hand.txt` | 右端项，头部为 `nedge: 682344`，后续含 id、实部、虚部等列。 | 当前使用第 2/3 列作为 rhs 实部/虚部。 |
| `slove.txt` | 对方给出的参考解。 | 第一列实部，第二列虚部；来自直接求解器。 |
| `jcol.txt` | 对方上三角 CSR 的行指针。 | 1 基，长度 682345，最后一个值 5932111。 |
| `irw.txt` | 对方上三角 CSR 的列编号。 | 1 基，值数量 5932110。 |
| `ssor_k1.txt` | 上三角矩阵值实部。 | 与 `irw` 一一对应。 |
| `ssor_k2.txt` | 上三角矩阵值虚部。 | 与 `irw` 一一对应。 |
| `CSR_I.txt` | 对方全量 CSR 行指针。 | 1 基。 |
| `CSR_J.txt` | 对方全量 CSR 列编号。 | 1 基。 |
| `csr_K1.txt` | 全量 CSR 值实部。 | 旧版本导出错误为全 0；新版本已替换。 |
| `csr_K2.txt` | 全量 CSR 值虚部。 | 旧版本导出错误为全 0；新版本已替换。 |

## 符号约定

当前代码中 `matshell_mult_a()` 实际实现的块矩阵乘法是：

```text
yr = C*xr - M*xi
yi = -M*xr - C*xi
```

所以内部块矩阵是：

```text
[  C   -M ]
[ -M   -C ]
```

对方物理复数系统可理解为 `A = Kr + i Ki`。为了让当前 cpp 内部解出等价物理系统，adapter 使用如下转换：

```text
C = Kr
M = -Ki
rhs_real_internal = br
rhs_imag_internal = bi
internal_xi = -physical_imag_solution
```

因此 `solve_eg1_fortran_upper_1based()` 在返回前会把内部虚部再取负，转换回调用者看到的物理虚部。

这个符号问题目前不是数学推导一定错误，而是当前代码实现和对方物理复数记号之间存在表示约定差异。迁移时必须保留这个事实，不要随意改 `matshell_mult_a()` 的符号，除非同步重新验证所有数据和参考解。

## 已验证证据

### 旧自带算例

当前数组化改造后，旧 `exmaple/eg_1` 输出和原结果完全一致：

```text
sha256(result.txt) = 342042dbcff12d00a13339367327a851158f0dc137f6af5caa1ee6779e9fd266
cmp_exit = 0
```

### 对方新全量 CSR 与上三角一致

对方后来重新给了 `csr_K1.txt/csr_K2.txt`，修复了旧文件全 0 的导出错误。已验证新全量 CSR 和上三角数据展开后是一致的：

```text
checked 11181876
missing 0
mismatch 0
maxdiff_k1 0.0
maxdiff_k2 9.999999881247492e-29
```

这说明“上三角展开成全量 CSR”这一步基本不是当前主要问题。

### `slove.txt` 与新 A/b 属于同一组方程

用新全量 CSR 和 rhs 代回 `slove.txt`，残差为：

```text
full_csr_reference_relative_residual 6.811145888676909e-07
full_csr_reference_max_abs_residual 8.649358490405724e-16 at 675215
```

上三角路径代回也得到同量级残差：

```text
reference_relative_residual 6.8111458886769823e-07
reference_max_abs_residual 8.6493584904057254e-16
```

这说明 `slove.txt`、新全量 A、上三角 A、rhs 基本是同一套数据。

### 原 cpp main 跑新全量 CSR 的结果

为了按“先让 cpp 跑他们给的新全量数据”验证，曾生成两个临时运行目录：

```text
run_full_csr_physical/
run_full_csr_internal/
```

这两个目录只包含转换后的旧格式输入、日志和运行结果，已经删除。

结果：

1. 原物理符号直接喂给旧 `FemAms`：

```text
result_nonzero_rows 0
relative_l2_vs_slove 1.0
max_abs_vs_slove 8.800695677823685e-11 at 339940
```

输出全 0，不匹配。

2. 按当前 cpp 内部符号转换后再喂给旧 `FemAms`：

运行进入正常长迭代，但按用户要求在不收敛时停止。停止前最新摘要：

```text
outer_A_entries 42
B_iter100_lines 60
A 32 rel 5.907183170157e-06
A 33 rel 5.826045509364e-06
A 34 rel 4.985762755140e-06
A 35 rel 4.907896118409e-06
A 36 rel 4.839645176220e-06
A 37 rel 4.833220361010e-06
A 38 rel 4.528659275827e-06
A 39 rel 4.474414029880e-06
A 40 rel 4.239496511327e-06
A 41 rel 4.133048511654e-06
```

目标 `rtol` 是 `1e-9`。停掉时还没有生成最终 `result.txt`。

此前完整上三角 adapter 运行到 `maxiter=100` 时也没有收敛：

```text
Linear A_ solve did not converge due to DIVERGED_ITS iterations 100
final true residual ratio about 1.893e-6
solve_relative_l2 4.8496141735955334e-01
solve_max_abs 8.2725335557340852e-11 at 339939
```

## 当前核心判断

证据较强：

1. 新全量 CSR 数据本身已经不像是问题。
2. 上三角展开逻辑已经不像是问题。
3. `slove.txt` 确实能代回新 A/b，说明参考解和矩阵/右端项基本匹配。
4. 当前 cpp AMS 路径无法复现 `slove.txt`，问题更可能在 AMS 所需辅助数据或求解配置。

仍是推测但优先级高：

1. A 矩阵第 i 行是否严格对应 `edges.txt` 中 edge id=i，仍是关键问题。AMS 的离散梯度矩阵 G 是从 `edges.txt/node.txt` 构造的，如果 A 的未知量顺序和 edge 顺序不一致，直接求解器仍可解 `A*x=b`，但 AMS 预条件器会失效或表现很差。
2. 对方数据中可能存在边界条件消元、自由度压缩、重排、方向约定等前处理；目前 cpp 假设 A/rhs/edge 都是同一套未重排顺序。
3. PETSc/Hypre 参数可能和对方个人电脑上“十分钟左右跑完”的 AMS 程序不同。
4. 直接解 `slove.txt` 的容差未知。它可以作为参考，但不能证明当前 AMS 设置应该达到完全相同的误差。

## 需要继续向对方确认的问题

优先级最高：

1. A 矩阵第 i 行/列是否就是 `edges.txt` 的 edge id=i？
2. 是否有对边未知量做过重排、压缩、删除边界自由度、合并自由度？
3. `right-hand.txt` 的行顺序是否也和 edge id 完全一致？
4. 他们个人电脑上十分钟跑的 cpp AMS 版本，具体命令、参数、PETSc/Hypre 版本是什么？
5. 他们原程序中传给 AMS 的 G/coords 是否就是由同一份 `edges.txt/node.txt` 构造？

次优先级：

1. `slove.txt` 直接求解器容差是多少？
2. 原直接求解器和 AMS 求解器是否使用完全相同的复数方程符号？
3. 是否有矩阵缩放、对角预处理、SSOR 预处理或额外归一化？
4. `ssor_k1/ssor_k2` 文件名中的 `ssor` 是否只是文件名，还是表示数据已经被 SSOR 相关流程处理过？

## 建议下一步

不要立刻做 Fortran 绑定，也不要先重写 Fortran 版本。建议顺序：

1. 固化一个正式的 full CSR 检查程序。
   - 读取 `CSR_I/CSR_J/csr_K1/csr_K2/right-hand/slove`。
   - 验证 `slove.txt` 残差。
   - 可选：调用当前 `solve_eg1()`，但要把输出比较写成机器可读摘要。

2. 固化一个 edge/order 检查程序。
   - 检查 `edges.txt` 是否严格 1 基连续。
   - 检查 node 引用范围。
   - 检查 A 行数、rhs 行数、edge 数是否一致。
   - 如果对方能提供 edge 重排表，加入重排验证。

3. 深查 `setup_ams()`。
   - 它是从 `edgesN/nodes` 构造 Hypre AMS 需要的离散梯度矩阵 G 和坐标。
   - 重点确认 1 基转 0 基、边方向、edge 顺序是否和 A 的未知量顺序一致。

4. 在确认 G/coords 和 A 顺序一致后，再看求解参数。
   - 当前外层 `rtol=1e-9`、`maxiter=100`。
   - 内层 B 求解频繁达到 100 次。
   - 需要对比对方原 AMS 程序参数，而不是盲目加大迭代数。

5. AMS 路径能复现后，再继续做接口迁移。
   - C/C++ 数组接口。
   - Fortran `iso_c_binding`。
   - 进一步决定是否保留 C++ 版本或重写 Fortran 版本。

## 给接手者的一句话摘要

新 full CSR 数据和上三角数据已经确认一致，`slove.txt` 也确实是这组 A/b 的直接解；但当前 cpp AMS 求解路径用原物理符号会输出全 0，用内部符号能迭代但停在 `1e-6` 量级，不能复现直接解。下一步应集中排查 AMS 的 G/coords、edge 顺序、边界条件/重排和 PETSc/Hypre 参数是否和对方原运行一致，而不是先做 Fortran 绑定。
