# 04 求解器架构：代码地图

> 前置：03。目标：拿着这篇能直接读 `PetscAMGInterface.cc`，知道每个函数在干什么、
> 两层迭代怎么嵌套、数据从对方文件到解向量怎么流动。

## 1. PETSc 速成（只讲本项目用到的）

PETSc 的对象模型，用你熟悉的概念类比：

| PETSc 对象 | 是什么 | 类比 |
| --- | --- | --- |
| `Mat` | 矩阵（多种存储格式） | 你的 CSR 数组三件套，带方法 |
| `MATAIJ` | 标准 CSR 稀疏矩阵 | 同上 |
| `MatShell` | **无矩阵**算子：只注册一个"乘向量"回调函数 | 你写的 `matvec` 子程序 |
| `Vec` | 向量 | 数组 |
| `KSP` | Krylov 迭代求解器（CG/GMRES/...统称） | 你的迭代法主循环 |
| `PC` | 预条件器，挂在 KSP 上 | SOR/ILU/多重网格那一层 |
| `PCSHELL` | 自定义预条件器：注册一个"应用"回调 | 你自己写的预条件子程序 |
| `-xxx_yyy` 命令行选项 | 运行时改参数不用重编译 | namelist |

选项有**前缀**机制：本代码外层求解器前缀 `A_`、内层前缀 `B_`，所以命令行
`-B_ksp_type gmres` 只改内层。`KSPSetFromOptions()` 调用时机决定命令行能否覆盖
代码里的硬编码（在 SetFromOptions **之后**再硬编码的参数覆盖不了——外层 maxiter
就栽在这，已通过 `-em_outer_max_it` 选项修复）。

## 2. 两层嵌套结构（核心）

```text
外层: FGMRES, 作用在 2N×2N 实数化系统上 (MatShell, 前缀 A_)
  │     算子 = matshell_mult_a():  yr = C·xr - M·xi;  yi = -M·xr - C·xi
  │     预条件器 = PCSHELL → pc_apply_b()
  │
  └─ pc_apply_b(): 把 2N 向量拆成实半/虚半，各做一次内层求解
        │
        内层: GMRES(原 FCG), 作用在 N×N 实矩阵 B = C + M 上 (前缀 B_)
              预条件器 = Hypre AMS (05 篇)
              容差 rtol=1e-2, 上限可调 (-B_ksp_max_it)
```

**为什么这样设计**：理想预条件器是 A⁻¹ 的近似。对这个 2×2 块系统，块对角近似
`diag(B⁻¹, B⁻¹)`（B = C+M）是文献标准做法——B 是对称正定的"实化代理"，可以用
AMS 高效近似求逆。外层用 **F**GMRES（Flexible GMRES）而不是普通 GMRES，因为内层
是"迭代解到 1e-2 就停"的**不精确/每次都略有不同**的预条件器，普通 GMRES 理论上
不允许预条件器变动，FGMRES 专门为此设计。

每步外层迭代成本 = 1 次 MatShell 乘（4 次稀疏矩阵乘向量）+ 2 次内层求解（各若干步
GMRES，每步 1 次 B 乘 + 1 次 AMS 应用）。**内层是绝对的成本大头**。

## 3. 函数地图（按调用顺序）

```text
UpperExampleCheck main()                    [upper_example_check.cpp]
 ├─ 读对方文件 (jcol/irw/ssor_k1/k2/edges/node/right-hand/slove)
 ├─ 用上三角数据+slove.txt 算"参考残差" (验证数据自洽，应得 ~6.8e-7)
 ├─ solve_eg1_fortran_upper_1based()        [PetscAMGInterface.cc]
 │   ├─ prepare_fortran_upper_1based_inputs()
 │   │    ├─ 边节点号 1基→0基
 │   │    ├─ expand_upper_triangle_1based_to_full_csr()  上三角→全量, 1基→0基
 │   │    ├─ data_imag 取负 (M = -Ki, 见03)
 │   │    └─ [修复1] -fortran_upper_unscale_edge_len: 值÷(len_r·len_c), rhs÷len_r
 │   ├─ solve_eg1()
 │   │    ├─ process_options()       读选项 (em_ctx.cc, 含 -em_outer_max_it/-em_outer_rtol)
 │   │    ├─ create_context()        分配 EMContext (全局状态包)
 │   │    ├─ load_context_from_arrays()  调用者数组 → ctx 内部 vector
 │   │    ├─ create_linear_system()  建 PETSc Vec/Mat 骨架
 │   │    ├─ assemble_matrix()       CSR 数据 → ctx->C (实部), ctx->M (虚部已取负)
 │   │    ├─ create_pc()             ★ 全文最重要的函数
 │   │    │    ├─ B = C + M (MatAXPY), 标记 SPD
 │   │    │    ├─ setup_ams():  从 edgesN/nodes 构造 ±1 的 G 和坐标数组
 │   │    │    ├─ B_ksp: FCG(可被 -B_ksp_type 覆盖) + PCHYPRE "ams"
 │   │    │    │    ├─ PCHYPRESetDiscreteGradient(G), PCSetCoordinates(坐标)
 │   │    │    │    └─ [修复2+3] -ams_beta_mass_poisson / -ams_beta_mass_shift:
 │   │    │    │         Aβ = Gᵀ·(M+shift·I)·G  → PCHYPRESetBetaPoissonMatrix
 │   │    │    ├─ A = MatShell(matshell_mult_a), 2N×2N
 │   │    │    └─ A_ksp: FGMRES + PCSHELL(pc_apply_b)
 │   │    ├─ assemble_rhs_csem()     rhs → PETSc Vec
 │   │    ├─ solve_linear_system()   KSPSetTolerances(外层) + KSPSolve  ← 真正求解
 │   │    └─ copy_result_to_arrays() 解 → 调用者数组
 │   └─ 后处理: 虚部取负还原; [修复1配对] 解÷len 还原对方约定
 ├─ 与 slove.txt 比对, 打印 solve_relative_l2 等
 └─ 落盘 ams_solution.txt (本轮新增, 供离线法证)
```

## 4. EMContext：全局状态包

`em_ctx.h` 定义的结构体，所有矩阵/向量/参数都挂在上面。和本次工作相关的成员：

| 成员 | 内容 |
| --- | --- |
| `rptr, cidx, data_real, data_imag` | 全量 0 基 CSR（虚部已按 M=-Ki 取负） |
| `b_real, b_imag` | 右端项 |
| `edgesN, nodes` | 边-节点表（0 基）与节点坐标，AMS 的原料 |
| `C, M, B, A, G` | PETSc 矩阵：实部、负虚部、内层 B=C+M、外层 MatShell、离散梯度 |
| `A_ksp, B_ksp` | 外层/内层求解器 |
| `K_max_it, dual_rtol` | 外层上限/容差（原硬编码 100/1e-9，现可选项覆盖） |

## 5. 日志去哪了

两个监视器（外层 A_、内层 B_）都把逐步残差写进一个 `*-group-000.log` 文件
（文件名前缀来自未初始化内存，是乱码——历史遗留小 bug，无害）。注意日志里
内外层行混杂、缩进层级因 PushTab/PopTab 动态变化，**不要用脚本盲目解析**，
判断收敛以程序最终打印的指标和"解代回残差"为准（10 篇有正确的解读方法）。
