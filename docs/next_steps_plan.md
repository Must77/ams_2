# 后续工作路线图：内存瘦身 → Fortran 混编 → 纯 Fortran 重写

> 面向压缩后的会话或零上下文接手者。
> 先读 `project_handoff.md`（项目背景 + 2026-06-13 节的根因链与验收结论），再读本文。

## 状态跟踪（活节，随工作推进更新）

**维护规则：完成或推翻下列任何一项的提交，必须在同一提交中更新本节。**
本节为唯一的进度真相源；正文各节描述"怎么做"，不描述"做没做"。

- [x] 1a 消除双重持有（ctx 数组组装后释放；另修复 C/M/s 重复创建泄漏，实测见 §1a）
- [x] 1b 外层 restart 降档实验与落地（restart 15 进标准命令而非代码默认，理由见 §1b）
- [x] ~~1c C/M 改 SBAIJ（可选深化）~~ 划掉：阶段一目标 2GB 已达成（实测 2.08GB），
      135MB 收益不抵 SBAIJ 改动面与 PtAP/hypre 兼容风险；若未来内存再紧再启
- [x] 2 dual-mode 验证关闭（2026-06-13：三个大算例日志各含 1 次该打印，重复求解确证
      不存在；误导文案已改为 "Solving linear system:"）
- [x] 3a C 接口层 ams_c_api（2026-06-13）
- [x] 3b Fortran 绑定模块 ams_solver_mod.f90（2026-06-13）
- [x] 3c Fortran 示例与数值对照（2026-06-13；**偏离原计划**：用内联 4 边小算例对照
      numpy 精确解（相对差 1.79e-9，秒级可常驻回归），而非读大算例文本——绑定层只是
      参数转发，大算例数值已由 UpperExampleCheck 验证；与对方真实 main.f90 的大算例
      集成测试留到对接时用对方自己的读文件代码做）
- [x] 3d CMake 静态库重构 AmsCore（2026-06-13）
- [ ] 4-探针 PETSc Fortran stub 覆盖度检查（决定纯 Fortran 路线）
- [ ] 4 纯 Fortran 移植（探针通过后细化）
- [x] 运维：push 到 origin/dev（2026-06-13 完成；凭据已存 credential.helper store，后续可直接 push）
- [ ] 运维：沟通包发给对方（5 点见 §0）

最后更新：2026-06-13，阶段一（3.2GB → 2.08GB）、阶段二、阶段三全部关闭。
下一步：阶段四探针（PETSc Fortran stub 覆盖度），以及等用户跟对方对接沟通包与混编集成。

阶段三使用说明（对方接入）：链接 `build_ubuntu/libAmsCore.a` + `ams_solver_mod.f90`，
Fortran 侧 `use ams_solver`，调用顺序 init→solve→finalize（见 fortran_example_check.f90 范例；
init 的推荐 options 字符串见 ams_c_api.h 注释，含 -A_ksp_gmres_restart 15）。
回归命令：`./build_ubuntu/FortranExampleCheck`（在 ~/ams_2 下运行，预期输出
fortran_check_relative_l2 ~1.8e-9 + PASSED）。

## 0. 现状快照

- **求解器已打通**：AMS 解代回原方程真实相对残差 7.1e-8（优于直接解 slove.txt 的 6.8e-7）。
  与 slove.txt 全局 L2 差 0.69 是空气层规范自由度假象（差异 100% 在低 σ 边、96% 为梯度模态），
  导电区吻合 3.2e-5。**勿把 0.69 当 bug 重新排查。**
- **标准运行命令**：
  ```bash
  cd ~/ams_2 && ./build_ubuntu/UpperExampleCheck ~/ams_example \
    -fortran_upper_unscale_edge_len -ams_beta_mass_poisson \
    -ams_beta_mass_shift 1e-6 -B_ksp_type gmres -B_ksp_max_it 15 \
    -em_outer_max_it 60 -em_outer_rtol 1e-7 -A_ksp_gmres_restart 15
  ```
  单核实测 ~8-9.5 分钟（2026-06-13；更早记录为 35 分钟），峰值 RSS 2.08GB
  （1a+1b 后实测；restart 15 比默认 30 多 ~19% 迭代/墙钟，换 312MB 内存）。
- **回归基线**（任何改动后必须全过）：
  1. `cd exmaple/eg_1 && ../../build_ubuntu/FemAms`，`result.txt` 的 sha256 必须是
     `342042dbcff12d00a13339367327a851158f0dc137f6af5caa1ee6779e9fd266`
     （跑完后 `git checkout -- "exmaple/eg_1/-group-000.log"` 还原日志噪声）；
  2. `./build_ubuntu/AdapterSelfCheck` 静默退出 0；
  3. 大算例标准命令：`reference_relative_residual` 必须仍为 6.8111458886769823e-07，
     `solve_relative_l2` ≈ 0.69（见上，预期值），解代回残差用
     `/tmp/solution_forensics.py`（远程机器）核 ~7e-8。
- **git**：分支 dev，user `magina <15947113700@163.com>`。**push 待办**——远端
  https://github.com/Must77/ams_2 是 HTTPS 且机器无凭据，等用户提供 PAT 或自行 gh auth login。
- **诊断工具**（已提交）：`UpperDumpSystem <data> <out> [-fortran_upper_unscale_edge_len]` 把
  B/M/G/coords/rhs 存 PETSc 二进制（已有 ~/ams_dump_whitney 与 ~/ams_dump_scaled）；
  `AmsLab <dump_dir> -lab_beta_mass [-lab_mass_shift x] [-lab_pc ams|jacobi|none]
  [-lab_random_rhs] [-lab_ksp_*]` 秒级试参。改 PC 相关代码前先在 AmsLab 验证。
- **待发沟通包**（等用户发给对方）：① 确认单位切向归一基函数；② AMS 残差 7.1e-8 已优于直接解；
  ③ 建议验收改为代回残差/导电区/接收点比对，并询问其原 AMS 程序的验收方式；
  ④ 若需空气区一致请告知规范约定（可加解后投影）；⑤ 请求 16 位精度重导数据
  （8 位精度是当前一切 β 矩阵 workaround 的根源，也是精度地板）。

## 1. 阶段一：内存瘦身（目标 3.2GB → ~2GB，约 4.9KB/边 → 3KB/边）

内存账（682344 边、上三角 nnz 5.93M、全量 11.18M、内部实自由度 1.36M）：

| 项 | 估算 | 处置 |
| --- | --- | --- |
| adapter 展开的全量 CSR（row_ptr/col_idx/data_real/data_imag）| ~230MB | 1a 释放 |
| `load_context_from_arrays` 拷进 ctx 的同样一份 | ~230MB | 1a 释放 |
| PETSc C、M（全量 AIJ）| ~270MB | 1c 可改 SBAIJ |
| B=C+M 的 PETSc 份 + hypre ParCSR 份 | ~270MB | 必须全量，不动 |
| AMS 层级（Π、粗网格）| ~300-400MB | 必要开销 |
| 外层 FGMRES V+Z 双基（restart 30 × 1.36M × 8B × 2）| ~650MB | 1b 降 restart |
| 内层 GMRES 基 + 杂项向量 | ~300MB | 顺带观察 |

### 1a 消除双重持有（2026-06-13 完成）
- 事实：`assemble_matrix()` 用 `MatCreateMPIAIJWithArrays`（**拷贝语义**，PETSc 文档确认 MPIAIJ
  变体复制数据）建 C/M。
- **实施时的两个新发现（比原计划多修的泄漏）**：
  1. `create_linear_system()` 原先用 NULL 值预建骨架 C + MatDuplicate 出 M，随后
     `assemble_matrix()` 重新 MatCreate 覆盖句柄而未销毁 → 骨架 C/M（约 230MB）纯泄漏；
     `s.re/s.im` 同样被 `assemble_rhs_csem` 的 VecCreateMPIWithArray 二次创建覆盖（约 22MB）。
     修复：create_linear_system 不再预建 C/M/s，dual_e/w 改由 VecCreateMPI 直接创建。
  2. **原计划"`ctx->b_real/b_imag` 在 assemble_rhs_csem 后同理释放"是错的**：
     `assemble_rhs_csem` 用 `VecCreateMPIWithArray`（**引用语义**，不拷贝），s.re/s.im 全程
     直接引用 ctx->b_real/b_imag 的内存，释放即悬垂指针。两数组必须保留到 destroy。
- 实际改法（全在 PetscAMGInterface.cc）：`assemble_matrix` 末尾释放 cidx/data_real/data_imag
  （约 224MB，`std::vector<T>().swap(v)`）；rptr 因 create_pc/assemble_rhs_csem 还要读尺寸，
  延至 `assemble_rhs_csem` 末尾释放；`setup_ams` 末尾释放 edgesN/nodes（G 已组装、坐标已转存
  v_coords，约 27MB）。v_coords 因 PCSetCoordinates 拷贝语义未查证，保守未释放（16MB）。
- 仍未做（可选深化）：adapter 层局部 vector 与 ctx 拷贝的生命周期重叠（~230MB），需把
  prepare 的输出 move 进 solve（改签名）。若 1b 后仍未达 2GB 目标再做。
- 实测（2026-06-13）：峰值 RSS 基线 ~3.2GB（run9，未留 time -v 精确值）→ **2.38GB**
  （2497024KB），墙钟 7m54s（此前文档估 35 分钟，差异可能来自内存压力降低或原估值偏保守）。
  四项回归全过：eg_1 黄金哈希、AdapterSelfCheck、大算例 reference_relative_residual
  精确不变 + solve_relative_l2=0.690、代回残差法证 7.140e-8（高 σ 区 3.1678e-05，与修改前一致）。

### 1b 外层 restart 降档（2026-06-13 完成，选项落地）
- 实验结果（大算例，1a 基础上加 `-A_ksp_gmres_restart 15`）：外层 44 步收敛到 9.94e-8
  （restart 30 时 ~37 步，+19%），峰值 RSS 2.38GB → **2.08GB**（省 312MB ≈ 预估），
  墙钟 7:54 → 9:24，reference_relative_residual 精确不变，代回残差 6.814e-8、
  高 σ 区 3.1601e-05 全部正常。
- **决策：不改代码默认值，restart 15 以选项形式进标准命令与阶段三推荐 options 字符串**。
  原因：eg_1 的外层要 25 步（>15），把 15 写成编译默认会改变 eg_1 的 FGMRES 收敛轨迹、
  破坏黄金哈希基线；而大算例生产配置本来就走选项字符串，收益一分不少。

### 1c C/M 改 SBAIJ 对称存储（2026-06-13 决定不做：阶段一目标已达成，收益不抵风险；保留原方案备查）
- `MatCreateSBAIJ`（bs=1）只存上三角；`MatMult` 兼容（matshell 不用改）。
- 两处下游需要全量：B=C+M 给 hypre（`MatConvert` SBAIJ→AIJ 后 AXPY，或直接从上三角数组
  另建 AIJ 的 B）；`MatPtAP(M,G)` 构造 Aβ（SBAIJ 不支持 PtAP，先 MatConvert 临时 AIJ，
  用完即毁——瞬态全量不可避免但常驻减半）。
- 该项做完后 expand 阶段也可改为只产上三角数组（省 1a 中展开数组的一半），属可选深化。

## 2. 阶段二：dual-mode 验证（大概率直接划掉）

- **2026-06-13 更正**：此前"同一系统解两遍"为误判，来源是对 636 字节日志 head/tail 重叠输出的
  误读。事实：`solve_eg1()` 只调用一次 `solve_linear_system()`（PetscAMGInterface.cc ~865），
  run9 完整迭代日志中 "Solving for dual mode:" 仅出现 1 次。
- 留 15 分钟动作：跑标准命令时 `grep -c "Solving for dual mode"` 完整日志确认一次性；顺手把
  这条误导性的打印文案改掉（如 "Solving linear system:"），然后关闭此阶段。

## 3. 阶段三：C++/Fortran 混合编程（iso_c_binding）

目标：对方 Fortran 程序像调普通子程序一样调用求解器，数组内存直传（项目原始目标）。

- **3a C 接口层**（新文件 `src_ams/ams_c_api.h/.cc`）：
  - `extern "C"` 包装三个函数：
    `ams_petsc_init(const char* options_string)`（内部 PetscInitialize + PetscOptionsInsertString，
    把 `-fortran_upper_unscale_edge_len -ams_beta_mass_poisson ...` 等作为字符串传入，避免
    Fortran 侧操心命令行）；`ams_petsc_finalize()`；
    `ams_solve_upper(...)` 直接转发 `solve_eg1_fortran_upper_1based`（参数已全是
    `PetscInt*/PetscReal*/double*` POD 指针，无需改造）。
  - 返回错误码 int（PetscErrorCode 透传），Fortran 侧检查非零。
- **3b Fortran 绑定模块**（`src_ams/ams_solver_mod.f90`）：
  - `bind(C)` interface；本机 PETSc 为 **32 位 PetscInt**（apt 默认）→ `integer(c_int)`；
    `PetscReal=double` → `real(c_double)`。若未来换 64 位 PETSc 需同步改 c_long——在模块注释里写明。
- **3c Fortran 示例**（`src_ams/fortran_example_check.f90`）：
  - 以对方 `~/ams_example/main.f90` 为骨架（其读文件代码可直接抄，注意它读的
    `right hand.txt` 文件名有空格、实际文件是 `right-hand.txt`，以实际为准），读上三角 →
    调 `ams_solve_upper` → 写解 → 与 C++ 版 `ams_solution.txt` 比对（相对差 < 1e-12）。
- **3d CMake**：项目已 `PROJECT(... Fortran)`；新增
  `ADD_EXECUTABLE(FortranExampleCheck fortran_example_check.f90 ams_solver_mod.f90 ...)`，
  `TARGET_LINK_LIBRARIES(... FemAmsCore? ${PETSC_LIBRARY} stdc++ m)`——当前源码是每个可执行
  重复编译 .cc，顺手抽一个 STATIC library target 避免三重编译（也加速 CI/构建）。
  链接器用 Fortran 驱动（CMake 自动），需显式补 `stdc++`。
- 验收：FortranExampleCheck 跑 ams_example 输出与 C++ 版数值一致；eg_1 回归不受影响。

## 4. 阶段四：纯 Fortran 重写

- **第一步是可行性探针（半天，先做再承诺）**：PETSc 原生 Fortran API 覆盖度检查——
  `PCHYPRESetDiscreteGradient` / `PCHYPRESetBetaPoissonMatrix` / `PCSetCoordinates` /
  `MatShellSetOperation` / `PCShellSetApply` / `MatCreateMPIAIJWithArrays` / `MatPtAP` /
  `VecNest` 系列在 `/usr/lib/petsc` 的 Fortran stub（ftn-auto/ftn-custom）中是否存在。
  历史上 PCHYPRE 一族 Fortran 绑定不全。**若缺**：保留 ~50 行 C shim（仅这几个调用），
  与用户确认"彻底 Fortran"是否接受这个例外；或自写 stub。
- 移植顺序（每步带数值对照测试）：
  1. 工具函数 `expand_upper_triangle_1based_to_full_csr` + `prepare_fortran_upper_1based_inputs`
     → 纯 Fortran，单测对照 AdapterSelfCheck 的用例数据；
  2. 矩阵/向量组装 + setup_ams（G、坐标）；
  3. MatShell/PCShell 回调（Fortran 回调注意 PETSc 的 `use petscksp` 模块接口与回调签名）；
  4. 求解驱动 + 选项处理；
  5. 双算例（eg_1 + ams_example）与 C++ 版数值一致后，C++ 版降级为参考实现保留在仓库。
- 风险登记：VecNest 的 Fortran 支持、回调里的 PushTab 日志（可砍）、PetscViewer 文件名
  乱码 bug（`create_context` 的 LS_log 文件名本来就有乱码——重写时顺手修）。

## 5. 环境与操作备忘（压缩保险）

- 远程：`sshpass -p 233719 ssh -p 9522 magina@10.170.108.78`（MaginaLaptop，Ubuntu 26.04，
  7.7GB RAM / 8 核）。仓库 `~/ams_2`（dev 分支），数据 `~/ams_example`（1.7G 文本）。
- 构建：`cmake -S src_ams -B build_ubuntu -DEIGEN_DIR=/usr -DMETIS_DIR=/usr
  -DNANOFLANN_DIR=/usr -DPETSC_DIR=/usr/lib/petsc -DPETSC_ARCH=` + `cmake --build build_ubuntu -j8`。
  PETSc 3.24.4（apt，带 hypre，32 位 PetscInt）。
- 已知陷阱：
  1. ssh 单引号外壳里写含 `'` 的 here-doc 会被截断——长代码一律 base64 过线；
  2. `pgrep/pkill -f` 会自匹配命令行——模式写成 `UpperExampleChec[k]` 且 pkill 仍可能杀到
     自身 shell（匹配 nohup 串），优先 `kill $(pgrep -f ...)`；
  3. 远程 python 输出有缓冲——`python3 -u`；
  4. eg_1 回归会弄脏被跟踪的 `-group-000.log`——提交前 `git checkout --` 它；
  5. **hypre AMS 的 G 必须是 ±1 关联矩阵**（缩放 G 的方案已实验证伪）；
  6. **内层禁用 FCG**（近奇异模态误报 INDEFINITE 第 2 步退出），用 GMRES；
  7. `-ams_beta_mass_shift` 只 shift `ctx->M` 的副本，**绝不能动 ctx->M 本体**（它在算子里）；
  8. 解的全局 L2 vs slove.txt ≈ 0.69 是预期值不是回归失败。
- 法证脚本（远程 /tmp，机器重启会丢，必要时从本文逻辑重写）：`solution_forensics.py`
  （代回残差 + 梯度投影 + σ 分区）、`g_orientation_probe.py`、`b_spd_probe.py` 等。
