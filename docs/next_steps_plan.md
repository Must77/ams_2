# 后续工作路线图：内存瘦身 → Fortran 混编 → 纯 Fortran 重写

> 面向压缩后的会话或零上下文接手者。
> 先读 `project_handoff.md`（项目背景 + 2026-06-13 节的根因链与验收结论），再读本文。

## 状态跟踪（活节，随工作推进更新）

**维护规则：完成或推翻下列任何一项的提交，必须在同一提交中更新本节。**
本节为唯一的进度真相源；正文各节描述"怎么做"，不描述"做没做"。

- [ ] 1a 消除双重持有（ctx 数组组装后释放）
- [ ] 1b 外层 restart 降档实验与落地
- [ ] 1c C/M 改 SBAIJ（可选深化）
- [ ] 2 dual-mode 验证关闭（含改打印文案）
- [ ] 3a C 接口层 ams_c_api
- [ ] 3b Fortran 绑定模块
- [ ] 3c Fortran 示例与数值对照
- [ ] 3d CMake 静态库重构
- [ ] 4-探针 PETSc Fortran stub 覆盖度检查（决定纯 Fortran 路线）
- [ ] 4 纯 Fortran 移植（探针通过后细化）
- [ ] 运维：push 到 origin/dev（等凭据）
- [ ] 运维：沟通包发给对方（5 点见 §0）

最后更新：2026-06-13，求解器打通后初始状态。

## 0. 现状快照

- **求解器已打通**：AMS 解代回原方程真实相对残差 7.1e-8（优于直接解 slove.txt 的 6.8e-7）。
  与 slove.txt 全局 L2 差 0.69 是空气层规范自由度假象（差异 100% 在低 σ 边、96% 为梯度模态），
  导电区吻合 3.2e-5。**勿把 0.69 当 bug 重新排查。**
- **标准运行命令**：
  ```bash
  cd ~/ams_2 && ./build_ubuntu/UpperExampleCheck ~/ams_example \
    -fortran_upper_unscale_edge_len -ams_beta_mass_poisson \
    -ams_beta_mass_shift 1e-6 -B_ksp_type gmres -B_ksp_max_it 15 \
    -em_outer_max_it 60 -em_outer_rtol 1e-7
  ```
  单核约 35 分钟（含 ~5 分钟文本解析），峰值 RSS ~3.2GB。
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

### 1a 消除双重持有（首做，收益 ~460MB 峰值）
- 事实：`assemble_matrix()` 用 `MatCreateMPIAIJWithArrays`（**拷贝语义**，PETSc 文档确认 MPIAIJ
  变体复制数据）建 C/M（PetscAMGInterface.cc ~387/390 行，骨架在 ~361）。
- 改法：C/M 组装完成（MatAssemblyEnd 后）即释放 `ctx->data_real/data_imag/cidx/rptr`
  （`std::vector<T>().swap(v)`）；`ctx->b_real/b_imag` 在 `assemble_rhs_csem` 后同理；
  adapter 层（`solve_eg1_fortran_upper_1based` 中的局部 vector）在 `solve_eg1` 返回后自动析构，
  但 `prepare_fortran_upper_1based_inputs` 输出与 `solve_eg1` 入参之间的生命周期重叠无法避免——
  可把 prepare 的输出 move 进 solve（需改签名为 vector&&或在 prepare 内部就地缩减 upper 数组）。
  第一刀先做 ctx 释放（不改签名、零风险）。
- 验证点：`matshell_mult_a` 只用 Mat C/M（已确认）；检查 `copy_result_to_arrays`、
  `destroy_*` 不再触碰已释放向量；`setup_ams` 用 `ctx->edgesN/nodes`（G 与坐标），
  这两个数组在 create_pc 之后才能释放。
- 计量：`/usr/bin/time -v ./build_ubuntu/UpperExampleCheck ...` 记录 Maximum resident set size，
  与基线 3.2GB 对比。

### 1b 外层 restart 降档（一行/纯选项，收益 ~325MB）
- 先用选项实验：标准命令加 `-A_ksp_gmres_restart 15`。注意外层在 run9 配置下 ~37 步收敛，
  restart 15 会引入一次重启，验证迭代数不显著增加（≤1.5×）后把默认写进 `create_pc()`
  （`KSPGMRESSetRestart(ctx->A_ksp, 15)`，放 `KSPSetFromOptions` 之前以便选项可覆盖）。
- 若迭代数恶化明显，试 20；权衡点写进提交信息。

### 1c C/M 改 SBAIJ 对称存储（最后做，收益 ~135MB 常驻，工作量中）
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
