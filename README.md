# rsksw2 — KSW2 的 Rust 移植

[![CI](https://github.com/keithyin/rsksw2/actions/workflows/ci.yml/badge.svg)](https://github.com/keithyin/rsksw2/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/keithyin/rsksw2)](https://github.com/keithyin/rsksw2/releases)
[![crates.io](https://img.shields.io/crates/v/rsksw2)](https://crates.io/crates/rsksw2)
[![docs.rs](https://docs.rs/rsksw2/badge.svg)](https://docs.rs/rsksw2)
[![License](https://img.shields.io/badge/license-MIT-blue)](#许可证)
[![Rust](https://img.shields.io/badge/rust-2024-dea584?logo=rust)](https://www.rust-lang.org/)
[![Platforms](https://img.shields.io/badge/platforms-x86__64%20linux%20%7C%20macos%20%7C%20windows-blue)](#构建与运行)

本仓库是 [KSW2](https://github.com/lh3/ksw2)（基于动态规划的双序列比对库）的 Rust 移植。
原始的 C 代码完整保留在 [`ksw2/`](ksw2) 目录中，既作为参考实现，也作为一致性测试的
oracle（`scripts/gen_expected.sh` 用 C 版本生成期望输出）。

移植目标：

1. **结果一致性** —— 对同一输入、同一参数，输出（score、max/max_t/max_q、mqe/mte、
   是否 zdropped、CIGAR）与 C 版本逐字节一致；
2. **性能一致性** —— 关键实现保持与 C 版本相同的算法结构与向量化方式，
   运行时间在同一量级（见[性能](#性能)）。

## 目录结构

| C 源文件 | Rust 模块 | 公开函数 | 说明 |
|:--|:--|:--|:--|
| `ksw2_gg.c` | `src/gg.rs` | `ksw_gg` | 全局比对，Green 经典形式 |
| `ksw2_gg2.c` | `src/gg2.rs` | `ksw_gg2` | 全局比对，Suzuki 对角形式 |
| `ksw2_gg2_sse.c` | `src/sse/gg2_sse.rs` | `ksw_gg2_sse` | 全局比对 + SSE4.1 |
| `ksw2_extz.c` | `src/extz.rs` | `ksw_extz` | 扩展比对，单段仿射罚分 |
| `ksw2_extz2_sse.c` | `src/sse/extz2_sse.rs` | `ksw_extz2_sse` | 扩展比对 + SSE4.1 |
| `ksw2_extd.c` | `src/extd.rs` | `ksw_extd` | 扩展比对，双段仿射罚分 |
| `ksw2_extd2_sse.c` | `src/sse/extd2_sse.rs` | `ksw_extd2_sse` | 双段仿射 + SSE4.1 |
| `ksw2_extf2_sse.c` | `src/sse/extf2_sse.rs` | `ksw_extf2_sse` | 线性罚分 + X-drop，仅打分 |
| `ksw2_exts2_sse.c` | `src/sse/exts2_sse.rs` | `ksw_exts2_sse` | 剪接感知扩展 + SSE4.1 |
| `cli.c` | `src/bin/ksw2-test.rs` | 可执行文件 `ksw2-test` | 与 C 版命令行完全兼容 |
| `ksw2.h`（结构体/工具函数） | `src/lib.rs`, `src/common.rs` | `KswExtz`, `ksw_*` | 标志位、CIGAR、回溯等 |

`ksw_*_sse` 系列是 C 源文件中 `__SSE4_1__` 分支的逐行移植（即 `-march=native`
时编译出来的那条路径），内部使用 `core::arch::x86_64` 的 SSE4.1 intrinsics，
因此与 C 版本在 score、坐标和 CIGAR 上完全一致。它们与 C 一样只在 x86_64 上编译，
并在运行时检查 `sse4.1`（2008 年之后的 x86 CPU 均支持）。

## 用法

```rust
use rsksw2::*;

// 5x5 打分矩阵：match +2、mismatch -4、wildcard 行/列为 0
let score_matrix: [i8; 25] = [ 2, -4, -4, -4, 0,
                             -4,  2, -4, -4, 0,
                             -4, -4,  2, -4, 0,
                             -4, -4, -4,  2, 0,
                              0,  0,  0,  0, 0];
let query:  &[u8] = &[0, 1, 2, 3];   // A C G T
let target: &[u8] = &[0, 1, 2, 3];

let mut result = KswExtz::new();                 // 可跨调用复用，内部 CIGAR 缓冲会被重复利用
ksw_extz2_sse(4, query, 4, target, 5, &score_matrix,
              4, 2,        // gap_open, gap_extension
              -1, -1,      // band_width = inf, z_drop = disabled
              0, 0,        // end_bonus, flags
              &mut result);
println!("score = {}, cigar = {}", result.score, result.cigar_string());
```

### 函数签名对照

`ksw_extz_t` → `KswExtz`：

| C | Rust |
|:--|:--|
| `uint32_t max:31` + `int zdropped:1` | `max: u32`（31 位，保留与 C 相同的截断语义）+ `zdropped: bool` |
| `max_q`, `max_t`, `mqe`, `mqe_t`, `mte`, `mte_q`, `score`, `reach_end` | 同名字段，类型 `i32` |
| `m_cigar` / `n_cigar` / `cigar` | `cigar: Vec<u32>`（`n_cigar()` 返回长度，容量自动复用） |
| `ksw_reset_extz()` | `KswExtz::reset()` / `KswExtz::new()` |

`ksw_gg` / `ksw_gg2` 在 C 中用三个指针表示「是否输出 CIGAR」，Rust 中改为
`Option<&mut Vec<u32>>`：`None` 等价于 `-s`（只算分数）；`ksw_gg2_sse` 在 C 中
总是回溯，Rust 版同样要求传入 `&mut Vec<u32>`。

## 一致性验证

三个层次，全部可复现：

1. **与 C 实现的逐字节对比**（`scripts/compare_with_c.sh`）
   对 10 个算法 × 14 组参数（`-r`、`-s`、`-g`、`-z`、`-w`、`-O/-E`、`-A/-B`、`-a` …）
   × 3 组数据（`t1/q1`、`t2/q2` 压缩输入、`MT-human/MT-orang` ≈16.5 kb）执行
   C 版和 Rust 版二进制，比较 stdout+stderr：

   ```
   $ ./scripts/compare_with_c.sh large
   compared 420 cases, 0 mismatches
   ```

2. **golden 回归测试**（`tests/consistency.rs` + `tests/expected/t1_q1.txt`）
   `scripts/gen_expected.sh` 用 **C 版** 生成期望输出并入库；
   `cargo test` 用 Rust 版重放同样的 140 个用例并要求完全一致。

3. **对 C 库的随机模糊对照**（`tests/fuzz_against_c.rs` + `tests/c/fuzz_ref.c`）
   `scripts/gen_fuzz_expected.sh` 用 C 库跑 2000 个随机用例（随机序列长度 1..60、
   随机打分/罚分、随机带宽、随机 z-drop、全部标志组合——包括 `cli.c` 无法设置的
   `KSW_EZ_EXTZ_ONLY`、`KSW_EZ_GENERIC_SC`、`KSW_EZ_EQX`、剪接 `junc` 数组等），
   覆盖全部 9 个函数共 2 万次调用；`cargo test` 用 Rust 重放同样用例并要求一致。

4. **实现间交叉验证**（`tests/impl_agreement.rs`）
   对随机序列比较 `ksw_gg` / `ksw_gg2` / `ksw_gg2_sse`，以及
   `ksw_extz` vs `ksw_extz2_sse`、`ksw_extd` vs `ksw_extd2_sse` 的
   score / CIGAR / max 坐标是否一致。

```
$ cargo test --release
test matches_c_reference_output ... ok            # 140 个 CLI 用例 vs C 二进制
test matches_c_library_on_randomised_cases ... ok # 2000 个随机用例 vs C 库
test affine_extension_agrees ... ok
test dual_affine_extension_agrees ... ok
test global_alignment_agrees ... ok
test scalar_ports_match_reference_example ... ok
```

> 注：`cli.c` 用 `"MID"[op]` 打印 CIGAR 操作符，操作符 3（`N`）会打印出 `'\0'`；
> Rust 版刻意保留这一行为以便逐字节对比（`=`/`X` 则打印为 `=`/`X`）。

## 性能

`scripts/bench.sh` 用 `ksw2/test/MT-human.fa` × `MT-orang.fa`（16.5 kb，KSW2
README 中的基准数据）对 C 与 Rust 两个二进制取 3 次运行的最小值：

| 算法 | 模式 | C (s) | Rust (s) | Rust/C |
|:--|:--|--:|--:|--:|
| `extz2_sse` | `-s` | 0.139 | 0.135 | 0.97 |
| `extz2_sse` | CIGAR | 0.430 | 0.433 | 1.01 |
| `extd2_sse` | `-s` | 0.171 | 0.178 | 1.04 |
| `extd2_sse` | CIGAR | 0.502 | 0.497 | 0.99 |
| `gg2_sse` | `-s` | 0.558 | 0.531 | 0.95 |
| `gg2_sse` | CIGAR | 0.558 | 0.521 | 0.93 |
| `extz` | `-s` | 1.130 | 0.961 | 0.85 |
| `extz` | CIGAR | 1.669 | 1.624 | 0.97 |
| `extd` | `-s` | 1.484 | 2.052 | 1.38 |
| `extd` | CIGAR | 2.700 | 2.712 | 1.00 |
| `gg` | `-s` | 0.948 | 0.684 | 0.72 |
| `gg` | CIGAR | 1.466 | 1.574 | 1.07 |
| `gg2` | `-s` | 1.054 | 1.429 | 1.35 |
| `gg2` | CIGAR | 2.104 | 1.985 | 0.94 |

结论：三个 SSE 实现（`extz2_sse`/`extd2_sse`/`gg2_sse`，也是 KSW2 推荐使用的接口）
与 C 持平（0.93–1.04，即 ±7% 内，差异在测量噪声量级）；标量 `extz`/`gg` 反而更快；
只有标量 `extd`/`gg2` 的**仅打分**路径约慢 1.35–1.38×（CIGAR 路径持平）——这两条是
纯标量、仅用于对照/教学意义的旧路径（KSW2 自身也不推荐），故未继续压榨。

为了得到上表，Rust 侧需与 C 侧同样的编译选项（脚本已经这样做）：

```
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

## 与 C 实现的差异（刻意为之的部分）

* **内存**：C 用 `kalloc` 内存池 + 裸指针；Rust 用 `Vec`，`KswExtz::cigar`
  会在多次调用间复用容量（等价于 C 的 `m_cigar`/`n_cigar`）。C 的 `km`/`-K`
  内存池参数不再需要，`cli` 的 `-K` 选项被接受但无作用。
* **未定义行为**：C 在越界输入（如 `query[i] >= m`、band 宽小于长度差、
  `ksw_extd2_sse` 的 `-s` 之外场景）下会读越界内存；Rust 版会 panic 或给出
  确定的（但可能与某次 C 运行的 UB 结果不同）结果。所有合法用法下两者一致。
* **C 中的隐藏缓冲越界**：SSE 版本打分循环按 16 字节为一组读写 `s[]`/`qr[]`，
  C 通过 `kcalloc(..., tlen_*6 + qlen_ + 1, 16)` 的 `+1` 个 16 字节块兜住了越界；
  Rust 版显式补上同样的 padding，避免 UB。
* **`uint8_t` 转换**：`extd2`/`exts2` 中的边界值在 C 里是 `_mm_cvtsi32_si128((uint8_t)x1)`
  （零扩展），而 `extz2`/`gg2` 是 `_mm_cvtsi32_si128(x1)`（符号扩展）。两者行为不同，
  Rust 版逐处保留了这个差异——这是早期调试中最容易写错的地方。
* **`extd2`/`exts2` 的近似 Z-drop**：C 把 `ksw_apply_zdrop()` 放在 `if (r > 0)`
  之外（`extz2` 放在里面），会影响 `max_t`/`max_q`；已照搬。
* **`extd2_sse` 的 `qe`**：C 在 q/q2 交换*之前*计算 `int qe = q + e`
  （ksw2_extd2_sse.c:68），r==0 时 `H[0] = v8[0] - qe` 用的是交换前的值；
  已照搬（交换发生时 r=0 的分值会与"教科书"定义差 `(q+e)-(q2+e2)`）。
* **`ksw_extd` 的平局处理**：`ksw2_extz.c` 的右对齐分支用 `max >= h`，
  而 `ksw2_extd.c` 三个分支都用 `max > h`（影响 `mte_q`）；已照搬。

### 复查中发现的 C 库自身缺陷（移植未照搬的部分）

模糊对照还顺带暴露了三处 C 库自身的缺陷，Rust 版行为如下：

1. 标量 `ksw_gg`/`ksw_extz`/`ksw_extd` 的回溯矩阵用 `kmalloc`（未初始化），
   且没有 SSE 版的 `off_end[]` 保护：当带宽不足以覆盖比对末端
   （`|qlen-tlen| > w`）时，C 会读越界内存（ASan 可检出的 UB），
   输出取决于堆上的随机内容；Rust 版对应场景会 panic。
   两者对合法输入（带宽覆盖比对）完全一致。
2. `ksw_gg2_sse` 的回溯矩阵同样未初始化：窄带下其 CIGAR 在 C 中
   *两次运行结果不同*（分数确定且与 Rust 一致，路径不确定）。
   Rust 版为确定的零填充语义。模糊测试的 C oracle 构建时把
   回溯矩阵分配改为 `kcalloc()` 以获得可复现的对照（见
   `scripts/gen_fuzz_expected.sh` 顶部的说明；这只影响 C 读到
   "从未写入的单元"时的行为，不影响 DP 本身）。
3. `ksw2.h` 的 `ksw_cigar2eqx()` 调用 `ksw_push_cigar()` 时丢弃了
   返回的（可能已 realloc 的）指针：EQX 展开超过原缓冲容量时，
   C 会通过悬垂指针读写（ASan 检出 heap-buffer-overflow）。
   Rust 版用 `Vec` 不存在此问题。模糊测试的 C oracle 通过预放大
   CIGAR 缓冲绕开该缺陷。

## 构建与运行

```
cargo build --release
target/release/ksw2-test -t extd2_sse ksw2/test/MT-human.fa ksw2/test/MT-orang.fa

# 与 C 版对比
./scripts/compare_with_c.sh          # 小/中数据集
./scripts/compare_with_c.sh large    # 追加 16.5 kb 数据

# 重新生成 golden 期望输出（需要先 make 出 C 版二进制）
(cd ksw2 && make) && ./scripts/gen_expected.sh

# 性能
./scripts/bench.sh
```

依赖：`flate2`（读取 `.gz` 输入，仅 CLI 使用）。只需要库时可用
`cargo build --release --no-default-features`，此时没有任何外部依赖。

## 许可证

移植自 KSW2（MIT，见 [`ksw2/LICENSE.txt`](ksw2/LICENSE.txt)），本项目同样以 MIT 发布。
如在使用本项目的工作中引用 KSW2，请同时引用 Suzuki & Kasahara (2018) 与 Li (2018)。
