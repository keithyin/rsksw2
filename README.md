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

| C 源文件 | Rust 模块（私有移植层） | 公开函数 | 说明 |
|:--|:--|:--|:--|
| `ksw2_gg.c` | `src/gg.rs` | `gg` | 全局比对，Green 经典形式 |
| `ksw2_gg2.c` | `src/gg2.rs` | `gg2` | 全局比对，Suzuki 对角形式 |
| `ksw2_gg2_sse.c` | `src/sse/gg2_sse.rs` | `gg2_sse` | 全局比对 + SSE4.1 |
| `ksw2_extz.c` | `src/extz.rs` | `extz` | 扩展比对，单段仿射罚分 |
| `ksw2_extz2_sse.c` | `src/sse/extz2_sse.rs` | `extz2_sse` | 扩展比对 + SSE4.1 |
| `ksw2_extd.c` | `src/extd.rs` | `extd` | 扩展比对，双段仿射罚分 |
| `ksw2_extd2_sse.c` | `src/sse/extd2_sse.rs` | `extd2_sse` | 双段仿射 + SSE4.1 |
| `ksw2_extf2_sse.c` | `src/sse/extf2_sse.rs` | `extf2_sse` | 线性罚分 + X-drop，仅打分 |
| `ksw2_exts2_sse.c` | `src/sse/exts2_sse.rs` | `exts2_sse` | 剪接感知扩展 + SSE4.1 |
| `cli.c` | `src/bin/ksw2-test.rs` | 可执行文件 `ksw2-test` | 与 C 版命令行完全兼容 |
| `ksw2.h`（结构体/工具函数） | `src/lib.rs`, `src/common.rs` | `Scoring`, `ExtzResult`, `Flags`, `CigarOp` 等 | 打分系统、结果结构、标志位、回溯等 |

crate 分为两层：**公开 API 收敛在 crate root**（`src/lib.rs`），签名是 Rust 风格的
（slice 自带长度、`Scoring`/`Flags`/`Option` 参数）；各模块内部则是与 C 逐行对应的
**私有移植层**——保留 C 形态的签名（长度参数、`-1` 哨兵、裸位标志），方便与 C 源码
逐行 diff 审查。两层语义完全一致。

`*_sse` 系列是 C 源文件中 `__SSE4_1__` 分支的逐行移植（即 `-march=native`
时编译出来的那条路径），内部使用 `core::arch::x86_64` 的 SSE4.1 intrinsics，
因此与 C 版本在 score、坐标和 CIGAR 上完全一致。它们与 C 一样只在 x86_64 上编译，
并在运行时检查 `sse4.1`（2008 年之后的 x86 CPU 均支持）。

## 用法

```rust
use rsksw2::*;

// 打分系统：5 种符号、match +2、mismatch −4、wildcard 行/列 0、gap 罚分 (4, 2)
let score_matrix = gen_simple_mat(5, 2, -4);
let scoring = Scoring::new(5, &score_matrix, 4, 2);
let query:  &[u8] = &[0, 1, 2, 3];   // A C G T
let target: &[u8] = &[0, 1, 2, 3];

let mut result = ExtzResult::new();  // 可跨调用复用，内部 CIGAR 缓冲会被重复利用
extz2_sse(query, target, &scoring,
          None,   // band_width: Option<i32>，None = 不限带宽
          None,   // z_drop: Option<i32>，None = 禁用
          0,      // end_bonus
          Flags::empty(), &mut result);
println!("score = {}, cigar = {}", result.score, result.cigar_string());
```

### 函数签名对照

C 长度参数（`qlen`/`tlen`）被 `&[u8]` slice 取代——slice 自带长度，不存在
「长度与 slice 不一致」这类 C 的未定义行为。

`ksw_extz_t` → `ExtzResult`：

| C | Rust |
|:--|:--|
| `uint32_t max:31` + `int zdropped:1` | `max: u32`（31 位，保留与 C 相同的截断语义）+ `zdropped: bool` |
| `max_q`, `max_t`, `mqe`, `mqe_t`, `mte`, `mte_q`, `score` | 同名字段，类型 `i32` |
| `int reach_end` | `reach_end: bool` |
| `m_cigar` / `n_cigar` / `cigar` | `cigar: Vec<u32>`（`n_cigar()` 返回长度，容量自动复用） |
| `ksw_reset_extz()` | `ExtzResult::reset()` / `ExtzResult::new()` |

其余 C 习惯的映射：

| C | Rust |
|:--|:--|
| `m` + `mat` + `gapo` + `gape` 四个散参数 | `Scoring` 结构体（`Scoring::new(m, &mat, gapo, gape)`） |
| `w` / `zdrop` / `xdrop`，负数表示禁用 | `Option<i32>`，`None` 表示禁用/不限 |
| `flag` 位掩码 + `KSW_EZ_*` 宏 | `Flags` 类型 + `Flags::SCORE_ONLY` 等，用 `\|` 组合 |
| `ksw_gg`/`ksw_gg2` 的「是否输出 CIGAR」三指针 | `Option<&mut Vec<u32>>`：`None` 等价 `-s`（只算分数） |
| `KSW_NEG_INF` | [`NEG_INF`](https://docs.rs/rsksw2) |
| `ksw_push_cigar` 操作码 `KSW_CIGAR_*` | `CigarOp` 枚举（`Match`/`Ins`/`Del`/`NSkip`/`Eq`/`X`） |

## API 参考

本节按功能分组介绍 crate 公开的全部 API。除特别说明外，语义与对应 C 函数一致；
与 C 的全部有意差异见[与 C 实现的差异](#与-c-实现的差异刻意为之的部分)，完整签名见
[docs.rs/rsksw2](https://docs.rs/rsksw2)。下文示例统一使用：

```rust
use rsksw2::*;

// A/C/G/T/N → 0/1/2/3/4（N 是 wildcard，对应矩阵的最后一行/列）
let enc = |s: &str| -> Vec<u8> {
    s.bytes().map(|b| match b {
        b'A' => 0, b'C' => 1, b'G' => 2, b'T' => 3, _ => 4,
    }).collect()
};
let mat = gen_simple_mat(5, 2, -4); // 匹配 +2，错配 −4，wildcard 0
let scoring = Scoring::new(5, &mat, 4, 2);
```

### 通用约定

* **序列编码**：`query`/`target` 是 `&[u8]`，每个字节是一个符号编号，必须
  `< scoring.alphabet_size`（剪接函数还假定 `0/1/2/3 = A/C/G/T`）。越界符号在 C 中是
  未定义行为，Rust 中会 panic。
* **打分系统**：`Scoring` 打包符号种数 `alphabet_size`、`alphabet_size × alphabet_size`
  的行主序打分矩阵与 gap 罚分；`score_matrix[i·m+j]` 为 target 符号 `i` 对 query 符号
  `j` 的得分。`Scoring::new` 会校验矩阵形状与 `alphabet_size ∈ 1..=127`。注意：
  **未设置 `Flags::GENERIC_SC` 时**，三个 SSE 扩展函数只用 `score_matrix[0]`（符号相同）
  与 `score_matrix[1]`（符号不同）构造快速路径，并把符号 `alphabet_size−1` 视为
  wildcard——其得分为 `mat[m·m−1]==0 ? −e : mat[m·m−1]`（C 原版行为）。标量函数永远读
  完整矩阵，不受该标志影响。
* **罚分**：`gap_open`/`gap_extension` 等按正数传入；长度 `l` 的 gap 花费
  `gap_open + l·gap_extension`。
* **带宽 `band_width: Option<i32>`**：`None` 表示不限制。要得到两端都到达的完整结果需
  `w ≥ |qlen−tlen|`；带宽过窄时 C 会读未初始化内存（UB），Rust 会 panic 或给出确定性结果。
* **`z_drop`/`x_drop: Option<i32>`**：`None` 禁用。
* **`flags: Flags`**：用 `\|` 组合，如 `Flags::EXTZ_ONLY | Flags::REV_CIGAR`；
  `Flags::from_bits(u32)` 接受 C 风格的原始位值（未知位原样保留）。
* **`*_sse` 函数**仅在 x86_64 上编译（与 C 的 `#ifdef __SSE2__` 一致），运行时缺 sse4.1 会 panic。
* **前提检查**：SSE 扩展函数在矩阵最低分 `< −2(gap_open+gap_extension)` 时直接返回空结果
  （保证内部 i8 差分不饱和，C 同）；`exts2_sse` 还要求 `gap_open2 > gap_open+gap_extension`。
* **结果复用**：扩展函数入口等价于 `ksw_reset_extz()`；`ExtzResult` 可 `new()` 一次后跨调用复用
  （CIGAR 缓冲容量保留，等价 C 的 `m_cigar` 复用）。
* **哨兵值**：`score == NEG_INF`（−0x40000000）表示未得到完整对齐分数。

### 全局比对：`gg` / `gg2` / `gg2_sse`

| 参数 | 说明 |
|:--|:--|
| `query`/`target` | 两条序列（长度即 slice 长度） |
| `scoring` | 打分系统（`Scoring::new`） |
| `band_width` | `Option<i32>`，`None` 不限制；须 `≥ |qlen−tlen|` |
| `cigar` | `gg`/`gg2` 为 `Option<&mut Vec<u32>>`：`None` 只算分数；`gg2_sse` 必须传 `&mut Vec<u32>`（总是回溯） |

返回全局比对分数；CIGAR 为 BAM 编码（`len<<4|op`）。

```rust
let (q, t) = (enc("ACGTACGT"), enc("ACGTACGT"));
assert_eq!(gg(&q, &t, &scoring, None, None), 16);

let mut cigar = Vec::new();
let s = gg2_sse(&q, &t, &scoring, None, &mut cigar);
assert_eq!((s, cigar), (16, vec![8 << 4])); // "8M"
```

### 扩展比对（单段仿射）：`extz` / `extz2_sse`

在全局比对参数基础上增加 `z_drop`、`flags`（`extz2_sse` 另有 `end_bonus`）：

| 参数 | 说明 |
|:--|:--|
| `z_drop` | `Option<i32>`，`None` 禁用；触发后 `zdropped=true` 并停止计算 |
| `end_bonus` | 仅 `extz2_sse`：配合 `Flags::EXTZ_ONLY` 使用，见下 |
| `flags` | `Flags` 组合，见[标志位](#标志位-flags) |

输出语义（字段详见[结果结构](#结果结构-extzresult)）：

* `score` 只有在计算**同时到达 query 末尾与 target 末尾**时才被赋值（该最终单元的分数）；
  Z-drop 提前终止或带宽不足时保持 `NEG_INF`。它不是路径最大值——最大值在 `max` 里。
* `mqe`/`mte` 分别记录到达 query 末端、target 末端时的最大分数及坐标。
* `Flags::EXTZ_ONLY`：CIGAR 回溯到 `(max_t, max_q)` 而非 `(0,0)`。`extz2_sse` 中若
  `mqe + end_bonus > max`，则改为回溯到 query 末端并把 `reach_end` 置 `true`。
* 标量 `extz` 无 `end_bonus`/`reach_end`，且忽略 `APPROX_*` 标志。

```rust
// query 比 target 多 3 个碱基 → 4M3I4M
let (q, t) = (enc("ACGTGGGACGT"), enc("ACGTACGT"));
let mut ez = ExtzResult::new();
extz(&q, &t, &scoring, None, None, Flags::empty(), &mut ez);
assert_eq!((ez.score, ez.cigar_string().as_str()), (6, "4M3I4M")); // 8×2 − (4+3×2)

// Z-drop：连续错配时提前终止，CIGAR 只回溯到最优点为止
let (q, t) = (enc("AAAAAGGGGG"), enc("AAAAACCCCC"));
let mut ez = ExtzResult::new();
extz(&q, &t, &scoring, None, Some(3), Flags::empty(), &mut ez);
assert!(ez.zdropped);
assert_eq!((ez.max, ez.max_t, ez.max_q), (10, 4, 4)); // 止步于 5 个匹配
assert_eq!(ez.cigar_string(), "5M");
assert_eq!(ez.score, NEG_INF);

// EXTZ_ONLY + end_bonus：延伸到 query 末端比内部最优更优时 reach_end=true
let (q, t) = (enc("ACGTGGGACGT"), enc("ACGTACGT"));
let mut ez = ExtzResult::new();
extz2_sse(&q, &t, &scoring, None, None, 5, Flags::EXTZ_ONLY, &mut ez);
assert!(ez.reach_end);                                  // mqe(6) + 5 > max(8)
assert_eq!((ez.score, ez.max, ez.max_t, ez.max_q), (6, 8, 3, 3));
```

### 扩展比对（两段仿射）：`extd` / `extd2_sse`

gap 费用为两段仿射取小：`cost(l) = min{gap_open + l·gap_extension, gap_open2 + l·gap_extension2}`。
第二段罚分（`gap_open2`/`gap_extension2`）作为独立参数传入。典型取法是
`gap_open+gap_extension < gap_open2+gap_extension2` 且 `gap_extension > gap_extension2`，
让短 gap 走第一段、长 gap 走第二段。Z-drop 以 `gap_extension2` 为斜率。
`extd2_sse` 是唯一支持 `Flags::EQX`（CIGAR 拆 `=`/`X`）的入口。

```rust
let mut ez = ExtzResult::new();
extd2_sse(&enc("ACGT"), &enc("AGGT"), &scoring, 13, 1,
          None, None, 0, Flags::EQX, &mut ez);
assert_eq!((ez.score, ez.cigar_string().as_str()), (2, "1=1X2=")); // 3×2 − 4
```

### 延伸打分（线性 gap + X-drop）：`extf2_sse`

* 参数是 `match_score`/`mismatch_score`/`gap_extension` 而非 `Scoring`；序列按字节相等比较，
  无 wildcard；`mismatch_score` 自动取负。
* gap 长度 `l` 花费 `l·gap_extension`（无 gap open）。
* `x_drop: Option<i32>`：`Some(x)` 且当前对角分比历史最大分低超过 `x` 时终止并置
  `zdropped=true`。
* 只打分：不产 CIGAR，`mqe`/`mte`/`reach_end` 恒为初值；`score` 只有跑完全部反对角线
  （两端到达）才有值。

```rust
let mut ez = ExtzResult::new();
extf2_sse(&enc("AAAAAGGGGG"), &enc("AAAAACCCCC"), 1, -1, 1, None, Some(3), &mut ez);
assert!(ez.zdropped);                                   // 第 4 个错配处 5−1 > 3
assert_eq!((ez.max, ez.max_t, ez.max_q), (5, 4, 4));
assert_eq!(ez.score, NEG_INF);

let mut ez = ExtzResult::new();
extf2_sse(&enc("AAAAAGGGGG"), &enc("AAAAACCCCC"), 1, -1, 1, None, None, &mut ez);
assert_eq!((ez.score, ez.zdropped), (0, false));        // 不设 X-drop：5×(+1)+5×(−1)
```

### 剪接延伸：`exts2_sse`

RNA 剪接感知的扩展：短 gap 按普通仿射罚分，超长 gap（内含子）的罚分**封顶**而非线性增长——
长度 `l ≤ long_thres ≈ (gap_open2−gap_open)/gap_extension − 1` 时按
`gap_open + l·gap_extension` 计费，超过后总罚分不再随长度增长（实测
`gap_open=2, gap_extension=3, gap_open2=10` 时 `l ≥ 3` 恒为 `gap_open2 = 10`）。
回溯时内含子输出 `N`（`CigarOp::NSkip`）。

剪接位点（源码假定编码 `0/1/2/3 = A/C/G/T`）：

* donor 位点看内含子**前两个**碱基：正链 `GT`、负链 `CT`；acceptor 位点看内含子
  **最后两个**碱基：正链 `AG`、负链 `AC`。由 `Flags::SPLICE_FOR`/`Flags::SPLICE_REV`
  启用（可同时）。
* 罚分：强信号 `GTr…`/`…yAG`（侧翼为嘌呤 r=A/G、嘧啶 y=C/T）罚 0；普通经典位点无
  `Flags::SPLICE_FLANK` 时罚 0、有则罚 `−non_canonical_penalty/2`（向零截断）；
  非经典罚 `−non_canonical_penalty`。
* `junctions`（可传空 `&[]`）按 target 每个碱基一个字节提供位点提示，命中位加
  `junction_bonus`（与位点是否经典无关）：

  | 位 | 含义 |
  |:--|:--|
  | bit0 (1) | donor（正链 GT-AG） |
  | bit1 (2) | acceptor（正链 GT-AG） |
  | bit2 (4) | acceptor（负链 CT-AC） |
  | bit3 (8) | donor（负链 CT-AC） |

  donor 位查询 `junctions[内含子首碱基下标]`，acceptor 位查询 `junctions[内含子末碱基下标]`；
  设置 `Flags::REV_CIGAR` 时（序列已先做反向互补）两者查询的位互换。
* `z_drop` 生效但斜率为 0（内含子长度不计入掉分）。

实测罚分构成（`gap_open=2, gap_extension=3, gap_open2=10, non_canonical=4`，
两侧 exon 全匹配共 8×2 分）：

| 场景 | 分值 |
|:--|:--|
| 经典内含子 `GTAAAAAG` | `16 − 10 = 6`，CIGAR `4M8N4M` |
| donor 非经典（内含子起 `CA…AG`） | `16 − 10 − 4 = 2` |
| 半经典 `GT…AG`（侧翼非 r/y） | 无 FLANK 与经典相同；有 FLANK 每个位点再 `−2` |
| 非经典 + `junction_bonus=3` 命中 | 位点罚分 `−4+3 = −1`（经典位点则为 `0+3 = +3`） |

```rust
let target = enc("ACGTGTAAAAAGGTAC"); // ACGT | GTAAAAAG | GTAC
let query  = enc("ACGTGTAC");
let splice_scoring = Scoring::new(5, &mat, 2, 3); // gapo=2, gape=3
let mut ez = ExtzResult::new();
exts2_sse(&query, &target, &splice_scoring, 10, 4, None, 0,
          Flags::SPLICE_FOR, &[], &mut ez);
assert_eq!(ez.cigar_string(), "4M8N4M");
assert_eq!(ez.score, 6);              // 8×2 − 10
```

### 结果结构 `ExtzResult`

对应 C 的 `ksw_extz_t`（C↔Rust 映射见[函数签名对照](#函数签名对照)）：

| 字段 | 初值 | 语义 |
|:--|:--|:--|
| `score` | `NEG_INF` | 两端都到达时的对齐分数（`extf2_sse` 为整条延伸的分数） |
| `max` | `0`（u32） | 全程最大单元分数（保留 C `uint32_t max:31` 位域的 31 位截断语义） |
| `max_t`、`max_q` | `-1` | `max` 所在的 (target, query) 坐标 |
| `mqe`、`mqe_t` | `NEG_INF`, `-1` | 到达 query 末端时的最大分数及其 target 坐标 |
| `mte`、`mte_q` | `NEG_INF`, `-1` | 到达 target 末端时的最大分数及其 query 坐标 |
| `zdropped` | `false` | 是否因 Z-drop / X-drop 提前终止 |
| `reach_end` | `false` | 仅 `extz2_sse`/`extd2_sse` + `EXTZ_ONLY`：`mqe+end_bonus > max` 时置 `true` |
| `cigar` | `[]` | BAM 编码 CIGAR；`n_cigar()` 返回长度；`cigar_string()` 渲染为 `10M2I5M` 形式 |

CIGAR 操作符用 `CigarOp` 枚举表示：`Match`(0)、`Ins`(1)、`Del`(2)、`NSkip`(3)、
`Eq`(7)、`X`(8)；`CigarOp::from_code(编码字低 4 位)` 解码，`cigar_string()` 已内置。

### 标志位 `Flags`

下表中「扩展函数」指 `extz` / `extz2_sse` / `extd` / `extd2_sse` /
`exts2_sse`（`extf2_sse` 没有 flags 参数）。

| 标志 | 值 | 生效范围 | 语义 |
|:--|:--|:--|:--|
| `Flags::SCORE_ONLY` | 0x01 | 扩展函数 | 不回溯、不产 CIGAR |
| `Flags::RIGHT` | 0x02 | 扩展函数 | gap 右对齐 |
| `Flags::GENERIC_SC` | 0x04 | 3 个 SSE 扩展函数 | 按完整矩阵查分（见[通用约定](#通用约定)）；标量版无效 |
| `Flags::APPROX_MAX` | 0x08 | SSE 扩展函数 | max 改用对角差分近似跟踪；此时只有再设 `APPROX_DROP` 才做 Z-drop |
| `Flags::APPROX_DROP` | 0x10 | SSE 扩展函数 | 近似模式下的 Z-drop（与 `APPROX_MAX` 成对使用，即 CLI 的 `-g`） |
| `Flags::EXTZ_ONLY` | 0x40 | 扩展函数 | 只延伸：CIGAR 回溯到最优/末端单元而非 (0,0) |
| `Flags::REV_CIGAR` | 0x80 | 扩展函数 | CIGAR 反序输出（配合「先反向互补再延伸」的工作流） |
| `Flags::SPLICE_FOR` / `Flags::SPLICE_REV` | 0x100 / 0x200 | `exts2_sse` | 启用正链 / 负链剪接信号 |
| `Flags::SPLICE_FLANK` | 0x400 | `exts2_sse` | 半经典位点罚 `−noncan/2` |
| `Flags::EQX` | 0x800 | 仅 `extd2_sse` | M 拆为 `=`/`X` |

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
   `Flags::EXTZ_ONLY`、`Flags::GENERIC_SC`、`Flags::EQX`、剪接 `junc` 数组等），
   覆盖全部 9 个函数共 2 万次调用；`cargo test` 用 Rust 重放同样用例并要求一致。

4. **实现间交叉验证**（`tests/impl_agreement.rs`）
   对随机序列比较 `gg` / `gg2` / `gg2_sse`，以及
   `extz` vs `extz2_sse`、`extd` vs `extd2_sse` 的
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

* **API 命名与参数形式**（0.2.0 起）：公开函数去掉 `ksw_` 前缀（`ksw_extz2_sse` →
  `extz2_sse`）；长度参数由 slice 自带长度取代；`m/mat/gapo/gape` 打包为
  [`Scoring`](https://docs.rs/rsksw2)；`w/zdrop/xdrop` 的负数哨兵改为 `Option<i32>`；
  `flag` 位掩码改为 `Flags`；`KswExtz` → `ExtzResult`（`reach_end` 为 `bool`）；
  `gen_simple_mat` 直接返回 `Vec<i8>`。私有移植层模块仍保留 C 形态签名，便于逐行
  对照 C 源码——**算法行为与 0.1.x 完全一致**（一致性/模糊测试逐字节通过）。
* **内存**：C 用 `kalloc` 内存池 + 裸指针；Rust 用 `Vec`，`ExtzResult::cigar`
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
