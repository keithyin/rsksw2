# 延伸比对（Extension Alignment）

一句话：**延伸比对是「种子-延伸」（seed-and-extend）比对策略的第二步——从一个已知的匹配锚点出发，
用动态规划沿序列方向延伸比对，分数掉得太狠就提前停手**。KSW2 的 `ext*` 系列函数就是为这个场景
设计的（minimap2 的核心比对引擎用的正是 KSW2）。

## 1. 为什么需要它：从读段比对（read mapping）说起

把一条测序 read 贴回参考基因组时，直接做全长动态规划有两个问题：

- **全局比对（Needleman-Wunsch，本仓库的 `gg` 系列）**要求两条序列从头到尾全部对上，两端必须
  到达。但 read 可能只来自基因组的某个局部区域，两端根本不该被强行比对；而且对 3 Gb 的基因组
  逐细胞填矩阵完全不可行。
- **局部比对（Smith-Waterman）**允许任意位置重新开始（分数与 0 取 max），能找到最佳子段，但它
  不知道「从哪里开始」——对长 read 来说搜索空间仍然太大。

实际的长读段比对器采用两阶段策略：

1. **找种子（seeding/chaining）**：先用 minimizer、MEM 等廉价方法找到 read 与参考之间的一串
   精确匹配短片段（锚点），并筛出共线的一组。这一步只有哈希查表，非常快。
2. **延伸（extension）**：锚点之间的空隙、以及 read 两端悬空的部分，交给动态规划去「补全」——
   从锚点出发向外延伸，允许错配和 gap，直到延伸不再划算为止。

`extz` 做的正是第 2 步。[`ksw2/ksw2.h`](../ksw2/ksw2.h#L45) 里对它的注释就是 "NW-like
extension"：**形式上像全局比对（起点固定在原点、没有 SW 那种清零重启），但终点自由、且可以
中途放弃**。

## 2. 三种模式的精确区别

| | 全局比对 `gg` | 局部比对 SW | 延伸比对 `extz` |
|:--|:--|:--|:--|
| 起点 | 固定 (0,0) | 任意（分数可与 0 取 max 重启） | 固定 (0,0)，即锚点位置 |
| 终点 | 必须到达两序列末端 | 任意最高分细胞 | 自由：停在历史最高分处即可 |
| 提前终止 | 无 | 无 | **Z-drop**：掉分超阈值即停 |
| 输出 | 一个最终分数 | 最佳子段 | `max`（全程最高分）+ 坐标 + 是否 zdropped |

对应到代码：[`src/extz.rs`](../src/extz.rs#L236) 只有在**同时**到达 query 末端与 target 末端时
才给 `result.score` 赋值，否则 `score` 保持 `NEG_INF`；真正有用的结果是 `max`/`max_t`/`max_q`
（延伸到最优点的分数和坐标），以及 `mqe`/`mte`（到达 query 末端 / target 末端时的最好分数——
用来判断延伸是否覆盖了整条 read）。

## 3. 具体例子：Z-drop 如何「及时止损」

用 [README](../README.md#L185) 中的例子，打分 match +2、mismatch −4，`z_drop = 3`：

```
query:  A A A A A G G G G G
target: A A A A A C C C C C
```

沿主对角线逐格填 DP：

| 位置 | 事件 | 对角分数 H | 历史最高 max |
|:--|:--|--:|--:|
| 0–4 | 5 个匹配 | 2→4→6→8→**10** | 10，坐标 (4,4) |
| 5 | 第一个错配 G/C | 10−4 = **6** | 仍为 10 |

在位置 5 处检查 Z-drop 条件（[`ksw2/ksw2.h`](../ksw2/ksw2.h#L201)）：

```
max − H > zdrop + l·e
10 − 6 = 4 > 3 + 0·2 = 3   → 成立，停止延伸
```

于是 `zdropped = true`，DP 循环 break（[`src/extz.rs`](../src/extz.rs#L225)），CIGAR 只回溯到
最优点 (4,4)，输出 `5M`，`score = NEG_INF`。生物学含义：前 5 个碱基的锚延伸是可信的，后面一片
错配说明已经伸出了真正的同源区，再算下去浪费时间还会产生垃圾比对。

公式里的 `l·e` 值得单独说：`l` 是当前对角线与历史最优点之间 **indel 长度之差**
（`l = |Δt − Δq|`），`e` 是 gap 延伸罚分。也就是说掉分阈值会随「合理的 gap 代价」放宽——一个
长 insertion 本来就要花 `l·e` 的罚分，不应该被误判为「同源区结束」。这是 Z-drop 比朴素
X-drop（`extf2_sse` 用的那种，纯分数差阈值）更精细的地方。

## 4. 双向延伸：从锚点向两侧延伸

单次调用只能沿一个方向延伸：DP 原点固定在 (0,0)，只向索引增大的方向填表
（[`src/extz.rs`](../src/extz.rs#L97)），没有「负坐标」的概念。向锚点两侧延伸是**两次调用
拼出来的**——KSW2 的 `Flags::REV_CIGAR` 就是为这个工作流准备的
（[`src/lib.rs`](../src/lib.rs#L178)）。

设 0-based 坐标，种子（锚点）落在 `target[200..200+k]` ↔ `query[100..100+k]`：

```
target:  [0 ......... 199][种子 200..200+k][... 到尾]
query:   [0 ..... 99][种子 100..100+k][... 到尾]
          ←── 左延伸（翻转后调用）         右延伸 ──→
```

**向右（朝尾部）**：后缀直接调用，此次调用的 (0,0) 就是锚点。

```rust
extz2_sse(&query[100 + k..], &target[200 + k..], &scoring,
          band, z_drop, 0, Flags::EXTZ_ONLY, &mut right);
```

**向左（朝头部）**：把锚点之前的**前缀做反向互补**，在翻转坐标系里同样「向右」延伸；
`REV_CIGAR` 让输出的 CIGAR 反序，直接拼回正向。

```rust
let q_flip = revcomp(&query[..100]);   // revcomp 由调用方提供
let t_flip = revcomp(&target[..200]);
extz2_sse(&q_flip, &t_flip, &scoring,
          band, z_drop, 0,
          Flags::EXTZ_ONLY | Flags::REV_CIGAR, &mut left);
```

最终 CIGAR = 左 CIGAR + 种子 `kM` + 右 CIGAR（相邻同类操作可合并）。翻转调用里的坐标
`(t', q')` 还原为原坐标 `(199−t', 99−q')`。

为什么是反向**互补**而不是简单反转：DNA 是双链，反向读左翼等于换到互补链，需 A↔T、C↔G
才能保持匹配语义（编码 0↔3、1↔2）。KSW2 与字母表无关，`revcomp` 由调用方提供；非 DNA
字母表纯反转即可。CIGAR 只需反序、不用交换 I/D——两条序列同时 revcomp 后，query 多出的
碱基仍是 I。

一个跑通的例子（种子 `GGCC`，左翼一个错配，match +2 / mismatch −4）：

```
target: ACGTAC GGCC TTAGCA
query:  ACGTTC GGCC TTAGCA
```

- 右调用：`6M`，score 12；左调用（revcomp 前缀 + `REV_CIGAR`）：`6M`，score 6
- 拼接：`6M + 4M + 6M`，总分 `6 + 8 + 12 = 26 = 15×2 − 4`

两个方向各自独立，各有自己的 `max`/`zdropped`/`score`。注意 `EXTZ_ONLY` 的语义：回溯从
最优点 `(max_t, max_q)` 出发，即延伸只做到还有利可图的地方为止——左翼若出现把分数拉低的
indel，CIGAR 会在最优点截断（与 Z-drop 一致的止损行为）；不加 `EXTZ_ONLY` 且两端都到达
时，则会强制对齐整个前缀。`extz2_sse` 的 `end_bonus` 还能在 `mqe + end_bonus > max` 时把
终点拉到 query 末端（`reach_end`）。

## 5. 函数命名对照

- **`extz` vs `extz2_sse`**：同一个算法的两条实现路径。`extz` 是标量、Green 经典 DP 形式；
  `extz2_sse` 用 Suzuki 差分形式 + SSE4.1 向量化，是 KSW2/minimap2 实际推荐使用的接口，两者
  结果逐字节一致（本仓库的一致性测试就在验证这件事）。
- **「单段仿射罚分」**：长度为 l 的 gap 花费 `gap_open + l·gap_extension`，一条直线。对应的
  `extd`（dual）用两段仿射取小：短 gap 罚得重、长 gap 斜率更缓，更符合长读段中真实 indel 的
  分布；`exts` 再加剪接感知（超长 gap 封顶为内含子，输出 `N`）；`extf` 是线性 gap + X-drop 的
  仅打分快速路径。
- **`Flags::EXTZ_ONLY`**：CIGAR 回溯从最优点 `(max_t, max_q)` 出发走回 (0,0)，而不是强行对齐
  到两序列末端——即「只要延伸段」。配合 `REV_CIGAR` 可做双向延伸，见第 4 节。

## 6. 小结

延伸比对 = **固定起点（锚点）+ 自由终点 + 仿射 gap + Z-drop 止损**的 NW 型动态规划，是把廉价
的精确种子变成完整比对的关键一步。

API 用法、参数与结果字段的完整说明见 [README 的延伸比对章节](../README.md#延伸比对单段仿射extz--extz2_sse)。