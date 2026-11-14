# analysis

下游分析与绘图的重构版本。输入是 `run_pipeline.R` 写出的网络结果。
可以用 `--network-dir` 把队列网络和 `--results` 分析输出分开；不指定时
退回旧布局（两者共用一棵目录树）。原目录
`/home/work/wenhai/meta/high_rare_net4` 未做任何修改。

## 范围

覆盖 `scripts/final_figures/README.md` 列出的最终图件脚本：

| 旧脚本 | 新位置 | 输出 |
|---|---|---|
| `stab_plot.R` | `R/01_stability.R` | `fig_net_1.pdf`、`stab/*` |
| `node_compare.R` | `R/02_node_overlap.R` | `fig_net_2.pdf`、`net_compare/node_overlap/*` |
| `degree_change_analysis_V2.R` | `R/03_degree_change.R` | `fig_net_3.pdf`、`mean_degree_change/*` |
| `net_compare2_V3.R` | `R/04_net_compare.R` | `fig_net_4.pdf`、`fig_net_4m.pdf`、`net_compare2_V3/*` |
| `rho_plot.R` | `R/05_rho.R`（`run_rho_global`） | `fig_net_5.pdf`、`rho_compare/*` |
| `rho_edge_plot.R` | `R/05_rho.R`（`run_rho_by_edge_type`） | `fig_net_6.pdf`、`rho_compare/*` |
| `meta_consensushub_abundance_count_V3.R` | `R/06_consensus_hub.R`（`run_meta_consensus_hub`） | `consensus_hub/*` |
| `plot_consensus_hub_rare_abundant.R` | `R/06_consensus_hub.R`（`run_hub_summary_plot`） | `ConsensusHub_RareAbundant_*` |
| `final_main_figure_from_rds_V1.R` | `R/07_final_figure.R` | `fig_net_main_2.pdf`、`main_split/*` |
| `02_utils.R`（部分函数） | `R/00_utils.R` | — |

`ERAS_plot_merge.R`、`prevance_curve.R`、`k_core_plot.R`、`species_info.R`、
`xishu.R`、`matrix_type.R`、旧的 `net_compare.R`/`plot_merge.R`/`07`–`09`
等脚本不在最终图件映射中，未纳入。

## 目录

```
run_analysis.R          入口；--step 选择步骤
config/analysis.yaml    结果根目录、后缀、跳过的疾病、hub 参数、路径
config/analysis.R       同一份配置的 R 格式示例（仅 --config 显式指定时使用）
R/00_utils.R            队列表、疾病名映射、输出路径、通用小工具
R/01_stability.R        稳定性 vs 样本量（fig_net_1）
R/02_node_overlap.R     节点重叠（fig_net_2）
R/03_degree_change.R    平均 degree 变化（fig_net_3）
R/04_net_compare.R      边保留/丢失（fig_net_4/4m）与 bundle RDS
R/05_rho.R              rho 总表与分边类型（fig_net_5/6）
R/06_consensus_hub.R    consensus hub 重算与汇总表出图
R/07_final_figure.R     组合主图与分面板
verification/           旧/新结果对照脚本
```

## 运行

```bash
export LD_LIBRARY_PATH=/home/wenhai/miniconda3/envs/r_4.3/lib
export R_LIBS=/home/wenhai/miniconda3/envs/r_4.3/lib/R/library:/home/wenhai/application/Qmatey/tools/R

# 全部步骤（不含 hub_summary）
Rscript analysis/run_analysis.R --step all

# 网络根目录与分析输出目录分开
Rscript analysis/run_analysis.R \
  --network-dir /path/to/network \
  --results /path/to/analysis \
  --suffix res_rc \
  --step net_compare,rho,rho_edge

# 单独的 hub 汇总表出图（需要外部表）
Rscript analysis/run_analysis.R --step hub_summary \
  --input /path/ConsensusHub_master_summary.csv \
  --outdir /path/hub_summary_out
```

步骤名：`stability`、`node`、`degree`、`net_compare`、`rho`、`rho_edge`、
`hub`、`hub_summary`、`final`。`--step all` 按顺序执行除 `hub_summary`
以外的全部步骤；`final` 依赖 `net_compare` 与 `hub` 的 RDS 输出。

也可以直接用仓库根目录的 `run.sh`（会先补齐 `scripts/` 下的映射文件和
输出目录，并记录日志）：

```bash
./run.sh analysis                 # 默认读 runs/full/network，写到 runs/full/analysis
./run.sh full                     # 先建 73 队列 + ERAS，再跑全部分析
ANALYSIS_STEP=net_compare,rho ./run.sh analysis
```

`--results` 是分析输出根目录，`--network-dir` 是队列网络根目录。不传
`--network-dir` 时退回旧行为（两者都取自 `results_dir`）。`--results` 不写时
使用配置里的旧结果根目录，会覆盖已有表格和图件；首次运行请显式指定新目录。

## 配置要点

默认配置是 `config/analysis.yaml`；`--config` 也接受 R 列表格式
（`config/analysis.R` 是同一份默认值的示例）。

- `results_dir`：分析输出根目录；默认旧目录。
- `network_dir`：队列网络根目录；`null` 表示与 `results_dir` 相同（旧布局）。
- `suffix`：结果后缀，可给多个（`--suffix a,b`）。旧脚本只处理 `res_rc`；
  依赖单后缀的步骤取第一个。
- `prevalence`：只用于拼 prevalence 文件名，保持旧行为。
- `skip_diseases`：每个步骤跳过的疾病，复现旧脚本各自的 `next` 条件。
- `hub`：consensus hub 的 top 比例、最小 degree、投票阈值、配色、柱间距。
- `paths`：输出子目录名，属于已发表图件布局，改名前需重新核对。
- `disease_map_file`：默认 `<results_dir>/scripts/disease_name_mapping.tsv`。

## 验证

`verification/compare_analysis.R` 对照旧/新两棵结果树。比较内容：

- `.csv`/`.tsv`：按表格值比较；
- `.rds`：数据框列名/行数/因子水平/数值（容差 1e-12）、列表结构；
- `.txt`：逐行；
- `ggsave()` 捕获的 ggplot/patchwork 对象：固定种子渲染成 PNG 后逐像素比较；
- `.png`/`.pdf`：跳过（设备/时间戳差异）。

合成 fixture 覆盖 7 个疾病 × 2 个队列（其中一个疾病名需要映射，例如
`Alzheimer → AD`）加 ERAS 5 个时间点；`net_compare` 的 ERAS 分支、
多队列 Control 比较、主图/附图切分、hub 主图/附图、最终组合图都会走到。
最终组合图额外用已发表的 5 个输入 RDS 验证，保证两边输入完全一致。

结果见 `REVIEW.md`，复现步骤见 `verification/README.md`。
