<p align="center">
  <strong>语言：</strong> <a href="README.md">English</a> | <a href="README.zh.md">中文</a>
</p>

# pantax-dbg-analysis

基于 `pantax-dbg` 产出的 profiling 数据，构建 74 个肠道菌群队列的微生物共现
网络（73 个疾病/对照队列 + ERAS 抗生素时间点队列），并分析网络稳定性、边保留、
hub 物种、稀有/丰富物种以及 ERAS 恢复等。另构建 20 个国家的核心物种健康
（Control）网络，并分析全球核心物种在各队列网络中的作用。

## 主要功能

- 使用 SPIEC-EASI（默认）、FlashWeave 或 SparCC 构建分组网络。
- 支持 read count 与相对丰度输入、流行率过滤、边权重/置信度过滤和合并队列模式。
- 输出每个队列的网络对象、节点表、中心性指标、边统计和丰度分档图。
- 使用核心物种矩阵为每个国家构建 Control-only 网络（`config/cohorts_core.tsv`）。
- 跨队列分析：稳定性 vs 样本量、节点重叠、degree 变化、边保留/丢失、rho 比较、
  consensus hub 和最终多面板图。
- 核心物种分析：健康网络边富集、疾病边丢失、定向移除稳健性、abundant–rare
  边中的核心稀有比例，以及 ERAS 抗生素恢复。
- 处理 ERAS 时间点（D0/D4/D8/D42/D180），其中 D0 重命名为 Control。

## 环境要求

- R 4.4（已在 4.4.1 测试）。
- 建网依赖：`SpiecEasi`、`igraph`、`RColorBrewer`、`data.table`、`ggplot2`、
  `reshape2`、`patchwork`、`tidyr`、`cowplot`、`dplyr`、`purrr`、`stringr`、
  `Matrix`、`ggsci`。
- 下游分析依赖：`igraph`、`dplyr`、`tidyr`、`ggplot2`、`patchwork`、
  `data.table`、`stringr`、`ggsci`、`scales`、`gridExtra`、`clinfun`、
  `RColorBrewer`；`pheatmap` 为可选项（缺失时跳过热图）。
- 核心物种分析：上述 `igraph`/`dplyr`/`tidyr`/`purrr` 加 `survival`，以及
  导出图片所需的 `ragg`/`svglite`。
- FlashWeave 和 SparCC 模式需要各自的外部工具；默认的 SPIEC-EASI 只需上述 R 包。

R 包路径可在 `config/pipeline.yaml` 的 `runtime.library_paths` 和
`runtime.spiec_easi_library` 中配置。集群上运行时按需设置
`LD_LIBRARY_PATH`/`R_LIBS`：

```bash
export LD_LIBRARY_PATH=/path/to/R/lib
export R_LIBS=/path/to/R/library
Rscript run_pipeline.R --help
```

## 目录结构

```
run_pipeline.R        建网入口（命令行与配置解析）
run.sh                常用场景和输出布局
R/                    建网与单队列流程
config/               默认配置、队列表和 ERAS 元数据
analysis/             下游分析与绘图入口
core_species/         核心物种健康网络与作用分析
verification/         重构代码的开发期对照检查
reproduction/         已发表运行的复现配方和代码快照
archive/              重构阶段的旧文档
```

## 输入数据

每个队列对应队列表（如 `config/cohorts.tsv`）中的一行，制表符分隔：

```text
disease  cohort  metadata_file  species_relative_abund  species_abund  maaslin3_res_file
```

- `cohort`：输出目录使用的队列名，例如 `PRJNA398089_598_CD`。
- `metadata_file`：样本元数据，需包含 `Sample_ID`、`Cohort` 和
  `Standard_Disease`（或 `Group`）列。
- `species_relative_abund`：相对丰度矩阵，行为样本、列为物种。
- `species_abund`：read count 矩阵，样本与相对丰度矩阵一致。
- `maaslin3_res_file`：可选的差异物种表；`NA` 表示没有。

也接受不含 `species_abund` 的 5 列格式。ERAS 队列使用单独的
`config/eras.tsv` 和元数据文件。

## 快速开始

构建全部非 ERAS 队列的网络：

```bash
Rscript run_pipeline.R \
  --config config/pipeline.yaml \
  --cohorts config/cohorts.tsv \
  --outdir results/network
```

把 ERAS 时间点网络写入同一棵树：

```bash
Rscript run_pipeline.R \
  --config config/eras.yaml \
  --outdir results/network
```

运行下游分析：

```bash
Rscript analysis/run_analysis.R \
  --network-dir results/network \
  --results results/analysis \
  --suffix res_rc \
  --step all
```

构建核心物种健康网络并分析核心物种作用：

```bash
# 20 个国家核心物种 Control 网络 -> results/all_country_core/
Rscript run_pipeline.R \
  --config config/core_network.yaml \
  --outdir results

# 对 results/network 中的 74 个队列网络做作用分析
Rscript core_species/run_core_species.R \
  --base-dir results/network \
  --step all
```

调试时只跑一个队列：

```bash
Rscript run_pipeline.R \
  --config config/pipeline.yaml \
  --cohorts config/cohorts.tsv \
  --only PRJNA398089_598_CD \
  --outdir results/network
```

等价的一键场景在 `run.sh` 中：

```bash
./run.sh spiec_rc          # read count + SPIEC-EASI + 10% 流行率
./run.sh eras              # ERAS 时间点
./run.sh core_network      # 20 个国家核心物种 Control 网络
./run.sh core_species      # 核心物种表格、图和汇总
./run.sh full              # 73 队列 + ERAS，再跑全部分析
./run.sh analysis          # 对已有结果树跑分析
DRY_RUN=1 ./run.sh full    # 只打印命令
```

`run.sh` 把队列网络写到 `<scenario>/network/`，下游结果写到
`<scenario>/analysis/`，日志和有效配置也保存在对应目录。需要自定义运行时，
直接给 `run_pipeline.R` 传 `--only`、`--outdir` 等参数即可。

## 配置

`config/pipeline.yaml` 包含以下默认值：

- 输入类型（`read_count` / `relative_abundance`）、流行率和过滤条件；
- 建网工具与 SPIEC-EASI 参数（含 StARS/pulsar 设置）；
- 丰度分档阈值、配色和图例文字；
- hub 定义；
- 输出目录/后缀和运行时库路径。

命令行参数优先于配置文件；任意嵌套字段可用 `--set`：

```bash
--set network/spiec_easi/pulsar/rep.num=50
--set abundance/thresholds/high=2e-3
```

`--print-config` 打印合并后的有效配置但不运行；`--save-config FILE` 先把配置
写成可再次读取的文件（按扩展名写 YAML 或 R 格式），然后继续运行。ERAS 使用
`config/eras.yaml` 和 `config/eras.tsv`。核心物种健康网络使用
`config/core_network.yaml` 和 `config/cohorts_core.tsv`，输出到
`<outdir>/all_country_core/`；核心物种分析有自己的配置
`core_species/config/core_species.yaml`，详见 `core_species/README.md`。

## 输出布局

```
results/
  network/
    PRJNA398089_598_CD/res_rc/
      g_Control_full_network.rds
      g_CD_full_network.rds
      Control_node_info.csv
      Control_centrality_metrics.csv
      ...
  analysis/
    stab/
    net_compare/node_overlap/
    mean_degree_change/
    net_compare2_V3/
    rho_compare/
    consensus_hub/
    scripts/final_figures/
```

后缀（默认 `res_rc`）由输入类型、过滤条件和建网工具推导，也可用 `--suffix`
显式指定。

## 下游分析

`analysis/run_analysis.R` 支持以下步骤：

```text
stability   稳定性 vs 样本量（fig_net_1）
node        节点重叠（fig_net_2）
degree      平均 degree 变化（fig_net_3）
net_compare 边保留/丢失（fig_net_4、fig_net_4m）
rho         rho 比较（fig_net_5）
rho_edge    分边类型的 rho（fig_net_6）
hub         consensus hub
hub_summary 单独的 hub 汇总图
final       组合主图与分面板
```

`--step all` 执行除 `hub_summary` 以外的全部步骤；`hub_summary` 需要外部汇总表
（`--input`/`--outdir`）。`--network-dir` 和 `--results` 可以指向同一棵树
（旧布局），也可以像上面一样分别指向网络和分析目录。

## 核心物种分析

`core_species/run_core_species.R` 读取相同的 `PR*/res_rc` 网络对象，输出：

```text
results/core_species/
  qc/ tables/ plot_data/          分析表格
  figures/                        主图和补充图
  manuscript_stats/               汇总表格和 markdown
```

使用固定的全球核心物种名单（`core_species/data/` 中 143 个稀有、35 个丰富），
保留每张网络的局部丰度分类；preserved 边要求边仍存在且两端丰度类别不变。
`--step` 支持 `analysis`、`figures`、`manuscript` 或 `all`。置换次数、随机种子
和 consensus hub 阈值在 `core_species/config/core_species.yaml` 中配置。
定义与完整输出列表见 [`core_species/README.md`](core_species/README.md)。

## 复现已发表运行

复现配方和代码快照位于 `reproduction/`。2026-11-14 的完整运行见
[`reproduction/2026-11-14-full-run/`](reproduction/2026-11-14-full-run/)；
其中的 `reproduce.sh` 会在启动 worker 前保存实际执行的代码和配置。

## License

MIT License，见 [LICENSE](LICENSE)。
