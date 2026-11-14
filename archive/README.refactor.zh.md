# pantax-dbg-analysis

**语言 / Language:** [中文](../README.zh.md) | [English](../README.md)

基于 `pantax-dbg` 产出的 profiling 数据，进行 74 个队列的肠道菌群共现网络
构建、下游分析与图件，以及 ERAS 时间点的网络恢复分析。代码重构自
`high_rare_net4`，原目录 `/home/work/wenhai/meta/high_rare_net4` 未做任何修改。

## 范围

建网部分覆盖从 `00_main.R` 可达的流程：

```
00_main.R  ->  run_pipeline.R
01_config.R ->  config/pipeline.yaml
02_utils.R  ->  R/01_abundance.R
03_functions_network.R -> R/02_network.R
04_functions_analysis.R -> R/03_hub_analysis.R
05_functions_plot.R -> R/04_edge_plots.R
                       R/05_centrality_plots.R
                       R/06_network_plot.R
06_process_one.R -> R/07_process_cohort.R
```

下游分析与绘图（`net_compare2_V3.R`、`rho_plot.R`、`stab_plot.R` 等）
在 `analysis/` 下单独重构，见 `analysis/README.md`。

## 目录

```
run_pipeline.R          入口；参数解析、加载包、按队列循环
run.sh                  可复现的场景配方（显式参数、日志、有效配置）
config/
  pipeline.yaml         所有运行参数的默认值（旧 01_config.R 的开关集中在此）
  eras.yaml             ERAS 专用配置（指向 eras.tsv，runs.eras: true）
  pipeline.R            同一份配置的 R 格式示例（仅 --config 显式指定时使用）
  read_config.R         配置加载器（YAML 默认，兼容 R 格式）
  cohorts.tsv           73 个非 ERAS 队列的输入表
  eras.tsv              单行 ERAS 队列输入表（含时间点数据路径）
  PRJEB20800_57_ERAS_metadata.csv  ERAS 样本分组元数据（随仓库保存）
R/
  00_config.R           队列表读取、输出路径、配置校验
  01_abundance.R        丰度分档、流行率过滤、矩阵合并
  02_network.R          建网（SPIEC-EASI/FlashWeave/SparCC）、拓扑指标、hub
  03_hub_analysis.R     hub 与差异物种重叠
  04_edge_plots.R       边数量/强度/正负/富集类型图
  05_centrality_plots.R Rich vs Rare 中心性比较
  06_network_plot.R     单组网络图与节点表
  07_process_cohort.R   单队列主流程
analysis/               下游分析与绘图（单独入口与验证，见 analysis/README.md）
REVIEW.md               重构差异与 code-humanizer 扫描结果
```

## 运行

参数全部通过命令行传入，`config/pipeline.yaml` 只是默认值；不需要复制粘贴
修改配置块，也不需要为每种模式保存一个配置副本。

```bash
export LD_LIBRARY_PATH=/home/wenhai/miniconda3/envs/r_4.3/lib
export R_LIBS=/home/wenhai/miniconda3/envs/r_4.3/lib/R/library:/home/wenhai/application/Qmatey/tools/R
Rscript run_pipeline.R --help
```

常用开关（均可与 `--config` 组合，命令行优先）：

| 参数 | 对应旧变量 |
|---|---|
| `--data-type read_count\|relative_abundance` | `is_rc` |
| `--tool spiec-easi\|flashweave\|sparcc` | `net_tool` |
| `--prevalence 0.1` / `--present 0` | `PREVALENCE_THRESHOLD` / `PRESENT_THRESHOLD` |
| `--skip-prevalence-filter TRUE` | `Notprefilter` |
| `--filter-weight TRUE` / `--filter-confidence TRUE` | `is_filter` / `is_edge_con_filter` |
| `--fill-zero-counts TRUE` / `--floor 1e-3` | `is_rc_fill` / `FilterLow` |
| `--merged TRUE` | `is_merge` |
| `--eras TRUE` | `ERAS` |
| `--network-only TRUE` / `--core-network-only TRUE` | `build_net_only` / `build_core_net_only` |
| `--spiec-easi-ra TRUE` | `spieceasi_ra` |
| `--suffix NAME\|auto` | `outdir_last` |
| `--cohorts FILE` / `--outdir DIR` / `--only A,B` | `config_file` / `all_out_dir` / 调试单队列 |
| `--hub-top-percent` / `--hub-min-degree` | hub 定义 |

任意嵌套字段可以用 `--set`，路径用 `/` 分隔：

```bash
--set network/spiec_easi/pulsar/rep.num=50
--set abundance/thresholds/high=2e-3
```

`--print-config` 以 YAML 打印合并后的完整配置和推导出的输出目录/后缀；
`--save-config FILE` 把同样的配置写成可再次传给 `--config` 的文件，然后继续运行。
保存格式由扩展名决定：`.yaml`/`.yml` 写 YAML，其他扩展名写 R 列表；
`--config` 两种格式都接受。

```bash
# 查看本次实际会用的配置（不运行）
Rscript run_pipeline.R --tool sparcc --outdir /tmp/out --print-config

# 保存有效配置并运行
Rscript run_pipeline.R --tool sparcc --outdir /tmp/out --save-config /tmp/out/config.yaml
```

### run.sh（可复现的配方）

`run.sh` 把常用组合写成带注释的场景，每条命令都显式列出参数，并把
实际命令、日志和有效配置保存到结果目录：

```bash
./run.sh spiec_rc                  # read count + SPIEC-EASI + 0.1
./run.sh eras                      # ERAS 表 + 时间顺序
./run.sh relative                  # 相对丰度 + SPIEC-EASI
./run.sh flashweave                # read count + FlashWeave
./run.sh sparcc                    # read count + SparCC
./run.sh filter_weight             # 边权重过滤
./run.sh filter_conf               # 边置信度过滤
./run.sh network_only              # 只建 Control 网络
./run.sh full                      # 73 队列 + ERAS，再跑全部分析
./run.sh analysis                  # 对已有结果树跑分析
DRY_RUN=1 ./run.sh full            # 只打印命令
```

默认输出到 `./runs/<scenario>/`，不会覆盖原项目：

```
runs/<scenario>/
  network/            每个队列一个目录（PRJ.../res_rc/...）
    configs/<scenario>.yaml   本次有效配置
  analysis/           下游分析表格与图件（full/analysis 场景）
logs/<scenario>.log   完整命令与输出
```

可用 `OUT_ROOT`、`COHORTS`、`ERAS_COHORTS`、`MERGED_TABLE`、
`ANALYSIS_RESULTS`、`ANALYSIS_NETWORK`、`ANALYSIS_SUFFIX`、
`ANALYSIS_STEP` 覆盖路径和步骤。额外的 `run_pipeline.R` 参数可以追加在
场景名后面，例如：

```bash
./run.sh spiec_rc --only PRJNA398089_598_CD
```

`--outdir` 不写时会落到旧目录 `/home/work/wenhai/meta/high_rare_net4/`，
会覆盖已有结果；`run.sh` 始终显式指定输出目录。

### ERAS 专用配置

`config/eras.yaml` 是 `pipeline.yaml` 的 ERAS 版本，只改了两处：

- `cohorts_file: config/eras.tsv`
- `runs.eras: true`（D0 重命名为 Control，时间点保持 D4/D8/D42/D180 顺序）

`config/eras.tsv` 是单行 ERAS 队列表，已经使用 `themis_abundance_matrix_0508`
下实际存在的矩阵路径；元数据 `config/PRJEB20800_57_ERAS_metadata.csv` 也已随仓库
保存，不再依赖 `high_rare_net3`。直接运行：

```bash
Rscript run_pipeline.R --config config/eras.yaml --outdir /path/to/results
```

`config/eras.yaml` 里的 `output.dir` 仍是旧结果根目录；不传 `--outdir` 会覆盖
旧目录中的 ERAS 结果。`./run.sh eras` 默认使用这份队列表，并写入
`$OUT_ROOT/eras/network/`。

## 配置要点

`config/pipeline.yaml` 保留了旧入口的全部模式，并去掉了大量互相冲突的注释块
（`config/pipeline.R` 是同一份配置的 R 格式示例，仅 `--config` 显式指定时读取）：

- `input$data_type`：`read_count` 或 `relative_abundance`
- `input$merged`：合并队列模式（旧 `is_merge`）
- `input$skip_prevalence_filter`：旧 `Notprefilter`
- `network$tool`：`spiec-easi`、`flashweave`、`sparcc`
- `network$filter_by_weight` / `filter_by_confidence`：旧边过滤开关
- `runs$network_only` / `core_network_only`：旧 `build_net_only`
- `runs$eras`：全局 ERAS 时间顺序开关
- `output$suffix`：为 `NULL` 时按旧规则自动推导（`res_rc`、`res_raclr_no_filter`、
  `res_flashweave_rc` 等）
- `abundance$thresholds` 与 `abundance$labels`：分档阈值和图例文字，改阈值时同步改标签

`input$prevalence` 同时用于流行率过滤和输出文件名，保持旧行为：即使
`skip_prevalence_filter = TRUE`，文件名仍带 `0.1`。

## 验证

用合成数据对旧代码和新代码做了逐文件对照，覆盖：

| 场景 | 比较文件数 | 结果 |
|---|---|---|
| read count + SPIEC-EASI + 10% 流行率 | 36 | 一致 |
| GML 固定网络输入 + 多组出图 | 40 | 一致 |
| 跳过流行率过滤 | 36 | 一致 |
| 相对丰度输入 | 40 | 一致 |
| 仅建 Control 网络 | 6 | 一致 |
| 合并队列 | 36 | 一致（旧分支已损坏，见 REVIEW） |
| 边权重过滤 | 4 | 一致 |
| 边置信度过滤 | 4 | 一致 |
| SparCC | 43 | 一致 |
| FlashWeave | 41 | 一致 |

比较内容：网络 RDS 的顶点/边属性、CSV/TSV、统计文本、ggplot 的 build data
与标签。PNG 因设备差异未做字节比较。详见 `REVIEW.md`。

## 复现

2026-11-14 的完整运行（73 个非 ERAS 队列 + ERAS，再跑全部分析）记录在
[`reproduction/2026-11-14-full-run/`](reproduction/2026-11-14-full-run/)。
其中的 `reproduce.sh` 会在启动 13 个 worker 和下游分析之前，把实际执行的
代码/配置快照保存到 `control/executed_files_<UTC>/`。

## License

MIT License. See [LICENSE](LICENSE).
