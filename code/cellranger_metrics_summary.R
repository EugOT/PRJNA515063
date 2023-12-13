library(tidyverse)
metrics_summary_SRR8441476 <- read_csv("/data/PRJNA515063/cellranger/SRR8441476/outs/metrics_summary.csv") %>% mutate(Run = "SRR8441476")
metrics_summary_SRR8441477 <- read_csv("/data/PRJNA515063/cellranger/SRR8441477/outs/metrics_summary.csv") %>% mutate(Run = "SRR8441477")
metrics_summary <-
    bind_rows(
        metrics_summary_SRR8441476,
        metrics_summary_SRR8441477)

metrics_summary |>
    select("Estimated Number of Cells", "Run")

write_tsv(metrics_summary, here::here("metrics_summary.tsv"))

