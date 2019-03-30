library(fastqcr)
library(tidyverse)

qc <- read_rds('qc.rds')
per_base <- lapply(qc, function(x) {
  df <- x[['per_base_sequence_quality']]
  df %>%
    select(Base, Mean) %>%
    transform(Base = strsplit(as.character(Base), '-')) %>%
    unnest(Base) %>%
    mutate(Base = as.numeric(Base))
}) %>%
  bind_rows(.id = 'run')

# find low quality runs
per_base <- per_base %>%
  group_by(run) %>%
  mutate(run_average = mean(Mean) > 30)

per_base %>%
  select(run, run_average) %>%
  unique() %>%
  with(table(run_average))

# plot average per base quality
per_base %>%
  ggplot(aes(x = Base, y = Mean, group = run, color = run_average)) +
  geom_line()

read.csv('runs.csv') %>%
  left_join(select(per_base, run, run_average) %>% unique()) %>%
  mutate(run_average = ifelse(is.na(run_average), TRUE, run_average)) %>%
  write_csv('runs_modified.csv')
