setwd("xxxx")
library(ggplot2)
library(ggdensity)

# -----------------------------------------------------------------------------
# Author: Guo Zhihao
# Created Date:2024/12
# -----------------------------------------------------------------------------
# This script is part of the project hosted at:
# https://github.com/JeremyQuo/Ecoli_004_ONT_DRS_scripts
# -----------------------------------------------------------------------------
data <- read.csv("feature.csv")
data$Sample <- factor(data$Sample, levels = c("ss&rd_004", "ss&rd_002"))

plot <- ggplot(data, aes(y = Read_length, x = Q_value ,fill = Sample)) +
  geom_hdr()+
  scale_fill_manual(values = c('#FF8884','#F8AC8C'))+
  facet_grid( Sample ~ .)+
  theme_bw()+
  xlim(0,30)+
  ylim(0,2000)+
  labs(x='Q score',y='Read length (bps)')+
  theme(
    strip.background = element_blank(),
    legend.position = 'bottom',
    text = element_text(size = 12),
    panel.grid.minor =element_blank(),
  )

ggsave("length_Q_002.pdf", plot, width = 3, height = 4)