library("ggplot2", warn.conflicts = FALSE)
library("dplyr", warn.conflicts = FALSE)
library("readr")
library("ggpmisc")
library("ggpubr")
require(scales)
library("MASS")
library("reshape2")


df <- read.table(file = "../results/strand_bias_results.tsv", sep = '\t', header = T)

df2 <- melt(df, id.vars = c("Tag")) 

df3 <- df2[df2$variable == "Perc_plus" | df2$variable == "Perc_minus",  ]

df3 %>%
  ggplot(aes(x = Tag , y = value, fill = variable)) + 
  geom_boxplot(alpha = 0,) +
  geom_point(size = 4.2,
             pch =21,
             position = position_dodge(width=0.75)) +
  scale_fill_manual(values = c("blue", "red"),
                    labels = c("L strand","H strand")) + 
  scale_y_continuous(limits = c(0,100),
                     breaks = seq(0,100, 10)) +
  labs(y = "Percentage of total reads",
       fill = "Strand") +
  theme_bw(base_size = 35) +
  theme(axis.title.x = element_blank(), 
        panel.grid = element_blank())

ggsave(filename = "strand_bias_nanopore.png",
       width = 15,
       height = 10,
       dpi = 300)
