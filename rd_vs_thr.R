library("ggplot2", warn.conflicts = FALSE)
library("dplyr", warn.conflicts = FALSE)
library("readr")
library("ggpmisc")
library("ggpubr")
require(scales)
library("MASS")

df_correct <- read.table(file = "../results/rd_results_cor_final.tsv", sep = '\t', header = F)
colnames(df_correct) <- c("sample", "rd")
df_correct$throughput <- "throughput"
df_correct$throughput <- c(1763967379, 1463877581, 745081558,
                            1191851646, 1748177382, 888874016, 1388724441, 947373947, 988759990,
                            554638990, 1221290459, 434993311, 1066234226, 452157430, 1053592851,
                            2534899858, 1666201957, 1585781094,
                            473399623, 742045001, 391434319, 907230105, 427184866, 915506169)
df_correct$tag <- "tag"
df_correct$tag <- c("Fragmentation","Fragmentation","Fragmentation",
                    "Fragmentation","BamHI","Fragmentation","BamHI","Fragmentation","BamHI",
                    "Fragmentation","BamHI","Fragmentation","BamHI","Fragmentation","BamHI",
                    "Fragmentation","Fragmentation","Fragmentation",
                    "Fragmentation","BamHI","Fragmentation","BamHI","Fragmentation","BamHI")


subset<- df_correct[(df_correct$throughput > 10^8.8),]
subset$rangeth <- "<10^9"
subset$rangeth[(subset$throughput >= 10^9)]<- ">=10^9~kbp"
subset$facets = factor(subset$rangeth, labels = c("'<'~10^{9}~kbp",
                                                  "'>='~10^{9}~kbp"))
subset %>% 
  ggplot(aes(y = rd)) +
  geom_boxplot(aes(x = subset$facets, fill = tag),
               alpha = 0,
               position = position_dodge(width = 1)) +
  geom_point(aes(x = subset$facets,fill = tag),
             size = 4.2,
             pch =21,
             position = position_dodge(width = 1)) +
  scale_fill_manual(values = c("dodgerblue", "red")) + 
  facet_grid(~subset$facets,
             labeller = "label_parsed",
             scales = "free") +
  labs(x = "Throughput range",
       y = "Average mitochondrial read depth",
       fill = "Library protocol") +
  scale_y_continuous(breaks = seq(0,290,20),
                     minor_breaks = F)+
  theme_bw(base_size = 35) +
  theme(panel.grid = element_blank(),
        axis.text.x.bottom = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        strip.background = element_rect(fill = F, color = "black"))

ggsave(filename = "rd_vs_through_box.png",
       width = 15,
       height = 10,
       dpi = 300)

df_correct %>% 
  ggplot(aes(x = throughput, y = rd, color = tag)) +
  geom_smooth(method = "lm", se = F, formula = y~x, size =.5) +
  stat_cor(method = "spearman", size = 9,  show.legend = FALSE) +
  geom_point(size = 5) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks() +
  labs(x = "Throughput",
       y = "Read depth chrMT",
       color = "Library preparation") +
  theme_bw(base_size = 15)

ggsave(filename = "rd_vs_through_point.png",
       width = 15,
       height = 10,
       dpi = 300)
