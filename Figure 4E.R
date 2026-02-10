library(readxl)
library(ggplot2)
library(dplyr)

### Bandplot all conditions ####
path <- "D:/Studies/REM-REST/Paper/analysis/Stim-Pre_NREM"
table_bl <- read_excel(paste(path, "baseline", "Table_baseline.xlsx", sep = "/"))
table_bl$cond <- rep("baseline", nrow(table_bl))
names(table_bl) <- c("Band", "mean_diff", "sd_diff", "p", "cond")

table_sham <- read_excel(paste(path, "sham", "Table_sham.xlsx", sep = "/"))
table_sham$cond <- rep("sham", nrow(table_sham))
names(table_sham) <- c("Band", "mean_diff", "sd_diff", "p", "cond")

table_lower <- read_excel(paste(path, "lower", "Table_lower.xlsx", sep = "/"))
table_lower$cond <- rep("lower", nrow(table_lower))
names(table_lower) <- c("Band", "mean_diff", "sd_diff", "p", "cond")

table_peak <- read_excel(paste(path, "peak", "Table_peak.xlsx", sep = "/"))
table_peak$cond <- rep("peak", nrow(table_lower))
names(table_peak) <- c("Band", "mean_diff", "sd_diff", "p", "cond")

d <- rbind(table_bl, table_peak, table_sham, table_lower)

d$Band <- factor(d$Band, levels = c("delta", "theta", "sigma", "beta", "gamma"), labels = c("Delta", "Theta", "Sigma", "Beta", "Gamma"))
d$cond <- factor(d$cond, levels = c("baseline", "sham", "peak", "lower"),labels = c("No Stimulation", "TES¹⁵ᵏᴴᶻ", "TES¹⁵ᵏᴴᶻ-TIᴾᵉᵃᵏ", "TES¹⁵ᵏᴴᶻ-TI¹⁰ᴴᶻ"))

d <- d %>%
  mutate(
    ymin = ifelse(mean_diff >= 0, mean_diff, mean_diff - sd_diff),
    ymax = ifelse(mean_diff >= 0, mean_diff + sd_diff, mean_diff)
  )

d <- d %>%
  mutate(
    band_num = as.numeric(factor(Band)),
    dodge_pos = case_when(
      cond == "No Stimulation" ~ -0.325,
      cond == "TES¹⁵ᵏᴴᶻ" ~ -0.125,
      cond == "TES¹⁵ᵏᴴᶻ-TIᴾᵉᵃᵏ" ~ 0.105,
      cond == "TES¹⁵ᵏᴴᶻ-TI¹⁰ᴴᶻ" ~ 0.325
    ),
    x_pos = band_num + dodge_pos,
    label_y = ifelse(mean_diff >= 0, ymax + 0.05, ymin - 0.15),
    label = ifelse(p < 0.05, "*", "")
  )


# Plot
b_plot <- ggplot(d, aes(x = Band, y = mean_diff, fill = cond)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_segment(aes(x = x_pos, xend = x_pos, y = ymin, yend = ymax),
               color = "black", linewidth = 1) +
  geom_text(aes(x = x_pos, y = label_y, label = label), size = 10, color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", lwd = 1.25) +
  scale_fill_manual(values = rev(c("#BBBBBB", "#999999", "#666666", "#333333"))) +
  labs(x = "", y = bquote(Delta * " Power"["STIM-PRE"] * " (log" [10] * mu * "V" ^ 2 * ")"), title = "") +
  guides(fill = guide_legend(title = NULL)) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    legend.text = element_text(size = 14) 
  )

ggsave("C:/Users/sbruno3/OneDrive - UW-Madison/Desktop/Spindle paper/band_plot.tiff", b_plot, device = "tiff", width = 12000, height = 6000, units = "px", dpi = 1000)


