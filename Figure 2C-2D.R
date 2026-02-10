library(ggplot2)
library(dplyr)

save_path <- "C:/Users/sbruno3/OneDrive - UW-Madison/Desktop/Spindle paper"
d_path <- "C:/Users/sbruno3/OneDrive - UW-Madison/Desktop/Spindle paper"
d <- read.table(paste(d_path, "Efields.csv", sep = "/"), sep = ",", dec = ".", header = T)

d$cond <- factor(d$cond, levels = c("peak", "10"), labels = c( "TES¹⁵ᵏᴴᶻ-TIᴾᵉᵃᵏ",  "TES¹⁵ᵏᴴᶻ-TI¹⁰ᴴᶻ"))

# Group by subject and stimulation protocol, then average columns with data on EF intensity and focality
d <- d %>%
  group_by(Sub, cond) %>%
  summarise(across(9:23, mean, na.rm = TRUE), .groups = "drop")

# perform test and plot average EF intensity at target
wilcox.test(d$M_Mean ~ d$cond)

tiff(paste(save_path, "Efield_mean.tiff", sep = "/"), width = 5000, height = 4000, res = 1000)
ggplot(d, aes(x = cond, y = M_Mean)) +
  geom_violin(fill = "grey", alpha = 0.3) +
  geom_jitter(width = 0.1) +
  geom_segment(aes(x = 1, xend = 2, y = 2.2, yend = 2.2)) +
  annotate("text", x = 1.5, y = 2.35, label = "ns", size = 5) +
  xlab("") +
  ylab("Mean Electric Field at target (V/m)") +
  theme(
    axis.title = element_text(size = 22, face = "bold"), # Axis labels
    axis.text = element_text(size = 20, face = "bold") # Tick labels
  ) +
  theme_classic()
dev.off()

# perform test and plot average EF focality
wilcox.test(d$M_Foc ~ d$cond)

tiff(paste(save_path, "Foc_mean.tiff", sep = "/"), width = 5000, height = 4000, res = 1000)
ggplot(d, aes(x = cond, y = M_Foc)) +
  geom_violin(fill = "grey", alpha = 0.3) +
  geom_jitter(width = 0.1) +
  geom_segment(aes(x = 1, xend = 2, y = 1.5, yend = 1.5)) +
  annotate("text", x = 1.5, y = 1.6, label = "ns", size = 5) +
  xlab("") +
  ylab("Focality (EF target/EF gray matter)") +
  theme(
    axis.title = element_text(size = 22, face = "bold"), # Axis labels
    axis.text = element_text(size = 20, face = "bold") # Tick labels
  ) +
  theme_classic()
dev.off()



