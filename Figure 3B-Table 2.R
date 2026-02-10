library(dplyr)
library(ggplot2)
library(lme4)

# custom function to plot p values of mixed-effect models
summary_lmer <- function(model) {
  a <- summary(model)
  a$coef <- cbind(a$coef, p = 2*pnorm(-abs(a$coef[,3])))
  a$coef
}

# load and prepare data
path <- "C:/Users/sbruno3/OneDrive - UW-Madison/Desktop/Spindle paper/YASA/int_spindle_activity"
d <- read.csv(paste(path, "Int_spin_act_summary.csv", sep = "/"))
freqs <- unique(d$cond)
d_model <- data.frame()

save_path <- paste(path, "plots", sep = "/")
if (!exists(save_path)) {
  dir.create(save_path)
}

# Compute STIM-PRE difference in ISA and average within participant
# prepare data for regression analysis
for (iFreq in 1:length(freqs)) {
  d_cond <- d[d$cond == freqs[iFreq], ]
  groups <- d_cond %>% pull(1)
  diff <- d_cond$mean_int_spin_act_STIM - d_cond$mean_int_spin_act_PRE
  diff <- data.frame(diff, groups) %>%
    group_by(groups) %>%
    summarise(across(everything(), mean))
  
  diff$cond <- rep(freqs[iFreq], nrow(diff)) 
  d_model[(nrow(d_model) + 1):(nrow(d_model) + nrow(diff)), 1:3] <- diff
}

# Plot
d_model$cond <- factor(d_model$cond, levels = c("sham", "peak", "10"), labels = c("TES¹⁵ᵏᴴᶻ", "TES¹⁵ᵏᴴᶻ-TIᴾᵉᵃᵏ", "TES¹⁵ᵏᴴᶻ-TI¹⁰ᴴᶻ"))

tiff(paste0(save_path, "/ISA.tiff"), width = 5000, height = 4000, res = 1000)
ggplot(d_model, aes(x = cond, y = diff)) +
  geom_violin(fill = "grey", alpha = 0.3) +
  geom_jitter(width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  xlab("") +
  ylab(bquote(Delta * " Integrated Spindle Activity"["STIM-PRE"] * " (" * mu * "V" * "/s)")) +
  annotate("text", x = 3, y = 60, label = "*", size = 8, color = "red") +
  #scale_x_discrete(labels = c(bquote("TES"^"15kHz" * "-TI"^"10Hz"), bquote("TES"^"15kHz" * "-TI"^"Peak"), bquote("TES"^"15kHz"))) +
  theme(
    axis.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 20, face = "bold")
  ) +
  theme_classic()
dev.off()

# Fit the model to compare STIM-PRE across stimulation protocols
m <- lmer(diff ~ cond + (1|groups), data = d_model, REML = F)
summary_lmer(m)
