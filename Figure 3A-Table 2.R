library(lme4)
library(readxl)
library(ggplot2)

# custom function to plot p values of mixed-effect models
summary_lmer <- function(model) {
  a <- summary(model)
  a$coef <- cbind(a$coef, p = 2*pnorm(-abs(a$coef[,3])))
  a$coef
}

#### Load and prepare power data ####
pow_path <- "D:/Studies/REM-REST/Paper/analysis/Stim-Pre_NREM"

d_10 <- read_excel(paste(pow_path, "lower", "Subj_pow_lower.xlsx", sep = "/"))
d_peak <- read_excel(paste(pow_path, "peak", "Subj_pow_peak.xlsx", sep = "/"))
d_sham <- read_excel(paste(pow_path, "sham", "Subj_pow_sham.xlsx", sep = "/"))

d <- rbind(d_peak, d_sham, d_10)
names(d)[5] <- "stim_cond"
d$stim_cond <- factor(d$stim_cond, levels = c("sham", "peak", "lower"), labels = c("TES¹⁵ᵏᴴᶻ", "TES¹⁵ᵏᴴᶻ-TIᴾᵉᵃᵏ", "TES¹⁵ᵏᴴᶻ-TI¹⁰ᴴᶻ"))

### Sigma plot difference across stimulation protocols ####
tiff("C:/Users/sbruno3/OneDrive - UW-Madison/Desktop/Spindle paper/Sigma_pow.tiff", width = 5000, height = 4000, res = 1000)
ggplot(d, aes(x = stim_cond, y = sigma_stim - sigma_pre)) +
  geom_violin(fill = "grey", alpha = 0.3) +
  geom_jitter(width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  xlab("") +
  ylab(bquote(Delta * " Sigma band power"["STIM-PRE"] * " (log" [10] * mu * "V" ^ 2 * ")")) +
  annotate("text", x = 3, y = 2.2, label = "*", size = 8, color = "red") +
  ylim(-1, 2.3) +
  theme(
    axis.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 20, face = "bold")
  ) +
  theme_classic()
dev.off()

#### Fit regression model on STIM-PRE sigma power ####
m <- lmer((sigma_stim - sigma_pre) ~ stim_cond + (1|sub), data = d, REML = F)
summary_lmer(m)
hist(residuals(m))