library(readxl)
library(tableone)
library(dplyr)

setwd("C:/Users/sbruno3/OneDrive - UW-Madison/Desktop/Spindle paper")

#### demographics data stratified by experiment condition ####
d_part <- read_excel("DemosAndNaps.xlsx", sheet = "ParticipantsR")
vars <- names(d_part)[c(-1,-ncol(d_part))]
table1a <- CreateTableOne(vars = vars, data = d_part, strata = "Condition", test = TRUE)
print(table1a)
print(table1a, nonnormal = "Age", exact = c("Gender", "Ethnicity", "Race"))

# write table 1a
write.table(print(table1a, nonnormal = vars), "descriptive_naps_stratified.csv", dec = ".", sep = ",")

#### nap characteristics stratified by experiment condition ####
d_naps <- read_excel("DemosAndNaps.xlsx", sheet = "NapsR")

# if more than 1 nap x subj x condition, average nap metrics
d_10 <- d_naps[d_naps$Condition == "10", ]
d_10 <- d_10[!duplicated(d_10$Sub), ]
d_10_metrics <- d_naps %>%
  filter(Condition == "10") %>%
  group_by(Sub) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))
d_10 <- data.frame(d_10[, 1:3], d_10_metrics[, 2:ncol(d_10_metrics)])
d_naps <- d_naps[d_naps$Condition != "10", ]
d_naps <- rbind(d_naps, d_10)

d_peak <- d_naps[d_naps$Condition == "peak", ]
d_peak <- d_peak[!duplicated(d_peak$Sub), ]
d_peak_metrics <- d_naps %>%
  filter(Condition == "peak") %>%
  group_by(Sub) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))
d_peak <- data.frame(d_peak[, 1:3], d_peak_metrics[, 2:ncol(d_peak_metrics)])
d_naps <- d_naps[d_naps$Condition != "peak", ]
d_naps <- rbind(d_naps, d_peak)


vars <- names(d_naps)[4:ncol(d_naps)]
table1b <- CreateTableOne(vars = vars, data = d_naps, strata = "Condition", test = TRUE)
print(table1b, nonnormal = vars)

# write table 1b
write.table(print(table1b, nonnormal = vars), "descriptive_naps_stratified.csv", dec = ".", sep = ",")
