# Title: Petauri Statistician Interview Task
# Author: Helen Zhang
# Date: 17/07/2025


# Install required packages and loading them
if(!require("dplyr")) install.packages("dplyr");
if(!require("survminer")) install.packages("survminer");
if(!require("ggplot2")) install.packages("ggplot2");
if(!require("survival")) install.packages("survival");
if(!require("bshazard")) install.packages("bshazard");
if(!require("haven")) install.packages("haven");
if(!require("survRM2")) install.packages("survRM2");

library(dplyr)
library(survminer)
library(ggplot2)
library(survival)
library(bshazard)
library(haven)
library(survRM2)


# ----- Question 1 -----
#### a. reading the dataset into R
survival_data <- read_xpt("adam_dat.xpt")

## Basic inspection of the dataset
class(survival_data)
# "tbl_df"     "tbl"        "data.frame"

dim(survival_data)
# 2199 rows, 19 columns

survival_data %>% group_by(TRT01PN, TRT01P) %>% summarise(n = n())
#       1 tablemab x 52 weeks                     551
#       2 vismab x 52 weeks                       555
#       3 tablemab x 12 week -> vismab 34 weeks   557
#       4 tablemab + vismab 52 weeks              536

colSums(is.na(survival_data))
# no missing data

sum(duplicated(survival_data))
# no duplicated rows


#### b. filter the dataset down to the two treatments of interest
## extract all tablemab and vismab patients
filtered_data <- survival_data %>% filter(TRT01P %in% c("vismab x 52 weeks", "tablemab x 52 weeks"))

##check if filtering missed any patients or included unwanted patients
nrow(filtered_data) 
# 1106 patients received vismab or tablemab

any(filtered_data$TRT01PN %in% c("tablemab x 12 week -> vismab 34 weeks", "tablemab + vismab 52 weeks"))
# FALSE, patients received other treatments are not included in the filtered dataset

# save the filtered dataset to the output folder
write.csv(filtered_data, "./output/FilteredDataset.csv")

#### c. make an "Event" variable
## create an Event variable and flip the coding in CNSR
filtered_data <- filtered_data %>% mutate(Event = dplyr::if_else(CNSR == 0, 1, 0))

## check number of events
table(filtered_data$CNSR, filtered_data$Event)
# 420 CNSR = 0, Event = 1
# 686 CNSR = 1, Event = 0

#### d. create a time-to-event variable in months
## create an additional column for converting follow-up time from years to months
filtered_data <- filtered_data %>% mutate(AVAL_M = 12*AVAL)

## create a time-to-event object in months
TTE <- Surv(filtered_data$AVAL_M, filtered_data$Event == 1)

# ---- Question 2 ----
### summary of characteristics for the overall cohort
## calculate the proportion of hormone receptor positive patients
overall_hormone_pos <- filtered_data %>% 
  group_by(STR01) %>%
  summarise(n = n()) %>%
  mutate(perc_hormone_pos = 100*n/sum(n)) %>%
  filter(STR01 %in% "Positive")

## calculate the proportion of prior radiotherapy patients
overall_radio_pos <- filtered_data %>% 
  group_by(STR02) %>%
  summarise(n = n()) %>%
  mutate(perc_radio_pos = 100*n/sum(n)) %>%
  filter(STR02 %in% "PRIOR USE")

## calculate the median age of the overall cohort
overall_median_age <- filtered_data %>% summarise(median_age = median(AGE))

## combined into one table
overall_summary <- overall_median_age %>% 
  mutate(perc_hormone_pos = overall_hormone_pos$perc_hormone_pos,
         perc_radio_pos = overall_radio_pos$perc_radio_pos)

## export the console output to the output folder
sink("./output/characterstics_overall_cohort.txt")
overall_summary
sink()

#### summary of characteristics by treatment
## calculate the proportion of hormone receptor positive patients
perc_hormone_pos <- filtered_data %>% 
  group_by(TRT01P, STR01) %>% 
  summarise(n = n()) %>% 
  mutate(perc_hormone_pos = 100*n/sum(n)) %>%
  filter(STR01 %in% "Positive")

## calculate the proportion of prior radiotherapy patients
perc_radiotherapy <- filtered_data %>% 
  group_by(TRT01P, STR02) %>% 
  summarise(n = n()) %>% 
  mutate(perc_radiotherapy = 100*n/sum(n)) %>%
  filter(STR02 %in% "PRIOR USE")

## calculate the median age by treatment
median_age <- filtered_data %>% group_by(TRT01P) %>% summarise(median_age = median(AGE))

treatment_summary <- median_age %>% 
  mutate(perc_hormone_pos = perc_hormone_pos$perc_hormone_pos,
         perc_radiotherapy = perc_radiotherapy$perc_radiotherapy)

## export the console output to the output folder
sink("./output/characterstics_by_treatment.txt")
treatment_summary
sink()


# ---- Question 3 ----
#### a. KM plot with three survival lines
## create a new variable group to indicate which treatment each patient received
## each unique patient ID should correspond to overall and one of the treatments of interest
km_data <- rbind(
  filtered_data %>% mutate(group = "overall"),
  filtered_data %>% mutate(group = TRT01P)
)

## create a survival obeject
km <- survfit(Surv(AVAL_M, Event) ~ group, data = km_data)

## create a KM plot based on the survival object
km_plot <- ggsurvplot(km, 
           conf.int = TRUE, 
           data = km_data, 
           break.time.by = 10,
           risk.table = TRUE,
           risk.table.y.text = FALSE,
           risk.table.title = "Number at risk",
           xlab = "Time (months)",
           legend.title = "Treatment",
           legend.labs = c("Overall", "Tablemab x 52 weeks", "Vismab x 52 weeks"),
           title = "Kaplan-Meier plot between treatment and survival") 

## save the generated plot to the output folder
ggsave("./output/KM_plot.pdf", print(km_plot))

#### b. median survival
## median survival could be directly obtained from the survival object created in part a
## export the results as a .txt file to the output folder
sink("./output/median_survival.txt")
km
sink()
# the median survival for the overall group is 33.7 months
# the median survival for the tablemab group is 33.1 months
# the median survival for the vismab group is 37.1 months


#### c. median survival for vismab patients with negative hormone receptor status
## check the number of patients who received vismab and had negative hormone receptor status
filtered_data %>% group_by(TRT01P, STR01) %>% summarise(n = n()) # 228 patients 

## filter the dataset to only include hormone receptor negative vismab patients
vismab_neg <- filtered_data %>% 
  filter(TRT01P %in% "vismab x 52 weeks") %>%
  filter(STR01L %in% "Hormone receptor negative")

## double check filtering
nrow(vismab_neg) #228 rows
any(vismab_neg$STR01N %in% 1) # FALSE, no hormone receptor positive patients

## create a survival object for the selected patientsx
km_vismab_neg <- survfit(Surv(AVAL_M, Event) ~ 1, data = vismab_neg)

## save the results as a .txt file to the output folder
sink("./output/median_survival_visNeg.txt")
km_vismab_neg
sink()
# the median survival for vismab patients with negative hormone receptor status is 26.2 months


# ---- Question 4 ----
#### Generate a smoothed hazard plot
## create the estimates of the hazard function by treatment
fit_trt1 <- bshazard(Surv(AVAL_M, Event) ~ 1, data = subset(filtered_data, TRT01PN == 1))
fit_trt2 <- bshazard(Surv(AVAL_M, Event) ~ 1, data = subset(filtered_data, TRT01PN == 2))

## combine the estimates into one dataframe 
hazard_data <- rbind(
  data.frame(time = fit_trt1$time,
             hazard = fit_trt1$hazard,
             lowerCI = fit_trt1$lower.ci,
             upperCI = fit_trt1$upper.ci,
             Treatment = "tablemab x 52 weeks"),
  data.frame(time = fit_trt2$time,
             hazard = fit_trt2$hazard,
             lowerCI = fit_trt2$lower.ci,
             upperCI = fit_trt2$upper.ci,
             Treatment = "vismab x 52 weeks")
)

## plot the hazard plot based on the estimates
ggplot(hazard_data, aes(time)) +
  geom_line(aes(y = hazard, col = Treatment)) + 
  geom_ribbon(aes(ymin = lowerCI, ymax = upperCI, fill = Treatment), alpha = 0.2) + 
  ggtitle("Smoothed hazard rate estimates (95% confidence intervals) by treatment over time") + 
  ylab("Hazard rate") + 
  xlab("Time (months)")
# the hazard curves cross over, suggesting non-proportional hazards

## save the hazard plot to the output folder
ggsave("./output/hazard_plot.png", width = 8, height = 6)

# ---- Question 5 ----
#### Fit a Cox regression model stratified by treatment arm
cox <- coxph(Surv(AVAL_M, Event) ~ AGE + STR01 + STR02 + strata(TRT01P), data = filtered_data)

## save the summary results of the model to the output folder
sink("./output/stratified_cox.txt")
summary(cox)
sink()

# test the proportional hazards assumption for the cox model
cox.zph(cox)

# ---- Question 6 ----
#### Calculate the RMST for each treatment arm at 48 months
## Add a new variable to indicate treatment by 0 (vismab) or 1 (tablemab)
rmean_data <- filtered_data %>% mutate(TRT01PB = if_else(TRT01P == "tablemab x 52 weeks", 1, 0))

## Compare RMST between tablemab and vismab patients
rmean_results <- rmst2(rmean_data$AVAL_M, rmean_data$Event, rmean_data$TRT01PB, tau = 48)

## save the results to the output folder
sink("./output/RMST_results.txt")
rmean_results
sink()
# the difference in RMST is -1.168 months
# this means that on average, patients on tablemab survive 1.168 months shorter than patients on vismab within the 48-month follow-up period
# while no statistical significance was observed (p = 0.298)

## RMST for the tablemab arm
sink("./output/tablemab_RMST_results.txt")
rmean_results$RMST.arm1
sink()

## RMST for the vismab arm
sink("./output/vismab.RMST_results.txt")
rmean_results$RMST.arm0
sink()


