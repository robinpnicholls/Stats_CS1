library(ggplot2)
library(dplyr)
library(viridis)
library(ggpubr)
library(fitdistrplus)
library(hrbrthemes)
library(reshape2)
library(car)
library(olsrr)
library(e1071)
library(caret)

# read in dataframe for covariates and biomarkers
df = read.csv("/covariates_biomarkers.csv", 
              na.strings=c("", "NA"))

# create df for patients that have zero-week data
zero_week_df = df[(!is.na(df$zero_weeks_IL.8) & (!is.na(df$Vas.12months))),]
VAS_df = df[(!is.na(df$zero_weeks_IL.8) & !is.na(df$Vas.12months) & 
               !is.na(df$twelve_months_IL.8)),]

#RANDOM VARIABLES
#create subsets of dataset based on binary covariates
male_smoker = zero_week_df[zero_week_df$Sex..1.male..2.female. == 1 & 
                             zero_week_df$Smoker..1.yes..2.no. == 1,]
female_smoker = zero_week_df[zero_week_df$Sex..1.male..2.female. == 2 & 
                               zero_week_df$Smoker..1.yes..2.no. == 1,]
male_nonsmoker = zero_week_df[zero_week_df$Sex..1.male..2.female. == 1 & 
                                zero_week_df$Smoker..1.yes..2.no. == 2,]
female_nonsmoker = zero_week_df[zero_week_df$Sex..1.male..2.female. == 2 & 
                                  zero_week_df$Smoker..1.yes..2.no. == 2,]

####### SEX ####################################################################

# get counts of patient sex and plot barchart
sex = data.frame(Sex=c('Male','Female'),Count=table(zero_week_df$Sex..1.male..2.female.))

mean(zero_week_df[zero_week_df$Sex..1.male..2.female. == 1,]$VAS.at.inclusion, na.rm = TRUE)
mean(zero_week_df[zero_week_df$Sex..1.male..2.female. == 2,]$VAS.at.inclusion, na.rm = TRUE)

# Barplot
ggplot(sex, aes(x=Sex, y=Count.Freq)) + 
  geom_bar(stat = "identity", fill = "#273E47") + 
  xlab("Sex") + 
  ylab("Count")

ggplot(zero_week_df, aes(x = VAS.at.inclusion, fill = as.factor(Sex..1.male..2.female.))) +
  geom_density(alpha = 0.8) +
  guides(fill = guide_legend(title = "Sex")) +
  scale_fill_manual(values = c("#273E47", "#D6E5E3"), labels = c("Male", "Female"))

####### AGE ####################################################################

zero_week_df %>%
  ggplot( aes(x=Age)) +
  geom_density(fill="#273E47", color="#D6E5E3", alpha=0.8) +
  ylab("Density")


# round ages up to next closest 10 and plot distribution of VAS at inclusion
ggplot(zero_week_df, aes(x = VAS.at.inclusion, fill = as.factor(round(Age+5,-1)))) +
  geom_density(alpha = 0.8) +
  guides(fill = guide_legend(title = "Age")) +
  scale_fill_manual(values = c("#E2C290", "#9DC0BC", "#D6E5E3", "#F6E2E3", "#E6E5E3"), 
                    labels = c("20", "30", "40", "50", "60"))

high_cas = zero_week_df$Age[zero_week_df$VAS.at.inclusion >= 5]
low_cas = zero_week_df$Age[zero_week_df$VAS.at.inclusion < 5]

mean(high_cas, na.rm = TRUE)
mean(low_cas, na.rm = TRUE)

# F.test at 5% significance to decide if variances are equal (only use f.test if normally distributed)
varience = (var.test(high_cas, low_cas, alternative = "two.sided")$p.value > 0.05) & 
  (shapiro.test(zero_week_df$Age)$p.value > 0.05)
t.test(x = high_cas, y = low_cas, var.equal = varience)

cov(VAS_df$Vas.12months, VAS_df$Age)
cor.test(VAS_df$Vas.12months, VAS_df$Age, method = "pearson")


####### SMOKER? ################################################################

# get counts of patient sex and plot barchart
smoker = data.frame(Smoker=c('Yes','No'),Count=table(zero_week_df$Smoker..1.yes..2.no.))

# Barplot
ggplot(smoker, aes(x=Smoker, y=Count.Freq)) + 
  geom_bar(stat = "identity", fill = "#273E47") + 
  xlab("Sex") + 
  ylab("Count")

ggplot(zero_week_df, aes(x = VAS.at.inclusion, fill = as.factor(Smoker..1.yes..2.no.))) +
  geom_density(alpha = 0.8) +
  guides(fill = guide_legend(title = "Smoker")) +
  scale_fill_manual(values = c("#273E47", "#D6E5E3"), labels = c("Yes", "No"))

################################################################################
############################# OUTLIERS #########################################
################################################################################
outlier_removed_df = zero_week_df

outlier_removed_df = outlier_removed_df %>%
  filter(!zero_weeks_IL.8 %in% boxplot.stats(outlier_removed_df$zero_weeks_IL.8)$out)

outlier_removed_df = outlier_removed_df %>%
  filter(!zero_weeks_VEGF.A %in% boxplot.stats(outlier_removed_df$zero_weeks_VEGF.A)$out)

outlier_removed_df = outlier_removed_df %>%
  filter(!zero_weeks_OPG %in% boxplot.stats(outlier_removed_df$zero_weeks_OPG)$out)

outlier_removed_df = outlier_removed_df %>%
  filter(!zero_weeks_TGF.beta.1 %in% boxplot.stats(outlier_removed_df$zero_weeks_TGF.beta.1)$out)

outlier_removed_df = outlier_removed_df %>%
  filter(!zero_weeks_IL.6 %in% boxplot.stats(outlier_removed_df$zero_weeks_IL.6)$out)

outlier_removed_df = outlier_removed_df %>%
  filter(!zero_weeks_CXCL9 %in% boxplot.stats(outlier_removed_df$zero_weeks_CXCL9)$out)

outlier_removed_df = outlier_removed_df %>%
  filter(!zero_weeks_CXCL1 %in% boxplot.stats(outlier_removed_df$zero_weeks_CXCL1)$out)

outlier_removed_df = outlier_removed_df %>%
  filter(!zero_weeks_IL.18 %in% boxplot.stats(outlier_removed_df$zero_weeks_IL.18)$out)

outlier_removed_df = outlier_removed_df %>%
  filter(!zero_weeks_CSF.1 %in% boxplot.stats(outlier_removed_df$zero_weeks_CSF.1)$out)

################################################################################

############################ VAS correlation ###################################

#male non smoker correlation and covariance with VAS 12-months
mean(VAS_df$VAS.at.inclusion, na.rm = TRUE)
mean(VAS_df$Vas.12months, na.rm = TRUE)
cov(VAS_df$VAS.at.inclusion, VAS_df$Vas.12months)
cor.test(VAS_df$VAS.at.inclusion, VAS_df$Vas.12months, method = "pearson")

################################################################################

# BIOMARKER DISTRIBUTIONS

####### IL-8? ##################################################################

outlier_removed_df %>%
  ggplot( aes(x=zero_weeks_IL.8)) +
  geom_density(fill="#273E47", color="#D6E5E3", alpha=0.8) +
  ylab("Density")

ggplot(outlier_removed_df, aes(x = zero_weeks_IL.8, fill = as.factor(Sex..1.male..2.female.))) +
  geom_density(alpha = 0.8) +
  guides(fill = guide_legend(title = "Sex")) +
  scale_fill_manual(values = c("#273E47", "#D6E5E3"), labels = c("Male", "Female"))

ggplot(outlier_removed_df[outlier_removed_df$Smoker..1.yes..2.no. == 2,], aes(x = zero_weeks_IL.8, fill = as.factor(Sex..1.male..2.female.))) +
  geom_density(alpha = 0.8) +
  guides(fill = guide_legend(title = "Sex")) +
  scale_fill_manual(values = c("#273E47", "#D6E5E3"), labels = c("Male", "Female"))

ggplot(outlier_removed_df, aes(x = zero_weeks_IL.8, fill = as.factor(Smoker..1.yes..2.no.))) +
  geom_density(alpha = 0.8) +
  guides(fill = guide_legend(title = "Sex")) +
  scale_fill_manual(values = c("#273E47", "#D6E5E3"), labels = c("Smoker", "Non-Smoker"))

shapiro.test(outlier_removed_df$zero_weeks_IL.8)

descdist(outlier_removed_df$zero_weeks_IL.8, discrete = FALSE, boot = 1000)

high_cas = zero_week_df$zero_weeks_IL.8[zero_week_df$VAS.at.inclusion >= 5]
low_cas = zero_week_df$zero_weeks_IL.8[zero_week_df$VAS.at.inclusion < 5]
# F.test at 5% significance to decide if variances are equal (only use f.test if normally distributed)
varience = (var.test(high_cas, low_cas, alternative = "two.sided")$p.value > 0.05) & 
  (shapiro.test(zero_week_df$zero_weeks_IL.8)$p.value > 0.05)
t.test(x = high_cas, y = low_cas, var.equal = varience)

ggplot(zero_week_df, aes(x = zero_weeks_IL.8, y = Vas.12months)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "green")

# Calculating covariance, correlations between zero_weeks_IL.8 and twelve_months_IL.8
mean(VAS_df$zero_weeks_IL.8, na.rm = TRUE)
mean(VAS_df$twelve_months_IL.8, na.rm = TRUE)
cov(VAS_df$twelve_months_IL.8, VAS_df$zero_weeks_IL.8)
cor.test(VAS_df$twelve_months_IL.8, VAS_df$zero_weeks_IL.8, method = "pearson")

# fit model and plot line of coefficients over data 
plot(VAS_df$twelve_months_IL.8, VAS_df$zero_weeks_IL.8)
M <- lm(twelve_months_IL.8 ~ zero_weeks_IL.8, VAS_df)
abline(coef(M), col = "red")

cov(VAS_df$Vas.12months, VAS_df$zero_weeks_IL.8)
cor.test(VAS_df$Vas.12months, VAS_df$zero_weeks_IL.8, method = "pearson")

#male non smoker correlation and covariance with VAS 12-months
male_nonsmoker = VAS_df[VAS_df$Sex..1.male..2.female. == 1 & VAS_df$Smoker..1.yes..2.no. ==2,]
mean(male_nonsmoker$zero_weeks_IL.8, na.rm = TRUE)
cov(male_nonsmoker$Vas.12months, male_nonsmoker$zero_weeks_IL.8)
cor.test(male_nonsmoker$Vas.12months, male_nonsmoker$zero_weeks_IL.8, method = "pearson")

####### VEGF-A? ################################################################

zero_week_df %>%
  ggplot( aes(x=zero_weeks_VEGF.A)) +
  geom_density(fill="#273E47", color="#D6E5E3", alpha=0.8) +
  ylab("Density")

ggplot(outlier_removed_df[outlier_removed_df$Smoker..1.yes..2.no. == 2,], aes(x = zero_weeks_VEGF.A, fill = as.factor(Sex..1.male..2.female.))) +
  geom_density(alpha = 0.8) +
  guides(fill = guide_legend(title = "Sex")) +
  scale_fill_manual(values = c("#273E47", "#D6E5E3"), labels = c("Male", "Female"))

shapiro.test(zero_week_df$zero_weeks_VEGF.A) 

descdist(zero_week_df$zero_weeks_VEGF.A, discrete = FALSE, boot = 1000)

high_cas = zero_week_df$zero_weeks_VEGF.A[zero_week_df$VAS.at.inclusion >= 5]
low_cas = zero_week_df$zero_weeks_VEGF.A[zero_week_df$VAS.at.inclusion < 5]
# F.test at 5% significance to decide if variances are equal (only use f.test if normally distributed)
varience = (var.test(high_cas, low_cas, alternative = "two.sided")$p.value > 0.05) & 
  (shapiro.test(zero_week_df$zero_weeks_VEGF.A)$p.value > 0.05)
t.test(x = high_cas, y = low_cas, var.equal = varience)

ggplot(zero_week_df, aes(x = zero_weeks_VEGF.A, y = Vas.12months)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "green")

# Calculating covariance, correlations between zero_weeks_VEGF.A and twelve_months_VEGF.A
mean(VAS_df$zero_weeks_VEGF.A, na.rm = TRUE)
mean(VAS_df$twelve_months_VEGF.A, na.rm = TRUE)
cov(VAS_df$twelve_months_VEGF.A, VAS_df$zero_weeks_VEGF.A)
cor.test(VAS_df$twelve_months_VEGF.A, VAS_df$zero_weeks_VEGF.A, method = "pearson")

# fit model and plot line of coefficients over data 
plot(VAS_df$twelve_months_VEGF.A, VAS_df$zero_weeks_VEGF.A)
M <- lm(twelve_months_VEGF.A ~ zero_weeks_VEGF.A, VAS_df)
abline(coef(M), col = "red")

cov(VAS_df$Vas.12months, VAS_df$zero_weeks_VEGF.A)
cor.test(VAS_df$Vas.12months, VAS_df$zero_weeks_VEGF.A, method = "pearson")

#male non smoker correlation and covariance with VAS 12-months
male_nonsmoker = VAS_df[VAS_df$Sex..1.male..2.female. == 1 & VAS_df$Smoker..1.yes..2.no. ==2,]
mean(male_nonsmoker$zero_weeks_VEGF.A, na.rm = TRUE)
cov(male_nonsmoker$Vas.12months, male_nonsmoker$zero_weeks_VEGF.A)
cor.test(male_nonsmoker$Vas.12months, male_nonsmoker$zero_weeks_VEGF.A, method = "pearson")

####### OPG? ###################################################################

zero_week_df %>%
  ggplot( aes(x=zero_weeks_OPG)) +
  geom_density(fill="#273E47", color="#D6E5E3", alpha=0.8) +
  ylab("Density")

ggplot(outlier_removed_df[outlier_removed_df$Smoker..1.yes..2.no. == 2,], aes(x = zero_weeks_OPG, fill = as.factor(Sex..1.male..2.female.))) +
  geom_density(alpha = 0.8) +
  guides(fill = guide_legend(title = "Sex")) +
  scale_fill_manual(values = c("#273E47", "#D6E5E3"), labels = c("Male", "Female"))

ggplot(zero_week_df, aes(x = zero_weeks_OPG, fill = as.factor(Sex..1.male..2.female.))) +
  geom_density(alpha = 0.8) +
  guides(fill = guide_legend(title = "Sex")) +
  scale_fill_manual(values = c("#273E47", "#D6E5E3"), labels = c("Male", "Female"))

ggplot(zero_week_df, aes(x = zero_weeks_OPG, fill = as.factor(Smoker..1.yes..2.no.))) +
  geom_density(alpha = 0.8) +
  guides(fill = guide_legend(title = "Sex")) +
  scale_fill_manual(values = c("#273E47", "#D6E5E3"), labels = c("Smoker", "Non-Smoker"))

shapiro.test(zero_week_df$zero_weeks_OPG)

descdist(zero_week_df$zero_weeks_OPG, discrete = FALSE, boot = 1000)

high_cas = zero_week_df$zero_weeks_OPG[zero_week_df$VAS.at.inclusion >= 5]
low_cas = zero_week_df$zero_weeks_OPG[zero_week_df$VAS.at.inclusion < 5]
# F.test at 5% significance to decide if variances are equal (only use f.test if normally distributed)
varience = (var.test(high_cas, low_cas, alternative = "two.sided")$p.value > 0.05) & 
  (shapiro.test(zero_week_df$zero_weeks_OPG)$p.value > 0.05)
t.test(x = high_cas, y = low_cas, var.equal = varience)

ggplot(zero_week_df, aes(x = zero_weeks_OPG, y = Vas.12months)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "green")

# Calculating covariance, correlations between zero_weeks_OPG and twelve_months_OPG
mean(VAS_df$zero_weeks_OPG, na.rm = TRUE)
mean(VAS_df$twelve_months_OPG, na.rm = TRUE)
cov(VAS_df$twelve_months_OPG, VAS_df$zero_weeks_OPG)
cor.test(VAS_df$twelve_months_OPG, VAS_df$zero_weeks_OPG, method = "pearson")

# fit model and plot line of coefficients over data 
plot(VAS_df$twelve_months_OPG, VAS_df$zero_weeks_OPG)
M <- lm(twelve_months_OPG ~ zero_weeks_OPG, VAS_df)
abline(coef(M), col = "red")

cov(VAS_df$Vas.12months, VAS_df$zero_weeks_OPG)
cor.test(VAS_df$Vas.12months, VAS_df$zero_weeks_OPG, method = "pearson")

#male non smoker correlation and covariance with VAS 12-months
male_nonsmoker = VAS_df[VAS_df$Sex..1.male..2.female. == 1 & VAS_df$Smoker..1.yes..2.no. ==2,]
mean(male_nonsmoker$zero_weeks_OPG, na.rm = TRUE)
cov(male_nonsmoker$Vas.12months, male_nonsmoker$zero_weeks_OPG)
cor.test(male_nonsmoker$Vas.12months, male_nonsmoker$zero_weeks_OPG, method = "pearson")

####### TGF-beta-1? ############################################################

zero_week_df %>%
  ggplot( aes(x=zero_weeks_TGF.beta.1)) +
  geom_density(fill="#273E47", color="#D6E5E3", alpha=0.8) +
  ylab("Density")

ggplot(outlier_removed_df[outlier_removed_df$Smoker..1.yes..2.no. == 2,], aes(x = zero_weeks_TGF.beta.1, fill = as.factor(Sex..1.male..2.female.))) +
  geom_density(alpha = 0.8) +
  guides(fill = guide_legend(title = "Sex")) +
  scale_fill_manual(values = c("#273E47", "#D6E5E3"), labels = c("Male", "Female"))

ggplot(zero_week_df, aes(x = zero_weeks_TGF.beta.1, fill = as.factor(Sex..1.male..2.female.))) +
  geom_density(alpha = 0.8) +
  guides(fill = guide_legend(title = "Sex")) +
  scale_fill_manual(values = c("#273E47", "#D6E5E3"), labels = c("Male", "Female"))

ggplot(zero_week_df, aes(x = zero_weeks_TGF.beta.1, fill = as.factor(Smoker..1.yes..2.no.))) +
  geom_density(alpha = 0.8) +
  guides(fill = guide_legend(title = "Sex")) +
  scale_fill_manual(values = c("#273E47", "#D6E5E3"), labels = c("Smoker", "Non-Smoker"))

shapiro.test(zero_week_df$zero_weeks_TGF.beta.1)

descdist(zero_week_df$zero_weeks_TGF.beta.1, discrete = FALSE, boot = 1000)

high_cas = zero_week_df$zero_weeks_TGF.beta.1[zero_week_df$VAS.at.inclusion >= 5]
low_cas = zero_week_df$zero_weeks_TGF.beta.1[zero_week_df$VAS.at.inclusion < 5]
# F.test at 5% significance to decide if variances are equal (only use f.test if normally distributed)
varience = (var.test(high_cas, low_cas, alternative = "two.sided")$p.value > 0.05) & 
  (shapiro.test(zero_week_df$zero_weeks_TGF.beta.1)$p.value > 0.05)
t.test(x = high_cas, y = low_cas, var.equal = varience)

ggplot(zero_week_df, aes(x = zero_weeks_TGF.beta.1, y = Vas.12months)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "green")

# Calculating covariance, correlations between zero_weeks_TGF.beta.1 and twelve_months_TGF.beta.1
mean(VAS_df$zero_weeks_TGF.beta.1, na.rm = TRUE)
mean(VAS_df$twelve_months_TGF.beta.1, na.rm = TRUE)
cov(VAS_df$twelve_months_TGF.beta.1, VAS_df$zero_weeks_TGF.beta.1)
cor.test(VAS_df$twelve_months_TGF.beta.1, VAS_df$zero_weeks_TGF.beta.1, method = "pearson")

# fit model and plot line of coefficients over data 
plot(VAS_df$twelve_months_TGF.beta.1, VAS_df$zero_weeks_TGF.beta.1)
M <- lm(twelve_months_TGF.beta.1 ~ zero_weeks_TGF.beta.1, VAS_df)
abline(coef(M), col = "red")

cov(VAS_df$Vas.12months, VAS_df$zero_weeks_TGF.beta.1)
cor.test(VAS_df$Vas.12months, VAS_df$zero_weeks_TGF.beta.1, method = "pearson")

#male non smoker correlation and covariance with VAS 12-months
male_nonsmoker = VAS_df[VAS_df$Sex..1.male..2.female. == 1 & VAS_df$Smoker..1.yes..2.no. ==2,]
mean(male_nonsmoker$zero_weeks_TGF.beta.1, na.rm = TRUE)
cov(male_nonsmoker$Vas.12months, male_nonsmoker$zero_weeks_TGF.beta.1)
cor.test(male_nonsmoker$Vas.12months, male_nonsmoker$zero_weeks_TGF.beta.1, method = "pearson")

####### IL-6? ##################################################################

zero_week_df %>%
  ggplot( aes(x=zero_weeks_IL.6)) +
  geom_density(fill="#273E47", color="#D6E5E3", alpha=0.8) +
  ylab("Density")

ggplot(outlier_removed_df[outlier_removed_df$Smoker..1.yes..2.no. == 2,], aes(x = zero_weeks_IL.6, fill = as.factor(Sex..1.male..2.female.))) +
  geom_density(alpha = 0.8) +
  guides(fill = guide_legend(title = "Sex")) +
  scale_fill_manual(values = c("#273E47", "#D6E5E3"), labels = c("Male", "Female"))

ggplot(zero_week_df, aes(x = zero_weeks_IL.6, fill = as.factor(Sex..1.male..2.female.))) +
  geom_density(alpha = 0.8) +
  guides(fill = guide_legend(title = "Sex")) +
  scale_fill_manual(values = c("#273E47", "#D6E5E3"), labels = c("Male", "Female"))

ggplot(zero_week_df, aes(x = zero_weeks_IL.6, fill = as.factor(Smoker..1.yes..2.no.))) +
  geom_density(alpha = 0.8) +
  guides(fill = guide_legend(title = "Sex")) +
  scale_fill_manual(values = c("#273E47", "#D6E5E3"), labels = c("Smoker", "Non-Smoker"))

shapiro.test(zero_week_df$zero_weeks_IL.6)

descdist(zero_week_df$zero_weeks_IL.6, discrete = FALSE, boot = 1000)

high_cas = zero_week_df$zero_weeks_IL.6[zero_week_df$VAS.at.inclusion >= 5]
low_cas = zero_week_df$zero_weeks_IL.6[zero_week_df$VAS.at.inclusion < 5]
# F.test at 5% significance to decide if variances are equal (only use f.test if normally distributed)
varience = (var.test(high_cas, low_cas, alternative = "two.sided")$p.value > 0.05) & 
  (shapiro.test(zero_week_df$zero_weeks_IL.6)$p.value > 0.05)
t.test(x = high_cas, y = low_cas, var.equal = varience)

ggplot(zero_week_df, aes(x = zero_weeks_IL.6, y = Vas.12months)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "green")

# Calculating covariance, correlations between zero_weeks_IL.6 and twelve_months_IL.6
mean(VAS_df$zero_weeks_IL.6, na.rm = TRUE)
mean(VAS_df$twelve_months_IL.6, na.rm = TRUE)
cov(VAS_df$twelve_months_IL.6, VAS_df$zero_weeks_IL.6)
cor.test(VAS_df$twelve_months_IL.6, VAS_df$zero_weeks_IL.6, method = "pearson")

# fit model and plot line of coefficients over data 
plot(VAS_df$twelve_months_IL.6, VAS_df$zero_weeks_IL.6)
M <- lm(twelve_months_IL.6 ~ zero_weeks_IL.6, VAS_df)
abline(coef(M), col = "red")

cov(VAS_df$Vas.12months, VAS_df$zero_weeks_IL.6)
cor.test(VAS_df$Vas.12months, VAS_df$zero_weeks_IL.6, method = "pearson")

#male non smoker correlation and covariance with VAS 12-months
male_nonsmoker = VAS_df[VAS_df$Sex..1.male..2.female. == 1 & VAS_df$Smoker..1.yes..2.no. ==2,]
mean(male_nonsmoker$zero_weeks_IL.6, na.rm = TRUE)
cov(male_nonsmoker$Vas.12months, male_nonsmoker$zero_weeks_IL.6)
cor.test(male_nonsmoker$Vas.12months, male_nonsmoker$zero_weeks_IL.6, method = "pearson")

####### CXCL9? #################################################################

zero_week_df %>%
  ggplot( aes(x=zero_weeks_CXCL9)) +
  geom_density(fill="#273E47", color="#D6E5E3", alpha=0.8) +
  ylab("Density")

ggplot(zero_week_df, aes(x = zero_weeks_CXCL9, fill = as.factor(Sex..1.male..2.female.))) +
  geom_density(alpha = 0.8) +
  guides(fill = guide_legend(title = "Sex")) +
  scale_fill_manual(values = c("#273E47", "#D6E5E3"), labels = c("Male", "Female"))

ggplot(zero_week_df, aes(x = zero_weeks_CXCL9, fill = as.factor(Smoker..1.yes..2.no.))) +
  geom_density(alpha = 0.8) +
  guides(fill = guide_legend(title = "Sex")) +
  scale_fill_manual(values = c("#273E47", "#D6E5E3"), labels = c("Smoker", "Non-Smoker"))

shapiro.test(zero_week_df$zero_weeks_CXCL9)

descdist(zero_week_df$zero_weeks_CXCL9, discrete = FALSE, boot = 1000)

high_cas = zero_week_df$zero_weeks_CXCL9[zero_week_df$VAS.at.inclusion >= 5]
low_cas = zero_week_df$zero_weeks_CXCL9[zero_week_df$VAS.at.inclusion < 5]
# F.test at 5% significance to decide if variances are equal (only use f.test if normally distributed)
varience = (var.test(high_cas, low_cas, alternative = "two.sided")$p.value > 0.05) & 
  (shapiro.test(zero_week_df$zero_weeks_CXCL9)$p.value > 0.05)
t.test(x = high_cas, y = low_cas, var.equal = varience)

ggplot(zero_week_df, aes(x = zero_weeks_CXCL9, y = Vas.12months)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "green")

# Calculating covariance, correlations between zero_weeks_CXCL9 and twelve_months_CXCL9
mean(VAS_df$zero_weeks_CXCL9, na.rm = TRUE)
mean(VAS_df$twelve_months_CXCL9, na.rm = TRUE)
cov(VAS_df$twelve_months_CXCL9, VAS_df$zero_weeks_CXCL9)
cor.test(VAS_df$twelve_months_CXCL9, VAS_df$zero_weeks_CXCL9, method = "pearson")

# fit model and plot line of coefficients over data 
plot(VAS_df$twelve_months_CXCL9, VAS_df$zero_weeks_CXCL9)
M <- lm(twelve_months_CXCL9 ~ zero_weeks_CXCL9, VAS_df)
abline(coef(M), col = "red")

cov(VAS_df$Vas.12months, VAS_df$zero_weeks_CXCL9)
cor.test(VAS_df$Vas.12months, VAS_df$zero_weeks_CXCL9, method = "pearson")

#male non smoker correlation and covariance with VAS 12-months
male_nonsmoker = VAS_df[VAS_df$Sex..1.male..2.female. == 1 & VAS_df$Smoker..1.yes..2.no. ==2,]
mean(male_nonsmoker$zero_weeks_CXCL9, na.rm = TRUE)
cov(male_nonsmoker$Vas.12months, male_nonsmoker$zero_weeks_CXCL9)
cor.test(male_nonsmoker$Vas.12months, male_nonsmoker$zero_weeks_CXCL9, method = "pearson")

####### CXCL1? #################################################################

zero_week_df %>%
  ggplot( aes(x=zero_weeks_CXCL1)) +
  geom_density(fill="#273E47", color="#D6E5E3", alpha=0.8) +
  ylab("Density")

ggplot(zero_week_df, aes(x = zero_weeks_CXCL1, fill = as.factor(Sex..1.male..2.female.))) +
  geom_density(alpha = 0.8) +
  guides(fill = guide_legend(title = "Sex")) +
  scale_fill_manual(values = c("#273E47", "#D6E5E3"), labels = c("Male", "Female"))

ggplot(zero_week_df, aes(x = zero_weeks_CXCL1, fill = as.factor(Smoker..1.yes..2.no.))) +
  geom_density(alpha = 0.8) +
  guides(fill = guide_legend(title = "Sex")) +
  scale_fill_manual(values = c("#273E47", "#D6E5E3"), labels = c("Smoker", "Non-Smoker"))

shapiro.test(zero_week_df$zero_weeks_CXCL1)

descdist(zero_week_df$zero_weeks_CXCL1, discrete = FALSE, boot = 1000)

high_cas = zero_week_df$zero_weeks_CXCL1[zero_week_df$VAS.at.inclusion >= 5]
low_cas = zero_week_df$zero_weeks_CXCL1[zero_week_df$VAS.at.inclusion < 5]
# F.test at 5% significance to decide if variances are equal (only use f.test if normally distributed)
varience = (var.test(high_cas, low_cas, alternative = "two.sided")$p.value > 0.05) & 
  (shapiro.test(zero_week_df$zero_weeks_CXCL1)$p.value > 0.05)
t.test(x = high_cas, y = low_cas, var.equal = varience)

ggplot(zero_week_df, aes(x = zero_weeks_CXCL1, y = Vas.12months)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "green")

# Calculating covariance, correlations between zero_weeks_CXCL1 and twelve_months_CXCL1
mean(VAS_df$zero_weeks_CXCL1, na.rm = TRUE)
mean(VAS_df$twelve_months_CXCL1, na.rm = TRUE)
cov(VAS_df$twelve_months_CXCL1, VAS_df$zero_weeks_CXCL1)
cor.test(VAS_df$twelve_months_CXCL1, VAS_df$zero_weeks_CXCL1, method = "pearson")

# fit model and plot line of coefficients over data 
plot(VAS_df$twelve_months_CXCL1, VAS_df$zero_weeks_CXCL1)
M <- lm(twelve_months_CXCL1 ~ zero_weeks_CXCL1, VAS_df)
abline(coef(M), col = "red")

cov(VAS_df$Vas.12months, VAS_df$zero_weeks_CXCL1)
cor.test(VAS_df$Vas.12months, VAS_df$zero_weeks_CXCL1, method = "pearson")

#male non smoker correlation and covariance with VAS 12-months
male_nonsmoker = VAS_df[VAS_df$Sex..1.male..2.female. == 1 & VAS_df$Smoker..1.yes..2.no. ==2,]
mean(male_nonsmoker$zero_weeks_CXCL1, na.rm = TRUE)
cov(male_nonsmoker$Vas.12months, male_nonsmoker$zero_weeks_CXCL1)
cor.test(male_nonsmoker$Vas.12months, male_nonsmoker$zero_weeks_CXCL1, method = "pearson")

####### IL-18? #################################################################

zero_week_df %>%
  ggplot( aes(x=zero_weeks_IL.18)) +
  geom_density(fill="#273E47", color="#D6E5E3", alpha=0.8) +
  ylab("Density")

shapiro.test(zero_week_df$zero_weeks_IL.18)

descdist(zero_week_df$zero_weeks_IL.18, discrete = FALSE, boot = 1000)

high_cas = zero_week_df$zero_weeks_IL.18[zero_week_df$VAS.at.inclusion >= 5]
low_cas = zero_week_df$zero_weeks_IL.18[zero_week_df$VAS.at.inclusion < 5]
# F.test at 5% significance to decide if variances are equal (only use f.test if normally distributed)
varience = (var.test(high_cas, low_cas, alternative = "two.sided")$p.value > 0.05) & 
  (shapiro.test(zero_week_df$zero_weeks_IL.18)$p.value > 0.05)
t.test(x = high_cas, y = low_cas, var.equal = varience)

ggplot(zero_week_df, aes(x = zero_weeks_IL.18, y = Vas.12months)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "green")

# Calculating covariance, correlations between zero_weeks_IL.18 and twelve_months_IL.18
mean(VAS_df$zero_weeks_IL.18, na.rm = TRUE)
mean(VAS_df$twelve_months_IL.18, na.rm = TRUE)
cov(VAS_df$twelve_months_IL.18, VAS_df$zero_weeks_IL.18)
cor.test(VAS_df$twelve_months_IL.18, VAS_df$zero_weeks_IL.18, method = "pearson")

# fit model and plot line of coefficients over data 
plot(VAS_df$twelve_months_IL.18, VAS_df$zero_weeks_IL.18)
M <- lm(twelve_months_IL.18 ~ zero_weeks_IL.18, VAS_df)
abline(coef(M), col = "red")

cov(VAS_df$Vas.12months, VAS_df$zero_weeks_IL.18)
cor.test(VAS_df$Vas.12months, VAS_df$zero_weeks_IL.18, method = "pearson")

#male non smoker correlation and covariance with VAS 12-months
male_nonsmoker = VAS_df[VAS_df$Sex..1.male..2.female. == 1 & VAS_df$Smoker..1.yes..2.no. ==2,]
mean(male_nonsmoker$zero_weeks_IL.18, na.rm = TRUE)
cov(male_nonsmoker$Vas.12months, male_nonsmoker$zero_weeks_IL.18)
cor.test(male_nonsmoker$Vas.12months, male_nonsmoker$zero_weeks_IL.18, method = "pearson")

####### CSF-1? #################################################################

zero_week_df %>%
  ggplot( aes(x=zero_weeks_CSF.1)) +
  geom_density(fill="#273E47", color="#D6E5E3", alpha=0.8) +
  ylab("Density")

shapiro.test(zero_week_df$zero_weeks_CSF.1)

descdist(zero_week_df$zero_weeks_CSF.1, discrete = FALSE, boot = 1000)

high_cas = zero_week_df$zero_weeks_CSF.1[zero_week_df$VAS.at.inclusion >= 5]
low_cas = zero_week_df$zero_weeks_CSF.1[zero_week_df$VAS.at.inclusion < 5]
# F.test at 5% significance to decide if variances are equal (only use f.test if normally distributed)
varience = (var.test(high_cas, low_cas, alternative = "two.sided")$p.value > 0.05) & 
  (shapiro.test(zero_week_df$zero_weeks_CSF.1)$p.value > 0.05)
t.test(x = high_cas, y = low_cas, var.equal = varience, conf.level = 0.95)

ggplot(zero_week_df, aes(x = zero_weeks_CSF.1, y = Vas.12months)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "green")

working = female_smoker
high_cas = working$zero_weeks_CSF.1[working$VAS.at.inclusion >= 5]
low_cas = working$zero_weeks_CSF.1[working$VAS.at.inclusion < 5]
# F.test at 5% significance to decide if variances are equal (only use f.test if normally distributed)
varience = (var.test(high_cas, low_cas, alternative = "two.sided")$p.value > 0.05) & 
  (shapiro.test(working$zero_weeks_CSF.1)$p.value > 0.05)
t.test(x = high_cas, y = low_cas, var.equal = varience)

ggplot(zero_week_df, aes(x = zero_weeks_CSF.1, y = Vas.12months)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "green")

# Calculating covariance, correlations between zero_weeks_CSF.1 and twelve_months_CSF.1
mean(VAS_df$zero_weeks_CSF.1, na.rm = TRUE)
mean(VAS_df$twelve_months_CSF.1, na.rm = TRUE)
cov(VAS_df$twelve_months_CSF.1, VAS_df$zero_weeks_CSF.1)
cor.test(VAS_df$twelve_months_CSF.1, VAS_df$zero_weeks_CSF.1, method = "pearson")

# fit model and plot line of coefficients over data 
plot(VAS_df$twelve_months_CSF.1, VAS_df$zero_weeks_CSF.1)
M <- lm(twelve_months_CSF.1 ~ zero_weeks_CSF.1, VAS_df)
abline(coef(M), col = "red")

cov(VAS_df$Vas.12months, VAS_df$zero_weeks_CSF.1)
cor.test(VAS_df$Vas.12months, VAS_df$zero_weeks_CSF.1, method = "pearson")

#male non smoker correlation and covariance with VAS 12-months
male_nonsmoker = VAS_df[VAS_df$Sex..1.male..2.female. == 1 & VAS_df$Smoker..1.yes..2.no. ==2,]
mean(male_nonsmoker$zero_weeks_CSF.1, na.rm = TRUE)
cov(male_nonsmoker$Vas.12months, male_nonsmoker$zero_weeks_CSF.1)
cor.test(male_nonsmoker$Vas.12months, male_nonsmoker$zero_weeks_CSF.1, method = "pearson")

################################################################################
####### Regression Modelling ###################################################
################################################################################

set.seed(88) # Set Seed to reproduce in future

# Split data into train:test 80:20
sample <- sample(seq_len(nrow(zero_week_df)), size =  floor(0.80 * nrow(zero_week_df)))

train <- zero_week_df[sample, ]
test <- zero_week_df[-sample, ]

################################################################################
######## lm - All Biomarkers and covariates ####################################
################################################################################

m <- lm(Vas.12months~VAS.at.inclusion+zero_weeks_CSF.1+zero_weeks_IL.18+zero_weeks_CXCL1+zero_weeks_CXCL9+
           zero_weeks_IL.6+zero_weeks_TGF.beta.1+zero_weeks_OPG+zero_weeks_VEGF.A+
           zero_weeks_IL.8+Sex..1.male..2.female.+Smoker..1.yes..2.no., train)

#find any outliers
ols_plot_resid_lev(m)

summary(m)

x_test = test[c("VAS.at.inclusion", "zero_weeks_CSF.1", "zero_weeks_IL.18", "zero_weeks_CXCL1", "zero_weeks_CXCL9", 
                  "zero_weeks_IL.6", "zero_weeks_TGF.beta.1", "zero_weeks_OPG", "zero_weeks_VEGF.A", 
                  "zero_weeks_IL.8", "Sex..1.male..2.female.", "Smoker..1.yes..2.no.")]
predictions = predict(m,newdata=x_test,interval="prediction", level=0.9)

plot_data <- data.frame(Predicted_value = predictions[,1],   
                        Observed_value = test$Vas.12months) 
# find number of prediction that are within 1 of actual value
sum((plot_data$Predicted_value >= plot_data$Observed_value - 1) & 
      (plot_data$Predicted_value <= plot_data$Observed_value + 1))

# plot predictions against actuals
ggplot(plot_data, aes(x = Predicted_value, y = Observed_value)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "green")

################################################################################
######## lm - All Biomarkers and covariates outliers removed ###################
################################################################################

sample <- sample(seq_len(nrow(outlier_removed_df)), size =  floor(0.80 * nrow(outlier_removed_df)))

o_train <- outlier_removed_df[sample, ]
o_test <- outlier_removed_df[-sample, ]

m <- lm(Vas.12months~VAS.at.inclusion+zero_weeks_CSF.1+zero_weeks_IL.18+zero_weeks_CXCL1+zero_weeks_CXCL9+
          zero_weeks_IL.6+zero_weeks_TGF.beta.1+zero_weeks_OPG+zero_weeks_VEGF.A+
          zero_weeks_IL.8+Sex..1.male..2.female.+Smoker..1.yes..2.no., o_train)

#find any outliers
ols_plot_resid_lev(m)

summary(m)

x_test = o_test[c("VAS.at.inclusion", "zero_weeks_CSF.1", "zero_weeks_IL.18", "zero_weeks_CXCL1", "zero_weeks_CXCL9", 
                "zero_weeks_IL.6", "zero_weeks_TGF.beta.1", "zero_weeks_OPG", "zero_weeks_VEGF.A", 
                "zero_weeks_IL.8", "Sex..1.male..2.female.", "Smoker..1.yes..2.no.")]
predictions = predict(m,newdata=x_test,interval="prediction", level=0.9)

plot_data <- data.frame(Predicted_value = predictions[,1],   
                        Observed_value = o_test$Vas.12months) 
# find number of prediction that are within 1 of actual value
sum((plot_data$Predicted_value >= plot_data$Observed_value - 1) & 
      (plot_data$Predicted_value <= plot_data$Observed_value + 1))

# plot predictions against actuals
ggplot(plot_data, aes(x = Predicted_value, y = Observed_value)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "green")

################################################################################
######## lm - change in pain ####################################
################################################################################
train_VAS_difference = train
test_VAS_difference = test
train_VAS_difference$VAS_difference <- (train$VAS.at.inclusion - train$Vas.12months)
test_VAS_difference$VAS_difference <- (test$VAS.at.inclusion - test$Vas.12months)

m <- lm(VAS_difference~VAS.at.inclusion+zero_weeks_CSF.1+zero_weeks_IL.18+zero_weeks_CXCL1+zero_weeks_CXCL9+
          zero_weeks_IL.6+zero_weeks_TGF.beta.1+zero_weeks_OPG+zero_weeks_VEGF.A+
          zero_weeks_IL.8+Sex..1.male..2.female.+Smoker..1.yes..2.no., train_VAS_difference)

#find any outliers
ols_plot_resid_lev(m)

summary(m)

x_test = test_VAS_difference[c("VAS.at.inclusion","VAS_difference", "zero_weeks_CSF.1", "zero_weeks_IL.18", "zero_weeks_CXCL1", "zero_weeks_CXCL9", 
                "zero_weeks_IL.6", "zero_weeks_TGF.beta.1", "zero_weeks_OPG", "zero_weeks_VEGF.A", 
                "zero_weeks_IL.8", "Sex..1.male..2.female.", "Smoker..1.yes..2.no.")]
predictions = predict(m,newdata=x_test,interval="prediction", level=0.9)

plot_data <- data.frame(Predicted_value = predictions[,1],   
                        Observed_value = test_VAS_difference$VAS_difference) 
# find number of prediction that are within 1 of actual value
sum((plot_data$Predicted_value >= plot_data$Observed_value - 1) & 
      (plot_data$Predicted_value <= plot_data$Observed_value + 1))

# plot predictions against actuals
ggplot(plot_data, aes(x = Predicted_value, y = Observed_value)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "green")

################################################################################
######## lm - Good Biomarkers and covariates + outliers removed ################
################################################################################

# Split data into train:test 80:20
sample <- sample(seq_len(nrow(outlier_removed_df)), size =  floor(0.80 * nrow(outlier_removed_df)))

train <- outlier_removed_df[sample, ]
test <- outlier_removed_df[-sample, ]

m <- lm(Vas.12months~VAS.at.inclusion+zero_weeks_CSF.1+zero_weeks_IL.6+zero_weeks_VEGF.A+
          zero_weeks_IL.8+Sex..1.male..2.female.+Smoker..1.yes..2.no., train)

#find any outliers
ols_plot_resid_lev(m)

summary(m)

x_test = test[c("VAS.at.inclusion","zero_weeks_CSF.1","zero_weeks_IL.6", "zero_weeks_VEGF.A", 
                "zero_weeks_IL.8", "Sex..1.male..2.female.", "Smoker..1.yes..2.no.")]
predictions = predict(m,newdata=x_test,interval="prediction", level=0.9)

plot_data <- data.frame(Predicted_value = predictions[,1],   
                        Observed_value = test$Vas.12months) 
# find number of prediction that are within 1 of actual value
sum((plot_data$Predicted_value >= plot_data$Observed_value - 1) & 
      (plot_data$Predicted_value <= plot_data$Observed_value + 1))

# plot predictions against actuals
ggplot(plot_data, aes(x = Predicted_value, y = Observed_value)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "green")

################################################################################
######## lm - IL18 ####################################
################################################################################

m <- lm(twelve_months_IL.18~VAS.at.inclusion+zero_weeks_IL.6+
          zero_weeks_IL.18+Sex..1.male..2.female.+Smoker..1.yes..2.no., train)

#find any outliers
ols_plot_resid_lev(m)

summary(m)

x_test = test[c("twelve_months_IL.6", "VAS.at.inclusion","zero_weeks_IL.6",
                "zero_weeks_IL.18", "Sex..1.male..2.female.", "Smoker..1.yes..2.no.")]
predictions = predict(m,newdata=x_test,interval="prediction", level=0.9)

plot_data <- data.frame(Predicted_value = predictions[,1],   
                        Observed_value = test$twelve_months_IL.18) 
# find number of prediction that are within 1 of actual value
sum((plot_data$Predicted_value >= plot_data$Observed_value - 0.5) & 
      (plot_data$Predicted_value <= plot_data$Observed_value + 0.5))

# plot predictions against actuals
ggplot(plot_data, aes(x = Predicted_value, y = Observed_value)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "green")

################################################################################
############################## SVM #############################################
################################################################################

train_sample = train[c("Vas.12months", "VAS.at.inclusion","zero_weeks_CSF.1","zero_weeks_IL.6", "zero_weeks_VEGF.A", 
                       "zero_weeks_IL.8", "Sex..1.male..2.female.", "Smoker..1.yes..2.no.")]

test_sample = test[c("Vas.12months", "VAS.at.inclusion","zero_weeks_CSF.1","zero_weeks_IL.6", "zero_weeks_VEGF.A", 
                     "zero_weeks_IL.8", "Sex..1.male..2.female.", "Smoker..1.yes..2.no.")]

train_sample[-1] = scale(train_sample[-1]) 
test_sample[-1] = scale(test_sample[-1]) 


classifier = svm(formula = as.integer(Vas.12months) ~ ., 
                 data = train_sample, 
                 type = 'C-classification', 
                 kernel = 'linear') 

y_pred = predict(classifier, newdata = test_sample[-1], decision.value=FALSE, probabilities=FALSE)

table(test_sample[, 1], y_pred)

################################################################################
############################# Random Forrest ###################################
################################################################################

trctrl <- trainControl(method = "none")

rfregFit <- train(Vas.12months~., 
                  data = train_sample, 
                  method = "ranger",
                  trControl=trctrl,
                  importance="permutation", 
                  tuneGrid = data.frame(mtry=7,
                                        min.node.size = 5,
                                        splitrule="variance")
)
# plot Observed vs OOB predicted values from the model
plot(train_sample$Vas.12months,rfregFit$finalModel$predictions,
     pch=19,xlab="actual Vas.12months",
     ylab="predicted Vas.12months")
mtext(paste("R-squared",
            format(rfregFit$finalModel$r.squared,digits=2)))

# plot residuals
plot(train_sample$Vas.12months,(rfregFit$finalModel$predictions-train_sample$Vas.12months),
     pch=18,ylab="residuals (predicted-actual)",
     xlab="observed Vas.12months",col="green")
abline(h=0,col="red",lty=2)







