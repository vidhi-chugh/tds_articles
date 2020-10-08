rm(list=ls())

library(data.table)
library(dplyr)
library(Hmisc)
library(bnlearn)
library(Rgraphviz)
library(qgraph)

set.seed(2)

df = fread("framingham.csv")
describe(df)

df$education[is.na(df$education)] = median(df$education, na.rm=TRUE)
df$cigsPerDay[is.na(df$cigsPerDay)] = median(df$cigsPerDay, na.rm=TRUE)
df$BPMeds[is.na(df$BPMeds)] = median(df$BPMeds, na.rm=TRUE)
df$totChol[is.na(df$totChol)] = median(df$totChol, na.rm=TRUE)
df$BMI[is.na(df$BMI)] = median(df$BMI, na.rm=TRUE)
df$heartRate[is.na(df$heartRate)] = median(df$heartRate, na.rm=TRUE)
df$glucose[is.na(df$glucose)] = median(df$glucose, na.rm=TRUE)

table(df$TenYearCHD)

df = df[, lapply(.SD, as.numeric)]
            
str(df)
describe(df$glucose)

######################################################################
# Recode the variables -- character to numeric and convert as factors
######################################################################
df$age_5lev = cut(df$age, breaks = c(-Inf, 40, 45, 53, 60, Inf), 
                  labels = c("age < 40", "age 40-45", "age 45-53", "age 53-60", "age 60+"))
df$cigsPerDay_4lev = cut(df$cigsPerDay, breaks = c(-Inf, 0, 10, 20, Inf), 
                         labels = c("Do not Smoke", "Smoke < 10", "Smoke 10-20", "Smoke > 20"))
df$totChol_5lev = cut(df$totChol, breaks = c(-Inf, 200, 220, 250, 280, Inf), 
                      labels = c("Total Cholestrol < 200", "Total Cholestrol 200-220", "Total Cholestrol 220-250", "Total Cholestrol 250-280", "Total Cholestrol 280+"))
df$sysBP_5lev = cut(df$sysBP, breaks = c(-Inf, 115, 125, 140, 160, Inf), 
                      labels = c("Systol BP < 115", "Systol BP 115-125", "Systol BP 125-140", "Systol BP 140-160", "Systol BP 160+"))
df$diaBP_4lev = cut(df$diaBP, breaks = c(-Inf, 75, 85, 90, Inf), 
                    labels = c("Diastol BP < 75", "Diastol BP 75-85", "Diastol BP 85-90", "Diastol BP 90+"))
df$BMI_5lev = cut(df$BMI, breaks = c(-Inf, 22, 24, 27, 30, Inf), 
                  labels = c("BMI < 22", "BMI 22-24", "BMI 24-27", "BMI 27-30", "BMI 30+"))
df$heartRate_4lev = cut(df$heartRate, breaks = c(-Inf, 70, 80, 90, Inf), 
                    labels = c("Heart Rate < 70", "Heart Rate 70-80", "Heart Rate 80-90", "Heart Rate 90+"))
df$glucose_4lev = cut(df$glucose, breaks = c(-Inf, 70, 75, 80, 90, Inf), 
                        labels = c("Glucose < 70", "Glucose 70-75", "Glucose 75-80", "Glucose 80-90", "Glucose 90+"))


df$male = factor(df$male, labels = c("Female", "Male"))
df$education = factor(df$education, labels = c("Some High School", "High School", "Vocational School", "College"))
df$currentSmoker = factor(df$currentSmoker, labels = c("Non Smoker", "Smoker"))
df$BPMeds = factor(df$BPMeds, labels = c("Not on BP Medication", "On BP Medication"))
df$prevalentStroke = factor(df$prevalentStroke, labels = c("No", "Yes"))
df$prevalentHyp = factor(df$prevalentHyp, labels = c("No", "Yes"))
df$diabetes = factor(df$diabetes, labels = c("No", "Yes"))
df$TenYearCHD = factor(df$TenYearCHD, labels = c("No CHD", "CHD"))

df = df[, !c("age", "cigsPerDay", "totChol",	"sysBP", "diaBP", "BMI", "heartRate", "glucose", "currentSmoker")]

describe(df)
str(df)

##################################################
# Specify the Prognostic Model Structure for PGM
##################################################
wl = matrix(c("age_5lev",	"education",
             "age_5lev",	"prevalentStroke",
             "age_5lev",	"prevalentHyp",
             "age_5lev",	"diabetes",
             "age_5lev",	"totChol_5lev",
             "male",	"BMI_5lev",
             "cigsPerDay_4lev",	"prevalentStroke",
             "cigsPerDay_4lev",	"prevalentHyp",
             "cigsPerDay_4lev",	"heartRate_4lev",
             "diabetes",	"heartRate_4lev",
             "diabetes",	"glucose_4lev",
             "BMI_5lev",	"heartRate_4lev",
             "sysBP_5lev",	"heartRate_4lev",
             "diaBP_4lev",	"heartRate_4lev",
             "totChol_5lev",	"diabetes",
             "totChol_5lev",	"glucose_4lev",
             "glucose_4lev",	"TenYearCHD",
             "heartRate_4lev",	"TenYearCHD"), ncol = 2, byrow = TRUE, dimnames = list(NULL, c("from", "to")))

dag = hc(df, whitelist = wl)
dag

# Estimate the strength of the nodes
bootstr = boot.strength(df, R = 500, algorithm = "hc", algorithm.args = list(score = "bde", iss = 100))
bootstr[(bootstr$strength > 0.75) & (bootstr$direction >= 0.5), ]

# run a cross validation on the network to make sure the loss function is as low as possible
cvmodel = bn.cv(dag, data = df, runs = 10, method = "k-fold", folds = 10, loss = "pred", loss.args = list(target = "TenYearCHD"))
cvmodel

# Fit the model
fit = bn.fit(dag, df, method = "bayes")

#####################
# Check some queries
#####################
prop.table(table(cpdist(fit, n = 10^6, nodes = c("TenYearCHD"), evidence = (cigsPerDay_4lev == "Do not Smoke"))))
prop.table(table(cpdist(fit, n = 10^6, nodes = c("TenYearCHD"), evidence = (cigsPerDay_4lev == "Smoke > 20"))))

prop.table(table(cpdist(fit, n = 10^6, nodes = c("diabetes"), evidence = (cigsPerDay_4lev == "Do not Smoke"))))
prop.table(table(cpdist(fit, n = 10^6, nodes = c("diabetes"), evidence = (cigsPerDay_4lev == "Smoke > 20"))))

prop.table(table(cpdist(fit, n = 10^6, nodes = c("TenYearCHD"), evidence = (cigsPerDay_4lev == "Do not Smoke" & heartRate_4lev == "Heart Rate 70-80"))))
prop.table(table(cpdist(fit, n = 10^6, nodes = c("TenYearCHD"), evidence = (cigsPerDay_4lev == "Smoke > 20" & heartRate_4lev == "Heart Rate 70-80"))))
prop.table(table(cpdist(fit, n = 10^6, nodes = c("TenYearCHD"), evidence = (cigsPerDay_4lev == "Smoke > 20" & heartRate_4lev == "Heart Rate 90+"))))

prop.table(table(cpdist(fit, n = 10^6, nodes = c("TenYearCHD"), evidence = (cigsPerDay_4lev == "Do not Smoke" & male == "Female"))))
prop.table(table(cpdist(fit, n = 10^6, nodes = c("TenYearCHD"), evidence = (cigsPerDay_4lev == "Do not Smoke" & male == "Male"))))
prop.table(table(cpdist(fit, n = 10^6, nodes = c("TenYearCHD"), evidence = (cigsPerDay_4lev == "Smoke > 20" & male == "Female"))))
prop.table(table(cpdist(fit, n = 10^6, nodes = c("TenYearCHD"), evidence = (cigsPerDay_4lev == "Smoke > 20" & male == "Male"))))

prop.table(table(cpdist(fit, n = 10^6, nodes = c("TenYearCHD"), evidence = (diabetes == "No"))))
prop.table(table(cpdist(fit, n = 10^6, nodes = c("TenYearCHD"), evidence = (diabetes == "Yes"))))

prop.table(table(cpdist(fit, n = 10^6, nodes = c("TenYearCHD"), evidence = (as.numeric(totChol_5lev) == 5 & diabetes == "Yes"))))
prop.table(table(cpdist(fit, n = 10^6, nodes = c("TenYearCHD"), evidence = (as.numeric(totChol_5lev) == 1 & diabetes == "No"))))
prop.table(table(cpdist(fit, n = 10^6, nodes = c("TenYearCHD"), evidence = (cigsPerDay_4lev == "Smoke > 20" & heartRate_4lev == "Heart Rate 70-80"))))
prop.table(table(cpdist(fit, n = 10^6, nodes = c("TenYearCHD"), evidence = (cigsPerDay_4lev == "Smoke > 20" & heartRate_4lev == "Heart Rate 90+"))))

prop.table(table(cpdist(fit, n = 10^6, nodes = c("heartRate_4lev"), evidence = (prevalentStroke == "No"))))
prop.table(table(cpdist(fit, n = 10^6, nodes = c("heartRate_4lev"), evidence = (prevalentStroke == "Yes"))))

prop.table(table(cpdist(fit, n = 10^6, nodes = c("heartRate_4lev"), evidence = (prevalentHyp == "No"))))
prop.table(table(cpdist(fit, n = 10^6, nodes = c("heartRate_4lev"), evidence = (prevalentHyp == "Yes"))))

sim = cpdist(fit, n = 10^6, nodes = c("totChol_5lev", "BMI_5lev", "age_5lev", "TenYearCHD"), evidence = (as.numeric(age_5lev) > 2))
nrow(sim[(as.numeric(sim$age_5lev) > 2) & (as.numeric(sim$totChol_5lev) > 3), ]) / nrow(sim[(as.numeric(sim$age_5lev) > 2), ])
nrow(sim[(as.numeric(sim$age_5lev) > 2) & (as.numeric(sim$BMI_5lev) > 2), ]) / nrow(sim[(as.numeric(sim$age_5lev) > 2), ])
nrow(sim[(as.numeric(sim$age_5lev) > 2) & (as.numeric(sim$totChol_5lev) > 3) & (as.numeric(sim$BMI_5lev) > 2) & (sim$TenYearCHD == "CHD"), ]) / nrow(sim[(as.numeric(sim$age_5lev) > 2) & (as.numeric(sim$totChol_5lev) > 3) & (as.numeric(sim$BMI_5lev) > 3), ])

nrow(sim[(as.numeric(sim$age_5lev) > 4) & (as.numeric(sim$totChol_5lev) > 3), ]) / nrow(sim[(as.numeric(sim$age_5lev) > 4), ])
nrow(sim[(as.numeric(sim$age_5lev) > 4) & (as.numeric(sim$BMI_5lev) > 2), ]) / nrow(sim[(as.numeric(sim$age_5lev) > 4), ])
nrow(sim[(as.numeric(sim$age_5lev) > 4) & (as.numeric(sim$totChol_5lev) > 3) & (as.numeric(sim$BMI_5lev) > 2) & (sim$TenYearCHD == "CHD"), ]) / nrow(sim[(as.numeric(sim$age_5lev) > 4) & (as.numeric(sim$totChol_5lev) > 3) & (as.numeric(sim$BMI_5lev) > 3), ])

simcurr = cpdist(fit, n = 10^6, nodes = c("diaBP_4lev", "sysBP_5lev", "cigsPerDay_4lev", "TenYearCHD"), evidence = (as.numeric(diaBP_4lev) == 2 & as.numeric(sysBP_5lev) == 2))
nrow(simcurr[(as.numeric(simcurr$diaBP_4lev) == 2) & (as.numeric(simcurr$sysBP_5lev) == 2) & (as.numeric(simcurr$cigsPerDay_4lev) == 2), ]) / nrow(simcurr[(as.numeric(simcurr$diaBP_4lev) == 2 & as.numeric(simcurr$sysBP_5lev) == 2), ])

simfuture = cpdist(fit, n = 10^6, nodes = c("diaBP_4lev", "sysBP_5lev", "cigsPerDay_4lev", "TenYearCHD"), evidence = (as.numeric(diaBP_4lev) == 3 & as.numeric(sysBP_5lev) == 3))
nrow(simfuture[(as.numeric(simfuture$diaBP_4lev) == 3) & (as.numeric(simfuture$sysBP_5lev) == 3) & (as.numeric(simfuture$cigsPerDay_4lev) > 2), ]) / nrow(simfuture[(as.numeric(simfuture$diaBP_4lev) == 3 & as.numeric(simfuture$sysBP_5lev) == 3), ])

simfuture = cpdist(fit, n = 10^6, nodes = c("diaBP_4lev", "sysBP_5lev", "cigsPerDay_4lev", "TenYearCHD"), evidence = (as.numeric(diaBP_4lev) >= 3 & as.numeric(sysBP_5lev) >= 3))
nrow(simfuture[(as.numeric(simfuture$diaBP_4lev) >= 3) & (as.numeric(simfuture$sysBP_5lev) >= 3) & (as.numeric(simfuture$cigsPerDay_4lev) > 2), ]) / nrow(simfuture[(as.numeric(simfuture$diaBP_4lev) >= 3 & as.numeric(simfuture$sysBP_5lev) >= 3), ])
