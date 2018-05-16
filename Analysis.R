library(tidyverse)
library(data.table)
library(scales)
library(lme4)
library(drc)


## Load files
file_path_4C = "../4C Toxicity Data for review/csv files/"
unionid_4C <- fread(paste0(file_path_4C,
                           "Unionid assessment 4C trial for review.csv"), na.strings = ".")

ZBM_lw_4C <- fread(paste0(file_path_4C,
                          "ZBM weights and length 4 C trial for review.csv"), na.strings = ".")

zebra_mussel_4C <- fread(paste0(file_path_4C,
                                "/Zebra Mussel assessment 4 C trial for review.csv"), na.strings = ".")

file_path_20C = "../20C Toxicity Data for review/csv files/"
list.files(file_path_20C)

unionid_20C <- fread(paste0(file_path_20C,
                            "Unionid Mussel Assessment 20C for review.csv"), na.strings = ".")

ZBM_lw_20C <- fread(paste0(file_path_20C,
                          "ZBM weights and length 20 C trial for review.csv"), na.strings = ".")

zebra_mussel_20C <- fread(paste0(file_path_20C,
                                "/Zebra Mussel assessment 20C trial for review.csv"), na.strings = ".")

## Check to make sure data formatted the same
## Fix what is not by renaming or deleting phantom columns (N.B., ":= NULL" is how data.table deletes columns)
## Add in temp as a new column 
## and then merge together across trials 

colnames(unionid_20C)
colnames(unionid_4C)

unionid_4C[ , V10 := NULL]
setnames(unionid_4C, "TrtGroup", "Trtgroup")

unionid_20C[ , Temp := 20]
unionid_4C[ , Temp := 4]
unionid <- rbind(unionid_20C, unionid_4C)

colnames(ZBM_lw_4C)
colnames(ZBM_lw_20C)

setnames(ZBM_lw_20C, "Expdur(h)", "ExpDur(h)")
colnames(ZBM_lw_4C) == colnames(ZBM_lw_20C)

ZBM_lw_20C[ , Temp := 20]
ZBM_lw_4C[  , Temp := 4]

ZBM_lw = rbind(ZBM_lw_20C, ZBM_lw_4C)
ZBM_lw

colnames(zebra_mussel_20C)
colnames(zebra_mussel_4C)

zebra_mussel_4C[ , c('V17', 'V18', 'V19') := NULL]
colnames(zebra_mussel_20C) == colnames(zebra_mussel_4C)

zebra_mussel_20C[ , temp := 20]
zebra_mussel_4C[ , temp := 4]

zebra_mussel <- rbind(zebra_mussel_20C, zebra_mussel_4C)

## Examine survival unionid data
summary(unionid)
unionid[ , Tank := factor(Tank)]

unionidPlot <- ggplot(unionid, aes(x = mpCO2, y = Alive)) + geom_jitter(width = 0, height = 0.05) +
    facet_grid( TimeCode ~ Temp,
               labeller = label_bquote( rows = .(TimeCode) * " hr",
                                       cols = .(Temp)*degree*"C"),
               scales = "free_x") +
    theme_minimal() +
    xlab(expression("mP"~CO[2])) + 
    scale_x_continuous(label = comma) +
    stat_smooth(method = 'glm', method.args = list(family = 'binomial') )

unionidPlot
ggsave("uninidPlotSurvival.pdf", unionidPlot)

## We can only run a glm on the 7hr, 20C data becaues it is the only one with a partial response 

unionid[ , pCO2 := mpCO2 / 1000.0]
summary(unionid[ , .(pCO2, mpCO2, Tank)])

## Run GLMER on results
unionidGlmerOut <- glmer( Alive ~ pCO2 + (1| Tank), unionid[ TimeCode == 7 & Temp == 20, ], family = 'binomial')
summary(unionidGlmerOut)


unoidByTank <- unionid[ , .(Alive = sum(Alive, na.rm = TRUE),
                            Byssus = sum(Byssus, na.rm = TRUE),
                            pCO2 = mean(pCO2), mpCO2 = mean(mpCO2),
                            pCOSD = sd(pCO2), .N),
                       by = .(Tank, TimeCode, Temp)]

unoidByTank[ , Dead := N - Alive]
unoidByTank[ , pAlive := Alive/N]
unoidByTank[ , pByssus := Byssus / N]

unoidByTank[ TimeCode == 7 & Temp == 20, ]

## Fit dose-response curve and estiamte LC end points 
unoidDRC <- drm(pAlive ~ pCO2, weights = N, data = unoidByTank[ TimeCode == 7 & Temp == 20, ],
                fct = LL.2(), type = 'binomial')
unoidDRC

plot(unoidDRC)
ED(unoidDRC, respLev = c(0.5, 0.01, 0.99), interval = 'delta')


## look at unionid Byssus
summary(unionid)

unionidPlotByssus <- ggplot(unionid, aes(x = mpCO2, y = Byssus)) +
    geom_jitter(width = 0, height = 0.05) +
    facet_grid( TimeCode ~ Temp,
               labeller = label_bquote( rows = .(TimeCode) * " hr",
                                       cols = .(Temp)*degree*"C"),
               scales = "free_x") +
    theme_minimal() +
    xlab(expression("mP"~CO[2])) + 
    scale_x_continuous(label = comma) +
    stat_smooth(method = 'glm', method.args = list(family = 'binomial') )

unionidPlotByssus
ggsave("uninidPlotByssus.pdf", unionidPlotByssus)

## Fit dose-response curve and estiamte LC end points 
unoidDRC_byssus_7_20 <- drm(pByssus ~ pCO2, weights = N,
                            data = unoidByTank[ TimeCode == 7 & Temp == 20, ],
                            fct = LL.2(), type = 'binomial')
unoidDRC_byssus_7_20

unoidDRC_byssus_96_20 <- drm(pByssus ~ pCO2, weights = N,
                            data = unoidByTank[ TimeCode == 96 & Temp == 20, ],
                            fct = LL.2(), type = 'binomial')
unoidDRC_byssus_96_20

plot(unoidDRC_byssus_96_20)
ED(unoidDRC_byssus_96_20, respLev = c(0.5, 0.01, 0.99), interval = 'delta')


unoidDRC_byssus_96_4 <- drm(pByssus ~ pCO2, weights = N,
                            data = unoidByTank[ TimeCode == 96 & Temp == 4, ],
                            fct = LL.2(), type = 'binomial')
unoidDRC_byssus_96_4

plot(unoidDRC_byssus_96_4)
ED(unoidDRC_byssus_96_4, respLev = c(0.5, 0.01, 0.99), interval = 'delta')


unoidDRC_byssus_7_4 <- drm(pByssus ~ pCO2, weights = N,
                            data = unoidByTank[ TimeCode == 7 & Temp == 4, ],
                            fct = LL.2(), type = 'binomial')
unoidDRC_byssus_7_4

plot(unoidDRC_byssus_7_4)
ED(unoidDRC_byssus_7_4, respLev = c(0.5, 0.01, 0.99), interval = 'delta')


### save session info
sink("SessionInfo.txt")
sessionInfo()
sink()
