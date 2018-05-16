# install.packages(c("tidyverse", "data.table", "scales", "drc", "lme4"))
library(tidyverse)
library(data.table)
library(scales)
library(drc)
library(lme4)

###########################################################################################################
###########################################################################################################
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

###########################################################################################################
###########################################################################################################
## Check and Format files 
## Fix what is not by renaming or deleting phantom columns
## (N.B., ":= NULL" is how data.table deletes columns)
## Add in Temp as a new column 
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

zebra_mussel_20C[ , Temp := 20]
zebra_mussel_4C[ , Temp := 4]

zebra_mussel <- rbind(zebra_mussel_20C, zebra_mussel_4C)

###########################################################################################################
###########################################################################################################
## Examine unionid survival data
summary(unionid)
unionid[ , Tank := factor(Tank)]

unionidPlot <- ggplot(unionid[ TimeCode == 7, ], aes(x = mpCO2, y = Alive)) + geom_jitter(width = 0, height = 0.05) +
    facet_grid( . ~ Temp,
               labeller = label_bquote( 
                                       cols = .(Temp)*degree*"C"),
               scales = "free_x") +
    theme_minimal() +
    xlab(expression(mu*"atm"~CO[2])) + 
    scale_x_continuous(label = comma) +
    stat_smooth(method = 'glm', method.args = list(family = 'binomial') )

unionidPlot
ggsave("uninidPlotSurvival.pdf", unionidPlot)

## We can only run a glm on the 7hr, 20C data becaues it is the only one with a partial response 

unionid[ , pCO2 := mpCO2 / 1000.0]
summary(unionid[ , .(pCO2, mpCO2, Tank)])

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
ED(unoidDRC, respLev = c( 50, 20, 1, 99), interval = 'delta')


## look at unionid Byssus
summary(unionid)

unionid[ , TimeCode2 := factor(TimeCode)]
levels(unionid$TimeCode2) <-  c("7 day", "96 hr")
unionid[ , TimeCode2 := factor(TimeCode2, levels = c( "96 hr", "7 day"))]

unionidPlotByssus <- ggplot(unionid, aes(x = mpCO2, y = Byssus)) +
    geom_jitter(width = 0, height = 0.05) +
    facet_grid( TimeCode2 ~ Temp,
               labeller = label_bquote( 
                                       cols = .(Temp)*degree*"C"),
               scales = "free_x") +
    theme_minimal() +
    xlab(expression("mP"~CO[2])) + 
    scale_x_continuous(label = comma) +
    stat_smooth(method = 'glm', method.args = list(family = 'binomial') )

unionidPlotByssus
ggsave("uninidPlotByssus.pdf", unionidPlotByssus)

## Fit dose-response curve and estiamte LC end points 
## unoidDRC_byssus_7_20 <- drm(pByssus ~ pCO2, weights = N,
##                             data = unoidByTank[ TimeCode == 7 & Temp == 20, ],
##                             fct = LL.2(), type = 'binomial')

unoidDRC_byssus_96_20 <- drm(pByssus ~ pCO2, weights = N,
                            data = unoidByTank[ TimeCode == 96 & Temp == 20, ],
                            fct = LL.2(), type = 'binomial')
unoidDRC_byssus_96_20

plot(unoidDRC_byssus_96_20)
ED(unoidDRC_byssus_96_20, respLev = c( 50, 20, 1, 99), interval = 'delta')


unoidDRC_byssus_96_4 <- drm(pByssus ~ pCO2, weights = N,
                            data = unoidByTank[ TimeCode == 96 & Temp == 4, ],
                            fct = LL.2(), type = 'binomial')
unoidDRC_byssus_96_4

plot(unoidDRC_byssus_96_4)
ED(unoidDRC_byssus_96_4, respLev = c( 50, 20, 1, 99), interval = 'delta')

unoidDRC_byssus_7_4 <- drm(pByssus ~ pCO2, weights = N,
                            data = unoidByTank[ TimeCode == 7 & Temp == 4, ],
                            fct = LL.2(), type = 'binomial')
unoidDRC_byssus_7_4

plot(unoidDRC_byssus_7_4)
##ED(unoidDRC_byssus_7_4, respLev = c( 50, 1, 99), interval = 'delta')


###########################################################################################################
###########################################################################################################
## Analyze Zebra Mussel length and weight data
setnames(ZBM_lw, gsub(" |\\(|\\)", "", colnames(ZBM_lw)))
ZBM_lw[ , TempFactor := factor(Temp)]

## Group B is not different
summary(ZBM_lw[ grep("Control", TrtGroup),  lm(Nettissdrywtg/Lengthmm ~ TrtGroup)])

ggplot(ZBM_lw, aes(x = ExpDurh, y = Lengthmm, fill = TrtGroup)) + geom_boxplot() +
    facet_grid( Size ~ Temp, scale = "free_y", labeller = label_bquote(cols = "Temperature"~ .(Temp) *degree*"C")) +
    scale_fill_manual(values = c("blue", "red", "orange")) + 
    theme_minimal() +
    ylab("Length (mm)") +
    xlab("Exposure duration")

summary(lm(Lengthmm ~ TrtGroup* Size + Temp, data = ZBM_lw[  ExpDurh == '72h', ]))


summary(ZBM_lw[ grep("Control", TrtGroup),  lm(Nettissdrywtg/Lengthmm ~ factor(Temp) + ExpDurh)])

ggplot(ZBM_lw[ grep("Control", TrtGroup), ], aes(x = factor(Temp), y = Nettissdrywtg/Lengthmm )) + 
	geom_boxplot(notch = TRUE) + theme_minimal() + xlab("Temperature") + ylab("Condition") +
	facet_grid(~ ExpDurh)

###########################################################################################################
###########################################################################################################
## Analyze Zebra Mussel mortality data
setnames(zebra_mussel, "ExpDur(h)", "ExpDur")

zebra_mussel[ , pCO2 := mpCO2/ 1000 ]

zebra_mussel[ , pAlive := Alive / Total ]

zebra_mussel

zmAlivePlot <-
    ggplot(data = zebra_mussel[ TimeCode == 7], aes(x = mpCO2, y = pAlive, weight = Total)) +
    geom_jitter(width = 0, height = 0.05) +
    facet_grid( Size ~ Temp,
               labeller = label_bquote(
                   rows = "Size: "~.(Size),
                   cols = .(Temp)*degree*"C"),
               scales = "free_x") +
    theme_minimal() +
    xlab(expression("mP"~CO[2])) + 
    scale_x_continuous(label = comma, trans = 'log10')  +    
    stat_smooth(method = 'glm', method.args = list(family = 'binomial'))
zmAlivePlot

zebra_mussel_DRM_7_A_4 <- drm(pAlive ~ pCO2 , weights = Total,
                              data = zebra_mussel[ TimeCode == 7 &  Size == "A" & Temp == 4, ] , 
                              fct = LL.2(), type = 'binomial')
zebra_mussel_DRM_7_A_4

plot(zebra_mussel_DRM_7_A_4)

ED(zebra_mussel_DRM_7_A_4, respLev = c( 50, 1, 99), interval = 'delta')


zebra_mussel_DRM_7_A_20 <- drm(pAlive ~ pCO2 , weights = Total,
                              data = zebra_mussel[ TimeCode == 7 &  Size == "A" & Temp == 20, ] , 
                              fct = LL.2(), type = 'binomial')
zebra_mussel_DRM_7_A_20

plot(zebra_mussel_DRM_7_A_20)

ED(zebra_mussel_DRM_7_A_20, respLev = c( 50, 1, 99), interval = 'delta')



zebra_mussel_DRM_7_S_4 <- drm(pAlive ~ pCO2 , weights = Total,
                              data = zebra_mussel[ TimeCode == 7 &  Size == "S" & Temp == 4, ] , 
                              fct = LL.2(), type = 'binomial')
zebra_mussel_DRM_7_S_4

plot(zebra_mussel_DRM_7_S_4)

ED(zebra_mussel_DRM_7_S_4, respLev = c( 50, 1, 99), interval = 'delta')


zebra_mussel_DRM_7_S_20 <- drm(pAlive ~ pCO2 , weights = Total,
                              data = zebra_mussel[ TimeCode == 7 &  Size == "S" & Temp == 20, ] , 
                              fct = LL.2(), type = 'binomial')
zebra_mussel_DRM_7_S_20

plot(zebra_mussel_DRM_7_S_20)

ED(zebra_mussel_DRM_7_S_20, respLev = c( 50, 1, 99), interval = 'delta')

###########################################################################################################
###########################################################################################################
## Attached Zebra Mussel mortality data
zebra_mussel[ , pAttached := Attached / Total ]

zmA <- zebra_mussel[ TimeCode != 0,]

zmAttachedPlot <-
    ggplot(data = zmA, aes(x = mpCO2, y = pAttached, weight = Total, color = factor(TimeCode) )) +
    geom_jitter(width = 0, height = 0.05) +
    facet_grid( Size ~ Temp,
               labeller = label_bquote(
                   rows = "Size: "~.(Size),
                   cols = .(Temp)*degree*"C"),
               scales = "free_x") +
    theme_minimal() +
    xlab(expression("mP"~CO[2])) + 
    scale_x_continuous(label = comma, trans = 'log10')  +    
    stat_smooth(method = 'glm', method.args = list(family = 'binomial')) +
    scale_color_manual("Observation day", values = c("red", "blue", "seagreen"))
zmAttachedPlot


zebra_mussel_DRM_7_A_4_attached <- drm(pAttached ~ pCO2 , weights = Total,
                                      data = zebra_mussel[ TimeCode == 7 &  Size == "A" & Temp == 4, ] , 
                                      fct = LL.2(), type = 'binomial')
zebra_mussel_DRM_7_A_4_attached
plot(zebra_mussel_DRM_7_A_4_attached)
ED(zebra_mussel_DRM_7_A_4_attached, respLev = c( 50, 1, 99), interval = 'delta')



zebra_mussel_DRM_7_S_4_attached <- drm(pAttached ~ pCO2 , weights = Total,
                                      data = zebra_mussel[ TimeCode == 7 &  Size == "S" & Temp == 4, ] , 
                                      fct = LL.2(), type = 'binomial')
zebra_mussel_DRM_7_S_4_attached
plot(zebra_mussel_DRM_7_S_4_attached)
ED(zebra_mussel_DRM_7_S_4_attached, respLev = c( 50, 1, 99), interval = 'delta')



zebra_mussel_DRM_7_A_20_attached <- drm(pAttached ~ pCO2 , weights = Total,
                                      data = zebra_mussel[ TimeCode == 7 &  Size == "A" & Temp == 20, ] , 
                                      fct = LL.2(), type = 'binomial')
zebra_mussel_DRM_7_A_20_attached
plot(zebra_mussel_DRM_7_A_20_attached)
ED(zebra_mussel_DRM_7_A_20_attached, respLev = c( 50, 1, 99), interval = 'delta')



zebra_mussel_DRM_7_S_20_attached <- drm(pAttached ~ pCO2 , weights = Total,
                                      data = zebra_mussel[ TimeCode == 7 &  Size == "S" & Temp == 20, ] , 
                                      fct = LL.2(), type = 'binomial')
zebra_mussel_DRM_7_S_20_attached
plot(zebra_mussel_DRM_7_S_20_attached)
ED(zebra_mussel_DRM_7_S_20_attached, respLev = c(50, 1, 99), interval = 'delta')

## Look at 1

zebra_mussel_DRM_1_A_4_attached <- drm(pAttached ~ pCO2 , weights = Total,
                                      data = zebra_mussel[ TimeCode == 1 &  Size == "A" & Temp == 4, ] , 
                                      fct = LL.2(), type = 'binomial')
zebra_mussel_DRM_1_A_4_attached
plot(zebra_mussel_DRM_1_A_4_attached)
ED(zebra_mussel_DRM_1_A_4_attached, respLev = c( 50, 1, 99), interval = 'delta')



zebra_mussel_DRM_1_S_4_attached <- drm(pAttached ~ pCO2 , weights = Total,
                                      data = zebra_mussel[ TimeCode == 1 &  Size == "S" & Temp == 4, ] , 
                                      fct = LL.2(), type = 'binomial')
zebra_mussel_DRM_1_S_4_attached
plot(zebra_mussel_DRM_1_S_4_attached)
ED(zebra_mussel_DRM_1_S_4_attached, respLev = c( 50, 1, 99), interval = 'delta')



zebra_mussel_DRM_1_A_20_attached <- drm(pAttached ~ pCO2 , weights = Total,
                                      data = zebra_mussel[ TimeCode == 1 &  Size == "A" & Temp == 20, ] , 
                                      fct = LL.2(), type = 'binomial')
zebra_mussel_DRM_1_A_20_attached
plot(zebra_mussel_DRM_1_A_20_attached)
ED(zebra_mussel_DRM_1_A_20_attached, respLev = c( 50, 1, 99), interval = 'delta')



zebra_mussel_DRM_1_S_20_attached <- drm(pAttached ~ pCO2 , weights = Total,
                                      data = zebra_mussel[ TimeCode == 1 &  Size == "S" & Temp == 20, ] , 
                                      fct = LL.2(), type = 'binomial')
zebra_mussel_DRM_1_S_20_attached
plot(zebra_mussel_DRM_1_S_20_attached)
ED(zebra_mussel_DRM_1_S_20_attached, respLev = c(50, 1, 99), interval = 'delta')



ED(zebra_mussel_DRM_1_S_4_attached, respLev = c(50, 1, 99), interval = 'delta')
ED(zebra_mussel_DRM_7_S_4_attached, respLev = c(50, 1, 99), interval = 'delta')


### Look at time to death
file_path_LT_4C_1 = "../LT100 4C Trial 1 Data for Review/csv files/"
ZM_LT_4C_1 <- fread(paste0(file_path_LT_4C_1, "Zebra Mussel Assessment LT100 4C trial 1 for review.csv"), na.strings = ".")

file_path_LT_4C_2 = "../LT100 4C Trial 2 Data for Review/csv files/"
ZM_LT_4C_2 <- fread(paste0(file_path_LT_4C_2, "Zebra Mussel Assessment LT100 4C trial 2 for review.csv"), na.strings = ".")

file_path_LT_12C = "../LT100 12C Data for Review/csv files/"
ZM_LT_12C  <- fread(paste0(file_path_LT_12C, "Zebra Mussel Assessment LT100 12C for review.csv"), na.strings = ".")

file_path_LT_20C = "../LT100 20C Data for Review/csv files/"
ZM_LT_20C <- fread(paste0(file_path_LT_20C, "Zebra Mussel Assessment LT100 20C for review.csv"), na.strings = ".")

colnames(ZM_LT_4C_1) == colnames(ZM_LT_4C_2 )
colnames(ZM_LT_12C ) == colnames(ZM_LT_20C )


ZM_LT_4C_1[ , Temp := 4]
ZM_LT_4C_2[ , Temp := 4]
ZM_LT_4C_1[ , TempB := '4a']
ZM_LT_4C_2[ , TempB := '4b']


ZM_LT_12C[ ,  Temp := 12]
ZM_LT_20C[ ,  Temp := 20]
ZM_LT_12C[ ,  TempB := '12']
ZM_LT_20C[ ,  TempB := '20']


ZM_TL <- rbind(ZM_LT_4C_1, ZM_LT_4C_2, ZM_LT_12C, ZM_LT_20C)

setnames(ZM_TL, "Expdur (h)", "Duration")

ZM_TL[ is.na(Alive), ]
ZM_TL[ , pAlive := Alive / Total]

ZM_TL[ , Treatment := factor(ifelse(Tank <= 4, "Control", "Elevated")) ]


ggplot(ZM_TL[ TimeCode != 24, ], aes(x = Duration, y = pAlive, weight = Total, group = Tank)) + 
	geom_jitter(width = 0, height = 0.05) +
    facet_grid( Treatment ~ Temp,
               labeller = label_bquote(
                          cols = .(Temp)*degree*"C")) +
    theme_minimal() +
    xlab("Time after start (hr)") +
    ylab("Survival") +
    stat_smooth(method = 'glm', method.args = list(family = 'binomial'), se = FALSE)

ZM_TL_20 <- drm(pAlive ~ Duration , weights = Total,
                                 data = ZM_TL[ TimeCode != 24 &  Temp == 20 & Treatment == "Elevated", ] , 
                                 fct = LL.2(), type = 'binomial')

ZM_TL_20
ED(ZM_TL_20, respLev = c(50, 1, 99), interval = 'delta')


ZM_TL_4 <- drm(pAlive ~ Duration , weights = Total,
                                 data = ZM_TL[ TimeCode != 24 &  Temp == 4 & Treatment == "Elevated", ] , 
                                 fct = LL.2(), type = 'binomial')

ZM_TL_4
ED(ZM_TL_4, respLev = c(50, 1, 99), interval = 'delta')


ZM_TL_12 <- drm(pAlive ~ Duration , weights = Total,
                                 data = ZM_TL[ TimeCode != 12 &  Temp == 4 & Treatment == "Elevated", ] , 
                                 fct = LL.2(), type = 'binomial')

ZM_TL_12
ED(ZM_TL_12, respLev = c(50, 1, 99), interval = 'delta')



### save session info
sink("SessionInfo.txt")
sessionInfo()
sink()
