# jonashaslbeck@gmail.com; May 9, 2022

# ------------------------------------------------------
# -------- What is happening here? ---------------------
# ------------------------------------------------------

# Load modality results & dataset characteristics
# Compute results & plot figures


# ------------------------------------------------------
# -------- Load Packages -------------------------------
# ------------------------------------------------------

library(xtable) # plotting latex table
library(RColorBrewer)
library(scales)
library(lme4)
library(lmerTest)
library(texreg)
library(performance)

# ------------------------------------------------------
# -------- Load Results --------------------------------
# ------------------------------------------------------

# M and density estimates (of DEN and GMM methods)
l_mB <- readRDS("Files/EMS_Modality.RDS")
l_density <- readRDS("Files/EMS_Densities.RDS")
l_skew <- readRDS("Files/EMS_Skewness.RDS")

# Load time series length business
l_TSlength <- readRDS("Files/EMS_TSlength.RDS")
l_TSlength_IDs <- readRDS("Files/EMS_TSlength_IDs.RDS")

# Load list of between-person variables
l_bet <- readRDS("Files/BetweenData.RDS")


# ------------------------------------------------------
# -------- Load Data Set Characteristics ---------------
# ------------------------------------------------------

# List of dataset characteristics
CD <- readRDS("Files/EMS_DS_info.RDS")
unique_items_ord <- CD$unique_items_ord

# ------------------------------------------------------------------
# -------- Explain Skewness with Sample Characteristics ------------
# ------------------------------------------------------------------

# ----- On level of distributions -----
# this allows me to add item-level predictors

l_reg <- list()

for(v in 1:7) {

  l_reg_v <- list()
  m_skew <- l_skew[[v]]
  m_mB <- l_mB[[v]]
  counter <- 1

  m_TSlength <- l_TSlength[[v]]

  for(i in 1:CD$unique_ptp_count[v]) {

    TS_length_v_i <- m_TSlength[i, 2]

    for(j in 1:CD$unique_item_count[v]) {

      # Get item type
      if(CD$items_per_DS[[v]][j] %in% CD$items_l_type[[1]]) val <- 1 # positive emotion
      if(CD$items_per_DS[[v]][j] %in% CD$items_l_type[[2]]) val <- 0 # negative emotion

      if(!is.na(m_mB[i,j]) & m_mB[i,j]==1) {

        # Oisín added - need unique ID numbers for the model with everyone in it
        if(v ==1) id <- i
        if(v != 1){
          idstart <- sum(CD$unique_ptp_count[1:v-1])
          id <- idstart + i}

        # Save M-estimate & item type
        l_reg_v[[counter]] <- c(m_skew[i,j], val,id,j,v, TS_length_v_i)
        counter <- counter + 1

      }

    }
  }

  # Add dataset characteristic
  m_reg_v <- do.call(rbind, l_reg_v)
  m_reg_v <- cbind(m_reg_v, CD$scale[v], CD$mDay[v], CD$popclin[v],
                   CD$item_retro[v])
  l_reg[[v]] <- m_reg_v

} # end for: datasets

m_fulldata_skew <- do.call(rbind, l_reg)
dim(m_fulldata_skew)

# --- Fit logistic regression ---

m_fulldata_skew <- as.data.frame(m_fulldata_skew)
colnames(m_fulldata_skew) <- c("Skew", "valence", "par","item","study",  "TSlength",
                          "scale", "MpD", "pop", "retro")
m_fulldata_skew$Skew <- abs(m_fulldata_skew$Skew)
head(m_fulldata_skew)




# ------------------------------------------------------------------
# -------- Multilevel plots ----------------------------------------
# ------------------------------------------------------------------
# Code below: Oisin Ryan

ids <- unique(m_fulldata_skew$par)
storage <- matrix(NA, length(ids),3, dimnames = list(NULL, c("neg","pos","study")))
for(i in 1:length(ids)){
  tmp <- subset(m_fulldata_skew, m_fulldata_skew$par == ids[i])
  storage[i,] <- c(mean(tmp[tmp$valence ==0, "Skew"]),mean(tmp[tmp$valence ==1, "Skew"]),tmp[1,"study"])
}


sc <- 0.8
pdf("Figures/Figure_Result_MeanSkew.pdf", width=8*sc, height=6*sc)

paper_names = c("Rowland et al. (2020)",
                "Bringmann et al. (2016)",
                "Vrijen et al. (2018)",
                "Fisher et al. (2017)",
                "Bringmann et al. (2013)",
                "Fried et al. (2021)",
                "Wendt et al. (2020)")

cols <- brewer.pal(8, "Set1")[-6]
plot.new()
plot.window(xlim = c(0,12), ylim = c(0,12))
axis(1); axis(2);
title(xlab = "Negative Valence",
      ylab = "Positive Valence")

for(i in c(7, 1:6)){
  tmp <- storage[storage[,"study"]==i,]
  points(tmp[,"neg"], tmp[,"pos"],pch=20, col = alpha(cols[i], 0.4), cex = 2)
}

for(i in c(7, 1:6)){
  tmp <- storage[storage[,"study"]==i,]
  abline(lm(tmp[,"pos"]~tmp[,"neg"]), col = cols[i], lwd = 2, lty = 2)
  print(cor(tmp[,"neg"], tmp[,"pos"], use = "complete.obs"))
}

legend("topright", legend=paper_names, text.col=cols, cex=.8, box.col = "white", bg = "white", box.lwd = 0)
dev.off()


# ------------------------------------------------------------------
# ----------------- Multilevel analysis ----------------------------
# ------------------------------------------------------------------

# ------------ Data prep and descriptives ----------------------

# center TSlength around the average TS_length across people
# need to do this here because unbalanaced multilevel dataset
TSl <- do.call(rbind,l_TSlength)
mean_TSlength <- mean(TSl[,2],na.rm = TRUE)

m_fulldata_skew$TSlength_c <- m_fulldata_skew$TSlength - mean_TSlength

# center MpD around study average
m_m_pd <- mean(sapply(1:7, function(s){
  m_fulldata_skew$MpD[m_fulldata_skew$study ==s][1]
}))

m_fulldata_skew$MpD_c <- m_fulldata_skew$MpD - m_m_pd

# dummy code scale
m_fulldata_skew$scale_d <- ifelse( m_fulldata_skew$scale == 100,1,0)

REML <- TRUE

# ---------------- NULL model ------------------------------------ #

outml_ips <- lmer(Skew ~ 1 + (1| par) + (1|study), data= m_fulldata_skew)

summary(outml_ips)
performance::icc(outml_ips, by_group = TRUE)
performance::icc(outml_ips)

#--------------  Model 1: Add valence with random slopes ----------#


# add valence as fixed predictor
outml1a <- lmer(Skew ~ 1 + valence + (1| study) + (1 | par),
               data= m_fulldata_skew,
               REML = REML)
summary(outml1a)

outml1b <- lmer(Skew ~ 1 + valence + (1| study) + (1 + valence | par),
                data= m_fulldata_skew,
                REML = REML)
summary(outml1b)

outml1c <- lmer(Skew ~ 1 + valence + (1 + valence| study) + (1 + valence| par),
               data= m_fulldata_skew,
               REML = REML)
summary(outml1c)

anova(outml_ips, outml1a, outml1b, outml1c)


# Table for intext model comparison
table_01 <- anova(outml_ips, outml1a, outml1c)
table_01

# Result: Accept random intercepts + random slopes of valence

# -------- Model 2: Add level-2 predictor time-series length --------

# add main effect
outml2a <- lmer(Skew ~ 1 + valence +  TSlength_c +
                 (1 + valence| study) + (1 + valence| par),
               data= m_fulldata_skew,
               REML = REML)

summary(outml2a)
anova(outml2a, outml1c)

# Add interaction effect of ts length with valence
outml2b <- lmer(Skew ~ 1 + valence +  TSlength_c + valence:TSlength_c +
                  (1 + valence| study) + (1 + valence| par),
                data= m_fulldata_skew)
summary(outml2b)


table_12 <- anova(outml2a, outml2b)

# Accept model with random intercepts, slopes valence, TSlength + interaction with valence
table_12
summary(outml2b)

# -------- Model 3: Add study-level predictors --------

# First - try adding model with all main effects
outml3a <- lmer(Skew ~ 1 + valence +  TSlength_c + valence:TSlength_c +
                   MpD_c +  retro  + pop + scale_d +
                  (1 + valence | study) + (1 + valence| par),
                data= m_fulldata_skew)
# Alternatively - add all main effects and interactions
outml3b <- lmer(Skew ~ 1 + valence +  TSlength_c + valence:TSlength_c +
                  MpD_c +  retro  + pop + scale_d +
                  valence:MpD_c +  valence:retro  + valence:pop + valence:scale_d +
                  (1 + valence | study) + (1 + valence| par),
                data= m_fulldata_skew)

summary(outml3a)
summary(outml3b)
# test both models against level-2 Model
anova(outml3a, outml2b)
anova(outml3b, outml2b)

# Result: Both models are worse than Model 2 according to BIC
# But both models are better than model 2 according to AIC and chi-square difference test

#---- Exploratory step: Find the best model which includes level 3 predictors!

# Start by trying out all individual level 3 predictors + their interactions with valence
outml3_1i <- update(outml2b, . ~ . + MpD_c)
outml3_2i <- update(outml2b, . ~ . + pop)
outml3_3i <- update(outml2b, . ~ . + retro)
outml3_4i <- update(outml2b, . ~ . + scale_d)
outml3_5i <- update(outml2b, . ~ . + MpD_c + MpD_c:valence)
outml3_6i <- update(outml2b, . ~ . + retro + retro:valence)
outml3_7i <- update(outml2b, . ~ . + pop + pop:valence)
outml3_8i <- update(outml2b, . ~ . + scale_d + scale_d:valence)

# check them
table3i <- anova(outml2b, outml3_1i, outml3_2i, outml3_3i, outml3_4i, outml3_5i,
                 outml3_6i,outml3_7i, outml3_8i)

par(mfrow = c(1,2))
plot(table3i$AIC, type = "b", col = "red")
plot(table3i$BIC, type = "b", col = "blue")

# all models worse in BIC, but "best" model according to AIC is the one with retro and interaction

summary(outml3_6i)
anova(outml3_6i, outml2b)
# best model in terms of AIC has retro + interaction in it

# ---- Step 2 exploratory: Try to improve our current best model by adding another level-3 predictor

outml3_6i_m1 <- update(outml3_6i, . ~ . + MpD_c)
outml3_6i_m2 <- update(outml3_6i, . ~ . + MpD_c + MpD_c:valence)
outml3_6i_p1 <- update(outml3_6i, . ~ . +  pop )
outml3_6i_p2 <- update(outml3_6i, . ~ . +  pop + pop:valence)
outml3_6i_s1 <- update(outml3_6i, . ~ . + scale_d )
outml3_6i_s2 <- update(outml3_6i, . ~ . + scale_d + scale_d:valence)

table3i_2 <- anova(outml2b, outml3_6i, outml3_6i_m1 ,outml3_6i_m2 ,outml3_6i_p1 ,outml3_6i_p2
                   ,outml3_6i_s1 ,outml3_6i_s2)

# this is the best one, including just the main effect of scale
summary(outml3_6i_s2)

par(mfrow = c(1,2))
plot(table3i_2$AIC, type = "b", col = "red")
plot(table3i_2$BIC, type = "b", col = "blue")

anova(outml3_6i_s2, outml3_6i)

# But the overall fit is still worse than the earlier model
anova(outml3_6i, outml3b)

# ok, so this model is best but scale is not really adding anything

# # AIC goes up everywhere if we add another predictor
# outml3_6i_s1_1 <- update(outml3_6i_s1, . ~ . + MpD_c)
# outml3_6i_s1_2  <- update(outml3_6i_s1, . ~ . + MpD_c + MpD_c:valence)
# outml3_6i_s1_3  <- update(outml3_6i_s1, . ~ . +  pop )
# outml3_6i_s1_4  <- update(outml3_6i_s1, . ~ . +  pop + pop:valence)
# outml3_6i_s1_5 <-  update(outml3_6i_s1, . ~ . +  pop + pop:valence + MpD_c + MpD_c:valence)
#
#
# table31_6_b <- anova(outml3_6i_s1,
#                      outml3_6i_s1_1,
#                      outml3_6i_s1_2,
#                      outml3_6i_s1_3,
#                      outml3_6i_s1_4,
#                      outml3_6i_s1_5
# )
#
# summary(outml3_6i_s1_5)
#
# par(mfrow = c(1,2))
# plot(table31_6_b$AIC, type = "b", col = "red")
# plot(table31_6_b$BIC, type = "b", col = "blue")
#

# ------------------------------------------------------------------
# ------------ Study-specific Analysis: Data prep ------------------
# ------------------------------------------------------------------


# # Load list of between-person variables
# l_bet <- readRDS(file = "files/l_between.RDS")

# create storage
ldat <- list()

# names of between person variables
nm <- list("group","neurotic","happybias", c("dep","anx"), "neurotic","DASS","neurotic")

# Process between person data in the following way
# 1. For all group level variables, create a dummy coding where 0 is "unhealthy" and 1 is "healthy"
# 2. For all continuous variables, standardize within the datraset
      # neuroticism shared across datasets 2,5, and 7. But measured on different scales.

# create new dataset with between person variables merged with Skew
for(s in 1:7){
  bet <- l_bet[[s]][,-1]
  if(s == 4){
    m <- apply(bet, 2, mean, na.rm = T)
    sds <- apply(bet,2, sd, na.rm = T)
    for(r in 1:nrow(bet)){
      bet[r,] <- c((bet[r,1] - m[1])/sds[1],(bet[r,2] - m[2])/sds[2])
    }
    tmp <- subset(m_fulldata_skew, m_fulldata_skew$study == s)

    ids <- unique(tmp$par)
    tmp <- cbind(tmp, NA, NA)
    colnames(tmp)[c(ncol(tmp) -1, ncol(tmp))] <- nm[[s]]
    for(i in 1:length(ids)){
      tmp[tmp$par == ids[i],c(ncol(tmp) -1, ncol(tmp))] <- bet[i,]
    }
    ldat[[s]] <- tmp
  }else{

    # standardization for continuous variables
    if(s %in% c(2,5,6,7)){

      m <- mean(bet ,na.rm = TRUE)
      sds <- sd(bet, na.rm = T)
      for(r in 1:length(bet)){
        bet[r] <- ifelse(is.na(bet[r]), NA, (bet[r] - m)/sds)
      }}

    # get study specific data from the skew dataset
    tmp <- subset(m_fulldata_skew, m_fulldata_skew$study == s)

    ids <- unique(tmp$par)
    tmp <- cbind(tmp, NA)
    colnames(tmp)[ncol(tmp)] <- nm[[s]]
    for(i in 1:length(ids)){
      tmp[tmp$par == ids[i],ncol(tmp)] <- bet[i]
    }
    ldat[[s]] <- tmp
  }
}

# ------------------------------------------------------------------
# --------------- Study-specific Analysis: Models ------------------
# ------------------------------------------------------------------

# make storage for tables
tablist <- list()

# ---------- Step 1: Check all neurotic models ------

# pick out
s = 2
  d <- ldat[[s]]
  nrow(d); length(unique(d$par))
  out <- lmer(Skew ~ 1 + valence +  (1 | par) , data=d)
  outb <- lmer(Skew ~ 1 + valence +  (1| par) + neurotic + neurotic:valence, data=d)
  # singular fit with random slopes - doesn't really seem to affect parameter estimates
  # out <- lmer(Skew ~ 1 + valence +  (1 + valence| par) , data=d)
  # outb <- lmer(Skew ~ 1 + valence +  (1 + valence| par) + neurotic + neurotic:valence, data=d)
  summary(out)
  summary(outb)

# model fit
anova(out, outb)
# save for table in appendix
tablist[[s]] <- outb


# dataset 5
s = 5
d <- ldat[[s]]
# drop rows with missing betwen person characteristic
d <- na.omit(d[,c("Skew","valence","par","neurotic")])
nrow(d); length(unique(d$par))

out <- lmer(Skew ~ 1 + valence +  (1 | par) , data=d)
outb <- lmer(Skew ~ 1 + valence +  (1| par) + neurotic + neurotic:valence, data=d)
# out <- lmer(Skew ~ 1 + valence +  (1 + valence| par) , data=d)
# outb <- lmer(Skew ~ 1 + valence +  (1 + valence| par) + neurotic + neurotic:valence, data=d)
summary(out)
summary(outb)

# save model fit
anova(out, outb)
tablist[[s]] <- outb

# dataset 7
s = 7
d <- ldat[[s]]
d <- na.omit(d[,c("Skew","valence","par","neurotic")])

nrow(d); length(unique(d$par))


out <- lmer(Skew ~ 1 + valence +  (1 | par) , data=d)
outb <- lmer(Skew ~ 1 + valence +  (1| par) + neurotic + neurotic:valence, data=d)
# out <- lmer(Skew ~ 1 + valence +  (1 + valence| par) , data=d)
# outb <- lmer(Skew ~ 1 + valence +  (1 + valence| par) + neurotic + neurotic:valence, data=d)
summary(out)
summary(outb)

# save model fit
anova(out, outb)
tablist[[s]] <- outb

# ---------- Step 2: Check all remaining models ------
s = 1
d <- ldat[[s]]
nrow(d); length(unique(d$par))


out <- lmer(Skew ~ 1 + valence +  (1 | par) , data=d)
outb <- lmer(Skew ~ 1 + valence +  (1| par) + group + valence:group, data=d)
# out <- lmer(Skew ~ 1 + valence +  (1 + valence| par) , data=d)
# outb <- lmer(Skew ~ 1 + valence +  (1 + valence| par) + group + valence:group, data=d)
summary(out)
summary(outb)

# save model fit
anova(out, outb)
tablist[[s]] <- outb


s = 3
d <- ldat[[s]]
nrow(d); length(unique(d$par))

# most people have NA for happybias, only 25 in each group
d <- d[!is.na(d$happybias),]

out <- lmer(Skew ~ 1 + valence +  (1 | par) , data=d)
outb <- lmer(Skew ~ 1 + valence +  (1| par) + happybias + valence:happybias, data=d)
summary(out)
summary(outb)

# save model fit
anova(out, outb)
tablist[[s]] <- outb

# dataset 6
s = 6
d <- ldat[[s]]
nrow(d); length(unique(d$par))

d <- na.omit(d[,c("Skew","valence","par","DASS")])

out <- lmer(Skew ~ 1 + valence +  (1 | par) , data=d)
outb <- lmer(Skew ~ 1 + valence +  (1| par) + DASS + valence:DASS, data=d)
summary(out)
summary(outb)

# save model fit
anova(out, outb)
tablist[[s]] <- outb

# dataset 4

s = 4
d <- ldat[[s]]
nrow(d); length(unique(d$par))
d <- na.omit(d[,c("Skew","valence","par","dep","anx")])
out <- lmer(Skew ~ 1 + valence +  (1 | par) , data=d)
outb <- lmer(Skew ~ 1 + valence +  (1| par) + dep + valence:dep + anx + valence:anx, data=d)
# try models with only one predictor
outb1 <- lmer(Skew ~ 1 + valence +  (1| par) + dep + valence:dep, data=d)
outb2 <- lmer(Skew ~ 1 + valence +  (1| par) + anx + valence:anx, data=d)

summary(out)
summary(outb)
summary(outb1)
summary(outb2)

# save model fit
anova(out, outb)
tablist[[s]] <- outb


# ---------- Output model results seperately for neuroticism and others ------
texreg(c(tablist[[2]], tablist[[5]], tablist[[7]]))
texreg(c(tablist[[1]], tablist[[3]], tablist[[4]], tablist[[6]]))
