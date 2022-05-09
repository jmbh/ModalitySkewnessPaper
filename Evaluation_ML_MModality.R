# jonashaslbeck@gmail.com; March 31, 2022

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
# -------- Explain Multimodality with Sample Characteristics -------
# ------------------------------------------------------------------

# ----- On level of distributions -----
# this allows me to add item-level predictors

l_reg <- list()

for(v in 1:7) {

  l_reg_v <- list()
  m_mB <- l_mB[[v]]
  v_TS_length <- rep(NA, length=nrow(m_mB))

  m_TSlength <- l_TSlength[[v]]

  counter <- 1
  for(i in 1:CD$unique_ptp_count[v]) {

    TS_length_v_i <- m_TSlength[i, 2]

    for(j in 1:CD$unique_item_count[v]) {

      # Get item type
      if(CD$items_per_DS[[v]][j] %in% CD$items_l_type[[1]]) val <- 1 # positive emotion
      if(CD$items_per_DS[[v]][j] %in% CD$items_l_type[[2]]) val <- 0 # negative emotion

      # Oisín added - need unique ID numbers for the model with everyone in it
      if(v ==1) id <- i
      if(v != 1){
        idstart <- sum(CD$unique_ptp_count[1:v-1])
        id <- idstart + i}

      # Save M-estimate & item type
      l_reg_v[[counter]] <- c(m_mB[i,j], val,id,j,v, TS_length_v_i)
      counter <- counter + 1

    }
  }

  # Add dataset characteristic
  m_reg_v <- do.call(rbind, l_reg_v)
  m_reg_v <- cbind(m_reg_v, CD$scale[v], CD$mDay[v], CD$popclin[v],
                   CD$item_retro[v])
  l_reg[[v]] <- m_reg_v

} # end for: datasets

m_fulldata <- do.call(rbind, l_reg)
dim(m_fulldata)


# --- Fit logistic regression ---

m_fulldata <- as.data.frame(m_fulldata)
colnames(m_fulldata) <- c("M","valence", "par","item","study",  "TSlength",
                               "scale", "MpD", "pop", "retro")
# make binary response
m_fulldata$M <- (m_fulldata$M > 1)*1
m_fulldata$scale <- ifelse(m_fulldata$scale == 100, 1, 0)

TSl <- do.call(rbind, l_TSlength)
mean_TSlength <- mean(TSl[,2],na.rm = TRUE)
m_fulldata$TSlength_c <- m_fulldata$TSlength - mean_TSlength

m_m_pd <- mean(sapply(1:7, function(s){
  m_fulldata$MpD[m_fulldata$study ==s][1]
}))

m_fulldata$MpD_c <- m_fulldata$MpD - m_m_pd


table(m_fulldata$M)
head(m_fulldata)

m_fulldata <- na.omit(m_fulldata)

write.csv(m_fulldata, "Files/ModalityML_Data.csv", row.names = FALSE)


# ------------------------------------------------------------------
# -------- Multilevel stuff ----------------------------------------
# ------------------------------------------------------------------

### 3-level ML model done in Julia; see file "julia_analysis.jl"


# ------------------------------------------------------------------
# ------------------ Study-specific Analysis------------------------
# ------------------------------------------------------------------

## ----- prep ------------

# storage
ldat <- list()

# names of between person variables (note: no info for dataset 4)
nm <- c("group","neurotic","happybias", "None", "neurotic","DASS","neurotic")


# neuroticism shared across datasets 2,5, and 7. But measured on different scales
# For ease of interpretation, within-person standardize all continuous variables like neuroticism

# create new dataset with between person variables merged with Skew
for(s in 1:7){
  # 4th dataset doesn't have anything
  if(s == 4){ ldat[[s]] <- matrix(NA,1,1) }else {
    bet <- l_bet[[s]][,2]
    # standardization for continuous variables
    if(s %in% c(2,5,6,7)){
      m <- mean(bet ,na.rm = TRUE)
      sds <- sd(bet, na.rm = T)
      for(r in 1:length(bet)){
        bet[r] <- ifelse(is.na(bet[r]), NA, (bet[r] - m)/sds)
      }
    }
    # get study specific data from the skew dataset
    tmp <- subset(m_fulldata, m_fulldata$study == s)

    ids <- unique(tmp$par)
    tmp <- cbind(tmp, NA)
    colnames(tmp)[ncol(tmp)] <- nm[s]
    for(i in 1:length(ids)){
      tmp[tmp$par == ids[i],ncol(tmp)] <- bet[i]
    }
    ldat[[s]] <- tmp
  }
}



# ---------- Step 1: Check all neurotic models ------

# pick out neurosis
s = 2
d <- ldat[[s]]
nrow(d); length(unique(d$par))
# singular fit with random slopes - doesn't really seem to affect parameter estimates
m0 <- glmer(M ~ 1 + valence + (1 | par) , data=d, family = 'binomial')
m1 <- glmer(M ~ 1 + valence + (1 + valence | par) , data=d, family = 'binomial')
m1a <- glmer(M ~ 1 + valence + (1 + valence | par) + neurotic, data=d, family = 'binomial')
m1b <- glmer(M ~ 1 + valence +  (1 + valence | par) + neurotic + neurotic:valence, data=d, family = 'binomial')
summary(m0)
summary(m1)
summary(m1a)
summary(m1b)

anova(m0, m1)
anova(m1, m1a)
anova(m1, m1b)



# dataset 5
s = 5
d <- ldat[[s]]
d <- na.omit(d[,c("M","valence","par","neurotic")])
nrow(d); length(unique(d$par))
m0 <- glmer(M ~ 1 + valence + (1 | par) , data=d, family = 'binomial')
m1 <- glmer(M ~ 1 + valence + (1 + valence | par) , data=d, family = 'binomial')
m1a <- glmer(M ~ 1 + valence + (1 + valence | par) + neurotic, data=d, family = 'binomial') # NON CONVERGENCE
m1b <- glmer(M ~ 1 + valence +  (1 + valence | par) + neurotic + neurotic:valence, data=d, family = 'binomial') # NON CONVERGENCE
summary(m0)
summary(m1)
summary(m1a)
summary(m1b)

anova(m0, m1)
anova(m1, m1a)
anova(m1, m1b)

# same pattern, but smaller parameter estimates, not significant!


# dataset 5
s = 7
d <- ldat[[s]]
d <- na.omit(d[,c("M","valence","par","neurotic")])

nrow(d); length(unique(d$par))
m0 <- glmer(M ~ 1 + valence + (1 | par) , data=d, family = 'binomial')
m1 <- glmer(M ~ 1 + valence + (1 + valence | par) , data=d, family = 'binomial')
m1a <- glmer(M ~ 1 + valence + (1 + valence | par) + neurotic, data=d, family = 'binomial') # NON CONVERGENCE
m1b <- glmer(M ~ 1 + valence +  (1 + valence | par) + neurotic + neurotic:valence, data=d, family = 'binomial') # NON CONVERGENCE
anova(m0, m1)
anova(m1, m1a)
anova(m1, m1b)


#------ Step 2: Check remaining models (1,3,6)
s = 1
d <- ldat[[s]]
nrow(d); length(unique(d$par))

m0 <- glmer(M ~ 1 + valence +  (1 | par) , data=d, family = 'binomial')
m1 <- glmer(M ~ 1 + valence +  (1 + valence | par) , data=d, family = 'binomial')
m1a <- glmer(M ~ 1 + valence +  (1 + valence| par) + group + valence:group, data=d, family = 'binomial')
anova(m0, m1)
anova(m1, m1a)
anova(m1, m1b)


# check happybias group 3
s = 3
d <- ldat[[s]]
nrow(d); length(unique(d$par))

# most people have NA for happybias, only 25 in each group
d <- d[!is.na(d$happybias),]

m0 <- glmer(M ~ 1 + valence + (1 | par) , data=d, family = 'binomial')
m1 <- glmer(M ~ 1 + valence + (1 | par) + happybias, data=d, family = 'binomial')
m1a <- glmer(M ~ 1 + valence + (1 | par) + happybias + valence:happybias, data=d, family = 'binomial')
anova(m0, m1)
anova(m1, m1a)
anova(m1, m1b)


s = 6
d <- ldat[[s]]
nrow(d); length(unique(d$par))

d <- na.omit(d[,c("M","valence","par","DASS")])
m0 <- glmer(M ~ 1 + valence +  (1 | par) , data=d, family = 'binomial')
m1 <- glmer(M ~ 1 + valence +  (1 + valence| par) , data=d, family = 'binomial')
m1a <- glmer(M ~ 1 + valence +  (1| par) + DASS + valence:DASS, data=d, family = 'binomial')
