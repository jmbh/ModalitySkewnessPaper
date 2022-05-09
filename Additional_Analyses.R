# jonashaslbeck@gmail.com; May 9, 2022

# ------------------------------------------------------
# -------- What is happening here? ---------------------
# ------------------------------------------------------

# (a) This produces the analysis and figures for Appendix F
# and (b) the analysis and figure for Appendix G

# ------------------------------------------------------
# -------- Load Packages -------------------------------
# ------------------------------------------------------

library(RColorBrewer)
library(scales)

library(mgm)
library(mlVAR)

source("aux_Functions.R")

# ------------------------------------------------------
# -------- Load Results --------------------------------
# ------------------------------------------------------

# M-estimates, density, skewness, quantile info
l_mB <- readRDS("files/EMS_Modality.RDS")
l_density <- readRDS("files/EMS_Densities.RDS")
l_skew <- readRDS("files/EMS_Skewness.RDS")
l_qntl <- readRDS("files/EMS_Quantiles.RDS")


# ------------------------------------------------------
# -------- Load Data Set Characteristics ---------------
# ------------------------------------------------------

# List of dataset characteristics
CD <- readRDS("files/EMS_DS_info.RDS")
unique_items_ord <- CD$unique_items_ord


# ------------------------------------------------------
# -------- Function to plot labels ---------------------
# ------------------------------------------------------

plotLab <- function(label, rot=0, cex=1, x=0, y=0) {
  par(mar=rep(0,4))
  plot.new()
  plot.window(xlim=c(-1, 1), ylim=c(-1, 1))
  text(x, y, label, srt=rot, cex=cex)
}


# ------------------------------------------------------------------
# -------- Make Proper Figure for Appendix -------------------------
# ------------------------------------------------------------------

# Plot continuous and Likert scale data separately

l_parts <- list(1:4, 5:7)

for(p in 1:2) {

  n_panels <- length(l_parts[[p]])

  sc <- 1.3
  v_width <- c(n_panels*2*sc-1, n_panels*2*sc)

  pdf(paste0("Figures/Relation_Qntl_Skew_Agg_AppP_", p, ".pdf"), width=v_width[p], height=4*sc*0.9)

  if(p==1) {
    lmat <- matrix(7:14, 2, 4, byrow = FALSE)
    lmat <- cbind(1:2, lmat)
    lmat <- rbind(c(0, 3:6), lmat)
    lo <- layout(lmat, widths = c(.1, 1, 1, 1), heights = c(.1, 1, 1))
    # layout.show(lo)
  } else {
    lmat <- matrix(6:11, 2, 3, byrow = FALSE)
    lmat <- cbind(1:2, lmat)
    lmat <- rbind(c(0, 3:5), lmat)
    lo <- layout(lmat, widths = c(.1, 1, 1), heights = c(.1, 1, 1))
    # layout.show(lo)
  }

  # --- Plot Labels ---
  cex=1.5
  plotLab("Positive/Neutral", rot=90, cex=cex)
  plotLab("Negative", rot=90, cex=cex)

  if(p==1) {
    plotLab(CD$datasetname_year[1], cex=cex, x=0.22)
    plotLab(CD$datasetname_year[2], cex=cex, x=0.22)
    plotLab(CD$datasetname_year[3], cex=cex, x=0.22)
    plotLab(CD$datasetname_year[4], cex=cex, x=0.22)
  } else {
    plotLab(CD$datasetname_year[5], cex=cex, x=0.22)
    plotLab(CD$datasetname_year[6], cex=cex, x=0.22)
    plotLab(CD$datasetname_year[7], cex=cex, x=0.22)
  }

  par(mar=c(3,3.5,0.75,1))

  # --- Plot Data ---

  for(v in l_parts[[p]]) {

    # Load dataset specific data
    m_qntl_v <- l_qntl[[v]]
    m_skew_v <- l_skew[[v]]

    # Get indicator for postive/negative valence
    items <- CD$items_per_DS[[v]]
    ind_PENeuE <- which(items %in% CD$items_l_type[[1]])
    ind_NE <- which(items %in% CD$items_l_type[[2]])
    l_ind_e <- list(ind_PENeuE, ind_NE)

    # Plot for each item
    for(j in 1:2) {

      ## Layout
      plot.new()
      plot.window(xlim = c(0, 1), ylim=c(-5, 15))
      axis(1)
      axis(2)
      title(xlab="Prop(low)", ylab="Skewness", line=2)

      ## Data
      # aggregate
      v_skew_v_type <- as.numeric(m_skew_v[, l_ind_e[[j]]])
      v_qntl_v_type <- as.numeric(m_qntl_v[, l_ind_e[[j]]])
      # draw
      points(v_qntl_v_type, v_skew_v_type, pch=20, col=alpha("black", alpha=0.2))

      # linear regression
      lm_obj <- lm(v_skew_v_type~ v_qntl_v_type)
      abline(lm_obj, col="tomato", lty=1, lwd=2)
      lm_obj_sum <- summary(lm_obj)
      text(0.3, 6.8, paste0("R2 = ", round(lm_obj_sum$adj.r.squared, 2)), col="tomato")

    } # end: i

  } # end: v

  dev.off()


} # part of figure



# ------------------------------------------------------
# -------- Load Data -----------------------------------
# ------------------------------------------------------

# Load first subject from Rowland data
data <- readRDS("Data/Rowland2020/data_Rowland2020.RDS")
u_subj <- unique(data$subj_id)

i <- 1
data_i <- data[data$subj_id==u_subj[i], ]
data_i_vars <- data_i[, 5:12]
data_i_noNA <- na.omit(data_i)
data_i_noNA_vars <- data_i_noNA[, 5:12]

# ------------------------------------------------------
# -------- Fit VAR & generate data ---------------------
# ------------------------------------------------------

# ----- Fit VAR -----
out <- mvar(data=data_i_noNA_vars,
            type=rep("g", 8),
            level=rep(1, 8),
            lags = 1,
            lambdaSel = "EBIC",
            lambdaSeq = 0,
            dayvar = data_i_noNA$dayno,
            beepvar = data_i_noNA$beep,
            threshold = "none",
            scale=FALSE)

# intercepts
ints <- unlist(out$intercepts)

# phi-matrix
out$signs[is.na(out$signs)] <- 1
phi <- out$wadj[, , 1] * out$signs[, , 1]

# residual variances
pred <- predict(out, data=data_i_noNA_vars,
                dayvar = data_i_noNA_vars$dayno,
                beepvar = data_i_noNA_vars$beep)
res <- pred$predicted - data_i_noNA_vars
cov_res <- cov(res)

# ----- Generate Data -----

set.seed(13)
data_sim <- simulateVAR(pars = phi, means = ints,
                        Nt = 240,
                        init = as.numeric(data_i_noNA_vars[1,]),
                        residuals = cov_res) # inital values
range(unlist(data_sim))

# Combine into array
l_data <- list(data_i_vars, data_sim)


# ------------------------------------------------------
# -------- Plotting ------------------------------------
# ------------------------------------------------------

sc <- 0.9
pdf("Figures/Illu_VAR_fit.pdf", width = 9*sc, height = 12*sc)

dmat1 <- matrix(11:42, ncol=4, byrow = TRUE)
dmat2 <- cbind(3:10, dmat1)
dmat3 <- rbind(c(0,1,0,2,0), dmat2)
lo <- layout(mat=dmat3, widths = c(0.15,1,.3,1,.3), heights = c(.3,rep(1,8)))
# layout.show(lo)

# ----- Labels -----
# Rows
cex=1.75
plotLab("Empirical Data", cex=cex)
plotLab("Data Generated by fitted VAR", cex=cex)

# Colums
cex=1.5
library(stringr)
varnames <- str_to_title(colnames(data_i_noNA_vars))
for(i in 1:8) plotLab(varnames[i], rot=90, cex=cex, x=-0.5, y=.2)

# ----- Plot Data -----

for(j in 1:8) {
  for(s in 1:2) {

    ylim <- c(-80, 160)

    # Time Series
    par(mar=c(2,0,0,0))
    plot.new()
    plot.window(xlim=c(1,240), ylim=ylim)
    lines(l_data[[s]][, j], col="black")

    # Marginals
    par(mar=c(2,0,0,.5))
    breaks <- seq(ylim[1], ylim[2], length=30)
    tb <- hist(l_data[[s]][, j], breaks=breaks, plot = FALSE)
    barplot(tb$counts, axes = FALSE, horiz = TRUE, ylab="")

  }
}

dev.off()

















