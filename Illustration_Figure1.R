# jonashaslbeck@gmail.com; March 14, 2022

# ------------------------------------------------------
# -------- What is happening here? ---------------------
# ------------------------------------------------------

# Get 2x4 examples of: unimodal, unimodal skewed, bimodal, >=3 modal
# .. separate for continuous/Likert scales and for the same item

# ------------------------------------------------------
# -------- Source & Load packages ----------------------
# ------------------------------------------------------

# Load auxiliary functions
source("aux_Functions.R")

# ------------------------------------------------------
# -------- Loading Data --------------------------------
# ------------------------------------------------------

baseDir <- "Data/"
data_Rowland <- readRDS(paste0(baseDir, "Rowland2020/data_Rowland2020.RDS"))
data_MM <- readRDS(paste0(baseDir, "Bringmann2013/data_Bringmann2013.RDS"))

# ------------------------------------------------------
# -------- Examples for paper --------------------------
# ------------------------------------------------------

l_ex <- list()

# Continuous (From Rowland)
u_ptp_1 <- unique(data_Rowland$subj_id)
l_ex[[1]] <- X_c_1 <- na.omit(data_Rowland$sad[data_Rowland$subj_id == u_ptp_1[15]]) # unimodal
l_ex[[2]] <- X_c_2 <- na.omit(data_Rowland$sad[data_Rowland$subj_id == u_ptp_1[31]]) # unimodal skew
l_ex[[3]] <- X_c_3 <- na.omit(data_Rowland$sad[data_Rowland$subj_id == u_ptp_1[1]]) # bimodal
l_ex[[4]] <- X_c_4 <- na.omit(data_Rowland$sad[data_Rowland$subj_id == u_ptp_1[104]]) # multimodal

# Likert-scale
u_ptp_5 <- unique(data_MM$subj_id)
l_ex[[5]] <- X_L_1 <- na.omit(data_MM$Sad[data_MM$subj_id==u_ptp_5[4]]) # unimodal
l_ex[[6]] <- X_L_2 <- na.omit(data_MM$Sad[data_MM$subj_id==u_ptp_5[8]]) # unimodal skew
l_ex[[7]] <- X_L_2 <- na.omit(data_MM$Sad[data_MM$subj_id==u_ptp_5[59]]) # bimodal
l_ex[[8]] <- X_L_4 <- na.omit(data_MM$Sad[data_MM$subj_id==u_ptp_5[65]]) # multimodal


# ------------------------------------------------------
# -------- Make nice illustration Example  -------------
# ------------------------------------------------------

set.seed(2)
sc <- 1
pdf("Figures/Fig1_8Examples_Sad.pdf", width=8*sc, height=4.5*sc)

# Make Layout
lmat <- rbind(c(0,1:4),
              c(5,7:10),
              c(6,11:14))
lo <- layout(lmat, widths = c(.2, 1, 1), heights = c(.25, 1, 1))
# layout.show(lo)

# Plotting Labels
k_names <- c("Unimodal", "Unimodal", "Bimodal", "Multimodal")
k_n_2 <- c("(symmetric)", "(skewed)")

# Line 1
for(type in 1:4) {
  par(mar=rep(0,4))
  plot.new()
  plot.window(xlim=c(0, 1), ylim=c(0, 1))
  text(0.6, 0.7, adj=0.5, labels = k_names[type], cex=1.5, font=1)
  if(type %in% 1:2) text(0.6, 0.2, adj=0.5, labels = k_n_2[type], cex=1.5, font=1)
}

# Side
PL("0-100 Scale", srt=90)
PL("1-7 Scale", srt=90)

# Plotting Data
par(mar=c(4, 4, 1, 1))

v_xlim_b <- c(rep(100, 4), rep(7, 4))
v_xlim_a <- c(rep(0, 4), rep(1, 4))
xlim <- cbind(v_xlim_a, v_xlim_b)
ylim <- c(rep(70, 4), rep(50, 4))

v_vec <- c(rep(1, 4), rep(5,4))

for(i in 1:8) {

  if(i %in% 1:4) {
    breaks <- seq(0, 100, length=20)
  } else {
    breaks <- seq(1,7, length=20)
  }

  hout <- hist(l_ex[[i]], axes=FALSE, xlab="", main="", breaks=breaks,
       xlim=xlim[i, ], ylim=c(0, ylim[i]), ylab="")

  # --- Density Estimation ---
  data_i_j <- l_ex[[i]]
  u_cat <- length(unique(data_i_j))
  max_th_range <- c(100, 100, 100, 100, 7, 5, 5)
  prop_modecat <- (max(table(data_i_j))/length(data_i_j))

  v <- v_vec[i]
  if(v %in% 1:4) { # continuous
    base_noise <- 0.035
    noise_adj0 <- base_noise + prop_modecat / 5
    noise <- noise_adj0 * max_th_range[v]
  } else { # ordinal
    noise <- prop_modecat
  }

  den_out <- DensMMdet(data_i_j, n=100, noise=noise, adjust=2)
  lines(den_out$den_x, (den_out$den_y/ max(den_out$den_y)) * max(hout$counts), col="black")

  if(i %in% c(1,5)) title(ylab="Frequency")
  if(i %in% c(5:8)) title(xlab="Responses")

  axis(1, cex.axis=.9)
  axis(2, las=2)

}


dev.off()










