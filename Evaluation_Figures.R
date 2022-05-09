# jonashaslbeck@gmail.com; May 8, 2022

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

# ------------------------------------------------------
# -------- Load Results --------------------------------
# ------------------------------------------------------

# List of dataset characteristics
CD <- readRDS("Files/EMS_DS_info.RDS")

# M and density estimates (of DEN and GMM methods)
l_mB <- readRDS("Files/EMS_Modality.RDS")
l_density <- readRDS("Files/EMS_Densities.RDS")
l_skew <- readRDS("Files/EMS_Skewness.RDS")


# ------------------------------------------------------
# -------- Load Data Set Characteristics ---------------
# ------------------------------------------------------

unique_items_ord <- CD$unique_items_ord

# ------------------------------------------------------
# -------- Evaluation: Exclusion Criteria --------------
# ------------------------------------------------------

# Create Table for Appendix
tab_exl <- round(m_excl / sum(v_u_pt * v_items), 2) * 100
tab_exl <- as.data.frame(tab_exl)
rownames(tab_exl) <- paper_names
colnames(tab_exl) <- c("No data", "SD<0.01", "Range=1", "Exclusion")

print(xtable(tab_exl, digits=0), include.rownames=TRUE)


# ------------------------------------------------------
# -------- Evaluation: Table with Skew -----------------
# ------------------------------------------------------

# --- New indicator matrix for 4 types ---

cutoff <- 2/3

l_iT <- list()
for(v in 1:7) {

  m_mB_v <- l_mB[[v]]
  m_skew_v <- l_skew[[v]]

  m_iT <- m_mB_v
  m_iT[m_mB_v == 1 & abs(m_skew_v) > cutoff] <- 0

  m_iT[m_mB_v > 3] <- 3

  l_iT[[v]] <- m_iT
}


# --- Make Tab ---
m_tab <- matrix(NA, 7, 4)
for(v in 1:7) {

  # load mode estimates & skewness
  m_mB_v <- l_mB[[v]]
  m_skew_v <- l_skew[[v]]

  # modality table
  tb <- table(m_mB_v)

  # make indicator for: unimodal skew/not
  m_ind_skew <- array(NA, dim=dim(m_mB_v))
  m_ind_skew[m_mB_v == 1 & abs(m_skew_v) < cutoff] <- 0
  m_ind_skew[m_mB_v == 1 & abs(m_skew_v) > cutoff] <- 1

  tb_skew <- table(m_ind_skew)

  m_tab[v, 1:2] <- tb_skew
  m_tab[v, 3] <- tb[2]
  m_tab[v, 4] <- sum(tb[3:length(tb)])

}

# Normalize
for(v in 1:7)   m_tab[v, ] <- m_tab[v, ]/sum(m_tab[v, ], na.rm=T)
m_tab <- round(m_tab, 3)
rownames(m_tab) <- CD$datasetname_scale
colnames(m_tab) <- c("Uni sym", "Uni skew", "Bimodal", "Multimodal")
m_tab[is.na(m_tab)] <- 0
m_tab
m_tab_DEN <- m_tab
m_tab_DEN

# ------------------------------------------------------------------
# -------- Evaluation: Key Figure with Density line graphs ---------
# ------------------------------------------------------------------

# ---- Some Settings for Plotting  -----

# Scale the grey-scale with the number of data points
v_u_pt_norm <- 1 / ( v_u_pt / max(v_u_pt))
s <- 1:3
v_ref <- c("Current", "Retro")

set.seed(1) # Reproducibility (for selecting example distributions)

sc <- .9
pdf(paste0("Figures/Fig2_MainRes_4patterns.pdf"), width=6.5*sc, height=7.5*sc)

# --- Create Layout ---
lmat_core <- matrix(12:(11+28), 7, 4, byrow=T)
lmat_core <- cbind(1:7, lmat_core)
lmat_core <- rbind(c(0,8:11), lmat_core)
lo <- layout(mat = lmat_core, heights = c(.5, 1,1,1,1,1,1,1), width=c(1, 1, 1, 1))
# layout.show(lo)

# --- Add Paper Labels ---

for(v in 1:7) {
  par(mar=c(.2, 0, .2, 0))
  plot.new()
  plot.window(xlim=c(0, 1), ylim=c(0, 1))

  adj_l <- 0.05
  cex <- 0.8
  v_l_seq <- seq(0.95, 0.1, length=8)

  # Paper
  text(0, v_l_seq[1], adj=0, labels = CD$datasetname_year[v], font=2, cex=cex+0.03)
  # Items
  text(0, v_l_seq[2], adj=0, labels = paste0("Items: ", CD$unique_item_count[v]), cex=cex)
  # Subjects
  text(0, v_l_seq[3], adj=0, labels = paste0("Subjects: ", CD$unique_ptp_count[v]), cex=cex)
  # Scale
  text(0, v_l_seq[4], adj=0, labels = paste0("Scale: ", CD$scale_range_char[v]), cex=cex)
  # Mean Time Series Length
  text(0, v_l_seq[5], adj=0, labels = paste0("Av Time Points: ", CD$TSlength[v]), cex=cex)
  # Measures per day
  text(0, v_l_seq[6], adj=0, labels = paste0("Measures/Day: ", CD$mDay[v]), cex=cex)
  # Reference: retrospective vs current
  text(0, v_l_seq[7], adj=0, labels = paste0("Reference: ", v_ref[CD$item_retro[v] + 1]), cex=cex)
  # Population
  text(0, v_l_seq[8], adj=0, labels = paste0("Population: ", c("Students", "Clinical")[CD$popclin[v]+1]), cex=cex)

}

# --- Add Modality labels ---

k_names <- c("Unimodal", "Unimodal", "Bimodal", "Multimodal")
k_n_2 <- c("(symmetric)", "(skewed)")

# Line 1
for(type in 1:4) {
  par(mar=rep(0,4))
  plot.new()
  plot.window(xlim=c(0, 1), ylim=c(0, 1))
  text(0.45, 0.7, adj=0.5, labels = k_names[type], cex=1.2, font=1)
  if(type %in% 1:2) text(0.45, 0.3, adj=0.5, labels = k_n_2[type], cex=1.2, font=1)
}


# --- Plot Data ---

for(v in 1:7) {

  m_iT <- l_iT[[v]]
  a_density <-  l_density[[v]]

  # Loop over Ks
  for(type in 1:4) {  # uni sym, uni skew, bimod, multimod

    type_val <- c(1,0,2,3)

    # Get Densities for given dataset and K
    n_p <- dim(a_density)[1]
    n_var <- dim(a_density)[2]
    l_collect_k_x <- list()
    l_collect_k_y <- list()

    # Just densities for K=2 case
    counter <- 0
    for(i in 1:n_p) {
      for(j in 1:n_var) {

        if(!is.na(m_iT[i,j]))

          if(m_iT[i,j] == type_val[type]) {
            l_collect_k_x[[counter+1]] <- a_density[i, j, 1, ]
            l_collect_k_y[[counter+1]] <- a_density[i, j, 2, ]
            counter <- counter + 1
          }

      }# end for: variables
    } # end for: individuals


    # PLOTTING
    title <- paste0(m_tab[v, type]*100, "%")

    if(!(length(l_collect_k_x)==0)) {

      m_v_Kx_x <- as.matrix(do.call(rbind, l_collect_k_x))
      m_v_Kx_y <- as.matrix(do.call(rbind, l_collect_k_y))
      n_v_Kx <- nrow(m_v_Kx_y) # number of distributions with K=k

      # --------- Alternative: a few sample graphs ----------
      par(mar=c(1,1.5,2,1.5))
      plot.new()
      plot.window(xlim=CD$scale_range_num[v, ], ylim=c(0,1)) # only plot original data range

      # pre-Select those that are most pronounced multi-modal
      nD_set <- 5
      nD_1 <- ifelse(round(nrow(m_v_Kx_y)/2) < nD_set*2, nrow(m_v_Kx_y), round(nrow(m_v_Kx_y)/4)) #
      sd_den <- apply(m_v_Kx_y, 1, max)
      # Take care of pathological n=1 case
      m_v_Kx_ss <- matrix(matrix(m_v_Kx_y[order(sd_den, decreasing = FALSE), ], n_v_Kx, 100)[1:nD_1,], nD_1, 100)
      m_v_Kx_ss_X <- matrix(matrix(m_v_Kx_x[order(sd_den, decreasing = FALSE), ], n_v_Kx, 100)[1:nD_1,], nD_1, 100)

      nD_2 <- min(5, nrow(m_v_Kx_ss))
      int_samp <- sample(1:nrow(m_v_Kx_ss), size=nD_2)
      m_v_Kx_ss_2 <- m_v_Kx_ss[int_samp, ]
      m_v_Kx_ss_X_2 <-m_v_Kx_ss_X[int_samp, ]

      # Take care of pathological n=1 case
      m_v_Kx_ss_2 <- matrix(m_v_Kx_ss_2, nD_2, 100)
      m_v_Kx_ss_X_2 <- matrix(m_v_Kx_ss_X_2, nD_2, 100)

      # fit into [0,1] interval
      m_v_Kx_ss_norm <- (m_v_Kx_ss_2 / max(m_v_Kx_ss_2)) * 0.9

      # Order left to right
      n_den <- ncol(m_v_Kx_ss_norm)
      v_w_loc <- apply(m_v_Kx_ss_norm, 1, function(x) mean(1:n_den * x))
      ord <- order(v_w_loc, decreasing = FALSE)
      m_v_Kx_ss_norm_ord <- m_v_Kx_ss_norm[ord, ]
      m_v_Kx_ss_norm_ord <- matrix(m_v_Kx_ss_norm_ord, nD_2, 100)
      m_v_Kx_ss_X_2 <- matrix(m_v_Kx_ss_X_2, nD_2, 100)
      m_v_Kx_ss_X_2_ord <- m_v_Kx_ss_X_2[ord, ]
      m_v_Kx_ss_X_2_ord <- matrix(m_v_Kx_ss_X_2_ord, nD_2, 100)

      # Plot
      nD_final <- nrow(m_v_Kx_ss_norm_ord)
      cols <- brewer.pal(nD_set, "Set1")
      for(nD_i in 1:nD_final) lines(m_v_Kx_ss_X_2[nD_i, ], m_v_Kx_ss_norm_ord[nD_i, ],
                                    col=alpha(cols[nD_i], alpha=.75),
                                    lwd=2)


    } else {
      plot.new()
      plot.window(xlim=c(0,1), ylim=c(0,1)) # leave empty
    }

    # Plot percentage
    mtext(text = title, side=3, cex=0.8, font=1, col="black")

  } # end for: Ks

} # end for: data sets (v)

dev.off()


# ------------------------------------------------------------------
# -------- Figure Results Item-level Analysis [Multi-modality] -----
# ------------------------------------------------------------------

N_items <- length(CD$unique_items)

# Choose colors
library(RColorBrewer)
cols <- brewer.pal(8, "Set1")[-6]

sc <- 0.8
pdf("Figures/Figure_Result_ItemLevel.pdf", width=8*sc, height=6*sc)

# --- Layout ---
layout(matrix(1:2, ncol=1), heights = c(1, .1))

# --- Layout top panel ---
par(mar=c(4,3.5,1,1))
plot.new()
plot.window(xlim=c(1,N_items), ylim=c(.9,2))
axis(1, labels = unique_items_ord, at=1:N_items, las=2, cex.axis=.55)
axis(2, seq(1, 2, length=3), las=2)
title(ylab="Average Number of Modes", line=2.5)
segments(x0 = 1:N_items,
         y0 = rep(0, N_items),
         x1 <-1:N_items,
         y1 = rep(2, N_items), col="grey")

# --- Plot Data ---
set.seed(1) # for jitter
for(v in 1:7) {
  m_BC <- l_mB[[v]]
  p <- ncol(m_BC)
  v_ind <- rep(NA, p)
  for(j in 1:p) v_ind[j] <- which(l_items[[v]][j]  == unique_items_ord)


  # Look at proportion
  prop_v <- apply(m_BC, 2, function(x) mean(x==1))

  # Means
  means_v <- colMeans(m_BC, na.rm = TRUE)
  points(jitter(v_ind, factor = .5), means_v,
         pch=20, col=cols[v], cex=1.5)

} # end for: within

# --- Legend ---
rect(38, 1.55, 56, 2, col="white", border=FALSE)
legend(37, 2, legend=paper_names, text.col=cols, cex=0.8, bg="white", bty="n")


# --- Bottom panel: indicating PE/nE/NE ---
par(mar=c(0,3.5,0,1))
plot.new()
plot.window(xlim=c(1,N_items), ylim=c(-1,1))
# Draw arrows
v_length_i_type <- unlist(lapply(CD$items_l_type, length))
a_head <- 0.09
a_space <- 0.15
arrows(1+a_space,.5,v_length_i_type[1]-a_space+0.5,0.5, code=3, length=a_head)
arrows(v_length_i_type[1]+a_space+0.5,.5,sum(v_length_i_type[1:2])-a_space+0.5,0.5, code=3, length=a_head)
# arrows(sum(v_length_i_type[1:2])+a_space+0.5,.5,sum(v_length_i_type[1:3])-a_space+0.5,0.5, code=3, length=a_head)
# Emotion label
cex.text <- 0.7
text(0.5+v_length_i_type[1]/2, -0.25, "Positive Emotions", cex=cex.text)
text(0.5+v_length_i_type[1] + v_length_i_type[2]/2, -0.25, "Negative Emotions", cex=cex.text)
# text(0.5+sum(v_length_i_type[1:2]) + v_length_i_type[3]/2, -0.25, "Negative Emotions", cex=cex.text)


dev.off()


# ------------------------------------------------------------------
# -------- Figure Results Item-level Analysis [Unimodal Skew] ------
# ------------------------------------------------------------------

N_items <- length(CD$unique_items)

# Choose colors
cols <- brewer.pal(8, "Set1")[-6]

sc <- 0.8
pdf("Figures/Figure_Result_ItemLevel_Skew.pdf", width=8*sc, height=6*sc)

# --- Layout ---
layout(matrix(1:2, ncol=1), heights = c(1, .1))

# --- Layout ---
par(mar=c(4,3.5,1,1))
plot.new()
# plot.window(xlim=c(1,55), ylim=c(-4,10))
plot.window(xlim=c(1,N_items), ylim=c(-1,6))
axis(1, labels = unique_items_ord, at=1:N_items, las=2, cex.axis=.55)
axis(2, las=2)
title(ylab="Skewness", line=2.5)
segments(x0 = 1:N_items,
         y0 = rep(-1, N_items),
         x1 = 1:N_items,
         y1 = rep(6, N_items), col="grey")
abline(h=0, col="grey")

# --- Plot Data ---
set.seed(1) # for jitter
for(v in 1:7) {

  # Get skewness & modality for v
  m_skew <- l_skew[[v]]
  m_BC <- l_mB[[v]]

  # Get indices of items of v
  p <- ncol(m_BC)
  v_ind <- rep(NA, p)
  for(j in 1:p) v_ind[j] <- which(CD$items_per_DS[[v]][j]  == CD$unique_items_ord)

  # Subset unimodal data
  l_skew_j <- list()
  for(j in 1:p) l_skew_j[[j]] <- m_skew[m_BC[, j]==1, j]

  # Plot means
  points(v_ind, lapply(l_skew_j, function(x) mean(x, na.rm = TRUE)),
         pch=20, col=alpha(cols[v], alpha=1), cex=1.5)

} # end for: within

# --- Legend ---
rect(0, 3.5, 16.9, 6, col="white", border=FALSE)
legend("topleft", legend=paper_names, text.col=cols, cex=0.8, bg="white", bty="n")


# --- Bottom panel: indicating PE/nE/NE ---
par(mar=c(0,3.5,0,1))
plot.new()
plot.window(xlim=c(1,N_items), ylim=c(-1,1))
# Draw arrows
v_length_i_type <- unlist(lapply(CD$items_l_type, length))
a_head <- 0.09
a_space <- 0.15
arrows(1+a_space,.5,v_length_i_type[1]-a_space+0.5,0.5, code=3, length=a_head)
arrows(v_length_i_type[1]+a_space+0.5,.5,sum(v_length_i_type[1:2])-a_space+0.5,0.5, code=3, length=a_head)
# arrows(sum(v_length_i_type[1:2])+a_space+0.5,.5,sum(v_length_i_type[1:3])-a_space+0.5,0.5, code=3, length=a_head)
# Emotion label
cex.text <- 0.7
text(0.5+v_length_i_type[1]/2, -0.25, "Positive Emotions", cex=cex.text)
text(0.5+v_length_i_type[1] + v_length_i_type[2]/2, -0.25, "Negative Emotions", cex=cex.text)
# text(0.5+sum(v_length_i_type[1:2]) + v_length_i_type[3]/2, -0.25, "Negative Emotions", cex=cex.text)


dev.off()


# ------------------------------------------------------------------
# -------- Figure Results Person-level Analysis [Modality] ---------
# ------------------------------------------------------------------
# New: May 5th

# Choose colors
library(RColorBrewer)
cols <- brewer.pal(8, "Set1")[-6]

sc <- 0.7
pdf("Figures/Fig_PersonLevel_Modality_May5.pdf", width=8*sc, height=4.7*sc)

# --- Layout ---
lmat <- matrix(1:4, nrow=2, byrow = TRUE)
lo <- layout(lmat, widths = c(1,.3), heights = c(1, .1))
# layout.show(lo)

# --- Data Canvas ---
par(mar=c(0,4,1,0))
plot.new()
plot.window(xlim=c(0.5,14), ylim=c(.95,3))
axis(2, seq(1, 3, length=5), las=2)
cex.lab <- 0.9
title(ylab="Average Number of Modes", line=2.5, cex.lab=cex.lab)

# -- Get PE/NE indicators for each dataset --

# get indicator: PE/NE
l_ind <- list()
for(v in 1:7) {
  m_BC <- l_mB[[v]]
  p <- ncol(m_BC)
  v_ind <- rep(NA, p)
  for(j in 1:p) {
    if(CD$items_per_DS[[v]][j] %in% CD$items_l_type$NE) v_ind[j] <- 0 else v_ind[j] <- 1
  }
  l_ind[[v]] <- v_ind
}

# --- Plot Data ---
v_mp <- c(-1, 1)
c1 <- 1
for(q in 1:2) {

  for(v in 1:7) {

    m_BC <- l_mB[[v]]

    # Plot Boxplots & data for PE
    Rm_PE <- rowMeans(m_BC[, l_ind[[v]]==(1:0)[q], drop=FALSE], na.rm = TRUE)
    boxplot(Rm_PE, add=TRUE, col="white", border=cols[v],
            axes=FALSE, at=c1, outline=FALSE, width=1.5)

    X_plot <- rep(c1, CD$unique_ptp_count[v])
    points(X_plot, Rm_PE, pch=20, col=alpha(cols[v], 0.075), cex=2)

    c1 <- c1 + 1

  } # end for: within

}


# --- Legend ---
par(mar=c(0,0,0,0))
plot.new()
plot.window(xlim=c(0,1), ylim=c(0,1))
legend(-0.2, .75, legend=paper_names, text.col=cols,
       cex=0.85, bg="white", bty="n", adj=0)

# --- PE/NE Axis ---
a_head <- 0.09
gap <- 0.05
arr_h <- .74
par(mar=c(0,4,0,0))
plot.new()
plot.window(xlim=c(1,14), ylim=c(-1,1))
arrows(1.25, arr_h, 7.5-gap, arr_h, code=3, length=a_head)
arrows(7.5+gap, arr_h, 14.25, arr_h, code=3, length=a_head)
text(2.35, 0, "Positive Emotions", adj=0, cex=0.87)
text(2.35+6.3, 0, "Negative Emotions", adj=0, cex=0.87)

dev.off()


# ------------------------------------------------------------------
# -------- Figure Results Person-level Analysis [Skewness] ---------
# ------------------------------------------------------------------
# New: May 5th

# Choose colors
library(RColorBrewer)
cols <- brewer.pal(8, "Set1")[-6]

sc <- 0.7
pdf("Figures/Fig_PersonLevel_skew_May5.pdf", width=8*sc, height=4.7*sc)

# --- Layout ---
lmat <- matrix(1:4, nrow=2, byrow = TRUE)
lo <- layout(lmat, widths = c(1,.3), heights = c(1, .1))
# layout.show(lo)

# --- Data Canvas ---
par(mar=c(0,4,1,0))
plot.new()
plot.window(xlim=c(0.5,14), ylim=c(-3,10))
axis(2, -3:10, las=2)
cex.lab <- 0.9
title(ylab="Skewness", line=2.5, cex.lab=cex.lab)

# Zero line
segments(0, 0, 14.5, 0, lty=1, col="grey")

# -- Get PE/NE indicators for each dataset --

# get indicator: PE/NE
l_ind <- list()
for(v in 1:7) {
  m_BC <- l_mB[[v]]
  p <- ncol(m_BC)
  v_ind <- rep(NA, p)
  for(j in 1:p) {
    if(CD$items_per_DS[[v]][j] %in% CD$items_l_type$NE) v_ind[j] <- 0 else v_ind[j] <- 1
  }
  l_ind[[v]] <- v_ind
}

# --- Plot Data ---
v_mp <- c(-1, 1)
c1 <- 1
for(q in 1:2) {

  for(v in 1:7) {

    m_skew <- l_skew[[v]]

    # Plot Boxplots & data for PE
    Rm_PE <- rowMeans(m_skew[, l_ind[[v]]==(1:0)[q], drop=FALSE], na.rm = TRUE)
    boxplot(Rm_PE, add=TRUE, col="white", border=cols[v],
            axes=FALSE, at=c1, outline=FALSE, width=1.5)

    X_plot <- rep(c1, CD$unique_ptp_count[v])
    points(X_plot, Rm_PE, pch=20, col=alpha(cols[v], 0.075), cex=2)

    c1 <- c1 + 1

  } # end for: within

}


# --- Legend ---
par(mar=c(0,0,0,0))
plot.new()
plot.window(xlim=c(0,1), ylim=c(0,1))
legend(-0.2, .75, legend=paper_names, text.col=cols,
       cex=0.85, bg="white", bty="n", adj=0)

# --- PE/NE Axis ---
a_head <- 0.09
gap <- 0.05
arr_h <- .74
par(mar=c(0,4,0,0))
plot.new()
plot.window(xlim=c(1,14), ylim=c(-1,1))
arrows(1.25, arr_h, 7.5-gap, arr_h, code=3, length=a_head)
arrows(7.5+gap, arr_h, 14.25, arr_h, code=3, length=a_head)
text(2.35, 0, "Positive Emotions", adj=0, cex=0.87)
text(2.35+6.3, 0, "Negative Emotions", adj=0, cex=0.87)

dev.off()




