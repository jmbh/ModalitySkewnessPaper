# jonashaslbeck@gmail.com; May 8, 2022

# ------------------------------------------------------
# -------- Load Packages -------------------------------
# ------------------------------------------------------

library(mclust)
library(mousetrap)
library(scales)
library(plyr)
library(moments)
library(dplyr)

# Load auxiliary functions
source("aux_Functions.R")

# ------------------------------------------------------
# -------- Specify Data Location & Variables -----------
# ------------------------------------------------------

baseDir <- "Data/"

data_files <- c("data_Rowland2020.RDS",
                "data_Bringmann2016.RDS",
                "data_Vrijen2018.RDS",
                "data_Fisher2017.RDS",
                "data_Bringmann2013.RDS",
                "data_Fried2021.RDS",
                "data_Wendt2020.RDS")

# Subsetting of emotion variables in each dataset
l_ind_vars <- list("Rowland2020" = 5:12, # 0-100 scale
                   "Bringmann2016" = 2:7, # 0-100
                   "Vrijen2018" = c(5,6,7,8), # 0-100
                   "Fisher2017" = c(2:16, 19:20), # 0-100
                   "Bringmann2013" = c(6,8:11), # 1-7
                   "Fried2021" = c(2:8,10:11), # 1-5
                   "Wendt2020" = 9:39) # 1-5 scale

names <- c("Rowland2020 (0-100)",
           "Bringmann2016 (0-100)",
           "Vrijen2018 (0-100)",
           "Fisher2017 (0-100)",
           "Bringmann2013 (1-7)",
           "Fried2021 (1-5)",
           "Wendt2020 (1-5)")

paper_names <- c("Rowland et al. (2020)",
                 "Bringmann et al. (2016)",
                 "Vrijen et al. (2018)",
                 "Fisher et al. (2017)",
                 "Bringmann et al. (2013)",
                 "Fried et al. (2021)",
                 "Wendt et al. (2020)")


# ------------------------------------------------------
# -------- Specify some Dataset characteristics --------
# ------------------------------------------------------

CD <- list()
CD$scale <- c(100, 100, 100, 100, 6, 4, 4)
CD$popclin <- c(0, 1, 0, 1, 1, 0, 1)
CD$mDay <- c(6, 10, 3, 4, 10, 4, 3.7)
CD$item_retro <- c(0, 0, 1, 1, 0, 1,0)

CD$datasetname_scale <- names
CD$datasetname_year <- paper_names

CD$scale_range_char <- c("0-100",
                         "0-100",
                         "0-100",
                         "0-100",
                         "1-7",
                         "1-5",
                         "1-5")

CD$scale_range_num <- rbind(c(0,100),
                            c(0,100),
                            c(0,100),
                            c(0,100),
                            c(1,7),
                            c(1,5),
                            c(1,5))

# ------------------------------------------------------
# -------- Descriptives I: Time series lengths ---------
# ------------------------------------------------------

# Get some summary infos from each data set
v_u_pt <- rep(NA, 7)
v_items <- rep(NA, 7)
v_TSlength <- rep(NA, 7)

for(v in 1:7) {

  # load
  data <- readRDS(paste0(baseDir, names(l_ind_vars)[v], "/", data_files[v]))
  data_var <- data[, l_ind_vars[[v]]]
  data_var <- cbind(data$subj_id, data_var)
  colnames(data_var)[1] <- "subj_id" # this is necessary, because bringmann2016 has weird neuroticism variable

  # number of participants and items
  v_u_pt[v] <- length(unique(data$subj_id))
  v_items[v] <- length(l_ind_vars[[v]])

  # mean time series length
  out <- ddply(data_var, .(subj_id), function(x) nrow(na.omit(x)) )[, 2]
  v_TSlength[v] <- mean(out)

}
round(v_TSlength)

CD$TSlength <- round(v_TSlength) # Mean Time Series lengths
CD$Ndist <- sum(v_u_pt * v_items) # Number of distributions
CD$unique_ptp_count <- v_u_pt # Unique subjects
CD$unique_item_count <- v_items # Unique items


# ------------------------------------------------------
# -------- Descriptives II: TS length per subject ------
# ------------------------------------------------------

l_TSlength <- list()
l_TSlength_IDs <- list()
for(v in 1:7) {

  data_v <- readRDS(paste0(baseDir, names(l_ind_vars)[v], "/", data_files[v]))

  # get relevant variables + subj_id
  data_v_aug <- cbind(data_v[, c(l_ind_vars[[v]])], data_v$subj_id)
  colnames(data_v_aug)[ncol(data_v_aug)] <- "subj_id"

  # Count rows
  TS_lengths_vs <- ddply(data_v_aug, .(subj_id), function(x)  {
    nrow(na.omit(x))
  } )
  l_TSlength[[v]] <- TS_lengths_vs

  # Save
  l_TSlength_IDs[[v]] <- TS_lengths_vs$subj_id

} # end loop: v


# ------------------------------------------------------
# -------- Descriptives III: Analyze total set of items-
# ------------------------------------------------------

l_items <- list()

for(v in 1:7) {
  data <- readRDS(paste0(baseDir, names(l_ind_vars)[v], "/", data_files[v]))
  var_set <- l_ind_vars[[v]]
  data_varset <- data[, var_set]
  l_items[[v]] <- tolower(colnames(data_varset))
}

# Manually recode some items:
l_items[[6]][1] <- "relaxed" # from "Relaxed"
l_items[[6]][3] <- "worried" # from "Worry"

items_col <- do.call(c, l_items)
items_unique <- unique(items_col)
n_items <- length(items_unique)

# Make table
m_items <- matrix(NA, nrow=n_items, ncol=7)
for(v in 1:7) {
  for(j in 1:n_items) {
    m_items[j,v] <- ifelse(items_unique[j] %in% l_items[[v]], 1, 0)
  }
}

# Save into dataset characteristic list
CD$items_per_DS <- l_items
CD$unique_items <- items_unique

# ----- Make unique item vector ordered by PA/NA -----

vPA <- c("happy", "excited", "relaxed", "satisfied", "joy",
         "energetic", "enthusiastic", "content", "positive",
         "determined", "inspired", "interested", "proud",
         "strong","active", "alert", "attentive")

vNA <- c("angry", "anxious", "depressed", "sad", "dysphoric",
         "irritated", "worried", "irritable", "restless",
         "guilty", "afraid", "anhedonia", "hopeless", "down",
         "fatigue", "tension", "ruminate",
         "nervous", "tired", "ashamed", "distressed",
         "hostile", "jittery", "scared", "upset", "frightened",
         "shaky", "scornful", "disgusted", "loathing", "blue",
         "downhearted", "lonely", "alone", "future", "threatened")

# As list
items_l_type <- list("PE"=vPA,
                     "NE"=vNA)
CD$items_l_type <- items_l_type

# Combined
items_unique_ord <- c(vPA, vNA)
CD$unique_items_ord <- items_unique_ord

# Sanity Check (sets should be the same)
items_unique %in% items_unique_ord
items_unique_ord %in% items_unique

# ------- Manually add retrospective or not ------
CD$item_ref <- c("current", "current", "retro", "retro", "current", "retro","current")


# ------------------------------------------------------
# -------- Estimate Modality ---------------------------
# ------------------------------------------------------

# Create storage
l_mB <- l_mB_Silver <- l_skew <- l_qntl <- list()
l_density <- l_density_Silver <- list()
l_bwvar <- list()

n <- 100 # no. of density-estimate x values

pdf("Figures/VisualCheck_Modality.pdf", 8, 8)

# Exclusion matrix
m_excl <- matrix(0, 7, 4) # 7 datasets, 3 reasons for exclusions, + cmb

par(mfrow=c(4,4))

# Loop over 7 data sets
set.seed(1)

for(v in 1:7) {

  print(v)

  # ----- Read data -----
  data <- readRDS(paste0(baseDir, names(l_ind_vars)[v], "/", data_files[v]))

  # ----- Fit univariate Mixtures -----

  # Prep indicators
  up <- unique(data$subj_id)
  n_p <- length(up)
  var_set <- l_ind_vars[[v]]
  n_vars <- length(var_set)

  # Create Storage
  m_BC <- m_skew <- m_qntl <- matrix(NA, n_p, n_vars)
  a_density <- array(NA, dim=c(n_p, n_vars, 2, n))

  # Loop and fit mixtures
  for(i in 1:n_p) {

    for(j in 1:n_vars) {

      # Get single time series
      data_i_j <- data[data$subj_id == up[i], var_set[j]]
      data_i_j <- na.omit(data_i_j)
      range_i_j <- max(data_i_j) - min(data_i_j)

      excl_ind <- 0

      # Common sense exclusion criteria
      if(length(data_i_j) == 0) {
        m_excl[v, 1] <- m_excl[v, 1] + 1
        excl_ind <- 1
      } else {
        if(sd(data_i_j) < 0.01) {
          m_excl[v, 2] <- m_excl[v, 2] + 1
          excl_ind <- 1
        }
        if(range_i_j==1) { # Avoid cases where range=1; here the noise adding would create multi-modality
          m_excl[v, 3] <- m_excl[v, 3] + 1
          excl_ind <- 1
        }
      } # end if: exclusion
      m_excl[v, 4] <- m_excl[v, 4] + excl_ind # to get total count (some of the exclusion criteria may be dependent)

      # --- My density-based Modality detection function ---

      if(excl_ind==0) {

        # noise hyper-parameter conditional on scale
        u_cat <- length(unique(data_i_j))
        max_th_range <- c(100, 100, 100, 100, 7, 5, 5)
        prop_modecat <- (max(table(data_i_j))/length(data_i_j))

        if(v %in% 1:4) { # continuous
          base_noise <- 0.035
          noise_adj0 <- base_noise + prop_modecat / 5
          noise <- noise_adj0 * max_th_range[v]
        } else { # ordinal
          noise <- prop_modecat
        }

        # RERUN 10 times to reduce arbitrariness of specific noise draws
        v_M <- rep(NA, 10)
        l_M <- list()
        for(r in 1:10) { # do 10 reruns
          l_M[[r]] <- DensMMdet(data_i_j, n=n, noise=noise, adjust=2)
          v_M[r] <- l_M[[r]]$M
        }

        tb <- table(v_M)
        m_BC[i,j] <- as.numeric(names(tb)[which.max(tb)]) # if equal, take smaller one
        den_out <- l_M[[which(v_M == m_BC[i,j])[1] ]] # save density of first one with selected #modes
        a_density[i, j, 1, ] <- den_out$den_x
        a_density[i, j, 2, ] <- den_out$den_y

        # ---- Quantile Information (for appendix) ----
        if(v %in% 1:4) { # continuous
          m_qntl[i,j] <- sum(data_i_j<=10)/length(data_i_j)
        } else { # ordinal
          m_qntl[i,j] <- sum(data_i_j==1)/length(data_i_j)
        }

      } # end if: ! exclusion criteria

      # --- PLOTTING ---

      if(length(data_i_j) != 0) {
        # Plot Histogram
        main <- paste0("v = ", v, " i=", i, " j=", j, " (", l_items[[v]][j], ")")
        hout <- hist(data_i_j, main = main, xlab="", breaks=20, xlim=CD$scale_range_num[v,])

        # Indicate whether excluded
        if(excl_ind==1) text((CD$scale[v]+1)/2, max(hout$counts)*.9, "EXCLUDED", col="red", adj=0.5)

        # If not excluded
        if(excl_ind==0) {

          # Add Modality & density est
          lines(den_out$den_x, (den_out$den_y/ max(den_out$den_y)) * max(hout$counts), col="red")
          text(max(data_i_j)*.9, max(hout$counts)*.9, paste0("Modes = ",den_out$M) , col="red", adj=1)

          # Add Skewness
          sk <- skewness(data_i_j)
          m_skew[i,j] <- sk
          sk_r <- round(skewness(data_i_j), 2)
          text(max(data_i_j)*0.9, max(hout$counts)*0.5, paste0("Skew = ", sk_r) , col="blue", adj=1)

        } # exclude: no cases at all

      } # end if: any data at all?

    } # end loop: variables

    # write between person variable for that person to a list

  } # end loop: individuals

  # Collect objects on v-level in listsw
  l_mB[[v]] <- m_BC
  l_density[[v]] <- a_density
  l_skew[[v]] <- m_skew
  l_qntl[[v]] <- m_qntl

} # end loop: datasets

dev.off()



# ------------------------------------------------------
# -------- Save Results --------------------------------
# ------------------------------------------------------

# List of dataset characteristics
saveRDS(CD, "Files/EMS_DS_info.RDS")

# M and density estimates (of DEN and GMM methods)
saveRDS(l_mB, "Files/EMS_Modality.RDS")
saveRDS(l_density, "Files/EMS_Densities.RDS")
saveRDS(l_skew, "Files/EMS_Skewness.RDS")

# save exclusion matrix
saveRDS(m_excl, "Files/EMS_Exclusion.RDS")

# save individual time series lengths
saveRDS(l_TSlength, "Files/EMS_TSlength.RDS")
saveRDS(l_TSlength_IDs, "Files/EMS_TSlength_IDs.RDS")

# save quantile business
saveRDS(l_qntl, "Files/EMS_Quantiles.RDS")










