# Use bootstrap to estimate uncertainty on the TPR/FDR/FPR estimates
# 
# Plot ROC (Receiver Operating Curve) and Precision-Recall curve, use them to
# visually select "good" thresholds, i.e. threshold providing the best compromise
# between True Positives and False Positives.




## Initializations ----
library(tidyverse)
library(boot)

library(wbData)

gids <- wb_load_gene_ids(274)



#~ Read ground truth ----

# Gold standard expression patterns (cf table S5)
hbox_truth <- readRDS("data/072820_homeobox_truth.rds")
metabo_truth <- as.data.frame(readRDS("data/072820_metabotropic_receptors.rds"))
inx_truth <- readRDS("data/072820_innexin_truth.rds")
pan_truth <- as.data.frame(readRDS("data/072820_pan-neuronal_non-neuronal.rds"))

datasets <- list(hbox_truth, metabo_truth, inx_truth, pan_truth)


#~ Read dataset ----

# Proportion of cells with at least 1 UMI (generated separately)
prop_by_type <- readRDS("data/072820_L4_proportion_of_cells_expressing_each_gene.rds")


ground_truth <- bind_rows(datasets)
ground_truth$type <- rep(c("hbox", "metabo", "inx", "pan"), times = map_int(datasets, nrow))
ground_truth <- ground_truth[rownames(ground_truth) %in% rownames(prop_by_type), ] # no expression data

# save type stratifying column with same reshaping
source_types <- ground_truth$type
ground_truth$type <- NULL

all.equal(dim(ground_truth), c(154, 128))



#~ Definitions ----

# hard thresholds fixed
low_hard <- 0.02
high_hard <- 0.01

get_thres <- function(vals, level = 0.04, hard_low = 0.013, hard_high = 0.013){
  if(sum(vals) == 0) return(1)
  if(max(vals) < hard_low) return(1)
  if(min(vals) > hard_high) return(0)
  return(level*max(vals))
}


#~~~~~~~~~~~~~~~~~~~~~~~ ----


# Bootstraps ----

# Idea: for each threshold level, compute the bin matrix and run bootstraps

binarize <- function(dyn_lev, props, truth){
  cbind(1L*(props >= apply(props, 1, get_thres, dyn_lev, low_hard, high_hard)),
        truth)
}

#~~ define threshold percentiles ----


levels <- c(0, seq(1e-3,1e-2, 1e-4),
            seq(1.1e-2, 1e-1, 1e-3),
            seq(1.1e-1,1, 1e-2))

all_bins <- map(levels,
                binarize,
                prop_by_type[rownames(ground_truth), colnames(ground_truth)],
                ground_truth)
names(all_bins) <- as.character(levels)

#~~ callback functions ----
boot_tpr <- function(df_boot, i){
  # boot provides the original df and a vector of indexes selected in that bootstrap
  boot_bin <- df_boot[i,1:128]
  boot_truth <- df_boot[i,129:256]
  TPR <- sum(boot_bin * boot_truth)/sum(boot_truth)
}
boot_fdr <- function(df_boot, i){
  # boot provides the original df and a vector of indexes selected in that bootstrap
  boot_bin <- df_boot[i,1:128]
  boot_truth <- df_boot[i,129:256]
  FDR <- sum(boot_bin * (!boot_truth))/(sum(boot_bin))
}
boot_fpr <- function(df_boot, i){
  # boot provides the original df and a vector of indexes selected in that bootstrap
  boot_bin <- df_boot[i,1:128]
  boot_truth <- df_boot[i,129:256]
  FPR <- sum(boot_bin * (!boot_truth))/sum(!(boot_truth))
}


#~~ Loop bootstraps ----
## Running on Windows with snow
res_tpr <- map_dfr(all_bins, ~{
  res_boot_tpr <- boot(., boot_tpr, R=5000, strata = as.integer(as.factor(source_types)),
                       parallel="snow", ncpus=6)
  set_names(c(res_boot_tpr$t0,boot.ci(res_boot_tpr, type = "bca")$bca[4:5]),
            c("t0", "lower", "upper"))
})
res_fdr <- map_dfr(all_bins, ~{
  res_boot_tpr <- boot(., boot_fdr, R=5000, strata = as.integer(as.factor(source_types)),
                       parallel="snow", ncpus=6)
  set_names(c(res_boot_tpr$t0,boot.ci(res_boot_tpr, type = "bca")$bca[4:5]),
            c("t0", "lower", "upper"))
})
res_fpr <- map_dfr(all_bins, ~{
  res_boot_tpr <- boot(., boot_fpr, R=5000, strata = as.integer(as.factor(source_types)),
                       parallel="snow", ncpus=6)
  set_names(c(res_boot_tpr$t0,boot.ci(res_boot_tpr, type = "perc")$perc[4:5]),
            c("t0", "lower", "upper"))
})

names(res_tpr) <- c("TPR", "TPR_lower", "TPR_upper")
names(res_fpr) <- c("FPR", "FPR_lower", "FPR_upper")
names(res_fdr) <- c("FDR", "FDR_lower", "FDR_upper")

res_boot <- bind_cols(percentile = names(all_bins), res_tpr, res_fpr, res_fdr)




#~~~~~~~~~~~~~~~~~----

# Plots results ----

# -> Figures S6C, D

chosen_levels <- c(0.02, 0.04, 0.09, 0.15)

ggplot(res_boot, aes(x=1-FDR, y = TPR)) +
  theme_classic() +
  geom_ribbon(aes(ymin=TPR_lower, ymax=TPR_upper), fill = "gray") +
  geom_ribbon(aes(xmin=1-FDR_lower, xmax=1-FDR_upper), fill = "gray") +
  # geom_errorbar(aes(ymin=TPR_lower, ymax=TPR_upper), color = "gray") +
  # geom_errorbarh(aes(xmin=1-FDR_lower, xmax=1-FDR_upper), color = "gray") +
  geom_point() +
  geom_line() +
  # scale_x_continuous(limits=c(0.5,1)) +
  # scale_y_continuous(limits=c(0,1)) +
  xlab("Precision") +
  ylab("Recall") +
  geom_point(data=filter(res_boot, percentile %in% chosen_levels),
             color = "red", size=2)

ggsave("output/PR_rib.png",
       width = 5, height = 4, unit = "in")


ggplot(res_boot, aes(x=FPR, y = TPR)) +
  theme_classic() +
  geom_ribbon(aes(ymin=TPR_lower, ymax=TPR_upper), fill = "gray") +
  geom_ribbon(aes(xmin=FPR_lower, xmax=FPR_upper), fill = "gray") +
  # geom_errorbar(aes(ymin=TPR_lower, ymax=TPR_upper), color = "gray") +
  # geom_errorbarh(aes(xmin=FPR_lower, xmax=FPR_upper), color = "gray") +
  geom_point() +
  geom_line() +
  # scale_x_continuous(limits=c(0,.5)) +
  # scale_y_continuous(limits=c(0,1)) +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  geom_point(data=filter(res_boot, percentile %in% chosen_levels),
             color = "red", size=2)

ggsave("output/ROC_rib.png",
       width = 5, height = 4, unit = "in")


#~~ Save csv ----

# Not included in paper
write_csv(res_boot, "output/bootstrap_levels.csv")

#~~~~~~~~~~~~~ ----
