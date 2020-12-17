# Compute TPR/FPR/FDR for each neuron type with bootstrap replicates.
# 
# We use the thresholds selected with the whole dataset bootstrapping.



## Initializations ----
library(tidyverse)
library(boot)

library(wbData)

gids <- wb_load_gene_ids(274)



#~ Read data ----

# Gold standard expression patterns (cf table S5)
hbox_truth <- readRDS("data/072820_homeobox_truth.rds")
metabo_truth <- as.data.frame(readRDS("data/072820_metabotropic_receptors.rds"))
inx_truth <- readRDS("data/072820_innexin_truth.rds")
pan_truth <- as.data.frame(readRDS("data/072820_pan-neuronal_non-neuronal.rds"))

datasets <- list(hbox_truth, metabo_truth, inx_truth, pan_truth)

# Proportion of cells with at least 1 UMI (generated separately)
prop_by_type <- readRDS("data/072820_L4_proportion_of_cells_expressing_each_gene.rds")
# Number of cells for each neuron type (generated separately)
nb_cells <- readRDS("data/081820_neuron_nb_cells_df.rds")


#~ Format data
ground_truth <- bind_rows(datasets)
ground_truth$type <- rep(c("hbox", "metabo", "inx", "pan"), times = map_int(datasets, nrow))
ground_truth <- ground_truth[rownames(ground_truth) %in% rownames(prop_by_type), ] # no expression data

# save type stratifying column with same reshaping
source_types <- ground_truth$type
ground_truth$type <- NULL

all.equal(dim(ground_truth), c(154, 128))






### Bootstraps ----

# -- Technical note --
# We want a df with 4 thresh, 128 cells at each, 3 diagnostics at each
# Pseudocode possibility 1:
# for (thresh in [1,2,3,4])
#   get thresholded matrix
#   for (cell in [128 neuron types])
#     do Bootstraps
#        get TPR, FDR, FPR
# 
# But actually more efficient to loop neurons inside bootstrap
# for each thresh
#   do Bootstrap
#     compute each combination neuron x diag




#~ Functions ----

get_thres_max_hard <- function(vals, level = 0.04, hard_low, hard_high){
  # threshold a vector by a percentile level
  if(sum(vals) == 0) return(1)
  if(max(vals) < hard_low) return(1)
  if(min(vals) > hard_high) return(0)
  return(level*max(vals))
}
  
boot_diags_per_neur <- function(df_boot, idx){
  # compute all 3 diagnostics for each neuron type
  boot_bin <- df_boot[idx,1:128]
  boot_truth <- df_boot[idx,129:256]
  tpr_types <- colSums(boot_bin * boot_truth)/colSums(boot_truth)
  fdr_types <- colSums(boot_bin * (!boot_truth))/colSums(boot_bin)
  fpr_types <- colSums(boot_bin * (!boot_truth))/colSums(!boot_truth)
  
  names(tpr_types) <- paste("TPR", names(tpr_types), sep="_")
  names(fdr_types) <- paste("FDR", names(fdr_types), sep="_")
  names(fpr_types) <- paste("FPR", names(fpr_types), sep="_")
  
  return(c(tpr_types,fdr_types,fpr_types))
}


#~ Run bootstrap ----
low_hard <- 0.02
high_hard <- 0.01

thres_levels <- c("level_1" = 0.02,
                  "level_2" = 0.04,
                  "level_3" = 0.09,
                  "level_4" = 0.15)

props <- prop_by_type[rownames(ground_truth), colnames(ground_truth)]


boot_res <- vector("list", length(thres_levels))
names(boot_res) <- thres_levels

for(thr in thres_levels){
  thr_mat <- 1L*(props >= apply(props,
                                1,
                                get_thres_max_hard, thr, low_hard, high_hard))
  
  res <- boot(cbind(thr_mat, ground_truth),
              boot_diags_per_neur,
              R=5000,
              strata = as.integer(as.factor(source_types)))
  
  boot_res[[as.character(thr)]] <- bind_cols(enframe(res$t0),
                               map_dfr(seq_along(res$t0),
                                       ~ tryCatch(boot.ci(res, type="bca", index=.x)$bca,
                                                  error = function(e) rep(NA_real_, 5)) %>%
                                         set_names(c("conf", "kl", "ku","lower","upper")))) %>%
    select(-c(conf, kl, ku)) %>%
    separate(name, into=c("diag", "neuron"), extra="merge") %>%
    left_join(nb_cells, by="neuron") %>%
    add_column(threshold = thr)
}

boot_res <- bind_rows(boot_res)





# ~~~~~~~~~~~~~~~~~ ----



# Plotting ----



boot_res %>%
  mutate(threshold = as.factor(threshold)) %>%
  ggplot(aes(x=nb_cells, y=value, color=threshold)) +
  theme_classic() +
  geom_errorbar(aes(ymin=lower, ymax=upper), color="grey") +
  geom_point() +
  facet_wrap(~diag) +
  scale_x_log10() +
  geom_smooth(method="lm") +
  ggtitle("All information")

boot_res %>%
  filter(threshold == 0.04) %>%
  ggplot(aes(x=nb_cells, y=value)) +
  theme_classic() +
  geom_errorbar(aes(ymin=lower, ymax=upper), color="grey") +
  geom_point() +
  facet_wrap(~diag) +
  scale_x_log10() +
  geom_smooth(method="lm", se = FALSE) +
  ggtitle("Medium threshold")

boot_res %>%
  filter(threshold == 0.04) %>%
  ggplot(aes(x=nb_cells, y=value)) +
  theme_classic() +
  geom_point() +
  facet_wrap(~diag) +
  scale_x_log10() +
  geom_smooth(method="lm") +
  ggtitle("Medium threshold, se of fit")

boot_res %>%
  filter(threshold == 0.04) %>%
  ggplot(aes(x=nb_cells, y=value)) +
  theme_classic() +
  geom_point() +
  facet_wrap(~diag) +
  scale_x_log10() +
  geom_smooth(method="lm", se=FALSE,
              data=filter(boot_res,
                          threshold == 0.04,
                          nb_cells < 100), color='red') +
  geom_smooth(method="lm", se=FALSE,
              data=filter(boot_res,
                          threshold == 0.04,
                          nb_cells >= 100), color='green') +
  ggtitle("Medium threshold, se of fit")

# => the TPR depends on the number of cells in the scRNA-Seq cluster, whereas the FPR and FDR are stable.
# In other words, more cells gives more statistical power, but we don't increase type I error when not enough cells

#~ Export table ----

# Table S6 (note, table was reformatted before inclusion)
boot_res %>%
  pivot_wider(id_cols= c(neuron, threshold, nb_cells),
              names_from = diag,
              names_glue="{diag}_{.value}",
              values_from=c(value, lower, upper)) %>%
  write_csv("output/table_S6_boot_per_type.csv")



