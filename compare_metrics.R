# Compare different metrics for thresholding
# 
# 
# We compare here 4 main metrics, all containing information for each gene in 
# each neuron type.  This includes:
#   1) the number of cells of each cell type expressing the gene,
#   2) the proportion of cells of each type expressing the gene,
#   3) the number of UMIs detected for each gene
#          in each neuron type (summed across all cells), and
#   4) the averaged TPM of each gene in each
#         neuron type (averaged across the cells in each neuron type).
# 
# In addition, several normalizations are also explored.
# Note: here we use a split of the ground truth into a training and testing dataset, later
# both are regrouped for the bootstrap.


## Init ----
library(tidyverse)




## Functions ----
get_tpr <- function(expression, truth, threshold, na.rm = TRUE){
  # True Positive Rate, aka sensitivity, aka recall
  # TPR = TP/(TP+FN) = TP/P
  bin <- expression >= threshold
  return(sum(bin * truth)/sum(truth))
}
get_fpr <- function(expression, truth, threshold, na.rm = TRUE){
  # False Positive Rate
  # FPR = FP/(FP+TN) = FP/N
  bin <- expression >= threshold
  return(sum(bin * (!truth))/sum(!(truth)))
}
get_fdr <- function(expression, truth, threshold, na.rm = TRUE){
  # False Discovery Rate
  # FDR = FP/(FP+TP) = 1 - PPV
  bin <- expression >= threshold
  fdr <- sum(bin * (!truth))/(sum(bin*(!truth)) + sum(bin*truth))
  if(is.nan(fdr))
    fdr <- 0
  return(fdr)
}


define_threshold <- function(metric = c("n_umi", "ncells",
                                        "proportion", "tpm",
                                        "umi_norm", "sct",
                                        "sf_monocle", "scater", "random"), max = 1){
  # generalized seq() suited for each case
  metric <- match.arg(metric)
  
  
  if(metric %in% c("tpm", "umi_norm", "sf_monocle", "scater")){
    # non-integer from 0 to max
    return(c(0,
             10**(seq(0, log10(max), 2*log10(max)/1000))/(10*max),
             10**(seq(0, log10(max), 2*log10(max)/1000))))
  }
  
  if(metric %in% c("proportion", "random"))
    # 0 to 1, smaller step for low values
    return(c(seq(0, 1e-2, 1e-5),
             seq(1.1e-2, 0.5, 1e-3),
             seq(0.51,1,0.1)))
  
  
  if(metric %in% c("n_umi", "ncells", "sct")){
    # Integer from 0 to max
    x <- 0:190
    for(i in 2:(floor(log10(max))-1)){
      x <- c(x, seq(from=2*10**i, to=10**(i+1), by=10**(i-1)))
    }
    return(c(x, seq(from=2*10**i, to=max, by=10**(i-1))))
  }
}

normalize_umi_mat <- function(mat){
  return(1e6*t(t(mat)/colSums(mat)))
}


# Read ground truth ----

# Golden standard expression patterns (cf table S5)
hbox_truth <- readRDS("github_summarized/data/072820_homeobox_truth.rds")
metabo_truth <- as.data.frame(readRDS("github_summarized/data/072820_metabotropic_receptors.rds"))
inx_truth <- readRDS("github_summarized/data/072820_innexin_truth.rds")
pan_truth <- as.data.frame(readRDS("github_summarized/data/072820_pan-neuronal_non-neuronal.rds"))

# Make training and test datasets by taking 70% of each dataset
set.seed(123)
datasets <- list(hbox_truth, metabo_truth, inx_truth, pan_truth)
training_genes <- datasets %>%
  map(~sample(rownames(.), 0.7*floor(nrow(.))))

testing_genes <- datasets %>%
  map2(training_genes, ~rownames(.x)[which(!rownames(.x) %in% .y)])

neurs <- datasets %>%
  map(~ colnames(.)) %>%
  reduce(intersect)

training_set <- datasets %>%
  map2(training_genes, ~ .x[.y, neurs]) %>%
  reduce(rbind)

testing_set <- datasets %>%
  map2(testing_genes, ~ .x[.y, neurs]) %>%
  reduce(rbind)




## Use more metrics ----
metrics <- c("n_umi", "ncells", "proportion", "tpm")

values <- map(metrics,
               ~ readRDS(paste0("raw/030420_L4_all_neurons_", .x,"_cell_type.rds"))) %>%
  set_names(metrics)


# Compute normalized UMI counts
values[["umi_norm"]] <- normalize_umi_mat(values[["n_umi"]])

# other normalizations (computed separately)
values[["sct"]] <- readRDS("data/L4_norm_sct.rds")
values[["sf_monocle"]] <- readRDS("data/L4_norm_sf.rds")
values[["scater"]] <- readRDS("data/L4_norm_scater.rds")

# randomize the proportion per cell input
values[["random"]] <-  matrix(sample(values[["proportion"]]),
                              nrow=nrow(values[["proportion"]]),
                              ncol=ncol(values[["proportion"]]),
                              dimnames = dimnames(values[["proportion"]]))

metrics <- c("n_umi", "ncells", "proportion", "tpm", "umi_norm", "sct", "sf_monocle", "scater", "random")




# ~~ Select ground truth to use ----
mytruth <- training_set

# ~~~ Compute diagnostics for each condition ----

# matrices need to all have same columns and rows
neurons <- intersect(Reduce(intersect, lapply(values, colnames)),
                      colnames(mytruth))

genes <- intersect(Reduce(intersect, lapply(values, rownames)),
                    rownames(mytruth))

length(neurons)
length(genes)


truth <- mytruth[genes,neurons]

values_restricted <- list()
for(mt in metrics){
  values_restricted[[mt]] <- values[[mt]][genes,neurons]
}

# get diagnostics
tib_diags <- list()
for(mt in metrics){
  thres <- define_threshold(mt, max(unlist(values_restricted[[mt]])))
  tib_diags[[mt]] <- tibble(threshold = thres,
                            TPR = map_dbl(thres, ~get_tpr(values_restricted[[mt]],
                                                          truth,
                                                          .x)),
                            FPR = map_dbl(thres, ~get_fpr(values_restricted[[mt]],
                                                          truth,
                                                          .x)),
                            FDR = map_dbl(thres, ~get_fdr(values_restricted[[mt]],
                                                          truth,
                                                          .x))
  )
}

# convert to tibble for plotting
tib <- tibble()
for(mt in metrics){
  tib <- bind_rows(tib,
                   tibble(metric = rep(mt, nrow(tib_diags[[mt]])),
                          TPR = tib_diags[[mt]][["TPR"]],
                          FPR = tib_diags[[mt]][["FPR"]],
                          FDR = tib_diags[[mt]][["FDR"]]
                   )
  ) %>% mutate(metric = fct_inorder(metric))
}

# ~~~ plot metrics ----
metrics_dictionary <- data.frame(full = c("1 Number of UMI (raw)",
                                          "2 Number of cells expressing (raw)",
                                          "3 Proportion of cells expressing",
                                          "4 Deconvolved size factors (Lun 2016, scater)",
                                          "5 scTransform (Seurat default)",
                                          "6 Monocle Size Factors + Packer TPM (recomputed)",
                                          "7 Monocle Size Factors + Packer TPM",
                                          "8 Number of UMI normalized by type",
                                          "9 Randomized proportions (control)"),
                                 row.names = metrics)

tib %>%
  ggplot(mapping = aes(x=FPR, y=TPR, color=metric)) +
  geom_step() +
  # geom_point(size=1) +
  scale_color_brewer(type="qual", palette = "Set3") +
  theme_classic() +
  ggtitle("ROC for all metrics")


tib %>%
  ggplot( mapping = aes(x=1-FDR, y=TPR, color=metric)) +
  geom_step() +
  # geom_point(size=1) +
  scale_color_brewer(type="qual", palette = "Set3") +
  theme_classic() +
  labs(title="PR curves for all metrics",
       x="Precision (1-FDR)", y="Recall (TPR)")


# partial AUC
tib %>%
  group_by(metric) %>%
  filter(FPR >= min(tib$FPR[tib$FPR > 0]),
         FPR <= max(tib$FPR[tib$FPR < 1])) %>%
  summarize(ROC_pAUC = DescTools::AUC(FPR, TPR)) %>%
  mutate(metric_descr = metrics_dictionary[metric, "full"]) %>%
  ggplot(aes(x=metric_descr, y=ROC_pAUC, fill=metric_descr)) +
  geom_col(position="dodge2") +
  scale_fill_brewer(type="qual", palette = "Set3") +
  scale_x_discrete(labels = 1:9) +
  theme_classic() +
  xlab("") +
  ylab("Area under ROC")


tib %>%
  group_by(metric) %>%
  filter(1-FDR >= min(1-tib$FDR[tib$TPR < 1 & tib$metric != "random"]),
         1-FDR <= max(1-tib$FDR[tib$TPR > 0 & tib$metric != "random"])) %>%
  summarize(PR_pAUC = DescTools::AUC(1-FDR, TPR)) %>%
  add_row(metric = "random", PR_pAUC = 0) %>%   # since all points filtered out
  mutate(metric_descr = metrics_dictionary[metric, "full"]) %>%
  ggplot(aes(x=metric_descr, y=PR_pAUC, fill=metric_descr)) +
  geom_col(position="dodge2") +
  scale_fill_brewer(type="qual", palette = "Set3") +
  scale_x_discrete(labels = 1:9) +
  theme_classic() +
  xlab("") +
  ylab("Area under PR")


# Full AUC
tib %>%
  mutate(metric_descr = metrics_dictionary[metric, "full"]) %>%
  group_by(metric_descr) %>%
  summarize(ROC_AUC = DescTools::AUC(FPR, TPR),
            PR_AUC = DescTools::AUC(1-FDR,TPR)) %>%
  pivot_longer(-metric_descr) %>%
  mutate(name = fct_inorder(name)) %>%
  ggplot(aes(x=metric_descr, y=value, fill=metric_descr)) +
  facet_wrap(~name) +
  geom_col(position="dodge2") +
  scale_fill_brewer(type="qual", palette = "Set3") +
  scale_x_discrete(labels = 1:9) +
  theme_classic() +
  xlab("") +
  ylab("Area under curve")


########################### End ################################