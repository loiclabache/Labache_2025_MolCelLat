################################################################################
# Written by Loïc Labache, Ph.D.                                               #
# Holmes Lab, Department of Psychiatry - Rutgers University                    #
# March 17, 2025                                                               #
################################################################################


# This script performs the CCA analysis described in the paper (DOI: XXXX/X.X.X)
# and generates the corresponding figures.


# Libraries.....................................................................
#...............................................................................
packages = c("data.table", "progress", "here", "factoextra", "PCAtools",
             "dplyr", "tidyr", "scrutiny", "acca", "corrplot", "pals",
             "reshape2", "viridis", "ggdensity", "ggnewscale", "forcats",
             "ggplot2", "ggdist", "prismatic", "ggpp")
lapply(packages, require, character.only = T)
set.seed(13)


# Data for CCA..................................................................
#...............................................................................
path_folder = "project_path"
neuro = read.csv(here(path_folder, "Data", 
                      "19neurotransmitters_163ROIs_AICHA_asym_zScored.csv"))[, -1]
asym = read.csv(here(path_folder, "Data", 
                     "4tasks_typicalStrong_185ROIs_AICHA_asym.csv"))
# Common Brain Regions: no sub-cortical and full signal.........................
common_roi = intersect(neuro$region, asym$region)
neuro = neuro[neuro$region %in% common_roi, ]
asym = asym[asym$region %in% common_roi, ]


# Principal Component Analysis on fMRI Data.....................................
#...............................................................................
# Reduce number of fMRI variables...............................................
pca_fMRI = prcomp(asym[, -c(1, 6:10)],
                  center = F, scale. = F)
# Select PCs which explain X amount of variance ................................
selectedPCs_fMRI = findElbowPoint(get_eig(pca_fMRI)$variance.percent)
print(get_eig(pca_fMRI))
fMRI_PCs = pca_fMRI$x[,1:selectedPCs_fMRI]
# Inspect loading's to see which scales load onto the selected PCs..............
loadings_fMRI = as.data.frame(pca_fMRI$rotation[,1:selectedPCs_fMRI] )
loadings_fMRI$task = rownames(loadings_fMRI)
# Plot the loadings to interpret the main components............................
lim_fMRI_ax = round_up(max(c(abs(max(loadings_fMRI[, 1:selectedPCs_fMRI])),
                             abs(min(loadings_fMRI[, 1:selectedPCs_fMRI])))))
loadings_fMRI_long = loadings_fMRI %>%
  pivot_longer(!task, names_to = "PCs", values_to = "loadings") %>% 
  group_by(PCs) %>% arrange(PCs, loadings)
loadings_fMRI_long_order = loadings_fMRI_long %>%
  do(tibble(al = levels(reorder(interaction(.$PCs, .$task, drop=TRUE), 
                                .$loadings)))) %>% 
  pull(al)
loadings_fMRI_long = loadings_fMRI_long %>%
  mutate(al = factor(interaction(PCs, task),
                     levels = loadings_fMRI_long_order))
# Extended Data Fig. 1-a:b.
ggplot(loadings_fMRI_long,
       aes(y = al, x = loadings)) +
  geom_vline(xintercept = 0, color = "darkgrey", linetype = "solid") +
  scale_x_continuous(limits = c(-lim_fMRI_ax, lim_fMRI_ax)) +
  geom_col(show.legend = FALSE,
           aes(fill = loadings)) +
  scale_fill_gradient2(low = coolwarm(100)[1],
                       mid = coolwarm(100)[50],
                       high = coolwarm(100)[100]) +
  facet_wrap( ~ PCs, scales = "free_y") +
  theme_classic()


# Principal Component Analysis on Neurotransmitters Data........................
#...............................................................................
# Reduce number of Neurotransmitters variables..................................
pca_neurotrans = prcomp(neuro[, -c(1)],
                        center = F, scale. = F)
# Select PCs which explain X amount of variance ................................
selectedPCs_neuro = findElbowPoint(get_eig(pca_neurotrans)$variance.percent)
print(get_eig(pca_neurotrans))
neurotrans_PCs = pca_neurotrans$x[,1:selectedPCs_neuro]
# Inspect loading's to see which scales load onto the selected PCs..............
loadings_neurotrans = as.data.frame(pca_neurotrans$rotation[,1:selectedPCs_neuro] )
loadings_neurotrans$neurotransmitters = rownames(loadings_neurotrans)
# Plot the loadings to interpret the main components............................
lim_neurotrans_ax = round_up(max(c(abs(max(loadings_neurotrans[,-ncol(loadings_neurotrans)])),
                                   abs(min(loadings_neurotrans[,-ncol(loadings_neurotrans)])))))
loadings_neurotrans_long = loadings_neurotrans %>%
  pivot_longer(!neurotransmitters, names_to = "PCs", values_to = "loadings") %>% 
  group_by(PCs) %>% arrange(PCs, loadings)
loadings_neurotrans_long_order = loadings_neurotrans_long %>%
  do(tibble(al = levels(reorder(interaction(.$PCs, .$neurotransmitters, drop=TRUE), 
                                .$loadings)))) %>% 
  pull(al)
loadings_neurotrans_long = loadings_neurotrans_long %>%
  mutate(al = factor(interaction(PCs, neurotransmitters),
                     levels = loadings_neurotrans_long_order))
# Extended Data Fig. 1-a.
ggplot(loadings_neurotrans_long,
       aes(y = al, x = loadings)) +
  geom_vline(xintercept = 0, color = "darkgrey", linetype = "solid") +
  scale_x_continuous(limits = c(-lim_neurotrans_ax, lim_neurotrans_ax)) +
  geom_col(show.legend = FALSE,
           aes(fill = loadings)) +
  scale_fill_gradient2(low = coolwarm(100)[1],
                       mid = coolwarm(100)[50],
                       high = coolwarm(100)[100]) +
  facet_wrap( ~ PCs, scales = "free_y") +
  scale_y_discrete(breaks = loadings_neurotrans_long$al,
                   labels = sub("^PC[0-7]+\\.", "",
                                loadings_neurotrans_long$al)) +
  theme_classic()


# Canonical Correlation Analysis................................................
#...............................................................................
mod = cc(X  = fMRI_PCs,
         Y  = neurotrans_PCs)
mod$cor #print the canonical correlation coef between each CCA mode
mod$prop_expl_var # Variance explained by each CCA mode
# Visualization CCA's modes.....................................................
#...............................................................................
loadings_mode_fMRI = as.data.frame(mod$corr$corr.X.xscores)
loadings_mode_fMRI$Comp = paste0("fMRI_",
                                 rownames(loadings_mode_fMRI))
loadings_mode_fMRI$Var = rep("fMRI", dim(loadings_mode_fMRI)[1])
loadings_mode_neurotrans = as.data.frame(mod$corr$corr.Y.yscores)
loadings_mode_neurotrans$Comp = paste0("neurotransmitters_", 
                                       rownames(loadings_mode_neurotrans))
loadings_mode_neurotrans$Var = rep("neurotransmitters",
                                   dim(loadings_mode_neurotrans)[1])
colnames(loadings_mode_fMRI)[1:2] = c("ModeOne", "ModeTwo")
colnames(loadings_mode_neurotrans)[1:2] = c("ModeOne", "ModeTwo")
# Bar Plot......................................................................
biplot_CCA_data = rbind(loadings_mode_fMRI, loadings_mode_neurotrans)
rownames(biplot_CCA_data) = c(1:dim(biplot_CCA_data)[1])
biplot_CCA_data = biplot_CCA_data[order(biplot_CCA_data$ModeOne), ]
biplot_CCA_data$Comp = factor(biplot_CCA_data$Comp,
                              levels = biplot_CCA_data$Comp,
                              ordered = TRUE)
nb_col = ncol(biplot_CCA_data)
lim_cca_ax = round_up(max(c(abs(max(biplot_CCA_data[, -c(nb_col-1, nb_col)])),
                            abs(min(biplot_CCA_data[, -c(nb_col-1, nb_col)])))))
# Extended Data Fig. 2-a:b.
ggplot(biplot_CCA_data,
       aes(y = Comp, x = ModeOne)) +
  geom_vline(xintercept = 0, color = "darkgrey", linetype = "solid") +
  scale_x_continuous(limits = c(-lim_cca_ax, lim_cca_ax)) +
  geom_col(show.legend = FALSE,
           aes(fill = ModeOne)) +
  scale_fill_gradient2(low = coolwarm(100)[1],
                       mid = coolwarm(100)[50],
                       high = coolwarm(100)[100]) +
  facet_wrap( ~ Var, scales = "free_y") +
  theme_classic() 
# Express CCA modes according to Original Variables.............................
fMRI_PCA_loadings = pca_fMRI$rotation[, 1:selectedPCs_fMRI]
neurotrans_PCA_loadings = pca_neurotrans$rotation[, 1:selectedPCs_neuro]
fMRI_CCA_loadings = mod$corr$corr.X.xscores
neurotrans_CCA_loadings = mod$corr$corr.Y.yscores
cor_asym = as.data.frame(fMRI_PCA_loadings %*% fMRI_CCA_loadings)
cor_asym$var = rep("fMRI", nrow(cor_asym))
cor_neurotrans = as.data.frame(neurotrans_PCA_loadings %*% neurotrans_CCA_loadings)
cor_neurotrans$var = rep("Neurotrans", nrow(cor_neurotrans))
colnames(cor_asym)[1:2] = c("ModeOne", "ModeTwo")
colnames(cor_neurotrans)[1:2] = c("ModeOne", "ModeTwo")
cor_var = rbind(cor_asym, cor_neurotrans)
cor_var$subVar = rownames(cor_var)
# Ordered variables to plot:
cor_var = cor_var[order(cor_var$ModeOne), ]
cor_var = cor_var[order(cor_var$var), ]
cor_var$subVar = factor(cor_var$subVar, 
                        levels = cor_var$subVar,
                        ordered = TRUE)
# For First Mode................................................................
# Plot correlations to initial data for first mode..............................
# For task Mode:
init_data_taskMode = melt(data.frame(prod = asym$sentMword_prod, 
                                     read = asym$sentMword_read,
                                     list = asym$sentMword_list,
                                     attention = asym$LBJ, 
                                     taskMode = mod$scores$xscores[, 1]),
                          id.vars = "taskMode")
init_cor_taskMode = init_data_taskMode %>%
  group_by(variable) %>%
  summarize(cor = round(cor(value, taskMode), 2))
init_cor_taskMode = init_cor_taskMode[order(init_cor_taskMode$cor), ]
init_cor_taskMode$variable = factor(init_cor_taskMode$variable, 
                                    levels = init_cor_taskMode$variable,
                                    ordered = TRUE)
lim_cor_ax = 1
# Fig. 1-b.
ggplot(init_cor_taskMode,
       aes(y = variable, x = cor)) +
  geom_vline(xintercept = 0, color = "darkgrey", linetype = "solid") +
  scale_x_continuous(limits = c(-lim_cor_ax, lim_cor_ax)) +
  geom_col(show.legend = FALSE,
           aes(fill = cor)) +
  scale_fill_gradient2(low = coolwarm(100)[1],
                       mid = coolwarm(100)[50],
                       high = coolwarm(100)[100]) +
  theme_classic()
# For neuro Mode:
init_data_neuroMode = melt(data.frame(neuro[, -1],
                                     neuroMode = mod$scores$yscores[, 1]),
                          id.vars = "neuroMode")
init_cor_neuroMode = init_data_neuroMode %>%
  group_by(variable) %>%
  summarize(cor = round(cor(value, neuroMode), 2)) %>%
  arrange(desc(cor)) %>%
  mutate(variable = fct_reorder(variable, cor))
lim_cor_ax = 1
# Fig. 1-a.
ggplot(init_cor_neuroMode,
       aes(y = variable, x = cor)) +
  geom_vline(xintercept = 0, color = "darkgrey", linetype = "solid") +
  scale_x_continuous(limits = c(-lim_cor_ax, lim_cor_ax)) +
  geom_col(show.legend = FALSE,
           aes(fill = cor)) +
  scale_fill_gradient2(low = coolwarm(100)[1],
                       mid = coolwarm(100)[50],
                       high = coolwarm(100)[100]) +
  theme_classic() 
# For Second Mode...............................................................
# Ordered variables to plot:
cor_var = cor_var[order(cor_var$ModeTwo), ]
cor_var = cor_var[order(cor_var$var), ]
cor_var$subVar = factor(cor_var$subVar, 
                        levels = cor_var$subVar,
                        ordered = TRUE)
# Extended Data Fig. 3-a.
ggplot(cor_var[cor_var$var == "Neurotrans", ],
       aes(y = subVar, x = ModeTwo)) +
  geom_vline(xintercept = 0, color = "darkgrey", linetype = "solid") +
  scale_x_continuous(limits = c(-lim_cor_ax, lim_cor_ax)) +
  geom_col(show.legend = FALSE,
           aes(fill = ModeTwo)) +
  scale_fill_gradient2(low = coolwarm(100)[1],
                       mid = coolwarm(100)[50],
                       high = coolwarm(100)[100]) +
  theme_classic() 
# Extended Data Fig. 3-b.
ggplot(cor_var[cor_var$var == "fMRI", ],
       aes(y = subVar, x = ModeTwo)) +
  geom_vline(xintercept = 0, color = "darkgrey", linetype = "solid") +
  scale_x_continuous(limits = c(-lim_cor_ax, lim_cor_ax)) +
  geom_col(show.legend = FALSE,
           aes(fill = ModeTwo)) +
  scale_fill_gradient2(low = coolwarm(100)[1],
                       mid = coolwarm(100)[50],
                       high = coolwarm(100)[100]) +
  theme_classic()
#...............................................................................
# Correlation between First Modes...............................................
main_shape_data = data.frame(fMRI = mod$scores$xscores[, 1],
                             Neurotransmitters = mod$scores$yscores[, 1], 
                             regions = asym$region)
# Plot the correlation between First Modes:
# Fig. 2-c.
ggplot(main_shape_data, aes(fMRI, Neurotransmitters)) +
  geom_vline(xintercept = 0, color = "darkgrey", linetype = "solid") +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = "solid") +
  geom_point(aes(fill = "black"), pch=21, size = 1) +
  scale_fill_identity() +
  geom_smooth(method='lm', color = "black", se = F, linewidth = 0.5) +
  theme_classic()
# Correlation between Second Modes..............................................
main_shape_data_modeTwo = data.frame(fMRI = mod$scores$xscores[, 2],
                                     Neurotransmitters = mod$scores$yscores[, 2], 
                                     regions = asym$region)
# Plot the correlation between Second Modes:
# Extended Data Fig. 3-c.
ggplot(main_shape_data_modeTwo, aes(fMRI, Neurotransmitters)) +
  geom_vline(xintercept = 0, color = "darkgrey", linetype = "solid") +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = "solid") +
  geom_point(aes(fill = "black"), pch=21, size = 1) +
  scale_fill_identity() +
  geom_smooth(method='lm', color = "black", se = F, linewidth = 0.5) +
  theme_classic() 


# Bootstrapping method to evaluate contribution of each variables to CCA modes..
#...............................................................................
set.seed(13)
B = 10000
n = nrow(fMRI_PCs)
task_df = data.frame(prod = asym$sentMword_prod,
                     read = asym$sentMword_read,
                     list = asym$sentMword_list,
                     attention = asym$LBJ)
loadings_boot = matrix(NA, nrow = B, ncol = ncol(task_df) + (ncol(neuro)-1))
pb = progress::progress_bar$new(format = "Processing [:bar] :percent in :elapsed",
                                total = B, clear = FALSE, width = 60)
for (b in seq_len(B)) {
  idx = sample.int(n, n, replace = TRUE)
  cca_res = cc(fMRI_PCs[idx, , drop = FALSE], neurotrans_PCs[idx, , drop = FALSE])
  
  # sign alignment to the full-sample direction
  sgn_x = sign(cor(cca_res$scores$xscores[, 1], mod$scores$xscores[idx, 1], use = "pairwise.complete.obs"))
  sgn_y = sign(cor(cca_res$scores$yscores[, 1], mod$scores$yscores[idx, 1], use = "pairwise.complete.obs"))
  sgn_x[is.na(sgn_x)] = 1
  sgn_y[is.na(sgn_y)] = 1
  
  loadX = cor(task_df[idx, ], sgn_x * cca_res$scores$xscores[, 1],
              use = "pairwise.complete.obs")
  loadY = cor(neuro[idx, -1], sgn_y * cca_res$scores$yscores[, 1],
              use = "pairwise.complete.obs")
  loadings_boot[b, ] = c(loadX, loadY)
  pb$tick()
}
SE = apply(loadings_boot, 2, sd, na.rm = TRUE)
names(SE) = c(colnames(task_df), colnames(neuro)[-1])
load_full = c(setNames(init_cor_taskMode$cor,  as.character(init_cor_taskMode$variable)),
              setNames(init_cor_neuroMode$cor, as.character(init_cor_neuroMode$variable)))
load_full = load_full[names(SE)]
z = load_full / SE
p = 2 * pnorm(-abs(z))
p_FDR_neurot = p.adjust(p[5:23], method = "holm")
p_FDR_neurot[which(p_FDR_neurot<0.05)]
p_FDR_task = p.adjust(p[1:4], method = "holm")
p_FDR_task[which(p_FDR_task<0.05)]


# Spin Test on Task Data........................................................
#...............................................................................
activation = read.csv(here(path_folder, "Data",
                           "4tasks_typicalStrong_370ROIs_AICHA.csv"))
activation = activation[activation$region %in% c(common_roi, common_roi + 1),
                        c(2:5)]
rownames(activation) = 1:nrow(activation)
spin_list = fread(here(path_folder, "Atlas/AICHA", "permutation_23k_AICHA_Seed1893.csv"),
                  header = FALSE)
spin_list = lapply(spin_list, function(spin) spin + 1)
result_perm_task = as.data.frame(matrix(data = NA, 
                                        nrow = 10000,
                                        ncol = 2))
colnames(result_perm_task) = c("Mode_One", "Mode_Two")
pb = progress_bar$new(format = "Progress: [:bar] :percent eta: :eta",
                      total = 10000)
for (p in 1:10000){
  pb$tick()
  
  # Permute brain region:
  col_permut = as.data.frame(spin_list)[, p]
  permuted_activation = activation[col_permut, ]
  # Compute permuted asymmetry:
  permuted_activation = permuted_activation[seq(1, nrow(permuted_activation), 2), ] - permuted_activation[seq(2, nrow(permuted_activation), 2), ]
  
  pca_fMRI = prcomp(permuted_activation,
                    center = F, scale. = F)$x[, 1:2]
  pca_neuro = prcomp(neuro[, -1],
                     center = F, scale. = F)$x[, 1:7]
  
  result_perm_task[p, ] = cc(X = pca_fMRI,
                             Y = pca_neuro)$cor
}
# Permutation Tasks:
res_task_One = table(result_perm_task$Mode_One >= mod$cor[1])
res_task_One/10000
res_task_Two = table(result_perm_task$Mode_Two >= mod$cor[2])
res_task_Two/10000
# Visualization of spin test results:
transformed_perm = result_perm_task %>%
  mutate(Row_Number = row_number()) %>%
  pivot_longer(cols = c(Mode_One, Mode_Two), 
               names_to = "Mode", 
               values_to = "Value")
transformed_perm$Mode = factor(transformed_perm$Mode, 
                               levels = c("Mode_One", "Mode_Two"))
# Extended Data Fig. 2-c.
ggplot(transformed_perm, aes(x = Mode, y = Value)) + 
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.25, color = "#000000") +
  ggdist::stat_halfeye(adjust = 0.5, width = 0.6, .width = 0,
                       point_colour = NA, aes(colour = Mode,
                                              fill = after_scale(clr_lighten(colour, space="combined")))) + 
  geom_point(shape = 21, size = 0.1, stroke = 0.25,
             aes(colour=Mode, fill=after_scale(clr_alpha(colour, 0.4))),
             position = position_jitternudge(width = 0.1,  x = -0.2,
                                             seed = 1, nudge.from = "jittered")) +
  geom_boxplot(width = .25, size = 0.5, alpha = 0.5,
               outlier.shape = NA,
               aes(colour = Mode, fill = after_scale(clr_lighten(colour, 
                                                                 space = "combined",
                                                                 shift = 0.755))),
               position = position_nudge(x = -0.2)) +
  geom_linerange(aes(xmin = 0.8 - 0.125, xmax = 0.8 + 0.125, y = 0.3943),
                 color = "#48006a", linewidth = 0.5) +
  geom_linerange(aes(xmin = 1.8 - 0.125, xmax = 1.8 + 0.125, 
                     y = 0.1749), color = "#084081", linewidth = 0.5) +
  scale_color_manual(values = c("Mode_One" = "#c51b8a", "Mode_Two" = "#43a2ca")) +
  annotate("text", label = "**", x = 0.8, y = 0.45) + 
  ylim(0, ceiling(max(transformed_perm$Value) * 100) / 100) +
  labs(y = "Correlation") +
  scale_x_discrete(labels = c("Mode_One" = "First Modes", "Mode_Two" = "Second Modes")) +
  theme_classic()
