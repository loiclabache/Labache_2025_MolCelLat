################################################################################
# Written by Loïc Labache, Ph.D.                                               #
# Holmes Lab, Department of Psychiatry - Rutgers University                    #
# August 16, 2024                                                              #
################################################################################


#===============================================================================
# Libraries.....................................................................
#...............................................................................
packages = c("here", "ggplot2", "pals", "viridis", "neurobase", "RNiftyReg",
             "progress", "data.table", "tidyr", "dplyr", "ggdist", "prismatic",
             "ggpp", "scrutiny", "factoextra", "PCAtools")
lapply(packages, require, character.only = T)


#===============================================================================
# Data..........................................................................
#...............................................................................
path_folder = "project_path"
mode_one = read.csv(here(path_folder, "Data",
                         "modeOne_CCA_163ROIs_AICHA.csv"))[, -1]
cog_atlas = read.csv(here(path_folder, "Atlas",
                          "Cognitive_Atlases.txt"))
cog_atlas = cog_atlas[order(cog_atlas$Index), ]
subtype_df = read.csv(here(path_folder, "Data", 
                           "24_CellularSubtypes_304ROIs_AICHA_zscored.csv"))[, -1]
spin = as.data.frame(fread(here(path_folder, "Atlas", "AICHA",
                                "permutation_23k_AICHA_Seed1942.csv"),
                           header = FALSE)) + 1


#===============================================================================
# Genetic Maps..................................................................
# Computing Asymmetries.........................................................
r_choice = 1 # 1 = Pearson | 2 = Spearman
mode_choice = 2 # 1 = Task | 2 = Neuro
r_met = c("pearson", "spearman")
# Finding missing couterpart region:
right_roi = subtype_df$region_idx[subtype_df$region_idx %% 2 == 0]
left_roi = subtype_df$region_idx[subtype_df$region_idx %% 2 == 1]
missing_left = right_roi[!(right_roi - 1) %in% left_roi]
missing_right = left_roi[!(left_roi + 1) %in% right_roi]
right_roi = right_roi[!(right_roi %in% missing_left)]
left_roi = left_roi[!(left_roi %in% missing_right)]
subtype_asym = subtype_df[subtype_df$region %in% left_roi, ] - subtype_df[subtype_df$region  %in% right_roi, ]
subtype_asym$region_idx = left_roi
common_roi = intersect(mode_one$regions, subtype_asym$region_idx)
subtype_asym = subtype_asym[subtype_asym$region %in% common_roi, ]
# Spin Test.....................................................................
cog_atlas_tmp = cog_atlas[cog_atlas$Index %in% c(mode_one$regions,
                                                 mode_one$regions + 1), ]
rownames(cog_atlas_tmp) = 1:nrow(cog_atlas_tmp)
cog_atlas_tmp$lane_nb = as.numeric(rownames(cog_atlas_tmp))
permuted_subtype = subtype_df[subtype_df$region_idx %in% c(subtype_asym$region_idx,
                                                           subtype_asym$region_idx + 1), ]
sel_roi = permuted_subtype$region_idx[permuted_subtype$region_idx %% 2 == 1]
mode_one_subtype = mode_one[mode_one$regions %in% sel_roi, ]
cog_atlas_tmp = cog_atlas_tmp[cog_atlas_tmp$Index %in% permuted_subtype$region_idx, ]
rownames(cog_atlas_tmp) = 1:nrow(cog_atlas_tmp)
#...............................................................................
# Using Stepwise Multilinear Model on Subtypes Data.............................
combined_df = merge(subtype_asym, mode_one_subtype[, -1],
                    by.x = "region_idx", by.y = "regions")
combined_df$region_idx = NULL
initial_model = lm(neuro_mode ~ 1, data = combined_df)
full_model = lm(as.formula(paste("neuro_mode ~",
                                 paste(colnames(combined_df[,-ncol(combined_df)]),
                                       collapse = " + "))),
                data = combined_df)
original_model = step(initial_model,
                      scope = list(lower = initial_model, upper = full_model),
                      direction = "both",
                      trace = FALSE)
summary(original_model)
selected_var = names(coef(original_model))[-1]
# Permute on t-value............................................................
original_t_values = summary(original_model)$coefficients[, "t value"][-1]
num_permutations = ncol(spin)
perm_results = as.data.frame(matrix(NA,
                                    nrow = num_permutations,
                                    ncol = length(selected_var)))
colnames(perm_results) = selected_var
pb = progress_bar$new(format = "Progress: [:bar] :percent eta: :eta",
                      total = length(original_t_values)*num_permutations)
asym_roi = permuted_subtype$region_idx[permuted_subtype$region_idx %% 2 == 1]
formula_str = paste("neuro_mode ~", paste(selected_var, collapse = " + "))
formula = as.formula(formula_str)
for (j in 1:length(original_t_values)) {
  for (i in 1:num_permutations) {
    pb$tick()
    
    spin_p = spin[, i]
    spin_p = spin_p[spin_p %in% cog_atlas_tmp$lane_nb]
    spin_roi = cog_atlas_tmp[match(spin_p, cog_atlas_tmp$lane_nb), ]$Index
    
    permut = permuted_subtype
    permut = permut[match(spin_roi, permut$region_idx), ]
    permut = permut[seq(1, nrow(permut), 2), ] - permut[seq(2, nrow(permut), 2), ]
    permut$region_idx = asym_roi
    
    combined_perm = merge(permut, mode_one_subtype[, -1],
                          by.x = "region_idx", by.y = "regions")
    
    permuted_model = lm(formula, data = combined_perm[, -1])
    perm_results[i, j] = summary(permuted_model)$coefficients[j + 1, "t value"]
  }
}


#===============================================================================
# Visualization Spin Cell Subtypes..............................................
#...............................................................................
transformed_perm = perm_results %>%
  mutate(Row_Number = row_number()) %>%
  pivot_longer(cols = c(L5.IT, Micro.PVM, L2.3.IT, L6.IT, L6.IT.Car3), 
               names_to = "Cell", 
               values_to = "Value")
transformed_perm$Cell = factor(transformed_perm$Cell, 
                               levels = c("L5.IT", "Micro.PVM", "L6.IT", "L2.3.IT", "L6.IT.Car3"))
# Fig. 4-c.
fig_spin_cell = ggplot(transformed_perm, 
                       aes(x = Cell,
                           y = Value)) + 
  geom_hline(yintercept = 0,
             linetype = 3,
             linewidth = 0.25,
             color = "#000000") +
  ggdist::stat_halfeye(adjust = 0.5,
                       width = 0.6,
                       .width = 0,
                       point_colour = NA,
                       aes(colour = Cell,
                           fill = after_scale(clr_lighten(colour, space="combined")))) + 
  geom_point(shape = 21,
             size = 0.1,
             stroke = 0.25,
             aes(colour=Cell, 
                 fill=after_scale(clr_alpha(colour, 0.4))),
             position = position_jitternudge(width = 0.1,
                                             x = -0.2,
                                             seed = 1,
                                             nudge.from = "jittered")) +
  geom_boxplot(width = .25, 
               size = 0.5,
               alpha = 0.5,
               outlier.shape = NA,
               aes(colour = Cell,
                   fill = after_scale(clr_lighten(colour, 
                                                  space = "combined",
                                                  shift = 0.755))),
               position = position_nudge(x = -0.2)) +
  geom_linerange(aes(xmin = 0.8 - 0.125, xmax = 0.8 + 0.125, 
                     y = -3.27),
                 color = "#1e4565",
                 linewidth = 0.5) +
  geom_linerange(aes(xmin = 1.8 - 0.125, xmax = 1.8 + 0.125, 
                     y = -2.36),
                 color = "#8c4600",
                 linewidth = 0.5) +
  geom_linerange(aes(xmin = 2.8 - 0.125, xmax = 2.8 + 0.125, 
                     y = -2.30),
                 color = "#542b5a",
                 linewidth = 0.5) +
  geom_linerange(aes(xmin = 3.8 - 0.125, xmax = 3.8 + 0.125, 
                     y = -2.17),
                 color = "#2a6029",
                 linewidth = 0.5) +
  geom_linerange(aes(xmin = 4.8 - 0.125, xmax = 4.8 + 0.125, 
                     y = -1.39),
                 color = "#8d094f",
                 linewidth = 0.5) +
  scale_color_manual(values = c("L5.IT" = "#377eb8", "Micro.PVM" = "#ff7f00",
                                "L6.IT" = "#984ea3", "L2.3.IT" = "#4daf4a",
                                "L6.IT.Car3" = "#F781BF")) + 
  annotate("text", label = c("**", "*", "*", "*"), 
           x = c(0.8, 1.8, 2.8, 3.8), 
           y = rep(4.345, 4)) + 
  ylim(-5, 5) +
  labs(y = "Correlation") +
  theme_classic() + 
  theme(strip.text = element_text(size = 16),
        strip.background = element_blank(),
        text = element_text(family = "Arial"), 
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.title.y = element_text(size = 10,colour = "black"),
        axis.text.x = element_text(size = 9, colour = "black"),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        legend.position = "none")
fig_spin_cell


#===============================================================================
# Visualization Cell Subtypes Relationship......................................
#...............................................................................
selected_df = combined_df[, c(25, 20, 13, 19, 16, 17)]
res_df = data.frame(neuro_mode = selected_df$neuro_mode)
res_df$res_L5IT = resid(lm(L5.IT ~ Micro.PVM + L6.IT + L2.3.IT + L6.IT.Car3, data = selected_df))
res_df$res_MicroPVM = resid(lm(Micro.PVM ~ L5.IT + L6.IT + L2.3.IT + L6.IT.Car3, data = selected_df))
res_df$res_L6IT = resid(lm(L6.IT ~ Micro.PVM + L5.IT + L2.3.IT + L6.IT.Car3, data = selected_df))
res_df$res_L23IT = resid(lm(L2.3.IT ~ Micro.PVM + L6.IT + L5.IT + L6.IT.Car3, data = selected_df))
res_df$res_L6ITCar3 = resid(lm(L6.IT.Car3 ~ Micro.PVM + L6.IT + L2.3.IT + L5.IT, data = selected_df))
cor(res_df$neuro_mode, res_df$res_L5IT)
cor(res_df$neuro_mode, res_df$res_MicroPVM)
cor(res_df$neuro_mode, res_df$res_L6IT)
cor(res_df$neuro_mode, res_df$res_L23IT)
cor(res_df$neuro_mode, res_df$res_L6ITCar3)
# Plot..........................................................................
# Fig. 4-b.
long_res_df = res_df %>%
  pivot_longer(cols = c(res_L5IT, res_MicroPVM, res_L6IT, res_L23IT, 
                        res_L6ITCar3),
               names_to = "Variable", values_to = "Value")
long_res_df$Variable = factor(long_res_df$Variable, 
                              levels = c("res_L5IT", "res_MicroPVM", "res_L6IT",
                                         "res_L23IT", "res_L6ITCar3"))
res_plot = ggplot(long_res_df, aes(x = neuro_mode , y = Value)) +
  geom_vline(xintercept = 0, color = "darkgrey", linetype = "solid") +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = "solid") +
  geom_point(aes(fill = Variable, color = Variable), pch=21, size = 0.5) +
  geom_smooth(method = 'lm', color = "black", se = F, linewidth = 0.5) +
  facet_wrap(~Variable, scales = "fixed", ncol = 5) +
  scale_fill_manual(values =   c("res_L5IT" = "#377eb8", "res_MicroPVM" = "#ff7f00",
                                 "res_L6IT" = "#984ea3", "res_L23IT" = "#4daf4a",
                                 "res_L6ITCar3" = "#F781BF")) +
  scale_color_manual(values = c("res_L5IT" = "#377eb800", "res_MicroPVM" = "#ff7f0000",
                                "res_L6IT" = "#984ea300", "res_L23IT" = "#4daf4a00",
                                "res_L6ITCar3" = "#F781BF00")) +
  theme_classic() + 
  theme(strip.text = element_blank(),
        strip.background = element_blank(),
        text = element_text(family = "Arial"), 
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.title.y = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(size = 9, colour = "black"),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        legend.position = "none")
res_plot

