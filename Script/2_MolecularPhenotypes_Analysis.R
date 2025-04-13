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
synapticDensity = read.csv(here(path_folder, "Data",
                                "SynapticDensityMap_384ROIs_AICHA_zscored.csv"))[, -1]
spin = as.data.frame(fread(here(path_folder, "Atlas", "AICHA",
                                "permutation_23k_AICHA_Seed1942.csv"),
                           header = FALSE))


#===============================================================================
# Synaptic Density Map..........................................................
#...............................................................................
# Computing Asymmetries.........................................................
r_choice = 1 # 1 = Pearson | 2 = Spearman
mode_choice = 2 # 1 = Task | 2 = Neuro
r_met = c("pearson", "spearman")
left_roi = synapticDensity$region[synapticDensity$region %% 2 == 1]
right_roi = synapticDensity$region[synapticDensity$region %% 2 == 0]
synapticDensity_asym = setNames(as.data.frame(matrix(nrow = length(left_roi),
                                                     ncol = ncol(synapticDensity))),
                                colnames(synapticDensity))
synapticDensity_asym$region = left_roi
synapticDensity_asym[, 2] = synapticDensity[synapticDensity$region %in% left_roi, 2] - synapticDensity[synapticDensity$region %in% right_roi, 2]
synapticDensity_asym = synapticDensity_asym[synapticDensity_asym$region %in% mode_one$regions, ]
# Correlation...................................................................
cor_synDens = cor(mode_one[, mode_choice],
                  synapticDensity_asym[, 2], 
                  method = r_met[r_choice])
df = data.frame(y = mode_one[, mode_choice],
                x = synapticDensity_asym[, 2])
plot(y ~ x, df)
abline(lm(y ~ x, df))
# Spin Test.....................................................................
selected_roi = c(mode_one$regions, mode_one$regions + 1)
permuted_synapticDensity = synapticDensity[synapticDensity$region %in% selected_roi, ]
result_permut = as.data.frame(matrix(data = NA, 
                                     nrow = ncol(spin),
                                     ncol = 1))
colnames(result_permut) = "Spin_Correlation"
pb = progress_bar$new(format = "Progress: [:bar] :percent eta: :eta",
                      total = 10000)
for (p in 1:ncol(spin)){
  pb$tick()
  permut = permuted_synapticDensity
  permut[, 2] = permut[spin[, p] + 1, 2]
  rownames(permut) = permut$region
  permut = permut[seq(1, nrow(permut), 2), ] - permut[seq(2, nrow(permut), 2), ]
  result_permut[p, 1] = cor(permut$Synaptic_Density,
                            mode_one[, mode_choice], 
                            method = r_met[r_choice])
}


#===============================================================================
# Mitochondrial Density Map.....................................................
#...............................................................................
# Computing Asymmetries.........................................................
choice = 1 # 1 = TRC | 2 = MitoD | 3 = MRC
r_choice = 1 # 1 = Pearson | 2 = Spearman
mode_choice = 2 # 1 = Task | 2 = Neuro
maps = c("TRC_384ROIs_AICHA_zscored.csv", 
         "MitoDMap_384ROIs_AICHA_zscored.csv",
         "MRC_384ROIs_AICHA_zscored.csv")
r_met = c("pearson", "spearman")
mitoDens = read.csv(here(path_folder, "Data", maps[choice]))[, -1]
left_roi = mitoDens$region[mitoDens$region %% 2 == 1]
right_roi = mitoDens$region[mitoDens$region %% 2 == 0]
mitoDens_asym = setNames(as.data.frame(matrix(nrow = length(left_roi),
                                              ncol = ncol(mitoDens))),
                         colnames(mitoDens))
mitoDens_asym$region = left_roi
mitoDens_asym[, 2] = mitoDens[mitoDens$region %in% left_roi, 2] - mitoDens[mitoDens$region %in% right_roi, 2]
mitoDens_asym = mitoDens_asym[mitoDens_asym$region %in% mode_one$regions, ]
# Correlation...................................................................
perm_res = data.frame(syn = result_permut$Spin_Correlation,
                      trc = rep(NA, 10000),
                      mitod = rep(NA, 10000),
                      mrc = rep(NA, 10000))
pb = progress_bar$new(format = "Progress: [:bar] :percent eta: :eta",
                      total = 30000)
for (m in 1:length(maps)){
  cor_mitoDens = cor(mode_one[, mode_choice],
                     mitoDens_asym[, 2],
                     method = r_met[r_choice])
  df = data.frame(y = mode_one[, mode_choice],
                  x = mitoDens_asym[, 2])
  plot(y ~ x, df)
  abline(lm(y ~ x, df))
  # Spin Test...................................................................
  selected_roi = c(mode_one$regions, mode_one$regions + 1)
  permuted_mitoDens = mitoDens[mitoDens$region %in% selected_roi, ]
  result_permut = as.data.frame(matrix(data = NA, 
                                       nrow = ncol(spin),
                                       ncol = 1))
  colnames(result_permut) = "Spin_Correlation"
  for (p in 1:ncol(spin)){
    pb$tick()
    permut = permuted_mitoDens
    permut[, 2] = permut[spin[, p] + 1, 2]
    rownames(permut) = permut$region
    permut = permut[seq(1, nrow(permut), 2), ] - permut[seq(2, nrow(permut), 2), ]
    result_permut[p, 1] = cor(permut$Density,
                              mode_one[, mode_choice],
                              method = r_met[r_choice])
  }
  perm_res[, m + 1] = result_permut$Spin_Correlation
}


#===============================================================================
# Visualization Spin Molecular Phenotypes.......................................
#...............................................................................
transformed_perm = perm_res %>%
  mutate(Row_Number = row_number()) %>%
  pivot_longer(cols = c(syn, trc, mitod, mrc), 
               names_to = "Cell", 
               values_to = "Value")
transformed_perm$Cell = factor(transformed_perm$Cell, 
                               levels = c("syn", "mrc", "mitod", "trc"))
# Fig. 3-c.
fig_spin = ggplot(transformed_perm, 
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
                     y = -0.0022),
                 color = "#7d0e0f",
                 linewidth = 0.5) +
  geom_linerange(aes(xmin = 1.8 - 0.125, xmax = 1.8 + 0.125, 
                     y = 0.2726),
                 color = "#7b7b00",
                 linewidth = 0.5) +
  geom_linerange(aes(xmin = 2.8 - 0.125, xmax = 2.8 + 0.125, 
                     y = 0.3389),
                 color = "#5b2f16",
                 linewidth = 0.5) +
  geom_linerange(aes(xmin = 3.8 - 0.125, xmax = 3.8 + 0.125, 
                     y = 0.3598),
                 color = "#545454",
                 linewidth = 0.5) +
  scale_color_manual(values = c("syn" = "#e41a1c", "mrc" = "#e0e000",
                                "mitod" = "#a65628", "trc" = "#999999")) +
  annotate("text", label = c("***", "***", "***"), 
           x = c(1.8, 2.8, 3.8), 
           y = rep(0.4, 3)) + 
  ylim(-0.4, 0.4) +
  labs(y = "Correlation") +
  theme_classic() 
fig_spin


#===============================================================================
# Visualization Molecular Relationship..........................................
#...............................................................................
# Synaptic Density..............................................................
left_roi = synapticDensity$region[synapticDensity$region %% 2 == 1]
right_roi = synapticDensity$region[synapticDensity$region %% 2 == 0]
synapticDensity_asym = setNames(as.data.frame(matrix(nrow = length(left_roi),
                                                     ncol = ncol(synapticDensity))),
                                colnames(synapticDensity))
synapticDensity_asym$region = left_roi
synapticDensity_asym[, 2] = synapticDensity[synapticDensity$region %in% left_roi, 2] - synapticDensity[synapticDensity$region %in% right_roi, 2]
synapticDensity_asym = synapticDensity_asym[synapticDensity_asym$region %in% mode_one$regions, ]
# Mitochondrial Phenotypes......................................................
maps = c("TRC_384ROIs_AICHA_zscored.csv", 
         "MitoDMap_384ROIs_AICHA_zscored.csv",
         "MRC_384ROIs_AICHA_zscored.csv")
mitoDens = data.frame(region = read.csv(here(path_folder, "Data", maps[1]))[, 2],
                      TRC = read.csv(here(path_folder, "Data", maps[1]))[, -c(1:2)],
                      MitoDMap = read.csv(here(path_folder, "Data", maps[2]))[, -c(1:2)],
                      MRC = read.csv(here(path_folder, "Data", maps[3]))[, -(1:2)])
left_roi = mitoDens$region[mitoDens$region %% 2 == 1]
right_roi = mitoDens$region[mitoDens$region %% 2 == 0]
mitoDens_asym = mitoDens[mitoDens$region %in% left_roi, ] - mitoDens[mitoDens$region %in% right_roi, ]
mitoDens_asym$region = left_roi
mitoDens_asym = mitoDens_asym[mitoDens_asym$region %in% mode_one$regions, ]
# Plot correlation between each map and the neuro mode..........................
to_plot_df = mitoDens_asym
to_plot_df$SynD = synapticDensity_asym$Synaptic_Density
to_plot_df$ModeN = mode_one$neuro_mode
long_df = to_plot_df %>%
  pivot_longer(cols = c(TRC, MitoDMap, MRC, SynD),
               names_to = "Variable", values_to = "Value")
long_df$Variable = factor(long_df$Variable, 
                          levels = c("SynD", "MRC", "MitoDMap", "TRC"))
# Fig. 3-b.
cor_plot = ggplot(long_df, aes(x = ModeN , y = Value)) +
  geom_vline(xintercept = 0, color = "darkgrey", linetype = "solid") +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = "solid") +
  geom_point(aes(fill = Variable, color = Variable), pch=21, size = 0.5) +
  geom_smooth(method = 'lm', color = "black", se = F, linewidth = 0.5) +
  facet_wrap(~Variable, scales = "fixed", ncol = 4) +
  scale_fill_manual(values = c("SynD" = "#e41a1c", "MRC" = "#e0e000",
                               "MitoDMap" = "#a65628", "TRC" = "#999999")) +
  scale_color_manual(values = c("SynD" = "#e41a1c00", "MRC" = "#4daf4a00",
                                "MitoDMap" = "#a6562800", "TRC" = "#99999900")) +
  theme_classic()
cor_plot

