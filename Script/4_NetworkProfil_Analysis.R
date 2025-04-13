################################################################################
# Written by Loïc Labache, Ph.D.                                               #
# Holmes Lab, Department of Psychiatry - Rutgers University                    #
# May 13, 2024                                                                 #
################################################################################


# Libraries.....................................................................
#...............................................................................
packages = c("here", "ggplot2", "tidyr", "RColorBrewer", "dplyr",
             "NbClust", "psych", "dendextend", "prismatic", "Barnard")
lapply(packages, require, character.only = T)


# Data..........................................................................
#...............................................................................
path_folder = "project_path"
mitho = read.csv(here(path_folder, "Data",
                      "MitoDMap_384ROIs_AICHA_zscored.csv"))[, -1]
syn = read.csv(here(path_folder, "Data",
                    "SynapticDensityMap_384ROIs_AICHA_zscored.csv"))[, -1]
mrc = read.csv(here(path_folder, "Data",
                    "MRC_384ROIs_AICHA_zscored.csv"))[, -1]
trc = read.csv(here(path_folder, "Data",
                    "TRC_384ROIs_AICHA_zscored.csv"))[, -1]
cellSub = read.csv(here(path_folder, "Data",
                        "24_CellularSubtypes_304ROIs_AICHA_zscored.csv"))[, -1]
# Networks Atlas:
cog_atlas = read.csv(here(path_folder, "Atlas",
                          "Cognitive_Atlases.txt"))
yan_atlas = read.csv(here(path_folder, "Atlas", "AICHA",
                          "comparison_AICHA_Yan2023_400_2mm.txt"))
atlas = merge(cog_atlas[, c(3, 6, 11, 13, 17)], yan_atlas[, c(1, 5, 11)],
              by.x = "Index", by.y = "AICHA_Index")
# To select the right regions:
modeOne = read.csv(here(path_folder, "Data",
                        "modeOne_CCA_163ROIs_AICHA.csv"))[, -1]
# Molecular data:
mol = data.frame(mitoDMap = mitho$Density,
                 syn = syn$Synaptic_Density,
                 mrc = mrc$Density,
                 trc = trc$Density,
                 region = mitho$region)
mol = mol[mol$region %in% 
            c(modeOne$regions, modeOne$regions + 1), ]
# Cellular data:
cell = cellSub[, c(1, 21, 14, 20, 17)]
right_roi = cell$region_idx[cell$region_idx %% 2 == 0]
left_roi = cell$region_idx[cell$region_idx %% 2 == 1]
missing_left = right_roi[!(right_roi - 1) %in% left_roi]
missing_right = left_roi[!(left_roi + 1) %in% right_roi]
right_roi = right_roi[!(right_roi %in% missing_left)]
left_roi = left_roi[!(left_roi %in% missing_right)]
cell = cell[cell$region_idx %in% 
              c(left_roi, right_roi), ]
cell = cell[cell$region_idx %in% 
              c(modeOne$regions, modeOne$regions + 1), ]


# Average Value by Yan Network..................................................
#...............................................................................
# For Molecular Biomarkers:
merged_mol = merge(mol, atlas, by.x = "region", by.y = "Index")
res_mol = merged_mol %>%
  group_by(AICHA_Hemisphere, Network_7_Yan400_Name) %>%
  summarise(wgth_TRC = mean(trc),
            wgth_MitoDMap = mean(mitoDMap),
            wgth_MRC = mean(mrc),
            wgth_SynD = mean(syn))
mol_diff = res_mol %>%
  pivot_wider(names_from = AICHA_Hemisphere,
              values_from = c(wgth_TRC, wgth_MitoDMap, wgth_MRC, wgth_SynD)) %>%
  mutate(diff_wgth_TRC = wgth_TRC_Left - wgth_TRC_Right,
         diff_wgth_MitoDMap = wgth_MitoDMap_Left - wgth_MitoDMap_Right,
         diff_wgth_MRC = wgth_MRC_Left - wgth_MRC_Right,
         diff_wgth_SynD = wgth_SynD_Left - wgth_SynD_Right) %>%
  select(Network_7_Yan400_Name, diff_wgth_TRC, diff_wgth_MitoDMap,
         diff_wgth_MRC, diff_wgth_SynD)
# For Cell Subtypes:
merged_cell = merge(cell, atlas, by.x = "region_idx", by.y = "Index")
res_cell = merged_cell %>%
  group_by(AICHA_Hemisphere, Network_7_Yan400_Name) %>%
  summarise(wgth_L5IT = mean(L5.IT),
            wgth_MicroPVM = mean(Micro.PVM),
            wgth_L6IT = mean(L6.IT),
            wgth_L23IT = mean(L2.3.IT))
cell_diff = res_cell %>%
  pivot_wider(names_from = AICHA_Hemisphere,
              values_from = c(wgth_L5IT, wgth_MicroPVM, wgth_L6IT, wgth_L23IT)) %>%
  mutate(diff_wgth_L5IT = wgth_L5IT_Left - wgth_L5IT_Right,
         diff_wgth_MicroPVM = wgth_MicroPVM_Left - wgth_MicroPVM_Right,
         diff_wgth_L6IT = wgth_L6IT_Left - wgth_L6IT_Right,
         diff_wgth_L23IT = wgth_L23IT_Left - wgth_L23IT_Right) %>%
  select(Network_7_Yan400_Name, diff_wgth_L5IT, diff_wgth_MicroPVM, 
         diff_wgth_L6IT, diff_wgth_L23IT)
# Merge Results: 
merged_res = as.data.frame(mol_diff %>%
                             inner_join(cell_diff, by = "Network_7_Yan400_Name"))
row.names(merged_res) = merged_res$Network_7_Yan400_Name
# Min-Max Normalization:
normalize = function(x) {
  return((2 * (x - min(x)) / (max(x) - min(x))) - 1)
  # return(x / max(abs(x)))
}
normalized_res = merged_res
normalized_res[, -1] = apply(normalized_res[, -1], 2, normalize)
# Visualization Profiles........................................................
#...............................................................................
# Lollipop Plot by Network:
# Fig. 5-a.
normalized_res$diff_wgth_SynD = NULL
colnames(normalized_res) = gsub("diff_", "", colnames(normalized_res))
long_df = normalized_res %>%
  pivot_longer(cols = wgth_TRC:wgth_L23IT,
               names_to = "Variable", values_to = "Value")
long_df$Network_7_Yan400_Name = factor(long_df$Network_7_Yan400_Name,
                                       levels = c("SalVentAttn", "Vis", "Default",
                                                  "Limbic", "DorsAttn", "Cont",
                                                  "SomMot"))
average_values = long_df %>%
  group_by(Variable) %>%
  summarise(Avg_Value = mean(Value)) %>%
  arrange(Avg_Value)
long_df$Variable = factor(long_df$Variable,
                          levels = average_values$Variable)
col_Yeo = c(Vis = "#781286", SomMot = "#4682b4", DorsAttn = "#00760e",
            SalVentAttn = "#c43afa", Limbic = "#dcf8a4", Cont = "#e69422",
            Default = "#cd3e4f")
lol_plot = ggplot(long_df, aes(y = Value, x = Variable)) +
  geom_hline(yintercept = 0, linetype = 3,
             linewidth = 0.25, color = "#000000") +
  geom_segment(aes(y = 0, yend = Value,
                   x = Variable, xend = Variable),
               color="grey",
               linewidth = 0.5) +
  geom_point(aes(color = Network_7_Yan400_Name), size = 1.5) +
  facet_wrap(~ Network_7_Yan400_Name, scales = "fixed", ncol = 7) +
  scale_color_manual(values = col_Yeo) + 
  coord_flip() +
  theme_classic()
lol_plot


# Networks Classification.......................................................
#...............................................................................
hclust_res = hclust(dist(normalized_res[, -c(1)],
                         method = "euclidean"),
                    method = "ward.D2")
plot(hclust_res)
col_grp = c("1" = "#40e0d0", "2" = "#ff00ff")
h_dend = as.dendrogram(hclust_res) %>%
  color_branches(k = 2) %>%
  set("branches_lwd", 0.5) %>%
  set("branches_k_color", value = col_grp, k = 2)
h_dend_gg = as.ggdend(h_dend)
# Fig. 5-c.
fig_CAH = ggplot(h_dend_gg, labels = FALSE, horiz = F) + 
  theme_void() + 
  geom_rect(aes(xmin = c(0.5, 3.5), xmax = c(3.5,7.5),
                ymin = c(-0.1, -0.1), ymax = c(4.5, 4.5), 
                fill=after_scale(clr_alpha(col_grp, 0.1)))) +
  coord_flip()
fig_CAH
# "NbClust" criteria, to select number of clusters:
nom_methode =  c("kl", "ball", "mcclain")
nbPartition = rep(NaN, length(nom_methode))
for (i in 1:length(nom_methode)){
  nbPartition[i] = NbClust(normalized_res[, -1],
                           diss = as.dist(((1 - cor(t(normalized_res[, -1]))) / 2),
                                          diag = FALSE),
                           distance = NULL,
                           min.nc = 1,
                           max.nc = 6,
                           method = "ward.D2",
                           index = nom_methode[i])$Best.nc[1]
  print(paste("Number of Clusters: ", nbPartition[i]," (criteria: ", 
              nom_methode[i], ")", sep=""))
}
nb_cluster_table = table(nbPartition)
nb_cluster = as.integer(names(nb_cluster_table[which.max(nb_cluster_table)]))
restingStateNetwork = data.frame(net = names(cutree(hclust_res, nb_cluster)),
                                 grp = cutree(hclust_res, nb_cluster))
# Average Groups Profile:
avg_normalized_res = merge(normalized_res, restingStateNetwork, 
                           by.x = "Network_7_Yan400_Name", by.y = "net")
avg_grp = avg_normalized_res %>%
  group_by(grp) %>%
  summarise(grp_TRC =  mean(wgth_TRC),
            grp_MitoDMap = mean(wgth_MitoDMap),
            grp_MRC = mean(wgth_MRC),
            grp_L5IT = mean(wgth_L5IT),
            grp_MicroPVM = mean(wgth_MicroPVM),
            grp_L6IT = mean(wgth_L6IT),
            grp_L23IT = mean(wgth_L23IT))
# Visualization Average Profile by Cluster:
# Fig. 5-d.
long_grp = avg_grp %>%
  pivot_longer(cols = grp_TRC:grp_L23IT,
               names_to = "Variable", values_to = "Value")
long_grp$grp = factor(long_grp$grp, levels = c("2", "1"))
long_grp$Variable = factor(long_grp$Variable,
                           levels = gsub("wgth", "grp", average_values$Variable))
lol_grp = ggplot(long_grp, aes(y = Value, x = Variable)) +
  geom_hline(yintercept = 0, linetype = 3,
             linewidth = 0.25, color = "#000000") +
  geom_segment(aes(y = 0, yend = Value,
                   x = Variable, xend = Variable),
               color="grey",
               linewidth = 0.5) +
  geom_point(aes(color = grp), size = 1.5) +
  ylim(-1, 1) +
  facet_wrap(~ grp, scales = "fixed", ncol = 7) +
  scale_color_manual(values = col_grp) + 
  coord_flip() +
  theme_classic()
lol_grp


# Sign Test by Group on Lateralization..........................................
#...............................................................................
long_data = avg_normalized_res[, -1] %>%
  pivot_longer(cols = c(wgth_TRC:wgth_L23IT),
               names_to = "Variable",
               values_to = "Value")
long_data = long_data %>%
  mutate(Sign = case_when(Value > 0 ~ 'Positive',
                          Value < 0 ~ 'Negative',
                          TRUE ~ 'Zero'))
long_data_no_zero = long_data %>% filter(Sign != 'Zero')
sign_counts = long_data_no_zero %>%
  group_by(grp, Sign) %>%
  summarise(Count = n()) %>%
  ungroup()
# Barnard's Exact Test..........................................................
barnard.test(14, 7, 8, 20)
# Odds Ratio:
(14*20)/(8*7)
