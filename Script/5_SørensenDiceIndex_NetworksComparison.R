################################################################################
# Written by Loïc Labache, Ph.D.                                               #
# Holmes Lab, Department of Psychiatry - Rutgers University                    #
# February 12, 2025                                                            #
################################################################################


# Libraries.....................................................................
#...............................................................................
packages = c("here", "neurobase", "tidyr")
lapply(packages, require, character.only = T)


# Data..........................................................................
#...............................................................................
path_folder = "project_path"
sensaas = readnii(here(path_folder,  "Atlas", "SENSAAS",
                       "SENSAAS_MNI_ICBM_152_2mm.nii"))
alans = readnii(here(path_folder,  "Atlas", "ALANs",
                     "ALANs_MNI_ICBM_152_2mm.nii"))
yan = readnii(here(path_folder,  "Atlas", "Yan2023_v0283",
                   "400Parcels_Yeo2011_7Networks_FSLMNI152_2mm.nii.gz"))
sensaas_dscrptn = read.csv(here(path_folder,  "Atlas", "SENSAAS",
                                "SENSAAS_description.txt"))
alans_dscrptn = read.csv(here(path_folder,  "Atlas", "ALANs",
                              "ALANs_description.txt"))
yan_dscrptn = read.table(here(path_folder,  "Atlas", "Yan2023_v0283",
                              "400Parcels_Yeo2011_7Networks_LUT.txt"),
                         header = F, stringsAsFactors = FALSE)[, -c(3:6)]
colnames(yan_dscrptn) = c("ParcelID", "ParcelName")
yan_dscrptn = separate(yan_dscrptn, col = ParcelName,
                       into = c("delete", "hemisphere", "network", "region", "subPart"),
                       sep = "_", extra = "merge")
yan_dscrptn$delete = NULL
cluster_molCel = data.frame(network = c("Cont", "Default", "DorsAttn",
                                        "Limbic", "SalVentAttn", "SomMot", "Vis"),
                            cluster = c(1, 2, 1, 2, 2, 1, 2))

# Yan Masks.....................................................................
#...............................................................................
yan_dscrptn_clust = merge(yan_dscrptn, cluster_molCel, 
                          by = "network", all.x = TRUE)
parcel_cluster_lookup = setNames(yan_dscrptn_clust$cluster, 
                                 yan_dscrptn_clust$ParcelID)
img_data = as.array(yan)
voxel_cluster = ifelse(img_data == 0, 0,
                       parcel_cluster_lookup[as.character(img_data)])
mask_cluster1_data = ifelse(voxel_cluster == 1, 1, 0)
mask_cluster2_data = ifelse(voxel_cluster == 2, 1, 0)

# Cognitive Masks...............................................................
#...............................................................................
# SENSAAS:
sensaas_data = as.array(sensaas)
networks_sensaas = unique(sensaas_dscrptn$Network)
sensaas_masks = list()
for(net in networks_sensaas) {
  parcel_cluster_lookup = setNames(rep(1, dim(sensaas_dscrptn[sensaas_dscrptn$Network == net,])[1]),
                                   sensaas_dscrptn[sensaas_dscrptn$Network == net, ]$Index)
  voxel_cluster = ifelse(sensaas_data == 0, 0,
                         parcel_cluster_lookup[as.character(sensaas_data)])
  mask_cluster = ifelse(voxel_cluster == 1, 1, 0)
  sensaas_masks[[net]] = mask_cluster
}
# ALANs:
alans_data = as.array(alans)
networks_alans = unique(alans_dscrptn$Network)
alans_masks = list()
for(net in networks_alans) {
  parcel_cluster_lookup = setNames(rep(1, dim(alans_dscrptn[alans_dscrptn$Network == net,])[1]),
                                   alans_dscrptn[alans_dscrptn$Network == net, ]$Index)
  voxel_cluster = ifelse(alans_data == 0, 0,
                         parcel_cluster_lookup[as.character(alans_data)])
  mask_cluster = ifelse(voxel_cluster == 1, 1, 0)
  alans_masks[[net]] = mask_cluster
}
parcel_cluster_lookup = setNames(rep(1, dim(alans_dscrptn[alans_dscrptn$Network %in% 
                                                            c("ParietoFrontal", "TemporoFrontal"),])[1]),
                                 alans_dscrptn[alans_dscrptn$Network %in% 
                                                 c("ParietoFrontal", "TemporoFrontal"),]$Index)
voxel_cluster = ifelse(alans_data == 0, 0,
                       parcel_cluster_lookup[as.character(alans_data)])
mask_cluster = ifelse(voxel_cluster == 1, 1, 0)
alans_masks[["alans"]] = mask_cluster

# Dice Coefficient..............................................................
#...............................................................................
dice_coef = function(mask1, mask2) {
  mask1 = as.integer(mask1 == 1)
  mask2 = as.integer(mask2 == 1)
  intersection = sum(mask1 & mask2, na.rm = TRUE)
  sum_mask1 = sum(mask1, na.rm = TRUE)
  sum_mask2 = sum(mask2, na.rm = TRUE)
  if ((sum_mask1 + sum_mask2) == 0) {
    return(NA)
  } else {
    return((2 * intersection) / (sum_mask1 + sum_mask2))
  }
}
dice_alans_parieto_cluster1 = dice_coef(alans_masks$ParietoFrontal, mask_cluster1_data)
dice_alans_parieto_cluster2 = dice_coef(alans_masks$ParietoFrontal, mask_cluster2_data)
dice_alans_temporo_cluster1 = dice_coef(alans_masks$TemporoFrontal, mask_cluster1_data)
dice_alans_temporo_cluster2 = dice_coef(alans_masks$TemporoFrontal, mask_cluster2_data)
dice_alans_cluster1 = dice_coef(alans_masks$alans, mask_cluster1_data)
dice_alans_cluster2 = dice_coef(alans_masks$alans, mask_cluster2_data)
dice_sensaas_sent_core_cluster1 = dice_coef(sensaas_masks$SENT_CORE, mask_cluster1_data)
dice_sensaas_sent_core_cluster2 = dice_coef(sensaas_masks$SENT_CORE, mask_cluster2_data)
results = data.frame(Mask = rep(c("alans_ParietoFrontal", "alans_TemporoFrontal", "sensaas_SENT_CORE",
                                  "alans"), each = 2),
                     Cluster_Mask = rep(c("mask_cluster1", "mask_cluster2"), times = 4),
                     Dice_Coefficient = c(dice_alans_parieto_cluster1, dice_alans_parieto_cluster2,
                                          dice_alans_temporo_cluster1, dice_alans_temporo_cluster2,
                                          dice_sensaas_sent_core_cluster1, dice_sensaas_sent_core_cluster2,
                                          dice_alans_cluster1, dice_alans_cluster2))
results$Jaccard = (results$Dice_Coefficient) / (2 - results$Dice_Coefficient)
print(results[c(5:8), ])

