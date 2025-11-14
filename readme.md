The molecular and cellular underpinnings of human brain lateralization
================

[![DOI](https://zenodo.org/badge/965483781.svg)](https://doi.org/10.5281/zenodo.15207453)

------------------------------------------------------------------------

## Contents

- [Background](#background)
- [Reference](#reference)
- [Code Release](#code-release)
- [Data](#data)
- [Atlas Used](#atlas-used)
- [Other related papers that might interest
  you](#other-related-papers-that-might-interest-you)
- [Questions](#questions)

------------------------------------------------------------------------

## Background

**Hemispheric specialization** is a *fundamental characteristic of human
brain organization*, where most individuals exhibit **left-hemisphere
dominance for language** and **right-hemisphere dominance for
visuospatial attention**. While some lateralized functions are evident
in other species, the human brain displays a strong, species-wide bias.
*Despite the evolutionary and functional significance of these
asymmetries, their molecular and cellular foundations remain poorly
understood*. Here, we identify key neurochemical and cellular
asymmetries that underpin cortical lateralization. Specifically, we
demonstrate lateralized gradients in neurotransmitter receptor
densities, particularly along the **acetylcholine-norepinephrine axis**,
as well as asymmetries in **mitochondrial distribution** and the spatial
prevalence of **microglia** and **glutamatergic excitatory neurons**.
Using a multimodal approach that integrates *in vivo* functional MRI,
PET imaging, and post-mortem transcriptomic and cellular data, we
delineate **two distinct cortical clusters**: a left-lateralized network
centered on language processing and a right-lateralized network
supporting visuospatial attention. These results highlight a
biologically embedded **substrate for lateralized cognition** that may
inform both evolutionary theory and our mechanistic understanding of
neuropsychiatric illnesses characterized by disrupted lateralization.

------------------------------------------------------------------------

## Reference

For usage of the ***manuscript***, please cite:

- **Labache, L.**, Chopra, S., Zhang (张喜寒), X-H., & Holmes, A. J.
  (2025). The molecular and cellular underpinnings of human brain
  lateralization. *BioRxiv*. DOI:
  [10.1101/2025.04.11.648388](https://doi.org/10.1101/2025.04.11.648388)

For usage of the associated ***code***, please also cite:

- **Labache, L.** (2025). The molecular and cellular underpinnings of
  human brain lateralization. loiclabache/Labache_2025_MolCelLat. DOI:
  [10.5281/zenodo.15207454](https://zenodo.org/doi/10.5281/zenodo.15207454)

------------------------------------------------------------------------

## Code Release

The `Script` folder includes five `R` scripts. These `R` scripts are
designed to facilitate the replication of results as detailed in the
`Method Section` of the **manuscript**.

- `1_CCA.R`: `R` script to conduct the Canonical Correlation Analysis
  (CCA) to assess the functional cortical lateralization –
  neurotransmitters density asymmetries associations.
- `2_MolecularPhenotypes_Analysis.R`: `R` script to identify the
  molecular biomarkers spatially correlated to the
  acetylcholine-norepinephrine axis.
- `3_CellularSubtypes_Analysis.R`: `R` script to explore the multimodal
  cellular underpinning of the acetylcholine-norepinephrine axis.
- `4_NetworkProfil_Analysis.R`: `R` script to identify the mitochondrial
  and cellular fingerprints of large-scale brain networks.
- `5_SørensenDiceIndex_NetworksComparison.R`: `R` script to compute
  similarities between the identified clusters’ topography and a set of
  task-defined networks using the Sørensen-Dice index (SDI).

Note that the `R` scripts also contain the code **to reproduce the
figures found in the manuscript**. The brain renderings in the paper
require [Surf Ice](https://www.nitrc.org/projects/surfice/).

------------------------------------------------------------------------

## Data

All the data necessary to reproduce the results are available in the
`Data` folder.

- `4tasks_typicalStrong_370ROIs_AICHA.csv` corresponds to the average
  brain activity of 125 strongly typical participants (see [Labache, L.,
  et al. 2020](https://doi.org/10.7554/eLife.58722)) during three
  language tasks and a visuospatial attention task in each of the 370
  regions of the AICHA atlas.
  `4tasks_typicalStrong_185ROIs_AICHA_asym.csv` corresponds to the
  asymmetry for each of the 185 regions of the AICHA atlas during these
  four tasks.

- `19neurotransmitters_163ROIs_AICHA_asym_zScored.csv` corresponds to
  the asymmetry of the 19 neurotransmitter systems maps (see [Hansen, J.
  Y., et al. 2022](https://doi.org/10.1038/s41593-022-01186-3)) for each
  of the 163 regions of the AICHA atlas.
  `19neurotransmitters_384ROIs_AICHA_zScored.csv` corresponds to the
  left and right values of the 19 neurotransmitter maps for each of the
  384 regions of the AICHA atlas.

- `24_CellularSubtypes_304ROIs_AICHA_zscored.csv` corresponds to 24
  cell-type maps (see [Zhang, X-H., et
  al. 2024](https://doi.org/10.1038/s41593-024-01812-2)) for each of the
  304 regions of the AICHA atlas.

- `MitoDMap_384ROIs_AICHA_zscored.csv`, `MRC_384ROIs_AICHA_zscored.csv`,
  and `TRC_384ROIs_AICHA_zscored.csv` corresponds to three mitochondrial
  phenotypes maps (see [Mosharov, E. V., et
  al. 2025](https://doi.org/10.1038/s41586-025-08740-6)) for each of the
  384 regions of the AICHA atlas.

- `modeOne_CCA_163ROIs_AICHA.csv` corresponds to the CCA weights
  (assymetries) for the neurotransmitter and task modes for each of the
  163 regions of the AICHA atlas.

- `SynapticDensityMap_384ROIs_AICHA_zscored.csv` corresponds to the
  synaptic density map (see [Johansen, A., et
  al. 2024](https://doi.org/10.1523/JNEUROSCI.1750-23.2024)) for each of
  the 384 regions of the AICHA atlas.

------------------------------------------------------------------------

## Atlas Used

The atlases used in the paper are available in the `Atlas` folder.

- `Atlas/AICHA` contains the **AICHA** atlas, a functional brain ROIs
  atlas (384 brain regions) based on resting-state fMRI data acquired in
  281 individuals. AICHA regions cover the whole brain, each regions are
  characterized by: 1) their homogeneity of its constituting voxels
  intrinsic activity, and 2) a unique homotopic contralateral
  counterpart with which it has maximal intrinsic connectivity. Full
  description of the atlas can be found there:
  [AICHA](https://www.gin.cnrs.fr/en/tools/aicha/), and the related
  paper there: [Joliot, M., et
  al. 2015](https://doi.org/10.1016/j.jneumeth.2015.07.013).
  - The version of AICHA used in the paper corresponds to the file
    `AICHA.nii` (MNI ICBM 152 space). `AICHA_description.txt` is a
    description of each atlas’ regions.
    `comparison_AICHA_Yan2023_400_2mm.txt` provides the correspondance
    between the AICHA regions and the seven canonical intrinsic network
    as defined by Yan and colleagues [Yan, X., et
    al. 2023](https://doi.org/10.1016/j.neuroimage.2023.120010). The
    compressed file `permutation_23k_AICHA.zip` contains the files
    `permutation_23k_AICHA_Seed1893.csv` and
    `permutation_23k_AICHA_Seed1942.csv`. Those are null maps used
    during the “*spin test*”. They have been generated with the
    [netneurotools](https://github.com/netneurolab/netneurotools)
    `Python` toolbox (release `0.2.5`).
- `Atlas/ALANs` contains **ALANs**, an atlas of the *Lateralized
  visuospatial Attentional Networks* in standardized MNI volume space
  (`ALANs_MNI_ICBM_152_2mm.nii`).
  - The related paper can be found there: [Labache, L., et
    al. 2024](https://doi.org/10.1162/imag_a_00208).
    `ALANs_description.txt` is a description of each atlas’ regions.
- `Atlas/SENSAAS` contains **SENSAAS**, an atlas of the *SENtence
  Supramodal Areas AtlaS* in standardized MNI volume space
  (`SENSAAS_MNI_ICBM_152_2mm.nii`).
  - The related paper can be found there: [Labache, L., et
    al. 2019](https://doi.org/10.1007/s00429-018-1810-2).
    `SENSAAS_description.txt` is a description of each atlas’ regions.
- `Atlas/Yan2023_v0283` contains the version of the atlas from Yan and
  colleagues used in this work
  (`400Parcels_Yeo2011_7Networks_FSLMNI152_2mm.nii.gz`, version
  `0.28.3`). See their [Github
  repository](https://github.com/ThomasYeoLab/CBIG/tree/834990110778fe649cb1b9dca3396d6b08b1d4d9/stable_projects/brain_parcellation/Yan2023_homotopic)
  for more relevant information.
  - The related paper can be found there: [Yan, X., et
    al. 2023](https://doi.org/10.1016/j.neuroimage.2023.120010).
    `400Parcels_Yeo2011_7Networks_LUT.txt` is a description of each
    atlas’ regions.
- The file `Atlas/Cognitive_Atlases.txt` corresponds to the association
  of each AICHA region with one of the four atlases of lateralized
  cognitive functions:
  [ALANs](https://github.com/loiclabache/ALANs_brainAtlas),
  [HAMOTA](https://github.com/loiclabache/HAMOTA_brainAtlas),
  [SENSAAS](https://github.com/loiclabache/SENSAAS_brainAtlas), and
  [WMCA](https://github.com/loiclabache/WMCA_brainAtlas).

------------------------------------------------------------------------

## Other related papers that might interest you

- Sentence Supramodal Areas Atlas; Labache, L., et al. 2019. DOI:
  [10.1007/s00429-018-1810-2](https://doi.org/10.1007/s00429-018-1810-2),
  and related GitHub repository:
  [SENSAAS](https://github.com/loiclabache/SENSAAS_brainAtlas)
- Lateralized visuospatial Attentional Networks; Labache, L., et
  al. 2024. DOI:
  [10.1162/imag_a_00208](https://doi.org/10.1162/imag_a_00208), and
  related GitHub repository:
  [ALANs](https://github.com/loiclabache/ALANs_brainAtlas)
- The cell-type underpinnings of the human functional cortical
  connectome; Zhang, X-H., et al. 2025. DOI:
  [10.1038/s41593-024-01812-2](https://doi.org/10.1038/s41593-024-01812-2),
  and related GitHub repository:
  [human-cellular-func-con](https://github.com/XihanZhang/human-cellular-func-con)
- A human brain map of mitochondrial respiratory capacity and diversity;
  Mosharov, E. V., et al. 2025. DOI:
  [10.1038/s41586-025-08740-6](https://doi.org/10.1038/s41586-025-08740-6)
- Influence of Language Lateralisation on Gradient Asymmetry: Labache,
  L., et al. 2023. DOI:
  [10.1038/s41467-023-39131-y](https://doi.org/10.1038/s41467-023-39131-y),
  and related GitHub repository:
  [Labache_2022_AO](https://github.com/loiclabache/Labache_2022_AO)

------------------------------------------------------------------------

## Questions

Please contact me (Loïc Labache) at: <loic.labache@rutgers.edu> and/or
<loic.labache@ensc.fr>
