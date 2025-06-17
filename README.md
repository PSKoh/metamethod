# Meta-Analysis in R — Traditional & Multilevel Approaches
This repository contains example datasets and R code to demonstrate how to conduct both traditional meta-analysis and multilevel meta-analysis in R. These materials are intended for educational and non-commercial purposes.

Each folder corresponds to a specific meta-analytic approach:

- DYSCRE contains files related to the traditional meta-analysis.

- SPC includes files for the multilevel meta-analysis using Hedge’s g as the effect size.

- NFCWB holds files for the multilevel meta-analysis using Pearson’s r as the effect size.

- Additional information holds the PRISMA diagram template, as well as the supplementary document accompanying the [tutorial](https://doi.org/10.31234/osf.io/fmw7t_v1) associated with this repository. 
  

If you use or adapt these materials, please provide proper credit by linking back to this repository.

# Repository Contents
"DYSCRE.csv": This dataset is from Majeed et al. (2021), which investigates the association between clinically diagnosed dyslexia and creativity. For more information, one may refer [here](https://doi.org/10.1002/dys.1677).

"SPC.csv": This dataset is from Hartanto et al. (2024), which investigates the effects of smartphone presence on cognitive functions. For more information, one may refer [here](https://doi.org/10.1037/tmb0000123).

"NFCWB.csv": This dataset is from Lua et al. (2023), which investigates the relationship between the need for cognition and well-being. For more information, one may refer [here](https://doi.org/10.1007/s11031-023-10047-w).

"PRISMA Example.pptx": A modifiable PRISMA template is available for researchers to adapt it to their own systematic review or meta-analysis. 

"Supplementary Document.pdf": This document contains information that was not covered in the [tutorial](https://doi.org/10.31234/osf.io/fmw7t_v1) such as data structure formats. 

# R Scripts and Notebooks
"meta_traditional.R": R script for conducting traditional meta-analysis, with Hedges' g as the effect size 

"meta_multilevel.R": R script for conducting multilevel meta-analysis, with Hedges' g as the effect size

"meta_multilevel (Pearson r).R": R script for conducting multilevel meta-analysis, with Pearson's r as the effect size 

Additionally, annotated R Notebooks are provided to walk you through each step of the analysis, including data preparation, model specification, and interpretation of results.
