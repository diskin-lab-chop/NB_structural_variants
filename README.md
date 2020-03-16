# NB_structural_variants

## Somatic structural variation targets neurodevelopmental genes and identifies SHANK2 as a tumor suppressor in neuroblastoma
Gonzalo Lopez<sup>1,2+</sup>, Karina L. Conkrite<sup>1,3+</sup>, Miriam Doepner<sup>1,3</sup>, Komal S. Rathi<sup>4</sup>, Apexa Modi<sup>1,3,11</sup>, 
Zalman Vaksman<sup>1,3,4</sup>, Lance M. Farra<sup>1,3</sup>, Eric Hyson<sup>1,3</sup>, Moataz Noureddine<sup>1,3</sup>, Jun S. Wei<sup>5</sup>, Malcolm A. Smith<sup>10</sup>, Shahab Asgharzadeh<sup>6,7</sup>, Robert C. Seeger<sup>6,7</sup>, Javed Khan<sup>5</sup>, Jaime Guidry Auvil<sup>9</sup>, Daniela S. Gerhard<sup>9</sup>, John M. Maris<sup>1,3,11,12</sup>, Sharon J. Diskin<sup>1,3,4,10,11*</sup>

1.	Division of Oncology, Children’s Hospital of Philadelphia, Philadelphia, PA, USA.
2.	Department of Genetics and Genomic Sciences and Icahn Institute for Data Science and Genomic Technology, Icahn School of Medicine at Mount Sinai, New York, NY 10029, USA
3.	Center for Childhood Cancer Research, Children’s Hospital of Philadelphia, Philadelphia, PA, USA.
4.	Department of Biomedical and Health Informatics, Children’s Hospital of Philadelphia, Philadelphia, PA, USA.
5.	Oncogenomics Section, Genetics Branch, Center for Cancer Research, National Cancer Institute, Bethesda, MD, USA
6.	Division of Hematology, Oncology and Blood and Marrow Transplantation, Keck School of Medicine of the University of Southern California, Los Angeles, CA, USA.
7.	The Saban Research Institute, Children’s Hospital of Los Angeles, Los Angeles, CA, USA.
8.	Department of Pediatrics, Perelman School of Medicine, University of Pennsylvania, Philadelphia, PA, USA.
9.	Office of Cancer Genomics, National Cancer Institute, Bethesda, MD, USA.
10.	Cancer Therapy Evaluation Program, National Cancer Institute, Bethesda, MD, USA.
11.	Genomics and Computational Biology, Biomedical Graduate Studies, Perelman School of Medicine, University of Pennsylvania, Philadelphia, PA, USA.
12.	Abramson Family Cancer Research Institute, Perelman School of Medicine at the University of Pennsylvania, Philadelphia, PA, USA.


This repository contains all code and processed data necessary to reproduce analysis and figures in the above cited manuscript.

## R code :
### R/sv_somatic_nbl-pantarget_analysis_V2.R
Structural variant analysis applied to all [TARGET](https://ocg.cancer.gov/programs/target) WGS cancer datasets; The code implements variant filtering (common variant and artifact removal) and variant annotations (gene-level recurrent alterations)
### R/sv_somatic_nbl-segment_data_V2.r
Analysis of CNV breakpoints across [TARGET](https://ocg.cancer.gov/programs/target) WGS cancer datasets and SNP array neuroblastoma datasets; The code implements identification of breakpoints and annotation of recurrently altered genes (amplifications, deletions, etc)
### R/sv_somatic_nbl-pantarget_figures_V2.R
Uses the output from SV and CNV analyses to generate additional analyses and figures from manuscript.
### R/sv_somatic_nbl-segment_data_V2.r, R/my_stat_functions.r and R/heatmap3.R 
Contain required functions used by the main scripts


