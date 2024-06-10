Conserving the beauty of the world's reef fish communities <img src="https://raw.githubusercontent.com/FRBCesab/templates/main/logos/compendium-sticker.png" align="right" style="float:right; height:120px;"/>
=========================================================

<!-- badges: start -->
[![License: GPL (>= 2)](https://img.shields.io/badge/License-GPL%20%28%3E%3D%202%29-blue.svg)](https://choosealicense.com/licenses/gpl-2.0/)
[![DOI](https://zenodo.org/badge/488270412.svg)](https://zenodo.org/doi/10.5281/zenodo.11551780)
<!-- badges: end -->



<p align="left">
  • <a href="#overview">Overview</a><br>
  • <a href="#data-sources">Data sources</a><br>
  • <a href="#content">Content</a><br>
  • <a href="#installation">Installation</a><br>
  • <a href="#usage">Usage</a><br>
  • <a href="#citation">Citation</a><br>
  • <a href="#contributing">Contributing</a><br>
  • <a href="#acknowledgments">Acknowledgments</a><br>
  • <a href="#references">References</a>
</p>



## Overview

This research compendium provides code and data used to reproduce analyses and figures of the paper: 

> McLean _et al._ (2024) Conserving the beauty of the world's reef fish communities. Submitted to **Nature Sustainability**.


Resources:

- Manuscript : [link](https://docs.google.com/document/d/1gE_Z4GQoFHjzEu0WmANyfg5s1eBWuYOTKFAO3kYfxmY/edit#)
- Supplementary materials: [link](https://docs.google.com/document/d/18sS0vOqeFM_fM3VgSM2i-P693Ey9xNt_jW1McoKyZew/edit)
- Figures: [link](https://docs.google.com/document/d/1r4qGwa2xJKDPo24HSP9H05pNCAT_VrS4sW-Qk_q5BZg/edit)



## Data sources

This project uses the following databases:

:warning: **List all databases used in this project**



## Content

This repository is structured as follow:

- :file_folder: [data/](https://github.com/mcleamj/aesthetic_value/tree/main/data) contains raw data files that should not be modified
- :file_folder: [scripts/](https://github.com/mcleamj/aesthetic_value/tree/main/scripts) contains R scripts for running analyses
- :file_folder: [outputs/](https://github.com/mcleamj/aesthetic_value/tree/main/outputs) contains outputs and results stored as `.Rdata`, `.rds`, etc.
- :file_folder: [figures_tables/](https://github.com/mcleamj/aesthetic_value/tree/main/figures_tables) contains figures and tables created during analyses
- :file_folder: [R/](https://github.com/mcleamj/aesthetic_value/tree/main/R) contains functions created during analyses
- :page_facing_up: [DESCRIPTION](https://github.com/mcleamj/aesthetic_value/tree/main/DESCRIPTION) contains project metadata



## Installation

To install this compendium:

- [Fork](https://docs.github.com/en/get-started/quickstart/contributing-to-projects) 
this repository using the GitHub interface.
- [Clone](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) 
your fork using `git clone fork-url` (replace `fork-url` by the URL of your fork). 
Alternatively, open [RStudio IDE](https://posit.co/products/open-source/rstudio/) 
and create a New Project from Version Control.



## Usage

Before running any analysis listed in [scripts/](https://github.com/mcleamj/aesthetic_value/tree/main/scripts), please install the required packages by running:

```r
# Install 'remotes' package
install.packages("remotes")

# Install required packages
remotes::install_deps()
```



## Citation

Please use the following citation: 

> McLean _et al._ (2024) Conserving the beauty of the world's reef fish communities. URL: <https://github.com/mcleamj/aesthetic_value/>.



## Contributing

All types of contributions are encouraged and valued. For more information, 
check out our [Contributor Guidelines](https://github.com/mcleamj/aesthetic_value/blob/main/CONTRIBUTING.md).

Please note that this project is released with a 
[Contributor Code of Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html). 
By contributing to this project, you agree to abide by its terms.
