# RLS_Aesthetic_value
Research compendium to reproduce analyses and figures of the following article: "XXX" Mouquet, McLean et al. Submitted to XX.

## content

:file_folder:[data](/data) contains raw data files that should not be modified

:file_folder:[scripts](/scripts) contains R scripts for running analyses

:file_folder:[outputs](/outputs) contains outputs and results stored as .Rdata, .rds, etc.

:file_folder:[figures_tables](/figures_tables) contains figures and table created during analyses

:file_folder:[R](/R) contains functions created during analyses

## Figures 
  The following Figures and Tables can be reproduced with the script indicated in brackets (all in `script/`):
  
  FIG_1.png is in `species_level/all_species.R`
  FIG_2.png is in `survey_level/survey_metrics.R`

## Usage

  Clone the repository and run this command in R/RStudio:

```r 
source("make.R")
```
All required packages will be installed (if necessary) and loaded.
> :boom: WARNING: running `make.R` calls all the scripts and takes days so if you want to work on one or a few scripts, you should run lines XX-XX of `make.R` and then go to the other script.

Enjoy!