# aesthetic_value
Repo for the global distribution and drivers of the aesthetic value of reef fish communities

## content

:file_folder:[data](/data) contains raw data files that should not be modified

:file_folder:[scripts](/scripts) contains R scripts for running analyses

:file_folder:[outputs](/outputs) contains outputs and results stored as .Rdata, .rds, etc.

:file_folder:[figures](/figures) contains figures created during analyses

:file_folder:[R](/R) contains functions created during analyses

## Usage

Clone the repository and run this command in R/RStudio:

```r 
source("make.R")
```
All required packages will be installed (if necessary) and loaded.
> :boom: WARNING: running `make.R` calls all the scripts and takes days so if you want to work on one or a few scripts, you should run lines XX-XX of `make.R` and then go to the other script.

Enjoy!