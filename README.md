# GENERAL INFORMATION

This README.txt file was updated on 2026-10-05 by Shailee S. Shah

## A. Paper associated with this archive 
Citation: Shah et al. In press. Lifetime fitness benefits of short-distance dispersal are associated with inbreeding tolerance despite multiple inbreeding avoidance mechanisms. The American Naturalist.

Brief abstract: It is commonly assumed that mechanisms that reduce inbreeding, such as sex-biased dispersal or kin avoidance, are always adaptive. However, theoretical models predict that inbreeding avoidance evolves only when the relative costs of inbreeding are higher than the costs of inbreeding avoidance. We empirically tested this prediction with 32 years of data on Florida scrub-jays (Aphelocoma coerulescens), birds with short dispersal distances and multiple mating opportunities over their lifetimes. Using simulations, we showed that limited dispersal increases the inbreeding risk by at least two-fold, with a stronger effect at later life stages. Sex-biased dispersal reduces, but does not eliminate, this risk. Additionally, we found that the fitness costs of inbreeding are lower than those of dispersal in terms of lifetime reproductive success (LRS) and survival. Inbreeding is associated with substantially lower LRS only when mates are first-order kin. Accordingly, we found that the population shows a degree of inbreeding tolerance but, within their dispersal constraints, Florida scrub-jays actively avoid mating with first-order kin. Overall, our results empirically demonstrate the context-specific nature of inbreeding avoidance, with selection favoring limited dispersal despite costs of inbreeding and inbreeding avoidance via mate choice only when the relative fitness costs of inbreeding are particularly high.

## B. Originators
Shailee S. Shah, University of Rochester and Cornell University
Jennifer Diamond, University of California, Davis
Sahas Barve, Archbold Biological Station
Elissa J. Cosgrove, Cornell University
Reed Bowman, Archbold Biological Station
John W. Fitzpatrick, Cornell University
Nancy Chen, University of Rochester and University of California, Los Angeles

## C. Contact information
Shailee S. Shah
sss253@cornell.edu

## D. Dates of data collection
1990 to 2021

## E. Geographic Location(s) of data collection
Archbold Biological Station, Venus, Florida, USA (27.10°N, 81.21°W)

## F. Funding Sources 
This work was supported by National Science Foundation (NSF) Grants DEB0855879 and DEB1257628 and by the Cornell Lab of Ornithology Athena Fund. S.S.S was supported by NIH grant 1R35GM133412 to N.C. and NSF Postdoctoral Research Fellowship in Biology 2305705.

# ACCESS INFORMATION

## 1. Licenses/restrictions placed on the data or code
CC0 1.0 Universal (CC0 1.0)
Public Domain Dedication

## 2. Data derived from other sources

none

## 3. Recommended citation for this data/code archive

xxx

# DATA & CODE FILE OVERVIEW

This data repository consist of 7 data files, 3 code scripts, and this README document, with the following data and code filenames and variables

## Repository Structure
├── data/ # folder with .csv or .csv.zip data files
├── src/  # Source code
│   ├── CustomFunctions/   # Two custom R functions for efficient analyses 
│   └── MateChoiceInbreedingToleranceAnalysis.R    # R script for all analyses
└── README.md            # This readme file

### Directory Descriptions

* **`data/`** – Contains data as seven .csv files
* **`src/`** – Contains all code files
  * **`CustomFunctions/`** – Contains two custom functions for efficient analyses of relatedness to potential mates with distance and observed vs. simulated kinship coefficients to mates

## Data files and variables

### 1. Metadata for "BreederData.csv"
This dataset contains one datapoint per breeding instance for all breeders in the population between 1990 - 2021. Associated data include data on the mate as well as year, sex, and natal year of the breeder.

|Column|Entry    |Value         |Explanation   |
|------|---------|--------------|-------|
|A     |FocalInd |1-983         |unique ID of the individual     |
|B     |Mate     |1-1080        |unique ID of the mate        |
|C     |Year     |1990-2021     |Year       |
|D     |Coeff    |0-0.25        |Kinship coefficient to mate |
|E     |Sex      |Male or Female|Sex of the individual   |
|F     |NatalYear|1978-2019     |Year that the individual was born   |
|G     |PairType |First or Later|whether this is the first pair the FocalInd formed in its lifetime or a later pairing after a mate switching event|
|H     |Divorce  |0 or 1        |Pairing ended in divorce (1) or not (0)   |
|I     |CloseRel |0 or 1        |whether the mate is a close relative (first cousin or higher) or not |
|J     |Dist     |0-19          |distance moved for pairing in number of territories crossed  |
|K     |FirstOrderRel|0 or 1    |whether the mate is a first order relative (parent-offspring or full sibling) or not |


### 2. Metadata for "AllPossPairsData.csv"
This dataset pairs each individual breeder per year with all available mates of the opposite sex and contains associated data on distance and relatedness to those possible mates.

|Column|Entry|Value|Explanation|
|------|--------|--------------|-----|
|A|FocalInd|1-983|unique ID of the individual|
|B|Mate|1-1080|unique ID of the mate |
|C|Sex|Male or Female|Sex of the individual|
|D|Year|1990-2021|Year|
|E|Coeff|0-0.34|Kinship coefficient to mate|
|F|Dist|0-19|distance moved for pairing in number of territories crossed|
|G|CloseRel|0 or 1|whether the mate is a close relative (first cousin or higher) or not|
|H|PairType|First or Later|whether this is the first pair the FocalInd formed in its lifetime or a later pairing after a mate switching event|
|I|PairID|1-24918|unique ID of the pairs|


### 3. Metadata for "SimData.csv"
This dataset contains output from 1000 iterations of random mating simulations without any dispersal constraints.

|Column | Entry   | Value         | Explanation |
| ------- | -------- | -------------- | -------|
| A | Mate     | 1-1080   | unique ID of the mate  |
| B       | FocalInd | 1-1027         | unique ID of the individual     |
| C       | Rep      | 1-1000         | number of rep in random mating     |
| D       | Year     | 1990-2021      | Year  |
| E       | Sex      | Male or Female | Sex of the individual |
| F       | Coeff    | 0-0.25         | kinship coefficient to simulated mate   |
| G       | Dist     | 0-19           | simulated distance moved for pairing in number of territories crossed |
| H       | CloseRel | 0 or 1         | whether the mate is a close relative (first cousin or higher) or not    |
| I       | PairType | First or Later | whether this is the first pair the FocalInd formed in its lifetime or a later pairing after a mate switching event |
| J       | PairID   | 1-24918|unique ID of the pairs|

### 4. Metadata for "SimData_displim.csv"
This dataset contains output from 1000 iterations of random mating simulations with distribution of dispersal distances limited to the observed, sex-biased dispersal distances observed in the population.

|Column | Entry   | Value         | Explanation |
| ------- | -------- | -------------- | -------|
| A | Mate     | 1-1080   | unique ID of the mate  |
| B       | FocalInd | 1-1027         | unique ID of the individual     |
| C       | Rep      | 1-1000         | number of rep in random mating     |
| D       | Year     | 1990-2021      | Year  |
| E       | Sex      | Male or Female | Sex of the individual |
| F       | Coeff    | 0-0.25         | kinship coefficient to simulated mate   |
| G       | Dist     | 0-19           | simulated distance moved for pairing in number of territories crossed |
| H       | CloseRel | 0 or 1         | whether the mate is a close relative (first cousin or higher) or not    |
| I       | PairType | First or Later | whether this is the first pair the FocalInd formed in its lifetime or a later pairing after a mate switching event |
| J       | datatype   | "Random mating with sex-biased limited dispersal" |This dataset has the results of a model of random mating within the observed sex-biased limited dispersal distance for all individuals|

### 5. Metadata for "SimData_displim_nosexbias.csv"
This dataset contains output from 1000 iterations of random mating simulations with distribution of dispersal distances limited to the observed dispersal distances observed in the population, but with the sex-biased removed.

|Column | Entry   | Value         | Explanation |
| ------- | -------- | -------------- | -------|
| A | Mate     | 1-1080   | unique ID of the mate  |
| B       | FocalInd | 1-1027         | unique ID of the individual     |
| C       | Rep      | 1-1000         | number of rep in random mating     |
| D       | Year     | 1990-2021      | Year  |
| E       | Sex      | Male or Female | Sex of the individual |
| F       | Coeff    | 0-0.25         | kinship coefficient to simulated mate   |
| G       | Dist     | 0-19           | simulated distance moved for pairing in number of territories crossed |
| H       | CloseRel | 0 or 1         | whether the mate is a close relative (first cousin or higher) or not    |
| I       | PairType | First or Later | whether this is the first pair the FocalInd formed in its lifetime or a later pairing after a mate switching event |
| J       | datatype   | "Random mating with limited dispersal" |This dataset has the results of a model of random mating within the observed limited dispersal distance across individuals but with the sex-bias removed|


### 6. Metadata for "FitnessData.csv"
This data contains lifetime reproductive success, relatedness to mate, and dispersal distance data for all breeders for whom we have complete lifetime data.

| Column | Entry | Value | Explanation |
| -------| ------| ----- | ------------|
|A|FocalInd | 1-889 | unique ID of the individual|
|B|TotalDistLaterPairing | 0-7, 10 |Total distance moved for later pairing in number of territories over the individual's lifetime|
|B|FirstYrBred|1990-2012| Year that the FocalInd first became a breeder|
|E|FitnessTotal| 0-18| Number of offspring that survived to 300 days old|
|F|YrLastObs| 1990-2019 | Year when the individual was last seen |
|G|BreedingLifespan|1-13|Number of years for which the individual bred|
|H|Sex      | Male or Female | Sex|
|I|DistFirstPairing|0-13|Distance moved for first pairing in number of territories|
|J|NatalYear| 1988-2009 | Year that the individual was born |
|K|WtCoeff |0-0.25| Sum of the relatedness of the focal individual to each of their mates multiplied by the number of years they were paired with that mate |
|L|MateSwitch| 1 or NA | Whether the individual ever switched mates in their lifetime|


### 7. Metadata for "SurvivalData.csv"
This data contains survival data (one entry per year from birth to death) for all breeders for whom we have complete lifetime data.


| Column | Entry | Value | Explanation |
| -------| ------| ----- | ------------|
|A|FocalInd | 1-889 | unique ID of the individual|
|B|Dist     | 0-13 | Distance moved for first pairing in number of territories|
|C|Sex      | Male or Female | Sex|
|D|EntryAge | 0-11| Individual's age at the beginning of the year|
|E|ExitAge  | 1-12| Individual's age at the end of the year|
|F|Dead     | 0 or 1 | Whether an individual died at the end of the year (1) or survived (0) |
|G|TerrLost | 0 or 1 | Whether an individual moved away from its breeding territory (1) or not (0) |
|H|NatalYear| 1988-2006, 2009 | Year that the individual was born |
|I|numMateSwitch| 0-6 | Number of times the FocalInd switched mates in its lifetime |


## Code scripts and workflow

1. MateChoiceInbreedingToleranceAnalysis.R

This script has code for running all analyses and making figures and tables in the Main Text. The script calls on scripts that define custom functions for certain analyses (described below)

2. Custom Function Scripts in the "CustomFunctions" directory

	1. ObsSimsFn.R

	This function compares the observed population means to the simulated distributions of means of dispersal distance and relatedness to mate and computes empirical p-values

	2. PotentialMateCoeffFn.R

	This function runs generalized linear mixed models with a binomial error structure to estimate the relationship between the likelihood of a potential mate being a close relative and dispersal distance. The models are run separately by sex and pairing type.

Workflow: Unzip the zipped .csv data files (all simulation results). The main analysis script (MateChoiceInbreedingToleranceAnalysis.R) can then be run as a standalone script. It will load data files from a "data" directory in the working directory and call custom functions from a "CustomFunctions" directory within the "src" directory in the working directory.

# SOFTWARE VERSIONS

R version 4.5.3 

R package versions:
ggeffects_2.3.2
patchwork_1.3.2
ggpattern_1.3.1
DHARMa_0.5.0
jtools_2.3.1
survival_3.8-6
scales_1.4.0
knitr_1.51
glmmTMB_1.1.14
MASS_7.3-65
lme4_2.0-1
Matrix_1.7-5
lubridate_1.9.5
forcats_1.0.1
stringr_1.6.0
dplyr_1.2.1
purrr_1.2.2
readr_2.2.0
tidyr_1.3.2
tibble_3.3.1
ggplot2_4.0.3
tidyverse_2.0.0

# REFERENCES

none


