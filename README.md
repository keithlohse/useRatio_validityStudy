# useRatio_validityStudy
 Repository for analysis code from our study on the concurrent validity of the Use Ratio compared to the ARAT and FM UE in adults with stroke.  Data are available from NICHD DASH; see Miller, A. E. et al. A Large Harmonized Upper and Lower Limb Accelerometry Dataset: A Resource for Rehabilitation Scientists. medRxiv 2024.08.15.24312066 (2024) doi:10.1101/2024.08.15.24312066.

Data are available from NICD DASH in two separate "studies". Data are indexed as two separate studies due to differences in data use limitations, in which use of data from Part1 is limited to research on neurological conditions and/or movement disorders and use of data from Part2 has no additional limitations:
 1. https://dash.nichd.nih.gov/study/426315
 2. https://dash.nichd.nih.gov/study/426433

After a user has downloaded these data from DASH and unzipped the files (the download will be compressed), there are three .csv files that need to be imported, as shown in lines 17-19 of the analysis script file:
1. UE_DATA <- read.csv("./UpperLimbAccelerometry_restricted.csv", header=TRUE, stringsAsFactors = TRUE)
2. LE_DATA <- read.csv("./LowerLimbAccelerometry_restricted.csv", header=TRUE, stringsAsFactors = TRUE)
3. CLINICAL <- read.csv("./DemographicClinicalData_restricted.csv", header=TRUE, stringsAsFactors = TRUE)

The rest of the analyses follow in script_accel_validity.R, and for the most part follow the order the analyses are presented in the paper. The script is divded into three major sections:
0. Data importing, filtering, and basic formatting
1. Analyses of the cross sectional cohort
2. Analyses of the longitudinal cohort
