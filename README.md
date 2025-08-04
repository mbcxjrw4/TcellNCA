# TcellNCA

**TcellNCA** is an R package for calculating non-compartmental analysis (NCA) metrics related to **T-cell persistence** (e.g., CAR-T or TCR-T cell therapies) and **serum cytokine** profiles in cell therapy studies.

It supports preprocessing, flexible metric selection, and quality control of pharmacodynamic biomarker data.

---

## 🚀 Features

- 📊 Calculates key NCA metrics: `Cmax`, `Tmax`, `AUC0-t`, `AUC0-inf`, `lambda_z`, `Halflife`, `Tlast`, `MRT`, etc.
- ✅ Handles preprocessing: removes missing values, sorts timepoints, inserts `time = 0` if missing
- 📦 Designed for pharmacokinetic/pharmacodynamic data from cell therapy trials
- 🔄 Flexible metric selection: calculate all or a custom subset

---

## 📦 Installation

You can install the development version from GitHub:
```r
# install.packages("devtools")
devtools::install_github("your-username/ncaCell")
```

## 🧬 Usage
```r
library(ncaCell)

# Example time and concentration data
time <- c(0, 1, 2, 4, 7, 14)
conc <- c(0, 200, 150, 90, 40, 10)

# Run NCA for T-cell persistence
nca_results <- calculate_nca_tcell(time, conc)

# Run NCA for serum cytokines, requesting specific metrics
cytokine_results <- calculate_nca_cytokine(time, conc, metrics = c("Cmax", "Tmax", "AUC0_t"))
```

## ⚙️ Functions Overview
| Function                   | Description                           |
| -------------------------- | ------------------------------------- |
| `calculate_nca_tcell()`    | NCA metrics for T-cell persistence    |
| `calculate_nca_cytokine()` | NCA metrics for serum cytokine levels |

## 📋 Available Metrics
- `NoSamples`:	This parameter reports the number of non-missing observations used in the analysis of the profile
- `Rsq`:	Goodness of fit statistic for the terminal elimination phase
- `RsqAdjusted`:	Goodness of fit statistic for the terminal elimination phase adjusted for the number of points used in the estimation of Lambda Z
- `CorrXY`:	Correlation between time (X) and log concentration (Y) for the points used in the estimation of Lambda Z
- `NoPointsLambdaz`:	Number of points used in computing Lambda Z
- `Lambdaz`:	First-order rate constant associated with the terminal (log-linear) portion of the curve. Estimated by linear regression of time vs log concentration
- `LambdazIntercept`:	Intercept on log scale estimated via linear regression of time vs. log concentration
- `LambdazLower`:	Lower limit on time for values to be included in the calculation of Lambda Z
- `LambdazUpper`:	Upper limit on time for values to be included in the calculation of Lambda Z
- `HalfLife`:	Terminal half-life
- `Span`:	(LambdazUpper – LambdazLower)/HalfLiveLambdaz
- `Tmax`:	Time of maximum observed concentration
- `Cmax`:	Maximum observed concentration
- `CmaxD`:	Cmax/Transduced Dose
- `Tlast`:	Time of last measurable (positive) observed concentration
- `Clast`:	Observed concentration corresponding to Tlast
- `ClastPred`:	Predicted concentration at Tlast
- `AucLast`:	Area under the curve from the time of dosing to the time of the last measurable (positive) concentration (Tlast)
- `AucAll`:	Area under the curve from the time of dosing to the time of the last observation. If the last concentration is positive AUClast=AUCall
- `AucInfObs`:	AUC from time of dosing extrapolated to infinity based on the last observed concentration. AUClast + (Clast/Lambda_z)
- `AucExtrapObs`:	Percentage of AucInfObs due to extrapolation from Tlast to infinity 100*[(AucInfObs– AucLast)/AucInfObs]
- `AucInfPred`:	AUC from time of dosing extrapolated to infinity based on the  last predicted concentration. AUClast + (ClastPred/Lambda_z)
- `AucExtrapPred`:	Percentage of AucInfPred due to extrapolation from Tlast to infinity 100*[(AucInfPred– AucLast)/AucInfPred]
- `AumcLast`:	Area under the moment curve from the time of dosing to the last measurable (positive) concentration
- `AumcInfObs`:	Area under the first moment curve (AUMC) extrapolated to infinity based on the last observed concentration (obs)
- `AumcExtrapObs`:	Percent of AumcInfObs that is extrapolated from Tlast to infinity
- `AumcInfPred`:	Area under the first moment curve (AUMC) extrapolated to infinity based on the  last predicted concentration (pred)
- `AumcExtrapPred`:	Percent of AumcInfPred) that is extrapolated from Tlast to infinity
- `MrtLast`:	Mean residence time from the time of dosing to the time of the last measurable concentration AumcLast / AucLast
- `MrtInfObs`:	Mean residence time (MRT) extrapolated to infinity based on AucInfObs
- `MrtInfPred`:	Mean residence time (MRT) extrapolated to infinity based on AucInfPred
- `AucD7`:	Area under the curve from day 1 (infusion day) to day 8.
- `AucD28`:	Area under the curve from day 1 (infusion day) to day 29.
- `AucM3`:	Area under the curve from day 1 (infusion day) to day 85.
- `AucM6`:	Area under the curve from day 1 (infusion day) to day 169.
- `VolumeOfDistributionObs`:	Volume of distribution based on the terminal phase Dose / (AucInfObs * LambdaZ)
- `ClearanceObs`:	Total body clearance for extravascular administration Dose/AucInfObs
- `VolumeOfDistributionPred`:	Volume of distribution based on the terminal phase Dose / (AucInfPred * LambdaZ)
- `ClearancePred`:	Total body clearance for extravascular administration Dose/AucInfPred
