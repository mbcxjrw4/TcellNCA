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

🧬 Usage
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

⚙️ Functions Overview
| Function                   | Description                           |
| -------------------------- | ------------------------------------- |
| `calculate_nca_tcell()`    | NCA metrics for T-cell persistence    |
| `calculate_nca_cytokine()` | NCA metrics for serum cytokine levels |

📋 Available Metrics

