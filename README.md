# HPLC van 't Hoff Data Analysis Workflow

This directory contains all R code scripts for the HPLC chromatographic data processing and van 't Hoff thermodynamic analysis.

---

## Directory Structure

```text
.
├── classicVantHoff/
│   └── classic_vantHoff_analyzer.R    # Core R implementation for van 't Hoff analysis
├── read_analyze_data/
│   ├── run_batch_chromReader.R        # Script for batch-reading raw HPLC chromatograms
│   ├── analyze_classic_vantHoff.R     # Script executing the van 't Hoff workflow
│   └── validate_classic_vantHoff.R    # Script executing the validation of error propagation
├── data/                              # Sample raw chromatogram files (.arw)
└── README.md                          # This documentation file
```

---

## File Descriptions

| Directory / File | Description |
| :--- | :--- |
| **`classicVantHoff/classic_vantHoff_analyzer.R`** | Contains the core implementation and mathematical functions for classic van 't Hoff modeling and thermodynamic parameter extraction ($\Delta H^\circ$, $\Delta S^\circ$). |
| **`read_analyze_data/run_batch_chromReader.R`** | Script designed to batch-process raw HPLC chromatograms located in the `data/` directory using the `autoChromReader` engine. |
| **`read_analyze_data/analyze_classic_vantHoff.R`** | Main analysis pipeline script that processes parsed data, fits thermodynamic models, and generates figures/results. |
| **`read_analyze_data/validate_classic_vantHoff.R`** | Validation analysis that gauges sensitivity of fit workflow. |
| **`data/`** | Folder containing raw sample chromatographic data files (e.g., `.arw` formats across multiple detection wavelengths). |

---

## Dependencies & External Tools

* **R Environment:** R version 4.0 or higher is recommended.
* **`autoChromReader` Engine:** The batch chromatographic reading depends on the `autoChromReader` function package. You can access the public source code and installation instructions on GitHub:
  * 🔗 **GitHub Repository:** [autoChromReader on GitHub](https://github.com/mirzahosseini-arash-semmelweis/autoChromReader)

---

## Usage / Workflow

To reproduce the analysis workflow from the raw data:

1. **Chromatogram Reading & Preprocessing:**
   Run `read_analyze_data/run_batch_chromReader.R` to parse raw sample files in the `data/` directory into structured peak and retention time data using `autoChromReader`.

2. **Thermodynamic Analysis:**
   Run `read_analyze_data/analyze_classic_vantHoff.R`, which sources `classicVantHoff/classic_vantHoff_analyzer.R` to compute retention parameters and perform linear regressions across temperature sets.
