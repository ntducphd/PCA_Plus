# 🌐 PCA & ANOVA Shiny App (ver 2.0)
https://ntducphd.shinyapps.io/pca_plus/

## Features

This Shiny App provides comprehensive statistical analysis tools for research data:

### PCA Analysis Module
- Upload a dataset (CSV format) of numeric traits (e.g. from plant breeding or phenotyping studies)
- Perform PCA (Principal Component Analysis)
- Visualize:
  - Scree Plot (Eigenvalues)
  - PCA Biplot
  - Variable Contributions (Vector Plot)
  - HCPC Dendrogram (Hierarchical Clustering on Principal Components)
- Export publication-quality plots (600 dpi)
- Adjust font and figure size for Nature/Science standards

### ANOVA Analysis Module (NEW)
- **Data Input**: Upload CSV or Excel files (.csv, .xlsx, .xls) or use sample data
- **ANOVA Models**: 
  - One-way ANOVA (single factor)
  - Two-way ANOVA (two factors with interaction)
- **Variable Selection**: Easy selection of categorical factors and numeric response variables
- **ANOVA Table**: Detailed analysis of variance table with F-statistics and p-values
- **Post-hoc Tests**: 
  - Tukey HSD test
  - Bonferroni correction
- **Visualization**: Nature journal style boxplots with:
  - White background and black borders
  - Clear data points overlay
  - No grid lines
  - High-quality font rendering
  - Downloadable in PNG or TIFF format (600 dpi)

## Required R Packages

```r
# Core packages
install.packages(c("shiny", "FactoMineR", "factoextra", "DT", "ggplot2", 
                   "plotly", "psych", "corrplot"))

# ANOVA packages
install.packages(c("car", "emmeans", "multcomp", "readxl"))
```

## Usage

1. **For PCA Analysis**: Select the "PCA Analysis" tab and follow the existing workflow
2. **For ANOVA Analysis**: 
   - Select the "ANOVA Analysis" tab
   - Upload your data or use sample data
   - Choose ANOVA type (one-way or two-way)
   - Select factor(s) and response variable
   - Click "Run ANOVA"
   - View results in the tabs: Data Preview, ANOVA Table, Post-hoc Tests, Boxplot

## Notes

- All existing PCA functionality remains unchanged
- ANOVA module is independent and can be extended in the future
- Both modules support high-quality export for publications


