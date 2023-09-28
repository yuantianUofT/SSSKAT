# SSSKAT
R function for Semi-Supervised SKAT (SSSKAT) Analysis for Rapid and Reliable Gene-based Association Test Using Phenotyping Algorithms Surrogates

# Installation
```{R, eval = FALSE}
devtools::install_github(repo = "https://github.com/yuantianUofT/SSSKAT")
```

# Example
```{R, eval = FALSE
library(SKAT)
library(SSSKAT)
data("example_data")
analysis_data <- example_data$my_dat
labeled_data_id <- example_data$id.t

# Standard SKAT with the labeled data only
labeled_SKAT_result <- labeled_SKAT(dat = analysis_data, id.t = labeled_data_id)

# Semi-supervised SKAT
SS_SKAT_results <- SS_SKAT(dat = analysis_data, id.t = labeled_data_id, weights = NULL, full_NR_evaluation = TRUE, NULL_nlog_like = NULL_nlog_like, prev_est = mean(analysis_data$Y[labeled_data_id]), nboot = 5)

# SKAT with thresholded S
Y_threshold <- get_Y_threshold(dat = analysis_data, prev_value = mean(analysis_data$Y[labeled_data_id]))
threshold_SKAT_results <- threshold_SKAT(dat = analysis_data, id.t = labeled_data_id, thre_value = Y_threshold)

# SKAT directly with S
naive_SKAT_results <- naive_SKAT(dat = analysis_data)

```
