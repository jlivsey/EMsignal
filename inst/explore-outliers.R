# Look at outliers in estimates of irregular covariance matrix
# Seen in Figure 3 and 4, boxplots

# Find index of these outliers
# Extract Sigma[1, 1] for Irregular component over all 200 reps
Sig11_irr = rep(NA, 200)
for(i in 1:200){
  Sig_irr = out[[i]]$Sig.save[52, 19:27] |>
    matrix(nrow = 3, ncol = 3)
  Sig11_irr[i] = Sig_irr[1, 1]
}

# index of max [1, 1] entry of Irregular Sigma estimate
idx = Sig11_irr |> which.max()

# Plot data
dataList[[idx]] |> plot()
dataList[[12]] |> plot()

# ACF/CCF of differenced data
dataList[[idx]] |> diff() |> diff(lag = 12) |> acf()
dataList[[12]] |> diff() |> diff(lag = 12) |> acf()

# ACF/CCF of original data
dataList[[idx]] |> spec()

# Print Sigma estimates for Trend, Seasonal and Irregular
out[[idx]]$Sig.save[52, 1:9] |>
  matrix(nrow = 3, ncol = 3)
out[[idx]]$Sig.save[52, 10:18] |>
  matrix(nrow = 3, ncol = 3)
out[[idx]]$Sig.save[52, 19:27] |>
  matrix(nrow = 3, ncol = 3)



# ---- Repeat but for the [3, 3] entry of the irregular

# Find index of these outliers
# Extract Sigma[3, 3] for Irregular component over all 200 reps
Sig33_irr = rep(NA, 200)
for(i in 1:200){
  Sig_irr = out[[i]]$Sig.save[52, 19:27] |>
    matrix(nrow = 3, ncol = 3)
  Sig33_irr[i] = Sig_irr[3, 3]
}

# index of max [1, 1] entry of Irregular Sigma estimate
idx = Sig33_irr |> which.max()

# Plot data
dataList[[idx]] |> plot()
dataList[[12]] |> plot()

# ACF/CCF of differenced data
dataList[[idx]] |> diff() |> diff(lag = 12) |> acf()
dataList[[12]] |> diff() |> diff(lag = 12) |> acf()

# ACF/CCF of original data
dataList[[idx]] |> spec()

# Print Sigma estimates for Trend, Seasonal and Irregular
out[[idx]]$Sig.save[52, 1:9] |>
  matrix(nrow = 3, ncol = 3)
out[[idx]]$Sig.save[52, 10:18] |>
  matrix(nrow = 3, ncol = 3)
out[[idx]]$Sig.save[52, 19:27] |>
  matrix(nrow = 3, ncol = 3)

# ---- Store min of Trend-sigma upper.tri
# Looking for spurious negatively correlated trend estimates
Sig_min_trend = rep(NA, 200)
for(i in 1:200){
  Sig_trend = out[[i]]$Sig.save[52, 1:9] |>
    matrix(nrow = 3, ncol = 3)
  Sig_min_trend[i] = min(Sig_trend[upper.tri(Sig_trend)])
}

Sig_max_irr = rep(NA, 200)
for(i in 1:200){
  Sig_irr = out[[i]]$Sig.save[52, 19:27] |>
    matrix(nrow = 3, ncol = 3)
  Sig_max_irr[i] = max(diag(Sig_irr))
}

plot(Sig_min_trend, Sig_max_irr); abline(h = 0); abline(v = 0)

# ---- Print Sigma for Trend, Seasonal and Irregular for any run
run_idx = 13
out[[run_idx]]$Sig.save[52, 1:9] |>
  matrix(nrow = 3, ncol = 3)
out[[run_idx]]$Sig.save[52, 10:18] |>
  matrix(nrow = 3, ncol = 3)
out[[run_idx]]$Sig.save[52, 19:27] |>
  matrix(nrow = 3, ncol = 3)
