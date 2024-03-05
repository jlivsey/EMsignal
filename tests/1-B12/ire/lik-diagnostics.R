# What was the difference in lik'd after 50 iters?

likVal = matrix(NA, nrow = 5, ncol = 200)
for(i in 1:200){
  likVal[, i] = out[[i]]$Sig.save[48:52, 28]
}

mindiff = rep(NA, 200)
for(i in 1:200){
  mindiff[i] = out[[i]]$Sig.save[, 28] |> diff() |> abs() |> min()
}
sum(mindiff < 10^-3)


myidx = 21
out_8_50[[myidx]]$Sig.save[50:300, 28] |> diff() |> abs() |> plot(); abline(h = 0)
out_8_50[[myidx]]$Sig.save[50:300, 28]
out_8_50[[myidx]]$run_time_secs / 60

out_8_50 |> length()

out_8_50[[1]]$Sig.save[, 28]


out = c(out_1_7,
        out_8_50,
        out_51_100,
        out_101_150)
