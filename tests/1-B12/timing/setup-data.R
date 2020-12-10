# ---- Input parameters ----
TT = c(50, 100, 200, 400, 800)
N  = c(2, 3, 4, 5)

# ---- Create data.frame of all param combos ----
paramDf = data.frame(expand.grid(TT, N))
colnames(paramDf) <- c("TT", "N")

# iniitalize output object
dataList <- vector(mode = "list", length = nrow(paramDf))

# Create data for each row of paramDf
for(i in 1:nrow(paramDf)){
  set.seed(i)
  # ---- sim data for ith row parameters ----
  data.ts <- gen_data(paramDf$N[i], paramDf$TT[i])
  # ---- input data to List ----
  dataList[[i]] <- data.ts
}

save(dataList, file = 'dataList.Rdata')
save(paramDf , file = 'paramDf.Rdata')
