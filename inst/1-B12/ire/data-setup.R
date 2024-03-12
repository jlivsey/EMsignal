# ---- Simulate Data ----------------------------------------------------------

nreps <- 4
dataList <- vector(mode = "list", length = nreps)

for(i in 1:nreps){

  set.seed(i)

  # ---- Input parameters ----
  N  <- 3
  T  <- 300
  TT <- T   # I use TT and Tucker uses T
  t  <- 1:T
  start.date <- c(1990, 1) # generic
  period     <- 12
  transform  <- "none"
  aggregate  <- FALSE
  subseries  <- 1:3

  Phi <- diag(N)
  Sig <- diag(N); Sig[1,2] <- Sig[2,1] <- .75

  # ---- generate structural components
  s1 = gen_trendComp(T, Phi, Sig)
  s2 = gen_seasComp(T, Phi, diag(N))
  s0 = rmvnorm(n = T, mean = rep(0,N), sigma = diag(N))

  data = s1+s2+s0
  data = demean(data)

  # ---- Load ----
  dataALL.ts <- sigex.load(data = data,
                           start.date = start.date,
                           period = period,
                           epithets = c("Dim1", "Dim2", "Dim3"),
                           plot = FALSE)

  # ---- Prep for sigex ----
  begin.date <- start.date
  end.date <- end(dataALL.ts)
  range <- list(begin.date,end.date)
  data.ts <- sigex.prep(data.ts = dataALL.ts,
                        transform = transform,
                        aggregate = aggregate,
                        subseries = subseries,
                        range = range,
                        plot = FALSE)

  # ---- input data to List ----
  dataList[[i]] <- data.ts

}

# save(dataList, file = 'dataList.Rdata')
