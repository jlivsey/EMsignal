dataList <- list()

for(i in 1:200){

set.seed(i)

# ---- Simulate Data ----------------------------------------------------------
N = 3
T <- TT <- 300
t = 1:T
Phi=diag(N)
Sig=diag(N); Sig[1,2] <- Sig[2,1] <- .75

s1 = gen_trendComp(T, Phi, Sig)
s2 = gen_seasComp(T, Phi, diag(N))
s0 = rmvnorm(n = T, mean = rep(0,N), sigma = diag(N))

data = s1+s2+s0
data = demean(data)

# plot(ts(data), main="simulated series")

# Load
start.date <- c(1990, 1)
period <- 12
dataALL.ts <- sigex.load(data = data,
                         start.date = start.date,
                         period = period,
                         epithets = c("Dim1","Dim2","Dim3"),
                         plot = TRUE)

# Prep for sigex
transform <- "none"
aggregate <- FALSE
subseries <- c(1,2,3)
begin.date <- start.date
end.date <- end(dataALL.ts)
range <- list(begin.date,end.date)
data.ts <- sigex.prep(data.ts = dataALL.ts,
                      transform = transform,
                      aggregate = aggregate,
                      subseries = subseries,
                      range = range,
                      plot = TRUE)

dataList[[i]] <- data.ts

}
