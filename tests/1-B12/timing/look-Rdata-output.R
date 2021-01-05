
setwd("~/GitHub/EMsignal/tests/1-B12/timing/Rdata-output")


# Load parameter values data.frame
load("~/GitHub/EMsignal/tests/1-B12/timing/paramDf.Rdata")
df <- paramDf

# Add empty columns to data.frame
df$mle_time <- NA
df$em_time  <- NA

# Add timing results to df if file exists
for(i in 1:20){
  tryCatch({
    filename <- paste0(df$TT[i], "_", df$N[i], ".Rdata")
    load(filename)
    df$mle_time[i] <- outList$mle_time[3]
    df$em_time[i]  <- outList$em_time[3]
  },
  error = function(e){},
  warning = function(w){})
}


library(dplyr)
# df2 <- filter(df, TT = 100)
df2 <- tidyr::drop_na(df)
df2

df3 <- reshape2::melt(df2, measure.vars = c("mle_time", "em_time"))


library(ggplot2)
ggplot(df3, aes(x = TT, y = value, group = variable)) +
  geom_line(aes(color = variable)) +
  geom_point() +
  facet_wrap(~N)


ggplot(df3, aes(x = N, y = value, group = variable)) +
  geom_line(aes(color = variable)) +
  geom_point() +
  facet_wrap(~TT)
