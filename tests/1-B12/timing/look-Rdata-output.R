
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
library(tidyr)
library(reshape2)
library(ggplot2)

# data.frame for plotting
df_ggplot <- df %>%
  rename(MLE_time = mle_time) %>%
  rename(EM_time  = em_time) %>%
  drop_na() %>%
  melt(measure.vars = c("MLE_time", "EM_time")) %>%
  mutate(N = as.factor(N)) %>%
  mutate(T = as.factor(TT))

# Change levels for plotting
levels(df_ggplot$N)
levels(df_ggplot$N) <- c("N=2", "N=3", "N=4", "N=5")

# Change levles for plotting
levels(df_ggplot$T)
levels(df_ggplot$T) <- c("T=50", "T=100", "T=200", "T=400", "T=800")


ggplot(df_ggplot, aes(x = TT, y = value, group = variable)) +
  geom_line(aes(color = variable)) +
  geom_point() +
  facet_wrap(~N)

ggplot(df_ggplot, aes(x = N, y = value, group = variable)) +
  geom_line(aes(color = variable)) +
  geom_point() +
  facet_wrap(~TT)
