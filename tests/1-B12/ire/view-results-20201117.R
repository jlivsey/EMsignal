# ---- Purpose ----
# Unpack results of a simulation run on the IRE.
# Ran 200 replications of 3-dimensional series.

library(EMsigex)
library(sigex)

load("out-IRE-20201117.Rdata")
load("dataList.Rdata")

# --- True Params ----

Sig1 <- diag(N)
Sig1[1,2] <- Sig1[2,1] <- .75
Sig.true <- list(Sig1, diag(N), diag(N))

# ---- MLE sigma ----
N <- dim(dataList[[1]])[2]
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"trend",c(1,-1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"seasonal", rep(1,12))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"irregular",1)
mdl <- sigex.meaninit(mdl, dataList[[1]], 0) # regressors

Sig.mle <- matrix(NA, nrow = 200, ncol = 27) # storage
for(i in 1:200){
  param = sigex.psi2par(out[[i]]$psi.mle, mdl, dataList[[i]])
  Sig.mle[i, ] = unlist(param2sig(param))
}

# pdf(file = "MLEestimates.pdf", width = 15, height = 10)
boxplot(Sig.mle, main = "MLE estimates", xaxt = 'n')
axis(1,
     at=1:27,
     labels=c('T11','T21','T31','T12','T22','T32','T13','T23','T33',
              'S11','S21','S31','S12','S22','S32','S13','S23','S33',
              'I11','I21','I31','I12','I22','I32','I13','I23','I33'))
# Add true value as red horizonal line
{
  n <- ncol(Sig.em)
  # width of each boxplot is 0.8
  x0s <- 1:n - 0.4
  x1s <- 1:n + 0.4
  # these are the y-coordinates for the horizontal lines
  # that you need to set to the desired values.
  y0s <- unlist(Sig.true)
  # add segments
  segments(x0 = x0s, x1 = x1s, y0 = y0s, col = "red")
}
# dev.off()

# ---- EM sigma ----
Sig.em <- matrix(NA, nrow = 200, ncol = 27) # storage
for(i in 1:200){
  Sig.em[i, ] = out[[i]]$Sig.save[52, 1:27]
}

# pdf(file = "EMestimates.pdf", width = 15, height = 10)
boxplot(Sig.em, main = "EM estimates", xaxt = 'n')
axis(1,
     at=1:27,
     labels=c('T11','T21','T31','T12','T22','T32','T13','T23','T33',
              'S11','S21','S31','S12','S22','S32','S13','S23','S33',
              'I11','I21','I31','I12','I22','I32','I13','I23','I33'))
# Add true value as red horizonal line
{
  n <- ncol(Sig.em)
  # width of each boxplot is 0.8
  x0s <- 1:n - 0.4
  x1s <- 1:n + 0.4
  # these are the y-coordinates for the horizontal lines
  # that you need to set to the desired values.
  y0s <- unlist(Sig.true)
  # add segments
  segments(x0 = x0s, x1 = x1s, y0 = y0s, col = "red")
}
# dev.off()

# ---- Likelihood compare ----
# storage
Lik <- data.frame(
    tru = rep(NA, 200),
    mom = rep(NA, 200),
    mle = rep(NA, 200),
    em  = rep(NA, 200)
)
# populate dataframe with results
for(i in 1:200){
  Lik$tru[i] <- out[[i]]$lik.true
  Lik$mom[i] <- out[[i]]$lik.mom
  Lik$mle[i] <- out[[i]]$lik.mle
  Lik$em[i]  <- out[[i]]$Sig.save[52, 28]
}
# scatter plot
plot(x = c(1, 200),
     y = range(Lik[, 3:4]),
     type = "n",
     xlab = "",
     ylab = "likelihood")
points(Lik$tru, pch = 19)
points(Lik$mom, pch = 19, col = 2)
points(Lik$mle, pch = 19, col = 3)
points(Lik$em,  pch = 19, col = 4)
# print MLE vs EM
attach(Lik)
cbind(mle, em, em - mle)
# boxplot
boxplot(Lik)


# ---- Plot parameter estimates over time -------------------------------------
library(ggplot2)
library(dplyr)
library(gridExtra)

plot_function <- function(Sig.save){
  lik.true <- Sig.save[1, 28]
  dat <- data.frame(Sig.save)
  dat$iter <- 1:dim(dat)[1]
  colnames(dat) <- c('t11','t21', 't31','t12', 't22', 't32', 't13', 't23', 't33',
                     's11','s21', 's31','s12', 's22', 's32', 's13', 's23', 's33',
                     'i11','i21', 'i31','i12', 'i22', 'i32', 'i13', 'i23', 'i33',
                     'lik', 'iter')
  # Plot Likelihood
  dat.gath <- dat %>%
    select(lik, iter) %>%
    tidyr::gather('variable', 'value', -iter)
  gg_lik <-
    ggplot(dat.gath, aes(x = iter, y = value, colour = variable)) +
    geom_line() +
    geom_hline(yintercept = lik.true)
  # Plot all parameters
  dat.gath <- dat %>%
    select(-lik) %>%
    tidyr::gather('variable', 'value', -iter)
  gg_param <-
    ggplot(dat.gath, aes(x = iter, y = value, colour = variable)) +
    geom_line()
  gridExtra::grid.arrange(gg_param, gg_lik, ncol = 2)
}



# ---- Box-plots ----
# Attempt to make nice looking boxplots for the paper
library(ggplot2)
library(dplyr)

Lik_names <-
  Lik %>%
  rename(True = tru) %>%
  rename(MOM  = mom) %>%
  rename(MLE  = mle) %>%
  rename(EM   = em)


# Reshape data to fit ggplot framework
df = reshape2::melt(Lik_names,
                    measure.vars = c('True', 'MOM', 'MLE', 'EM'),
                    variable.name = "method",
                    value.name = "likelihood")

# remove any MOM estimates over 3000 (makes plot more condensed)
df <-
  df %>%
  filter(likelihood < 3000)

## Plot it in ggplot
ggplot() +
  geom_boxplot(data = df, aes(x = method, y = likelihood)) +
  theme_gray(base_size = 14) +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=14))

ggsave(filename = "estimMethodBoxplot.pdf")


# ---- Sig parameters comparison ----

trendVar <- data.frame(
  mle = rep(NA, 200),
  em  = rep(NA, 200)
)
























