# ---- Purpose ----
# Unpack results of a simulation run on the IRE.
# Ran 200 replications of 3-dimensional series.

load("out-IRE-20201117.Rdata")

# --- True Params ----

Sig1 <- diag(N)
Sig1[1,2] <- Sig1[2,1] <- .75
Sig.true <- list(Sig1, diag(N), diag(N))

# ---- MLE sigma ----
par.mle <- sigex.psi2par(psi.mle, mdl, data.ts)
Sig.mle = param2sig(param)


# ---- Likelihood compare ----

Lik <- data.frame(
    tru = rep(NA, 200),
    mom = rep(NA, 200),
    mle = rep(NA, 200),
    em  = rep(NA, 200)
)

for(i in 1:200){
  Lik$tru[i] <- out[[i]]$lik.true
  Lik$mom[i] <- out[[i]]$lik.mom
  Lik$mle[i] <- out[[i]]$lik.mle
  Lik$em[i]  <- out[[i]]$Sig.save[52, 28]
}

plot(x = c(1, 200),
     y = range(Lik[, 3:4]),
     type = "n",
     xlab = "",
     ylab = "likelihood")
points(Lik$tru, pch = 19)
points(Lik$mom, pch = 19, col = 2)
points(Lik$mle, pch = 19, col = 3)
points(Lik$em,  pch = 19, col = 4)

cbind(mleLik, emLik, emLik - mleLik)

plot_function(out[[189]]$Sig.save)

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
  geom_boxplot(data = df, aes(x = method, y = likelihood, fill = method ))



