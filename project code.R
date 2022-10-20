library("coda")
library("survival")
library("fields")
library("BayesX")
library("R2BayesX")
library(spBayesSurv)
data(LeukSurv)
d <- LeukSurv[order(LeukSurv$district), ]
summary(d)
nwengland <- read.bnd(system.file("otherdata/nwengland.bnd", package = "spBayesSurv"))
adj.mat <- bnd2gra(nwengland)
E <- diag(diag(adj.mat)) - as.matrix(adj.mat)
set.seed(1)
mcmc <- list(nburn = 5000, nsave = 2000, nskip = 4, ndisplay = 1000)
prior <- list(maxL = 4, a0 = 5, b0 = 1)
ptm <- proc.time()
####################################################
#PO model with ICAR frailties
########################################################################
set.seed(1)
mcmc <- list(nburn = 5000, nsave = 2000, nskip = 4, ndisplay = 1000)
prior <- list(maxL = 15)
ptm <- proc.time()
res1 <- survregbayes(formula = Surv(time, cens) ~ age + sex + wbc + tpi + frailtyprior("car", district), data = d, survmodel = "PO", dist = "loglogistic", mcmc = mcmc, prior = prior, Proximity = E)
proc.time() - ptm
prior <- list(maxL = 15, beta0 = rep(0, 3), S0 = diag(10, 3),theta0 = rep(0, 2), V0 = diag(10, 2), a0 = 1, b0 = 1, taua0 = 1, taub0 = 1)
(sfit1 <- summary(res1))
###################################################
#PO model with GRF fralties
##################################
set.seed(1)
mcmc <- list(nburn = 5000, nsave = 2000, nskip = 4, ndisplay = 1000)
prior <- list(maxL = 15, nu = 1, nknots = 100, nblock = 1043)
d$ID <- 1:nrow(d)
locations <- cbind(d$xcoord, d$ycoord);
ptm <- proc.time()
res2 <- survregbayes(formula = Surv(time, cens) ~ age + sex + wbc + tpi +frailtyprior("grf", ID), data = d, survmodel = "PO",dist = "loglogistic", mcmc = mcmc, prior = prior,Coordinates = locations)
proc.time() - ptm
(sfit11 <- summary(res11))
################################################################
#PO model with ICAR frailties with parametric distributio
#############################
set.seed(1)
prior <- list(maxL = 15, a0 = -1, thete0 = rep(0, 2), V0 = diag(1e10, 2))
state <- list(alpha = Inf)
ptm <- proc.time()
res11 <- survregbayes(formula = Surv(time, cens) ~ age + sex + wbc + tpi +
                           + frailtyprior("car", district), data = d, survmodel = "PO",
                         + dist = "loglogistic", mcmc = mcmc, prior = prior, state = state,
                         + Proximity = E, InitParamMCMC = FALSE)
proc.time() - ptm
(sfit11 <- summary(res11))

############################################################


res1 <- survregbayes(formula = Surv(time, cens) ~ age + sex + wbc + tpi +
                       baseline(age, sex, wbc, tpi),
                     data = d, mcmc = mcmc, prior = prior, survmodel= "PO", Proximity = E)
prior <- list(maxL = 15, beta0 = rep(0, 3), S0 = diag(10, 3),theta0 = rep(0, 2), V0 = diag(10, 2), a0 = 1, b0 = 1,taua0 = 1, taub0 = 1)


#res1 <- survregbayes(formula = Surv(time, cens) ~ age + sex + wbc + tpi +frailtyprior("car", district), data = d, mcmc = mcmc, prior = prior, Proximity = E)

(sfit1 <- summary(res1))

##################################################################

#WEISBULL MODEL PH


############################################################


#res1 <- survregbayes(formula = Surv(time, cens) ~ age + sex + wbc + tpi +
                       #baseline(age, sex, wbc, tpi)
                     #$data = d, mcmc = mcmc, prior = prior, survmodel= "PH",dist="loglogistic", Proximity = E)
#prior <- list(maxL = 15, beta0 = rep(0, 3), S0 = diag(10, 3),theta0 = rep(0, 2), V0 = diag(10, 2), a0 = 1, b0 = 1,taua0 = 1, taub0 = 1)


#res1 <- survregbayes(formula = Surv(time, cens) ~ age + sex + wbc + tpi +frailtyprior("car", district), data = d, mcmc = mcmc, prior = prior, Proximity = E)

#(sfit1 <- summary(res1))





res1 <- survregbayes(formula = Surv(time, cens) ~ age + sex + wbc + tpi +
                       frailtyprior("car", district), data = d, survmodel = "PO",
                     dist = "loglogistic", mcmc = mcmc, prior = prior, Proximity = E)
prior <- list(maxL = 15, beta0 = rep(0, 3), S0 = diag(10, 3),theta0 = rep(0, 2), V0 = diag(10, 2), a0 = 1, b0 = 1, taua0 = 1, taub0 = 1)
(sfit1 <- summary(res1))
proc.time() - ptm

##################################################################

#WEISBULL MODEL AFT


############################################################


res1 <- survregbayes(formula = Surv(time, cens) ~ age + sex + wbc + tpi +
                       baseline(age, sex, wbc, tpi) + frailtyprior("car", district),
                     data = d, mcmc = mcmc, prior = prior, survmodel= "PH", Proximity = E)
prior <- list(maxL = 15, beta0 = rep(0, 3), S0 = diag(10, 3),theta0 = rep(0, 2), V0 = diag(10, 2), a0 = 1, b0 = 1,taua0 = 1, taub0 = 1)


#res1 <- survregbayes(formula = Surv(time, cens) ~ age + sex + wbc + tpi +frailtyprior("car", district), data = d, mcmc = mcmc, prior = prior, Proximity = E)

(sfit1 <- summary(res1))

##################################################################


##################################################################

#WEISBULL MODEL 


############################################################


res1 <- survreg(formula = Surv(time, cens) ~ age + sex + wbc + tpi +
                       baseline(age, sex, wbc, tpi),
                     data = d,  dist = "weibull")
prior <- list(maxL = 15, beta0 = rep(0, 3), S0 = diag(10, 3),theta0 = rep(0, 2), V0 = diag(10, 2), a0 = 1, b0 = 1,taua0 = 1, taub0 = 1)


#res1 <- survregbayes(formula = Surv(time, cens) ~ age + sex + wbc + tpi +frailtyprior("car", district), data = d, mcmc = mcmc, prior = prior, Proximity = E)

(sfit1 <- summary(res1))

##################################################################
frail0 <- apply(res1$v, 1, mean)
frail <- frail0[as.integer(names(nwengland))]
values <- cbind(as.integer(names(nwengland)), frail)
op <- par(no.readonly = TRUE)
par(mar = c(3, 0, 0, 0))
plotmap(nwengland, x = values, col = (gray.colors(10, 0.3, 1))[10:1],
        pos = "bottomleft", width = 0.5, height = 0.04)

##############################################################
par(mfrow = c(3, 2))
par(cex = 1, mar = c(2.5, 4.1,1, 1))
traceplot(mcmc(res1$beta[1,]),xlab = "", main = "age")
traceplot(mcmc(res1$beta[2,]),xlab = "", main = "sex")
traceplot(mcmc(res1$beta[3,]),xlab = "", main = "wbc")
traceplot(mcmc(res1$beta[4,]), xlab = "", main = "tpi")
traceplot(mcmc(res1$tau2), xlab = "", main = "tau^2")

