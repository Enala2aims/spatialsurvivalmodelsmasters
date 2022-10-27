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
(sfit11 <- summary(res2))
################################################################
#PO model with ICAR frailties with parametric distributio
#############################
set.seed(1)
prior <- list(maxL = 15, a0 = -1, thete0 = rep(0, 2), V0 = diag(1e10, 2))
state <- list(alpha = Inf)
ptm <- proc.time()
res11 <- survregbayes(formula = Surv(time, cens) ~ age + sex + wbc + tpi +frailtyprior("car", district), data = d, survmodel = "PO",dist = "loglogistic", mcmc = mcmc, prior = prior, state = state,Proximity = E, InitParamMCMC = FALSE)
proc.time() - ptm
(sfit11 <- summary(res11))
#################################################
#PH model
################################################
set.seed(1)
mcmc <- list(nburn = 5000, nsave = 2000, nskip = 4, ndisplay = 1000)
prior <- list(maxL = 15)
ptm <- proc.time()
res12 <- survregbayes(formula = Surv(time, cens) ~ age + sex + wbc + tpi + frailtyprior("car", district), data = d, survmodel = "PH", dist = "loglogistic", mcmc = mcmc, prior = prior, Proximity = E)
proc.time() - ptm
prior <- list(maxL = 15, beta0 = rep(0, 3), S0 = diag(10, 3),theta0 = rep(0, 2), V0 = diag(10, 2), a0 = 1, b0 = 1, taua0 = 1, taub0 = 1)
(sfit11 <- summary(res12))
############################################################
#variable selection
#############################################################
set.seed(1)
res3 <- survregbayes(formula = Surv(time, cens) ~ age + sex + wbc + tpi +
                          frailtyprior("car", district), data = d, survmodel = "PO",
                        dist = "loglogistic", mcmc = mcmc, prior = prior, Proximity = E,
                        selection = TRUE)
(sfit3 <- summary(res3))
##################################################################
#PLOTTING MEAN FRAILTIES WITH VALUES
###################################################################
frail <- round(apply(res2$v, 1, mean), 3)
nclust <- 5
frail.cluster <- cut(frail, breaks = nclust)
frail.names <- names(table(frail.cluster))
rbPal <- colorRampPalette(c('blue', 'red'))
frail.colors <- rbPal(nclust)[as.numeric(frail.cluster)]
par(mar = c(3, 0, 0, 0))
plot(nwengland)
points(cbind(d$xcoord,d$ycoord), col = frail.colors)
legend("topright", title = "frailty values", legend = frail.names,
       col = rbPal(nclust), pch = 20, cex = 1.7)

##############################################################
#PLOTTING SURVIVAL FUNCTIONS
##############################################################
summary(res2)
tgrid <- seq(0.1, 180, length.out = 300);
xpred <- data.frame(age = c(50, 50,50), sex = c(1, 1,1), wbc = c(38.59, 38.59,38.59), tpi = c(-3, 0.38, 8),row.names = c("tpi=-3", "tpi=0.38", "tpi= 8"))
plot(res2, xnewdata = xpred, tgrid = tgrid, cex = 1, )