## =============================================================================
## LAB CODE: Analysis of Krackhardt Advice Network
## =============================================================================

library(ergm)
setwd('~/Tresorit/Teaching/Networks\ Class/Chapter_3_Basic_ERGM/Lab\ Files/')
load("Krackhardt.RData")
Krack <- Krackhardt
rm(Net)


## -----------------------------------------------------------------------------
## Descriptive statistics

Krack
network.density(Krack)
summary(Krack ~ idegree(1:10) + odegree(1:10) + istar(2:3)
        + ostar(2:3) + triangles)


## -----------------------------------------------------------------------------
## Obligatory plots

## Basic plot
set.seed(510)
plot(Krack, edge.col="gray", label.col="black", vertex.cex=1.5,
     vertex.col="red",
     main="Krackhardt Manager Advice Network")

## Plot by department
set.seed(510)
dept <- get.vertex.attribute(Krack, "Department")
cols <- c("red", "green", "blue", "black", "orange")
col.dept <- cols[match(dept, unique(dept))]
plot(Krack, edge.col="gray", label.col="black", vertex.cex=1.5,
     vertex.col=col.dept, label=dept,
     main="Krackhardt Manager Advice Network: Department")

## Plot by Level
set.seed(510)
level <- get.vertex.attribute(Krack, "Level")
cols <- c("black", "red", "blue")
col.level <- cols[match(level, unique(level))]
plot(Krack, edge.col="gray", label.col="black", vertex.cex=1.5,
     vertex.col=col.level, label=level,
     main="Krackhardt Manager Advice Network: Level")

## Plot by tenure
set.seed(510)
Pal <- colorRampPalette(c("blue", "red"))
tenure <- get.vertex.attribute(Krack, "Tenure")
col.tenure <- Pal(10)[as.numeric(cut(tenure, breaks = 10))]
plot(Krack, edge.col="gray", label.col="black", vertex.cex=1.5,
     vertex.col=col.tenure, label=tenure,
     main="Krackhardt Manager Advice Network: Tenure")

## Plot by age
set.seed(510)
Pal <- colorRampPalette(c("blue", "red"))
age <- get.vertex.attribute(Krack, "Age")
col.age <- Pal(10)[as.numeric(cut(age, breaks = 10))]
plot(Krack, edge.col="gray", label.col="black", vertex.cex=1.5,
     vertex.col=col.tenure, label=age,
     main="Krackhardt Manager Advice Network: Age")


## -----------------------------------------------------------------------------
## Exponential Random Graph Model explaining network formation

?ergm.terms

## Covariate model without network structural parameters.
mod0 <- ergm(Krack ~ edges
             + edgecov("reportsto")
             + nodeicov("Tenure") + nodeocov("Tenure") + absdiff("Tenure")
             + nodeicov("Age") + nodeocov("Age") + absdiff("Age"),
             verbose=TRUE)
gof0 <- gof(mod0, verbose=TRUE)
par(mfrow=c(2,2)); plot(gof0)
summary(mod0)

## Model with reciprocity and covariates. NOTE: Burnin and sample sizes are too
## low to get a good est. Increase by a factor of 10 to get reliable estimates.
set.seed(5)
mod1 <- ergm(Krack ~ edges + mutual
             + edgecov("reportsto")
             + nodeicov("Tenure") + nodeocov("Tenure") + absdiff("Tenure")
             + nodeicov("Age") + nodeocov("Age") + absdiff("Age"),
             control=control.ergm(
                 MCMC.samplesize=5000,
                 MCMC.burnin=10000,
                 MCMLE.maxit=10),
             verbose=TRUE)
mcmc.diagnostics(mod1, vars.per.page=5)
plot(mod1$sample[[1]])                       # using coda instead
gof1 <- gof(mod1, verbose=TRUE)
par(mfrow=c(2,2)); plot(gof1)
summary(mod1)

## Add out-stars
set.seed(10)
mod2 <- ergm(Krack ~ edges + mutual + ostar(2:3)
             + edgecov("reportsto")
             + nodeicov("Tenure") + nodeocov("Tenure") + absdiff("Tenure")
             + nodeicov("Age") + nodeocov("Age") + absdiff("Age"),
             control=control.ergm(
                 MCMC.samplesize=5000,
                 MCMC.burnin=10000,
                 MCMLE.maxit=10),
             verbose=TRUE)
mcmc.diagnostics(mod2)
plot(mod1$sample[[1]])                       # using coda instead
gof2 <- gof(mod2, verbose=TRUE)
par(mfrow=c(2,2)); plot(gof2)
summary(mod2)

## Add transitivity ; Degenerate
set.seed(15)
mod3 <- ergm(Krack ~ edges + mutual + ostar(2:3) + transitive
             + edgecov("reportsto")
             + nodeicov("Tenure") + nodeocov("Tenure") + absdiff("Tenure")
             + nodeicov("Age") + nodeocov("Age") + absdiff("Age"),
             control=control.ergm(
                 MCMC.samplesize=5000,
                 MCMC.burnin=10000,
                 MCMLE.maxit=10),
             verbose=TRUE)
mcmc.diagnostics(mod3)
gof3 <- gof(mod3, verbose=TRUE)
par(mfrow=c(2,2)); plot(gof3)
summary(mod3)

## Use GWESP
set.seed(20)
mod4 <- ergm(Krack ~ edges + mutual + ostar(2:3) + gwesp(0, fixed=TRUE)
             + edgecov("reportsto")
             + nodeicov("Tenure") + nodeocov("Tenure") + absdiff("Tenure")
             + nodeicov("Age") + nodeocov("Age") + absdiff("Age"),
             control=control.ergm(
                 MCMC.samplesize=5000,
                 MCMC.burnin=10000,
                 MCMLE.maxit=10),
             verbose=TRUE)
mcmc.diagnostics(mod4)
gof4 <- gof(mod4, verbose=TRUE)
par(mfrow=c(2,2)); plot(gof4)
summary(mod4)

## With GWEST, estimate decay parameter (curved exponential family)
set.seed(25)
mod5 <- ergm(Krack ~ edges + mutual + ostar(2:3) + gwesp(0, fixed=FALSE)
             + edgecov("reportsto")
             + nodeicov("Tenure") + nodeocov("Tenure") + absdiff("Tenure")
             + nodeicov("Age") + nodeocov("Age") + absdiff("Age"),
             control=control.ergm(
                 MCMC.samplesize=50000,
                 MCMC.burnin=100000,
                 MCMLE.maxit=10),
             verbose=TRUE)
mcmc.diagnostics(mod5)
gof5 <- gof(mod5, verbose=TRUE)
par(mfrow=c(2,2)); plot(gof5)
summary(mod5)
