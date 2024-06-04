# TERGM analysis of "Golden Age of Hollywood" Data

# Data referece: "Eigenvector-Based Centrality Measures for Temporal Networks" Dane Taylor, Sean A. Myers, Aaron Clauset, Mason A. Porter, Peter J. Mucha. Preprint, arXiv:1507.01266 (2015)


# Load up
library(statnet)
#setwd("~/Dropbox/writings/Networks Book/Data/Golden Age Of Hollywood/") # SJC

t0 <- as.matrix(read.table("HollywoodGoldenAge_matrix_1909-1919_s2.txt", header=FALSE))
t0[t0 > 0] <- 1 # No real structure here... too sparce... won't use


attributes <- read.table("HollywoodGoldenAge_names.txt", header=TRUE)

t1 <- as.matrix(read.table("HollywoodGoldenAge_matrix_1920-1929_s2.txt", header=FALSE)) # read in the adjacency matrix
t1[t1 > 0] <- 1 # replace edge weights with binary values
n1 <- network(t1) # create object of class "network"
network.vertex.names(n1) <- as.character(attributes$name) # add vertex names
set.vertex.attribute(n1, "female", attributes$female) # add sex attribute of vertices

t2 <- as.matrix(read.table("HollywoodGoldenAge_matrix_1930-1939_s2.txt", header=FALSE))
t2[t2 > 0] <- 1
n2 <- network(t2)
network.vertex.names(n2) <- as.character(attributes$name)
set.vertex.attribute(n2, "female", attributes$female)

t3 <- as.matrix(read.table("HollywoodGoldenAge_matrix_1940-1949_s2.txt", header=FALSE))
t3[t3 > 0] <- 1
n3 <- network(t3)
network.vertex.names(n3) <- as.character(attributes$name)
set.vertex.attribute(n3, "female", attributes$female)

t4 <- as.matrix(read.table("HollywoodGoldenAge_matrix_1950-1959_s2.txt", header=FALSE))
t4[t4 > 0] <- 1
n4 <- network(t4)
network.vertex.names(n4) <- as.character(attributes$name)
set.vertex.attribute(n4, "female", attributes$female)

t5 <- as.matrix(read.table("HollywoodGoldenAge_matrix_1960-1969_s2.txt", header=FALSE))
t5[t5 > 0] <- 1
n5 <- network(t5)
network.vertex.names(n5) <- as.character(attributes$name)
set.vertex.attribute(n5, "female", attributes$female)

t6 <- as.matrix(read.table("HollywoodGoldenAge_matrix_1970-1979_s2.txt", header=FALSE))
t6[t6 > 0] <- 1
n6 <- network(t6)
network.vertex.names(n6) <- as.character(attributes$name)
set.vertex.attribute(n6, "female", attributes$female)



hga <- list(n1, n2, n3, n4, n5, n6) # create list from oldest to newest. This will be the outcome object for TERGM analysis


t7 <- as.matrix(read.table("HollywoodGoldenAge_matrix_1980-1989_s2.txt", header=FALSE)) # not quite enough structure here, won't use
t7[t7 > 0] <- 1
n7 <- network(t7)
network.vertex.names(n7) <- as.character(attributes$name)
set.vertex.attribute(n7, "female", attributes$female)


save(list=c("attributes","n1","n2","n3","n4","n5","n6","n7"),file = "gaoh.RData")

library(ina)
data(gaoh)

hga <- list(n1, n2, n3, n4, n5, n6) # create list from oldest to newest. This will be the outcome object for TERGM analysis


## Plot the network
set.seed(5)
par(mfrow = c(3,2))
hgat <- c("20's", "30's", "40's", "50's", "60's", "70's")  
for (i in 1:length(hga)){
  plot(hga[[i]],displaylabels=F,label=network.vertex.names(n1),vertex.cex=2,label.cex=1,edge.col=rgb(150,150,150,100,maxColorValue=255),label.pos=5,vertex.col=c("lightblue", "pink")[get.vertex.attribute(n1,"female")+1], main=hgat[[i]]    )
}

# TERGM Analysis

library(bergm)

# Basic model
set.seed(12345)
m1 <- btergm(hga ~ edges 
             + mutual 
             #+ ttriple 
             #+ transitiveties 
             #+ ctriple 
             + gwesp(0.5, fixed=TRUE)
             #+ cycle(4)
             + idegreepopularity
             #+ odegreepopularity
             #+ istar(2)
             #+ ostar(2)
             #+ nodematch("female", diff=TRUE) 
             + absdiff("female")
             #+ nodemix("female", base=c(1,4))
             + nodefactor("female")
             + delrecip
             + memory(type="stability")
             , R=100)
summary(m1)
gof1 <- gof(m1, statistics = c(esp, dsp, geodesic,deg, triad.undirected, walktrap.modularity))
plot(gof1)

# Same model, but condensed code
set.seed(12345)
m1 <- btergm(hga ~ edges + mutual + gwesp(0.5, fixed=TRUE) + idegreepopularity + absdiff("female") + nodefactor("female") + delrecip + memory(type="stability"), R=100)
summary(m1)
gof1 <- gof(m1, statistics = c(esp, dsp, geodesic,deg, triad.undirected, walktrap.modularity))
plot(gof1)



# Same model estimated by MCMC
set.seed(5)
m1.5 <- m1 <- btergm(hga ~ edges + mutual + gwesp(0.5, fixed=TRUE) + idegreepopularity + absdiff("female") + nodefactor("female") + delrecip + memory(type="stability"))
summary(m1.5)


# Out of sample prediction
set.seed(12345)
m1oos <- btergm(hga[2:4] ~ edges + mutual + gwesp(0.5, fixed=TRUE) + idegreepopularity + absdiff("female") + nodefactor("female") + delrecip + memory(type="stability"), R=100)
gof1oos <- gof(m1oos, target=hga[[5]], statistics = c(esp, dsp, geodesic,deg, triad.undirected, walktrap.modularity))
plot(gof1oos)

interpret(m1, type="dyad", i=1, j=2, t=3)


dyads <- edgeprob(m1) # also errors
checkdegeneracy(m1, nsim=1000) # produces error "object 'mat' not found"


