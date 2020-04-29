library(nnet) # For multinom
library(car) # For logit
library(sn) # For rmst
library(tidyverse)

set.seed(12) 

# Example configurations
n1 <- 1200 # sample size of treatment group 1
n2 <- 2400
n3 <- 4800
N <- n1 + n2 + n3

b <- 0.50 # initial covariate bias between treatment groups
eta <- 3.5
clustnum <- 5 # Number of clusters

# Function for calculating MaxMax2SB
match3.results <- function(SD, x1, x2, x3){
  
  max(abs((mean(x1)-mean(x2))/SD), abs((mean(x1)-mean(x3))/SD), abs((mean(x2)-mean(x3))/SD))
}

# Generate 18 covariates
mu1 <- as.matrix(c(b,0,0,b,0,0,b,0,0,b,0,0,b,0,0,b,0,0))
mu2 <- as.matrix(c(0,b,0,0,b,0,0,b,0,0,b,0,0,b,0,0,b,0))
mu3 <- as.matrix(c(0,0,b,0,0,b,0,0,b,0,0,b,0,0,b,0,0,b))
cov1 <- matrix(0, nrow = 18, ncol = 18); diag(cov1) <- 1
cov2 <- matrix(0, nrow = 18, ncol = 18); diag(cov2) <- 1
cov3 <- matrix(0, nrow = 18, ncol = 18); diag(cov3) <- 1

C1 <- rmst(n1, as.vector(mu1), cov1, rep(eta, 18), nu=7)
C2 <- rmst(n2, as.vector(mu2), cov2, rep(eta, 18), nu=7)
C3 <- rmst(n3, as.vector(mu3), cov3, rep(eta, 18), nu=7)
X <- rbind(C1, C2, C3)
covariates <- data.frame(X)
colnames(covariates) <- c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10",
                          "X11", "X12", "X13", "X14", "X15", "X16", "X17", "X18")

# Treatment indicator
treat <- as.matrix(rep(c(1,2,3), c(n1,n2,n3)))
treat <- factor(treat, levels = c(1, 2, 3), 
                labels = c("Treatment 1", "Treatment 2", "Treatment 3"))

data <- data.frame(treat, covariates)

data$T1 <- data$treat == "Treatment 1"
data$T2 <- data$treat == "Treatment 2"
data$T3 <- data$treat == "Treatment 3"

# Generate binary potential outcomes using logistic link
prob <- with(data, 
             1/(1 + exp(-(2*X1 + 4*X2 + 6*X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + 
                            X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18))))
Ypot1 <- rbinom(n1+n2+n3, 1, prob)
Ypot2 <- rbinom(n1+n2+n3, 1, prob)
Ypot3 <- rbinom(n1+n2+n3, 1, prob)

data <- cbind(Ypot1, Ypot2, Ypot3, data)

# Get observed outcomes
data$Yobs <- rep(NA, nrow(data))
data[data$treat == "Treatment 1", ]$Yobs <- data$Ypot1[1:n1]
data[data$treat == "Treatment 2", ]$Yobs <- data$Ypot2[(n1+1):(n1+n2)]
data[data$treat == "Treatment 3", ]$Yobs <- data$Ypot3[(n1+n2+1):(n1+n2+n3)]

# Fit GPS model
fit <- multinom(treat ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 +
                  X12 + X13 + X14 + X15 + X16 + X17 + X18, 
                data = data, trace = F)
Rx <- fitted(fit)
colnames(Rx) <- c("p1", "p2", "p3")
data <- cbind(data.frame(data), Rx)

# Determine eligibility
min.max.Ps <- data %>%
  group_by(treat) %>%
  summarise(min1 = min(p1), max1 = max(p1), 
            min2 = min(p2), max2 = max(p2), 
            min3 = min(p3), max3 = max(p3))

data$Eligible <- 
  data$p1 >= max(min.max.Ps$min1) & data$p1 <= min(min.max.Ps$max1) &
  data$p2 >= max(min.max.Ps$min2) & data$p2 <= min(min.max.Ps$max2) &
  data$p3 >= max(min.max.Ps$min3) & data$p3 <= min(min.max.Ps$max3) 

data <- filter(data, Eligible)
data <- data.frame(data)

# Re-fit GPS model
fit.E <- multinom(treat ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 +
                    X12 + X13 + X14 + X15 + X16 + X17 + X18, 
                  data = data, trace = F)
Rx.E <- fitted(fit.E) 
colnames(Rx.E) <- c("p1", "p2", "p3")
data <- data %>% 
  dplyr::select(-p1, -p2, -p3)
data <- cbind(data.frame(data), Rx.E)

covariates <- data %>%
  select(X1, X2, X3, X4, X5, X6, X7, X8, X9, 
         X10, X11, X12, X13, X14, X15, X16, X17, X18)

# MaxMax2SB before ABB
sds <- c()
for (i in 1:ncol(covariates)){
  sds[i] <- sd(covariates[data$treat == "Treatment 1", i], na.rm = T)
}

msbs.pre <- c()
for (i in 1:ncol(covariates)){
  msbs.pre[i] <- match3.results(sds[i],
                                covariates[,i][which(data$treat == "Treatment 1")], 
                                covariates[,i][which(data$treat == "Treatment 2")], 
                                covariates[,i][which(data$treat == "Treatment 3")])
}
maxmax2sb.pre <- max(msbs.pre, na.rm = T)

log <- data.frame(cbind(logit(data$p1), logit(data$p2), logit(data$p3)))

temp12 <- filter(data, treat == "Treatment 1" | treat == "Treatment 2")
temp13 <- filter(data, treat == "Treatment 1" | treat == "Treatment 3")
temp23 <- filter(data, treat == "Treatment 2" | treat == "Treatment 3")

#########################
#########################

# Sample sizes of eligible units in each treatment group
n1 <- sum(data$treat == "Treatment 1")
n2 <- sum(data$treat == "Treatment 2")
n3 <- sum(data$treat == "Treatment 3")

# k means clustering with Q = 5
data$k.clust <- kmeans(cbind(logit(data$p1), logit(data$p2), logit(data$p3)), clustnum)$cluster

##### ABB Multiple Imputation #####

M <- 25 #number of imputations
data$id <- rownames(data)
d1 <- filter(data, treat == "Treatment 1")
data.imp <- d1 %>% select(treat, k.clust)
data.imp <- data.imp %>% mutate(Y1 = NA, Y2 = NA, Y3 = NA, id1 = NA, id2 = NA, id3 = NA)
imp.list <- list(data.imp, data.imp, data.imp, data.imp, data.imp,
                 data.imp, data.imp, data.imp, data.imp, data.imp,
                 data.imp, data.imp, data.imp, data.imp, data.imp,
                 data.imp, data.imp, data.imp, data.imp, data.imp,
                 data.imp, data.imp, data.imp, data.imp, data.imp)

weights <- with(data, table(treat, k.clust))[1,]/n1

maxmax2sb <- c()

for (i in 1:M){
  set.seed(i)
  
  msbs <- matrix(nrow = clustnum, ncol = 18)
  for (k in 1:clustnum){ 
    tmp <- filter(data, k.clust == k)
    n1.tmp <- sum(tmp$treat == "Treatment 1")
    n2.tmp <- sum(tmp$treat == "Treatment 2")
    n3.tmp <- sum(tmp$treat == "Treatment 3")
    
    donors2 <- as.numeric(sample(tmp$id[tmp$treat == "Treatment 2"], n2.tmp, replace = T))
    donors3 <- as.numeric(sample(tmp$id[tmp$treat == "Treatment 3"], n3.tmp, replace = T))
    
    imp.list[[i]]$id1[imp.list[[i]]$k.clust == k] <- as.numeric(tmp$id[tmp$treat == "Treatment 1"])
    imp.list[[i]]$id2[imp.list[[i]]$k.clust == k] <- as.numeric(sample(donors2, n1.tmp, replace = T))
    imp.list[[i]]$id3[imp.list[[i]]$k.clust == k] <- as.numeric(sample(donors3, n1.tmp, replace = T))
    
    imp.list[[i]]$Y2[imp.list[[i]]$k.clust == k] <- 
      data$Yobs[imp.list[[i]]$id2[imp.list[[i]]$k.clust == k]]
    imp.list[[i]]$Y3[imp.list[[i]]$k.clust == k] <- 
      data$Yobs[imp.list[[i]]$id3[imp.list[[i]]$k.clust == k]]
    
    triplets <- imp.list[[i]][imp.list[[i]]$k.clust == k,c("id1", "id2", "id3")]
    
    for (ii in 1:ncol(covariates)){
      msbs[k,ii] <- match3.results(sds[ii], 
                                   covariates[,ii][triplets[,1]], 
                                   covariates[,ii][triplets[,2]], 
                                   covariates[,ii][triplets[,3]])
    }
    msbs[k, ][msbs[k, ] == -Inf] <- NA
  }
  imp.list[[i]]$Y1 <- data$Yobs[1:n1]
  
  imp.list[[i]]$tau12.hat <- with(imp.list[[i]], Y1 - Y2)
  imp.list[[i]]$tau13.hat <- with(imp.list[[i]], Y1 - Y3)
  
  rm(tmp, n2.tmp, n3.tmp)
  
  msbs.post <- apply(msbs, 2, function(x) sum(x*weights, na.rm = T))
  
  maxmax2sb[i] <- max(msbs.post, na.rm = T)
}


tau12.hat <- c()
tau13.hat <- c()
Y1.hat.imp.vec <- c()
Y2.hat.imp.vec <- c()
Y3.hat.imp.vec <- c()
W12 <- c()
W13 <- c()
for (i in 1:length(imp.list)){
  imp.list[[i]] <- filter(imp.list[[i]], treat == "Treatment 1")
  tau12.hat[i] <- mean(imp.list[[i]]$tau12.hat)
  tau13.hat[i] <- mean(imp.list[[i]]$tau13.hat)
  
  Y1.hat.imp.vec[i] <- mean(data.frame(imp.list[[i]])$Y1)
  Y2.hat.imp.vec[i] <- mean(data.frame(imp.list[[i]])$Y2)
  Y3.hat.imp.vec[i] <- mean(data.frame(imp.list[[i]])$Y3)
  
  W12[i] <- (Y1.hat.imp.vec[i]*(1-Y1.hat.imp.vec[i]))/nrow(imp.list[[i]]) + 
    (Y2.hat.imp.vec[i]*(1-Y2.hat.imp.vec[i]))/nrow(imp.list[[i]])
  W13[i] <- (Y1.hat.imp.vec[i]*(1-Y1.hat.imp.vec[i]))/nrow(imp.list[[i]]) + 
    (Y3.hat.imp.vec[i]*(1-Y3.hat.imp.vec[i]))/nrow(imp.list[[i]])
  
}
tau1.hat.imp <- c(mean(tau12.hat), mean(tau13.hat))
W12.hat <- mean(W12)
W13.hat <- mean(W13)

# Between-imputation variance
K <- 25
B12 <- var(tau12.hat)
B13 <- var(tau13.hat)

# Total variance
V12 <- W12.hat + ((K+1)/K)*B12
V13 <- W13.hat + ((K+1)/K)*B13

# Confidence intervals
LB1 <- tau1.hat.imp[1] - 1.96*sqrt(V12)
UB1 <- tau1.hat.imp[1] + 1.96*sqrt(V12)
LB2 <- tau1.hat.imp[2] - 1.96*sqrt(V13)
UB2 <- tau1.hat.imp[2] + 1.96*sqrt(V13)

imp.summary <- data.frame(rbind(c(tau1.hat.imp[1], LB1, UB1, sqrt(V12)), 
                                c(tau1.hat.imp[2], LB2, UB2, sqrt(V13))))
names(imp.summary) <- c("Estimate", "95% Lower", "95% Upper", "SE")

# Summary stats
imp.summary

# Comparing MaxMax2SB before/after ABB
c(maxmax2sb.pre, mean(maxmax2sb))

