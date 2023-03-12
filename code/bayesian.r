#### load libraries ####
library(nimble, warn.conflicts = FALSE, quietly = TRUE)
library(dplyr)
library(ggplot2)
library(mcmcplots)
library(coda)

#### Data Preparation ####
df = as.data.frame(matrix(nrow=6,ncol=2))

colnames(df) = c("reactive", "tested")
rownames(df) = c("Montreal-Laval-prev", "surrounding-prev", "other-prev", "Montreal-Laval", "surrounding", "other")
region = rownames(df)

df[1,] = c(90, 3061)
df[2,] = c(48, 1925)
df[3,] = c(35, 2705)
df[4,] = c(215, 1578)
df[5,] = c(128, 1422)
df[6,] = c(372, 4304)

N <- dim(df)[1]

prevData <- list(Xi = df$reactive)
prevConst <- list(N = N,ni = df$tested)

#### starting values ####
prevInits   <- list(Se = 0.88, Sp = 0.98, Sr = 0.29,
                    Tpi = rep(0.5, prevConst$N))

##########  Non-Informative Model  ##########
prevFlat  <- nimbleCode({
  Se ~ dbeta(205, 29)
  Sp ~ dbeta(288, 2)
  Sr ~ dbeta(32, 109-32)
  
  for (i in 1:(N/2)){
    Tpi[i] ~ dbeta(1, 1)
    Opi[i] <- Tpi[i]*Se + (1 - Tpi[i])*(1-Sp)
    Xi[i] ~ dbin(Opi[i], ni[i])
  }
  
  for (i in (N/2+1):N){
    Tpi[i] ~ T(dbeta(1, 1), Sr * Tpi[i-(N/2)], 1)
    Opi[i] <- (Tpi[i] - Sr * Tpi[i-(N/2)]) * Se + (1 - (Tpi[i] - Sr * Tpi[i-(N/2)])) * (1-Sp)
    Xi[i] ~ dbin(Opi[i], ni[i])
  }
})

burnin  <- 100000
modelF  <- nimbleMCMC(code = prevFlat,
                      constants = prevConst,
                      data = prevData,
                      inits = prevInits,
                      monitors = c("Se", "Sp", "Sr",
                                   paste('Tpi[',1:N,']', sep = '')),
                      nburnin = burnin, niter = 2*burnin, nchains = 4)

### (1) True Sero-prevalence:

## a. Convergence Diagnostics.
# Visual inspection
par(mfrow = c(2,3))
for(i in 1:N){
  Tp1 <- modelF$chain1[,paste('Tpi[',i,']', sep = '')] 
  Tp2 <- modelF$chain2[,paste('Tpi[',i,']', sep = '')]
  Tp3 <- modelF$chain3[,paste('Tpi[',i,']', sep = '')] 
  Tp4 <- modelF$chain4[,paste('Tpi[',i,']', sep = '')]
  plot(1:burnin, Tp1, type = 'l', xlab = 'b', ylab = region[i],
       xlim = c(0, 4*burnin))
  lines(burnin + 1:burnin, Tp2, col = 2)
  lines(2*burnin + 1:burnin, Tp3, col = 3)
  lines(3*burnin + 1:burnin, Tp4, col = 4)
}

# Gelman-Rubin diagnostic
TpF <- vector(mode="list", length=N)
GRF <- vector(mode="list", length=N)
for (i in 1:N){
  TpF[[i]] <- c(modelF$chain1[,paste('Tpi[',i,']', sep = '')], 
                 modelF$chain2[,paste('Tpi[',i,']', sep = '')],
                 modelF$chain3[,paste('Tpi[',i,']', sep = '')],
                 modelF$chain4[,paste('Tpi[',i,']', sep = '')])
  GRF[[i]] <- mcmc.list(as.mcmc(modelF$chain1[,paste('Tpi[',i,']', sep = '')]), 
                         as.mcmc(modelF$chain2[,paste('Tpi[',i,']', sep = '')]), 
                         as.mcmc(modelF$chain3[,paste('Tpi[',i,']', sep = '')]), 
                         as.mcmc(modelF$chain4[,paste('Tpi[',i,']', sep = '')]))
  print(gelman.diag(GRF[[i]]))
}

## b.Credible Intervals for Each Region
Tp_q_F <- vector(mode="list", length=N)

for (i in 1:N){
  Tp_q_F[[i]] <- quantile(TpF[[i]], probs = c(0.5, 0.025, 0.975))
}

tp_q_F <- data.frame(matrix(unlist(Tp_q_F), nrow=N, byrow=T), stringsAsFactors = FALSE)
tp_q_F$regions <- region
new_tp_q_F <- tp_q_F %>% 
  select(regions, everything())
names(new_tp_q_F) <- c("regions","50%", "2.5%", "97.5%")
new_tp_q_F

#             regions         50%         2.5%      97.5%
# 1 Montreal-Laval-prev 0.027866366 0.0171102175 0.03744966
# 2    surrounding-prev 0.022850311 0.0114875743 0.03362774
# 3          other-prev 0.009083538 0.0008942495 0.01697984
# 4      Montreal-Laval 0.159175601 0.1368821371 0.18297204
# 5         surrounding 0.104568252 0.0847029309 0.12536096
# 6               other 0.096026050 0.0819250807 0.10977501

### (2) Sensitivity,Specificity, and Seroreversion

SeF <- c(modelF$chain1[,'Se'], modelF$chain2[,'Se'],
          modelF$chain3[,'Se'], modelF$chain4[,'Se'])
SpF <- c(modelF$chain1[,'Sp'], modelF$chain2[,'Sp'],
          modelF$chain3[,'Sp'], modelF$chain4[,'Sp'])
SrF <- c(modelF$chain1[,'Sr'], modelF$chain2[,'Sr'],
          modelF$chain3[,'Sr'], modelF$chain4[,'Sr'])

## Convergence Diagnostics
#  Visual inspection
par(mfrow = c(1,3))

for(i in 1:3){
  if (i == 1){
    Tp1 <- modelF$chain1[,paste('Se', sep = '')]
    Tp2 <- modelF$chain2[,paste('Se', sep = '')]
    Tp3 <- modelF$chain3[,paste('Se', sep = '')] 
    Tp4 <- modelF$chain4[,paste('Se', sep = '')]
    plot(1:burnin, Tp1, type = 'l', xlab = 'b', ylab = "Sensitivity",
         xlim = c(0, 4*burnin))
    lines(burnin + 1:burnin, Tp2, col = 2)
    lines(2*burnin + 1:burnin, Tp3, col = 3)
    lines(3*burnin + 1:burnin, Tp4, col = 4) 
  } else if (i == 2){
    Tp1 <- modelF$chain1[,paste('Sp', sep = '')]
    Tp2 <- modelF$chain2[,paste('Sp', sep = '')]
    Tp3 <- modelF$chain3[,paste('Sp', sep = '')] 
    Tp4 <- modelF$chain4[,paste('Sp', sep = '')]
    plot(1:burnin, Tp1, type = 'l', xlab = 'b', ylab = "Specificity",
         xlim = c(0, 4*burnin))
    lines(burnin + 1:burnin, Tp2, col = 2)
    lines(2*burnin + 1:burnin, Tp3, col = 3)
    lines(3*burnin + 1:burnin, Tp4, col = 4)
  } else {
    Tp1 <- modelF$chain1[,paste('Sr', sep = '')]
    Tp2 <- modelF$chain2[,paste('Sr', sep = '')]
    Tp3 <- modelF$chain3[,paste('Sr', sep = '')] 
    Tp4 <- modelF$chain4[,paste('Sr', sep = '')]
    plot(1:burnin, Tp1, type = 'l', xlab = 'b', ylab = "Seroreversion",
         xlim = c(0, 4*burnin))
    lines(burnin + 1:burnin, Tp2, col = 2)
    lines(2*burnin + 1:burnin, Tp3, col = 3)
    lines(3*burnin + 1:burnin, Tp4, col = 4)
  }
}

# Gelman-Rubin diagnostic
GRF_Se <- mcmc.list(as.mcmc(modelF$chain1[,'Se']), as.mcmc(modelF$chain2[,'Se']), 
                     as.mcmc(modelF$chain3[,'Se']), as.mcmc(modelF$chain4[,'Se']))
print(gelman.diag(GRF_Se))

GRF_Sp <- mcmc.list(as.mcmc(modelF$chain1[,'Sp']), as.mcmc(modelF$chain2[,'Sp']), 
                     as.mcmc(modelF$chain3[,'Sp']), as.mcmc(modelF$chain4[,'Sp']))
print(gelman.diag(GRF_Sp))

GRF_Sr <- mcmc.list(as.mcmc(modelF$chain1[,'Sr']), as.mcmc(modelF$chain2[,'Sr']), 
                     as.mcmc(modelF$chain3[,'Sr']), as.mcmc(modelF$chain4[,'Sr']))
print(gelman.diag(GRF_Sr))

## Credible Interval
quantile(SeF, probs = c(0.5, 0.025, 0.975)) # 0.8738667 0.8272994 0.9128040 
quantile(SpF, probs = c(0.5, 0.025, 0.975)) # 0.9947740 0.9870746 0.9991958 
quantile(SrF, probs = c(0.5, 0.025, 0.975)) # 0.2926801 0.2127489 0.3824631 