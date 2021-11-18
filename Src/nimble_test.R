library(nimble)

nbCode <- nimbleCode({
   ###top-level regression
   B1 ~ dnorm(0, 100)
   B0 ~ dnorm(0, 100)
   sigB1 ~ dunif(1e-10, 100)
   sigB0 ~ dunif(1e-10, 100)
   for (k in 1:K) {
      ###low-level regression
      b1[k] ~ dnorm(mean = B1, sd = sigB1)
      b0[k] ~ dnorm(mean = B0, sd = sigB0)
      for (n in 1:N){
        p[n, k] ~ dunif(1e-10, 1)
        r[n, k] <- exp(b0[k] + X[n, k] * b1[k])
        Y[n, k] ~ dnegbin(size = r[n, k], prob = p[n, k])
      }
   }
})

span_index <- which(simdata$new_times %in% times)
n <- dim(rece[span_index, ])[1]
K = 50

nbData <- list(Y = Y[span_index, 1:50],
                X = X[span_index, 1:50])

nbConsts <- list(N = n,
                K = K)

nbInits <- list(B0 = 0,
                B1 = 0,
                b0 = rep(0, K),
                b1 = rep(0, K),
                sigB0 = 0.0001,
                sigB1 = 0.0001)

nbModel <- nimbleModel(code = nbCode,
                        data = nbData,
                        inits = nbInits,
                        constants = nbConsts)

#compile nimble model to C++ code—much faster runtime
C_nbModel <- compileNimble(nbModel, showCompilerOutput = FALSE)

#configure the MCMC
nbModel_conf <- configureMCMC(nbModel)

nbModel_conf$monitors <- c("B0", "B1", "sigB0", "sigB1")
nbModel_conf$addMonitors2(c("b0", "b1"))

#samplers
nbModel_conf$removeSamplers(c("B0", "B1", "sigB0", "sigB1", "b0", "b1"))
nbModel_conf$addSampler(target = c("B0", "B1", "sigB0", "sigB1"),
                    type = "AF_slice")
for(k in 1:K){
   nbModel_conf$addSampler(target = c(paste("b0[", k, "]", sep = ""),
                                    paste("b1[", k, "]", sep = "")),
                    type = "AF_slice")
}

#thinning to conserve memory when the samples are saved below
nbModel_conf$setThin(1)
nbModel_conf$setThin2(1)

#build MCMC
nbModelMCMC <- buildMCMC(nbModel_conf)

#compile MCMC to C++—much faster
C_nbModelMCMC <- compileNimble(nbModelMCMC, project = nbModel)

#number of MCMC iterations
niter=1000000

#set seed for replicability
set.seed(1)

#save the MCMC chain (monitored variables) as a matrix
samples <- runMCMC(C_nbModelMCMC, niter = niter)

#save samples
save(samples, file = "mcmc_samples_Nimble.RData")
