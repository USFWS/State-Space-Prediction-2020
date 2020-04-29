#Predict EMGO observation from two state space models
#devtools::install_github("USFWS/AKaerial")
library(AKaerial)
library(jagsUI)
library(dplyr)
library(ggplot2)

data <-  tibble::as_tibble(YKGHistoric$combined) %>% 
  mutate(Species = as.character(Species)) %>% 
  filter(Species == "EMGO") %>% 
  select(years=Year, N = itotal, SE = itotal.se)

#estimate average CV of index and mse around cv expectation
fit <- summary(lm(data$SE~data$N-1))
fit
#define jags model from AKaerial::print(StateSpace)
cat("
    model{
    # Priors
    logN.est[1] ~ dunif(log(pN1min), log(pN1max))  # Prior for initial population size log scale
    mean.r ~ dunif(prmin, prmax)                   # Prior for mean growth rate
    sigma.proc ~ dunif(psigmamin, psigmamax)       # Prior for sd of state process log scale
    tau.proc <- pow(sigma.proc, -2)  
    # Likelihood
    # State process
    for (t in 1:(T-1) ) {
    r[t] ~ dnorm(mean.r, tau.proc)
    logN.est[t+1] <- logN.est[t] + r[t]
    }
    # Observation process
    for (i in 1:(T-1)) {
    tau.obs[i] <- pow(sigma.obs[i], -2)
    y[i] ~ dnorm(exp(logN.est[i]), tau.obs[i])
    }
    # predict new observation
    beta.se ~ dnorm(BETA, 1/SE^2)           #from linear model
    s ~ dchisq(DF) 
    sigma.new ~ dnorm(beta.se*exp(logN.est[T]), s/((SIGMA^2)*DF) ) #from linear model
    tau.new <- pow(sigma.new, -2)
    y.new ~ dnorm(exp(logN.est[T]), tau.new)
    }", file = "ssm1.jags", fill = TRUE)

jags.data <- list(T = length(data$years)+1, y = data$N, sigma.obs = data$SE, 
                  prmin = -0.3, prmax = 0.3, psigmamin = 0, psigmamax = 0.3, 
                  pN1min = 18000, pN1max = 19000,
                  DF = fit$df[2], SIGMA = fit$sigma, 
                  BETA=fit$coefficients[1], SE=fit$coefficients[2])

parameters <- c("logN.est", "mean.r", "sigma.proc", "y.new", "beta.se", "sigma.new")
inits <- function() {
  list(logN.est = c(runif(1, log(18500), log(18900)), rep(NA, 35)), 
       mean.r = runif(1, -0.01, -0.01), 
       sigma.proc = runif(1, 0.01, 0.011), 
       r = runif(35, -0.01, 0.01),
       beta.se = runif(1, 0.04, 0.05), 
       s = runif(0.9, 1.1))}

out <- jags(jags.data, inits, parameters, model.file = "ssm1.jags", 
            n.chains = 4, n.thin = 1, n.burnin = 1000, n.iter = 5000)

summary(out)
print(out)

#posterior probability index is < threshold
sum(out$sims.list$y.new < 23000)/length(out$sims.list$y.new)

#Plot index and prediction
plotData <- data.frame(Year=data$years, Index=data$N, 
                      lower=data$N - 1.96*data$SE, 
                      upper = data$N + 1.96*data$SE, 
                      Observed = "Yes")
newData <- data.frame(Year=2020, Index=out$mean$y.new, 
                      lower=out$q2.5$y.new, 
                      upper = out$q97.5$y.new, 
                      Observed ="No")
plotData <- rbind(plotData, newData)
plotData <- data.frame(plotData, N=apply(exp(out$sims.list$logN.est), 2, mean), 
                       Model="State Space Mean")

gplot <- ggplot() +
  geom_pointrange(data=plotData, aes(x=Year, y=Index, ymin=lower, ymax=upper, 
                                     color=Observed, group = Observed), size=1) +
  geom_point(data=plotData, aes(x=Year, y=N, shape=Model)) +
  geom_hline(yintercept=23000) + 
  annotate(geom="text", x=1989, y=24000, label="Closure Threshold") + 
  labs(title = "Emperor goose index estimate and state space model prediction for 2020", 
       x="Year", y = "Index") + 
  theme(legend.position=c(0.2, 0.75)) + 
  guides(color=guide_legend("Index estimate")) 
print(gplot)
ggsave("ssm1.png")

##################################################
# define ssm2: different trend line after opening of hunting
# H define a 'knot' position where the trend line changes
#   not a parameter, must be specified
# Note for EMGO beta.r is determined almost entirely by the prior, only 3 observations of r after hunting
# cat("
#     model{
#     # Priors
#     logN.est[1] ~ dunif(log(pN1min), log(pN1max))  # Prior for initial population size log scale
#     mean.r ~ dunif(prmin, prmax)                   # Prior for mean growth rate, pre-hunting
#     beta.r ~dnorm(0, 10000) T(,0)                  #effect of hunting on trend
#     sigma.proc ~ dunif(psigmamin, psigmamax)       # Prior for sd of state process log scale
#     tau.proc <- pow(sigma.proc, -2)  
#     # Likelihood
#     # State process
#     for (t in 1:(H-1) ) {
#     r[t] ~ dnorm(mean.r, tau.proc)
#     logN.est[t+1] <- logN.est[t] + r[t]
#     }
#     for (t in H:(T-1) ) {
#     r[t] ~ dnorm(mean.r+beta.r, tau.proc)
#     logN.est[t+1] <- logN.est[t] + r[t]
#     }
#     # Observation process
#     for (i in 1:(T-1)) {
#     tau.obs[i] <- pow(sigma.obs[i], -2)
#     y[i] ~ dnorm(exp(logN.est[i]), tau.obs[i])
#     }
#     # predict new observation
#     beta.se ~ dnorm(BETA, 1/SE^2)           #from linear model
#     s ~ dchisq(DF) 
#     sigma.new ~ dnorm(beta.se*exp(logN.est[T]), s/((SIGMA^2)*DF) ) #from linear model
#     tau.new <- pow(sigma.new, -2)
#     y.new ~ dnorm(exp(logN.est[T]), tau.new)
#     }", file = "ssm2.jags", fill = TRUE)
# 
# jags.data$H <- 32
# 
# parameters <- c("logN.est", "mean.r", "sigma.proc", "y.new", "beta.se", "sigma.new", "beta.r")
# out <- jags(jags.data, inits, parameters, model.file = "ssm2.jags", 
#             n.chains = 4, n.thin = 1, n.burnin = 1000, n.iter = 5000)
# 
# summary(out)
# print(out)
# #posterior probability index is < threshold
# sum(out$sims.list$y.new < 23000)/length(out$sims.list$y.new)
# 
# plotData <- data.frame(Year=data$years, Index=data$N, 
#                        lower=data$N - 1.96*data$SE, 
#                        upper = data$N + 1.96*data$SE, 
#                        Observed = "Yes")
# newData <- data.frame(Year=2020, Index=out$mean$y.new, 
#                       lower=out$q2.5$y.new, 
#                       upper = out$q97.5$y.new, 
#                       Observed ="No")
# plotData <- rbind(plotData, newData)
# 
# plotData <- data.frame(plotData, N=apply(exp(out$sims.list$logN.est), 2, mean))
# 
# gplot <- ggplot() +
#   geom_pointrange(data=plotData, aes(x=Year, y=Index, ymin=lower, ymax=upper, 
#                                      color=Observed, group = Observed), size=1) +
#   geom_point(data=plotData, aes(x=Year, y=N)) +
#   geom_hline(yintercept=23000) + 
#   annotate(geom="text", x=1989, y=24000, label="Closure Threshold") + 
#   labs(title = "Emperor goose index estimate and state space model prediction for 2020", 
#        x="Year", y = "Index") + 
#   theme(legend.position=c(0.2, 0.75)) + 
#   guides(color=guide_legend("Index estimate"))
# print(gplot)

# n=100
# plotData2 <- data.frame(x1=rep(2019, n), x2=rep(2020, n), 
#                         y1=rep(data$N[35], n), y2=out$sims.list$y.new[1:n])
# gplot2 <- gplot
# for(i in 1:n){
#   df = data.frame(X=c(2019,2020), Y=unlist(plotData2[i,c("y1","y2")]))
#   gplot2 <- gplot2 + geom_line(data=df, aes(x=X, y=Y), color="#00CCCC")
# }
# print(gplot2)
