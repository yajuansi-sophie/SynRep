# Application to the ACS data


library(rstan)
library(Matrix)
library(survival)
library(polyapost)
library(StanHeaders)
library(doParallel)
library(foreach)

library(tidyverse)
library(polyapost)
library(kableExtra)
library(LaplacesDemon)
library(sampling)
library(pps)
library(survey)
library(PracTools)
library(haven)
library(gridExtra)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

wlm_posterior_acs <- stan_model('wgtlikelihood_y.stan')

L=1
M = 10 #10
R=10
c <- 50
run_num <- 1

runs <- list(qboot = matrix(NA, run_num, 3), #naive estimate from sample not accounting for weights
             qbootvar = matrix(NA, run_num, 3), #variance estimate for naive estimate from sample not accounting for weights
             htboot = matrix(NA, run_num, 3), #design-based estimate from sample
             htbootvar = matrix(NA, run_num, 3), #design-based estimate from sample
             qsyn = matrix(NA, run_num, 3), #estimate from synthetic data generated from original sample
             qsynvar = matrix(NA, run_num, 3), #variance estimate for qsyn
             wtreg = matrix(NA, run_num, 3), #estimate from weighted likelihood
             wtregvar = matrix(NA, run_num, 3), #variance estimate for weighted lkl estimate
             Qbar = matrix(NA, run_num, 3), #pop level Q est
             qbar = matrix(NA, run_num, 3), #synth level Q est
             qbar_rep = matrix(NA, run_num, 3), #rep level Q est
             qbar_rep_1 = matrix(NA, run_num, 3), # R=1 rep  
             Bm = matrix(NA, run_num, 3), #pop level Bm est
             bm = matrix(NA, run_num, 3), #synth level bm est
             Bm_rep = matrix(NA, run_num, 3), #rep level Bm est
             Bm_rep_1 = matrix(NA, run_num, 3), #R=1 rep level Bm est
             vbar = matrix(NA, run_num, 3), #synth level vbar
             vbar_rep = matrix(NA, run_num, 3), #rep level vbar
             vbar_rep_1 = matrix(NA, run_num, 3), #R=1 rep vbar
             bm_rep = matrix(NA, run_num, 3), #rep level bm (var between qlbar's and qbar)
             bm_rep_1 = matrix(NA, run_num, 3), #R=1 rep level bm (var between qlbar's and qbar)
             wbar = matrix(NA, run_num, 3), #mean of wl's
             Qvar = matrix(NA, run_num, 3), #pop level Qvar est
             Tm = matrix(NA, run_num, 3),  #synth level Qvar est
             Tm_rep = matrix(NA, run_num, 3), #rep level Qvar est
             Tm_rep_1 = matrix(NA, run_num, 3), #R=1 rep level Qvar est
             Tm_0 = matrix(NA, run_num, 3),  #possible negative var synthe level Qvar est
             Tm_rep_0 = matrix(NA, run_num, 3), #rep level Qvar est
             Tm_rep_1_0 = matrix(NA, run_num, 3), #R=1 rep level Qvar est
             Tm_wtreg_rep = matrix(NA, run_num, 3), #rep weighted Qvar est
             result_pop = matrix(NA, run_num, 3),
             result_synth = matrix(NA, run_num, 3),
             result_rep = matrix(NA, run_num, 3),
             result_rep_1 = matrix(NA, run_num, 3),
             results_srs = matrix(NA, run_num, 3),
             result_ht = matrix(NA, run_num, 3),
             result_srssyn = matrix(NA, run_num, 3),
             result_wtreg = matrix(NA, run_num, 3),
             result_wtreg_rep =  matrix(NA, run_num, 3),
             wtreg_rep = matrix(NA, run_num, 3),
             wtregvar_rep = matrix(NA, run_num, 3)
             
) 


# load('/Users/yajuan/Dropbox (University of Michigan)/SyntheticSurveys/data/ACS2021_1yr/psam_pus_combined.rda')
# 
# acs_mi = psam_pus_combined %>% filter( ST == 26) %>% filter ( !is.na(PINCP)) %>% select(PWGTP, AGEP, PINCP)
# 
# n = dim(acs_mi)[1]
# 
# 
# acs_use = acs_mi %>% mutate(
#   X = ifelse(AGEP >= 65, 1, 0),
#   wts = PWGTP,
#   Y = (abs(PINCP))**(1/3) * sign(PINCP)
# ) %>% select(X, wts, Y,PINCP)
# 
# 
# saveRDS(acs_use, file ="acs_data.RDS")

orig_sample <- readRDS("acs_data.RDS")

N = sum(orig_sample$wts)
n = dim(orig_sample)[1]



#record sample non-design-based estimate
lm_ori_reg = summary(lm(Y ~ X, data = orig_sample))

j=1

runs$qboot[j,] = c(mean(orig_sample$X), mean(orig_sample$PINCP), lm_ori_reg$coef[2,1])
runs$qbootvar[j,] = c(mean(orig_sample$X)*(1-mean(orig_sample$X))/n, var(orig_sample$PINCP)/n, (lm_ori_reg$coef[2,2])^2)

#record design-based estimate
dsgn <- svydesign(ids=~1, data = orig_sample, weights =~wts)

htmean_y <- svymean(~PINCP, dsgn)
htmean_x <- svymean(~X, dsgn)
wt_svyglm = summary(svyglm(Y ~ X, dsgn))

runs$htboot[j,] = c(coef(htmean_x), coef(htmean_y), wt_svyglm$coef[2,1])
runs$htbootvar[j,] = c(SE(htmean_x)^2, SE(htmean_y)^2, (wt_svyglm$coef[2,2])^2)

#simulate from original sample and take mean
simple_synth_x <- rbinom(n, 1, mean(orig_sample$X))
simple_synth_y <- rnorm(n) * lm_ori_reg$sigma + lm_ori_reg$coef[1,1] + lm_ori_reg$coef[2,1] * simple_synth_x

simple_synth_y3 = (simple_synth_y)** 3 * sign(simple_synth_y)
runs$qsyn[j,] = c(mean(simple_synth_x), mean(simple_synth_y3), lm_ori_reg$coef[2,1])
runs$qsynvar[j,] = c(mean(simple_synth_x) *(1-mean(simple_synth_x))/n, var(simple_synth_y3)/n, (lm_ori_reg$coef[2,2])^2)

### pseudo-likelihood

stan_data_LW <- list(
  n = n,
  n_pred = n,
  y = as.vector(orig_sample$Y), 
  x = as.vector(orig_sample$X), 
  weights = as.vector(orig_sample$wts)
)

sim_est_LW <- sampling(
  object = wlm_posterior_acs, 
  data = stan_data_LW,
  chains = 4,
  iter = 4000, warmup = 2000, thin = 5,
  verbose = FALSE, refresh=0
)

posterior_samples <- rstan::extract(sim_est_LW, permuted = TRUE)

runs$wtreg[j,] = c(mean(posterior_samples$q_x)/n, mean(posterior_samples$q_y),mean(posterior_samples$beta2))
runs$wtregvar[j,] = c(var(posterior_samples$q_x)/n/n, var(posterior_samples$q_y), var(posterior_samples$beta2))


wt_max = posterior_samples$max_y3

### Bootstrap
#### 
synth_samples <- list(Q_y = rep(NA, M),
                      synth_q_y = rep(NA, M),
                      synth_v_y = rep(NA,M),
                      Q_x = rep(NA, M),
                      synth_q_x = rep(NA, M),
                      synth_v_x = rep(NA,M),
                      Q_beta = rep(NA, M),
                      synth_q_beta = rep(NA, M),
                      synth_v_beta = rep(NA,M),
                      wtreg_boot = matrix(NA, nrow = M, ncol=3), 
                      wtregvar_boot = matrix(NA, nrow = M, ncol=3), 
                      rep_y = matrix(data = NA, nrow = M, ncol = R),
                      rep_v_y = matrix(data = NA, nrow = M, ncol = R),
                      rep_x = matrix(data = NA, nrow = M, ncol = R),
                      rep_v_x = matrix(data = NA, nrow = M, ncol = R),
                      rep_beta = matrix(data = NA, nrow = M, ncol = R),
                      rep_v_beta = matrix(data = NA, nrow = M, ncol = R)
                      
)

wtboot_max = matrix(data = NA, nrow = M, ncol = 1600)

syn_max_y3 = matrix(data = NA, nrow = M, ncol = R)
####      
dsgn.rw <- as.svrepdesign(design=dsgn, type="subbootstrap",replicates=M)
repwt <- as.matrix(dsgn.rw$repweights)
########

for (i in seq(1:M)) {
  #create bootstrap sample for obtaining weights and then create population
  n_synth <- n
  
  wts <- repwt[,i]*orig_sample$wts
  wtind <- wts!=0
  wts <- wts[wtind]
  
  polya_ind <- (1:n)[wtind]
  
  wts <- (c*n)*(wts/sum(wts))
  bootsamp_size <- length(wts)
  
  ### pseudo-likelihood
  
  stan_data_boot <- list(
    n = bootsamp_size,
    n_pred = n,
    y = unlist(subset(orig_sample[polya_ind,], select = "Y")),
    x = unlist(subset(orig_sample[polya_ind,], select = "X")),
    weights = wts
  )
  
  sim_est_boot <- sampling(
    object = wlm_posterior_acs,
    data = stan_data_boot,
    chains = 4,
    iter = 4000, warmup = 2000, thin = 5,
    verbose = FALSE, refresh=0
  )
  
  posterior_samples_boot <- rstan::extract(sim_est_boot, permuted = TRUE)
  
  synth_samples$wtreg_boot[i,] = c(mean(posterior_samples_boot$q_x)/n, mean(posterior_samples_boot$q_y),mean(posterior_samples_boot$beta2))
  synth_samples$wtregvar_boot[i,] = c(var(posterior_samples_boot$q_x)/n/n, var(posterior_samples_boot$q_y), var(posterior_samples_boot$beta2))
  
  wtboot_max[i, ] = posterior_samples_boot$max_y3
  
  ###
  
  
  synpop <- integer()
  
  for (l in seq(1:L)) {
    pop_synth <- c(polyapost::wtpolyap(polya_ind, wts, c*n - bootsamp_size))
    synpop <- c(synpop, pop_synth)
  }
  
  sample_synth <- sample(synpop, n_synth) #srs from psedo-pop
  
  #calculate relevant quantities
  
  synth_samples$Q_y[[i]] <- apply(orig_sample[synpop,] %>% select(PINCP), 2, mean)
  
  synth_samples$synth_q_y[[i]] <- q_l <- apply(orig_sample[sample_synth,] %>% select(PINCP), 2, mean)
  
  synth_samples$synth_v_y[[i]] <- v_l <- (1-n_synth/(c*n*L))*apply(orig_sample[sample_synth,] %>% select(PINCP), 2, var)/n_synth
  
  
  synth_samples$Q_x[[i]] <- apply(orig_sample[synpop,] %>% select(X), 2, mean)
  
  synth_samples$synth_q_x[[i]] <- q_l_x<- apply(orig_sample[sample_synth,] %>% select(X), 2, mean)
  
  synth_samples$synth_v_x[[i]] <- v_l_x <- (1-n_synth/(c*n*L))*q_l_x * (1 - q_l_x)/n_synth
  
  
  synth_samples$Q_beta[[i]] <- summary(lm(Y ~ X, orig_sample[synpop,]))$coef[2,1]
  
  synth_samples$V_beta[[i]] <- (summary(lm(Y ~ X, orig_sample[synpop,]))$coef[2,2]) ^ 2
  
  lm_srs_reg = summary(lm(Y ~ X, orig_sample[sample_synth,]))
  
  synth_samples$synth_q_beta[[i]] <- lm_srs_reg$coef[2,1]
  
  synth_samples$synth_v_beta[[i]]  <- (lm_srs_reg$coef[2,2]) ^ 2
  

  
  #model synthetic dataset
  for (k in seq(1:R)) {

    ## Plug in without PPD
    
    rep_x = rbinom(n_synth, 1, q_l_x)
    
    rep_y = rnorm(n_synth) * lm_srs_reg$sigma + lm_srs_reg$coef[1,1] + lm_srs_reg$coef[2,1] * rep_x
    
    rep_y3 = (abs(rep_y)) ** 3 * sign(rep_y)
    
    syn_max_y3[i,k] = max(rep_y3)
    
    synth_samples$rep_y[i, k] <- mean(rep_y3)
    synth_samples$rep_v_y[i, k] <- (1-n_synth/(c*n*L))*var(rep_y3)/n_synth 
    
    synth_samples$rep_x[i, k] <- mean(rep_x)
    synth_samples$rep_v_x[i, k] <- (1-n_synth/(c*n*L))*mean(rep_x)*(1-mean(rep_x))/n_synth
    
    synth_samples$rep_beta[i, k]  = summary(lm(rep_y~rep_x))$coef[2,1]
    synth_samples$rep_v_beta[i, k]  = (summary(lm(rep_y~rep_x))$coef[2,2]) ^ 2
    
    
  }
}


#store information
Q <- cbind(synth_samples$Q_x,synth_samples$Q_y,synth_samples$Q_beta)
q <- cbind(synth_samples$synth_q_x, synth_samples$synth_q_y, synth_samples$synth_q_beta)
v <- cbind(synth_samples$synth_v_x, synth_samples$synth_v_y,synth_samples$synth_v_beta)


##
runs$wtreg_rep[j,] = apply(synth_samples$wtreg_boot, 2, mean)
runs$wtregvar_rep[j,] = apply(synth_samples$wtregvar_boot, 2, mean)


runs$Qbar[j,] <- Qbar <- apply(Q,2,mean)
runs$qbar[j,] <- qbar <- apply(q, 2, mean)
runs$vbar[j,] <- vbar <- apply(v, 2, mean)

runs$qbar_rep[j,] <- qbar_rep <- c(mean(synth_samples$rep_x), mean(synth_samples$rep_y), mean(synth_samples$rep_beta))
runs$vbar_rep[j,] <- vbar_rep <- c(mean(synth_samples$rep_v_x), mean(synth_samples$rep_v_y), mean(synth_samples$rep_v_beta))


ql_bar <- cbind(apply(synth_samples$rep_x,1,mean), apply(synth_samples$rep_y,1,mean), apply(synth_samples$rep_beta,1,mean))
wbar <- cbind(mean(apply(synth_samples$rep_x,1,var)), mean(apply(synth_samples$rep_y,1,var)), mean(apply(synth_samples$rep_beta,1,var)))


#R=1 
q_rep_1 <- cbind(synth_samples$rep_x[,R], synth_samples$rep_y[,R], synth_samples$rep_beta[,R])
v_rep_1 <- cbind(synth_samples$rep_v_x[,R], synth_samples$rep_v_y[,R], synth_samples$rep_v_beta[,R])

runs$qbar_rep_1[j,] <- qbar_rep_1 <- c(mean(synth_samples$rep_x[,R]), mean(synth_samples$rep_y[,R]), mean(synth_samples$rep_beta[,R]))
runs$vbar_rep_1[j,] <- vbar_rep_1 <- c(mean(synth_samples$rep_v_x[,R]), mean(synth_samples$rep_v_y[,R]), mean(synth_samples$rep_v_beta[,R]))

#compute parameters for posterior from pop level estimates
for (tt in 1:3){
  
  runs$Bm[j,tt] <- Bm <- (1/(M-1))*sum((Qbar[tt] - Q[,tt])^2)
  runs$Qvar[j, tt] <- (1 + (1/M))*Bm
  
  #compute parameters for approx posterior from sample level estimates
  
  runs$bm[j,tt] <- bm <- sum((q[, tt]-qbar[tt])^2)/(M-1)
  runs$Tm[j,tt] <- tm <- if_else((1 + (1/M))*bm - vbar[tt] <= 0, 
                                 (1 + (2/M))*vbar[tt], 
                                 (1 + (1/M))*bm - vbar[tt])
  runs$Tm_0[j,tt]  = (1 + (1/M))*bm - vbar[tt]
  
  
  runs$Tm_wtreg_rep[j,tt] = if_else((1 + (1/M))*var(runs$wtreg_rep[j,tt] - synth_samples$wtreg_boot[,tt])  - runs$wtregvar_rep[j,tt] <=0, 
                                    (1 + (2/M))*runs$wtregvar_rep[j,tt],
                                    (1 + (1/M))*var(runs$wtreg_rep[j,tt] - synth_samples$wtreg_boot[,tt])  - runs$wtregvar_rep[j,tt])
  #compute parameters for approx posterior from replicate sample estimates
  runs$bm_rep[j,tt] <- bm_rep <- sum((ql_bar[,tt] - qbar_rep[tt])^2)/(M-1)
  
  runs$wbar[j,tt] <- wbar[tt]
  
  runs$Bm_rep[j,tt] <- Bm_rep <- bm_rep - vbar_rep[tt] - wbar[tt]/R
  
  runs$Tm_rep[j,tt] <- tm_rep <- if_else(((1/M) + 1)*bm_rep - vbar_rep[tt] - wbar[tt]/R <= 0, 
                                         (1 + (2/M))*vbar_rep[tt] + (1/(R*M))*wbar[tt], 
                                         ((1/M) + 1)*bm_rep - vbar_rep[tt] - wbar[tt]/R)
  
  runs$Tm_rep_0[j,tt] <-  ((1/M) + 1)*bm_rep - vbar_rep[tt] - wbar[tt]/R
  
  #R=1
  runs$bm_rep_1[j,tt] <- bm_rep_1 <- var(q_rep_1[,tt])
  runs$Bm_rep_1[j,tt] <- Bm_rep_1 <- bm_rep_1 - 2*vbar_rep_1[tt]
  runs$Tm_rep_1_0[j,tt] <- ((1/M) + 1)*bm_rep_1 - 2*vbar_rep_1[tt]
  runs$Tm_rep_1[j,tt] <- tm_rep_1 <- if_else(((1/M) + 1)*bm_rep_1 - 2*vbar_rep_1[tt] <= 0, 
                                             (1 + (3/M))*vbar_rep_1[tt], 
                                             ((1/M) + 1)*bm_rep_1 - 2*vbar_rep_1[tt])
  
}

ci_pop <- cbind(qt(0.025, M-1) *sqrt(runs$Qvar[j,]) + runs$Qbar[j,], qt(0.975, M-1) *sqrt(runs$Qvar[j,]) + runs$Qbar[j,])
ci_synth <- cbind(qt(0.025, M-1) *sqrt(runs$Tm[j,]) + runs$qbar[j,], qt(0.975, M-1) *sqrt(runs$Tm[j,]) + runs$qbar[j,])
ci_rep <- cbind(qt(0.025, M-1) *sqrt(runs$Tm_rep[j,]) + runs$qbar_rep[j,], qt(0.975, M-1) *sqrt(runs$Tm_rep[j,]) + runs$qbar_rep[j,]) 
ci_rep_1 <- cbind(qt(0.025, M-1) *sqrt(runs$Tm_rep_1[j,]) + runs$qbar_rep_1[j,], qt(0.975, M-1) *sqrt(runs$Tm_rep_1[j,]) + runs$qbar_rep_1[j,]) 
ci_srs <- cbind(qt(0.025, M-1) *sqrt(runs$qbootvar[j,]) + runs$qboot[j,], qt(0.975, M-1) *sqrt(runs$qbootvar[j,]) + runs$qboot[j,]) 
ci_ht <- cbind(qt(0.025, M-1) *sqrt(runs$htbootvar[j,]) + runs$htboot[j,], qt(0.975, M-1) *sqrt(runs$htbootvar[j,]) + runs$htboot[j,])
# ci_ht <- confint(htmean)
ci_srssyn <- cbind(qt(0.025, M-1) *sqrt(runs$qsynvar[j,]) + runs$qsyn[j,], qt(0.975, M-1) *sqrt(runs$qsynvar[j,]) + runs$qsyn[j,]) 
ci_wtreg <- cbind(qt(0.025, M-1) *sqrt(runs$wtregvar[j,]) + runs$wtreg[j,], qt(0.975, M-1) *sqrt(runs$wtregvar[j,]) + runs$wtreg[j,])
ci_wtreg_rep <- cbind(qt(0.025, M-1) *sqrt(runs$Tm_wtreg_rep[j,]) + runs$wtreg_rep[j,], qt(0.975, M-1) *sqrt(runs$Tm_wtreg_rep[j,]) + runs$wtreg_rep[j,])



save.image("application-0216.RData")


#load("application.RData")


# ap_output =  data.frame(rbind(cbind(rep("Direct",3), c("x", "y", "beta"), runs$qboot[1,], ci_srs),
#       cbind(rep("HT",3), c("x", "y", "beta"), runs$htboot[1,], ci_ht),
#       cbind(rep("SRSsyn",3), c("x", "y", "beta"), runs$qsyn[1,], ci_srssyn),
#       cbind(rep("Pseudo-Pop",3), c("x", "y", "beta"), runs$Qbar[1,], ci_pop),
#       cbind(rep("Pseudo-SRS",3), c("x", "y", "beta"), runs$qbar[1,], ci_synth),
#       cbind(rep("SynRep-R",3), c("x", "y", "beta"), runs$qbar_rep[1,], ci_rep),
#       cbind(rep("SynRep-1",3), c("x", "y", "beta"), runs$qbar_rep_1[1,], ci_rep_1),
#       cbind(rep("Wtreg",3), c("x", "y", "beta"), runs$wtreg[1,], ci_wtreg),
#       cbind(rep("Wtreg-Boot",3), c("x", "y", "beta"), runs$wtreg_rep[1,], ci_wtreg_rep)))
# 
# names(ap_output) = c("method","term","est","low","up")
# 
# method_list = c("Direct", "SRSsyn", "HT", "Pseudo-Pop", "Pseudo-SRS", "SynRep-R", "SynRep-1", "Wtreg","Wtreg-Boot")
# 
# 
# 
# ap_output = ap_output%>% mutate(
#   est = as.numeric(as.character(est)),
#   low = as.numeric(as.character(low)),
#   up = as.numeric(as.character(up)),
#   method = ordered(method, levels=method_list),
#   term = ordered(term, levels=c("x","y","beta"))
# )
# 
# 
# 
# # p1=ggplot(ap_output %>% filter(term == "x"), aes(x = method, y=est, ymin=low, ymax=up)) +
# #   geom_pointrange() +
# #   coord_flip() + #facet_grid(rows = vars(term)) +
# #   facet_wrap(~term, scales = "free")
# 
#  p1= ggplot(ap_output, aes(x = method, y=est, ymin=low, ymax=up)) +
#   geom_pointrange(fatten = 0.5, size = 1) +
#   facet_grid(rows = vars(term),scales = "free",labeller = as_labeller(c("beta" = "scriptstyle(beta)", "x"="scriptstyle(bar(Y)[1])", "y"="scriptstyle(bar(Y)[2])"),label_parsed))+
#   labs(title="", y="Estimate", x="")+
#   theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size =5, angle=10), axis.title.y = element_text(size =7), legend.title = element_blank(),legend.position = "bottom")
# 
# 
#  ggsave(p1,
#  width = 5.5, height =7, units = "in", device="pdf", file=paste0(path,"fig/acs-ap-0216.pdf"))
#  
#  #disclosure risk
#  
#  table(orig_sample$X, rep_x)
#  
#  4407/n
#  49439/n
#  
#  max(orig_sample$PINCP) - syn_max_y3
#  
#  max(orig_sample$PINCP) - wt_max
#  
#  max(orig_sample$PINCP) - wtboot_max
#  
#  
#  mean(order(orig_sample$Y) == order(rep_y))
#  mean(orig_sample$Y == rep_y)
#  
#  
#  orig_sample$Y[sum(order(orig_sample$Y) == order(rep_y))]
# 
#  max(orig_sample$Y) 
#  max(rep_y)
