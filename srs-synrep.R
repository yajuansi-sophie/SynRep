## Simulation from PPS
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
#wlm_posterior <- stan_model('wgtlikelihood.stan')
#lm_posterior <- stan_model('likelihood.stan')

acs_w = readRDS("ACS2021_w.RDS")

syn = function(M, R){
 
#M = 50; R = 10;  n=1000
#M=10
#R=10
n=500

L = 1
c <- 50
run_num <- 1000 #5000

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
             qbar_ppd_rep = matrix(NA, run_num, 3),  #ppd        
             Bm = matrix(NA, run_num, 3), #pop level Bm est
             bm = matrix(NA, run_num, 3), #synth level bm est
             Bm_rep = matrix(NA, run_num, 3), #rep level Bm est
             Bm_ppd_rep = matrix(NA, run_num, 3), #rep ppd level Bm est
             Bm_rep_1 = matrix(NA, run_num, 3), #R=1 rep level Bm est
             vbar = matrix(NA, run_num, 3), #synth level vbar
             vbar_rep = matrix(NA, run_num, 3), #rep level vbar
             vbar_rep_1 = matrix(NA, run_num, 3), #R=1 rep vbar
             vbar_ppd_rep = matrix(NA, run_num, 3), #rep ppd level vbar
             bm_rep = matrix(NA, run_num, 3), #rep level bm (var between qlbar's and qbar)
             bm_ppd_rep = matrix(NA, run_num, 3), #rep ppd level bm (var between qlbar's and qbar)
             bm_rep_1 = matrix(NA, run_num, 3), #R=1 rep level bm (var between qlbar's and qbar)
             wbar = matrix(NA, run_num, 3), #mean of wl's
             vbar = matrix(NA, run_num, 3),
             wbar_ppd = matrix(NA, run_num, 3), #mean of ppd wl's
             Qvar = matrix(NA, run_num, 3), #pop level Qvar est
             Tm = matrix(NA, run_num, 3),  #synth level Qvar est
             Tm_rep = matrix(NA, run_num, 3), #rep level Qvar est
             Tm_ppd_rep = matrix(NA, run_num, 3), #rep level ppd Qvar est
             Tm_rep_1 = matrix(NA, run_num, 3), #R=1 rep level Qvar est
             Tm_0 = matrix(NA, run_num, 3),  #possible negative var synthe level Qvar est
             Tm_rep_0 = matrix(NA, run_num, 3), #rep level Qvar est
             Tm_rep_02 = matrix(NA, run_num, 3), #rep level Qvar est
             Tm_ppd_rep_0 = matrix(NA, run_num, 3),#rep level PPD Qvar est
             Tm_rep_1_0 = matrix(NA, run_num, 3), #R=1 rep level Qvar est
             Tm_wtreg_rep = matrix(NA, run_num, 3), #rep weighted Qvar est
             result_pop = matrix(NA, run_num, 3),
             result_synth = matrix(NA, run_num, 3),
             result_rep = matrix(NA, run_num, 3),
             result_ppd_rep = matrix(NA, run_num, 3),
             result_rep_1 = matrix(NA, run_num, 3),
             results_srs = matrix(NA, run_num, 3),
             result_ht = matrix(NA, run_num, 3),
             result_srssyn = matrix(NA, run_num, 3),
             result_wtreg = matrix(NA, run_num, 3),
             result_wtreg_rep =  matrix(NA, run_num, 3),
             wtreg_rep = matrix(NA, run_num, 3),
             wtregvar_rep = matrix(NA, run_num, 3)
             
) 

  
  N = length(acs_w)
  population <- data.frame(index = c(1:N),
                           W = acs_w,
                           X = rep(NA, N),
                           Y = rep(NA, N),
                           p = rep(NA, N))
 
  beta1 = 20; beta2=50; sigma=50 #10
  
  population$X <- rbinom(N, 1, exp(-7+2*log(population$W))/(1+exp(-7+2*log(population$W)))) #depend on W
  
  population$Y <- rnorm(N, beta1 + beta2*population$X, sigma)
  
  population$p <- sampling::inclusionprobabilities(population$W, n)
  population$wts <- 1/sampling::inclusionprobabilities(population$W, n)

  # cor(population$Y, population$wts) #-0.1554094
  # 
  # summary(lm(Y ~ X, data = population))
  # summary(lm(Y ~ X + W, data = population))
  #  
  # summary(lm(Y ~ X, data = orig_sample))
  # summary(lm(Y ~ X + W, data = orig_sample)) 
  
  true_val <- c(mean(population$X), mean(population$Y), beta2)
  
  
  for (j in seq(1:run_num)) {

    #take sample pps
    orig_sample <- population[pps::ppss(population$W, n),]
    orig_sample$wts <- N*orig_sample$wts/sum(orig_sample$wts)
    
    #record sample non-design-based estimate
    lm_ori_reg = summary(lm(Y ~ X, data = orig_sample))
    
    runs$qboot[j,] = c(mean(orig_sample$X), mean(orig_sample$Y), lm_ori_reg$coef[2,1])
    runs$qbootvar[j,] = c(mean(orig_sample$X)*(1-mean(orig_sample$X))/n, var(orig_sample$Y)/n, (lm_ori_reg$coef[2,2])^2)
    
    #record design-based estimate
    dsgn <- svydesign(ids=~1, data = orig_sample, weights =~wts)
    
    htmean_y <- svymean(~Y, dsgn)
    htmean_x <- svymean(~X, dsgn)
    wt_svyglm = summary(svyglm(Y ~ X, dsgn))
    
    runs$htboot[j,] = c(coef(htmean_x), coef(htmean_y), wt_svyglm$coef[2,1])
    runs$htbootvar[j,] = c(SE(htmean_x)^2, SE(htmean_y)^2, (wt_svyglm$coef[2,2])^2)
    
    #simulate from original sample and take mean
    simple_synth_x <- rbinom(n, 1, mean(orig_sample$X))
    simple_synth_y <- rnorm(n) * lm_ori_reg$sigma + lm_ori_reg$coef[1,1] + lm_ori_reg$coef[2,1] * simple_synth_x
    
    runs$qsyn[j,] = c(mean(simple_synth_x), mean(simple_synth_y), lm_ori_reg$coef[2,1])
    runs$qsynvar[j,] = c(mean(simple_synth_x) *(1-mean(simple_synth_x))/n, var(simple_synth_y)/n, (lm_ori_reg$coef[2,2])^2)
 
    # ### pseudo-likelihood
    # 
    # stan_data_LW <- list(
    #   n = n,
    #   n_pred = n,
    #   y = as.vector(orig_sample$Y),  #orig_pop$Y[sample_index_mat[,i]]
    #   x = as.vector(orig_sample$X), #orig_pop$X[sample_index_mat[,i]]
    #   weights = as.vector(orig_sample$wts)
    # )
    # 
    # sim_est_LW <- sampling(
    #   object = wlm_posterior, 
    #   data = stan_data_LW,
    #   chains = 4,
    #   iter = 4000, warmup = 2000, thin = 5,
    #   #control = list(adapt_delta = 0.8, max_treedepth = 10),
    #   verbose = FALSE, refresh=0
    # )
    # 
    # posterior_samples <- rstan::extract(sim_est_LW, permuted = TRUE)
    # 
    # # plot(posterior_samples$beta2)
    # # 
    # # length(apply(posterior_samples$syn_x, 2, mean)) #samples
    # #mean(posterior_samples$q_x)/n *(1-mean(posterior_samples$q_x)/n) /n
    # 
    # runs$wtreg[j,] = c(mean(posterior_samples$q_x)/n, mean(posterior_samples$q_y),mean(posterior_samples$beta2))
    # runs$wtregvar[j,] = c(var(posterior_samples$q_x)/n/n, var(posterior_samples$q_y), var(posterior_samples$beta2))
    # 
    
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
                          rep_ppd_y = matrix(data = NA, nrow = M, ncol = R),
                          rep_ppd_v_y = matrix(data = NA, nrow = M, ncol = R),
                          rep_ppd_x = matrix(data = NA, nrow = M, ncol = R),
                          rep_ppd_v_x = matrix(data = NA, nrow = M, ncol = R),
                          rep_ppd_beta = matrix(data = NA, nrow = M, ncol = R),
                          rep_ppd_v_beta = matrix(data = NA, nrow = M, ncol = R),   
                          rep_y = matrix(data = NA, nrow = M, ncol = R),
                          rep_v_y = matrix(data = NA, nrow = M, ncol = R),
                          rep_x = matrix(data = NA, nrow = M, ncol = R),
                          rep_v_x = matrix(data = NA, nrow = M, ncol = R),
                          rep_beta = matrix(data = NA, nrow = M, ncol = R),
                          rep_v_beta = matrix(data = NA, nrow = M, ncol = R)  
                          
    )
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

      # ### pseudo-likelihood
      # 
      # stan_data_boot <- list(
      #   n = bootsamp_size,
      #   n_pred = n,
      #   y = unlist(subset(orig_sample[polya_ind,], select = "Y")),
      #   x = unlist(subset(orig_sample[polya_ind,], select = "X")),
      #   weights = wts
      # )
      # 
      # sim_est_boot <- sampling(
      #   object = wlm_posterior,
      #   data = stan_data_boot,
      #   chains = 2,
      #   iter = 3000, warmup = 1000, thin = 5,
      #   #control = list(adapt_delta = 0.8, max_treedepth = 10),
      #   verbose = FALSE, refresh=0
      # )
      # 
      # posterior_samples_boot <- rstan::extract(sim_est_boot, permuted = TRUE)
      # 
      # # plot(posterior_samples$beta2)
      # 
      # synth_samples$wtreg_boot[i,] = c(mean(posterior_samples_boot$q_x)/n, mean(posterior_samples_boot$q_y),mean(posterior_samples_boot$beta2))
      # synth_samples$wtregvar_boot[i,] = c(var(posterior_samples_boot$q_x)/n/n, var(posterior_samples_boot$q_y), var(posterior_samples_boot$beta2))
      # 

      ###


      synpop <- integer()
      
      for (l in seq(1:L)) {
        pop_synth <- c(polyapost::wtpolyap(polya_ind, wts, c*n - bootsamp_size))
        synpop <- c(synpop, pop_synth)
      }
      
      sample_synth <- sample(synpop, n_synth) #srs from psedo-pop
      
      #calculate relevant quantities
      
      synth_samples$Q_y[[i]] <- apply(orig_sample[synpop,] %>% select(Y), 2, mean)
      
      synth_samples$synth_q_y[[i]] <- q_l <- apply(orig_sample[sample_synth,] %>% select(Y), 2, mean)
      
      synth_samples$synth_v_y[[i]] <- v_l <- (1-n_synth/(c*n*L))*apply(orig_sample[sample_synth,] %>% select(Y), 2, var)/n_synth
      
   
      synth_samples$Q_x[[i]] <- apply(orig_sample[synpop,] %>% select(X), 2, mean)
      
      synth_samples$synth_q_x[[i]] <- q_l_x<- apply(orig_sample[sample_synth,] %>% select(X), 2, mean)
      
      synth_samples$synth_v_x[[i]] <- v_l_x <- (1-n_synth/(c*n*L))*q_l_x * (1 - q_l_x)/n_synth
      
      
      synth_samples$Q_beta[[i]] <- summary(lm(Y ~ X, orig_sample[synpop,]))$coef[2,1]
      
      synth_samples$V_beta[[i]] <- (summary(lm(Y ~ X, orig_sample[synpop,]))$coef[2,2]) ^ 2
        
      lm_srs_reg = summary(lm(Y ~ X, orig_sample[sample_synth,]))
      
      synth_samples$synth_q_beta[[i]] <- lm_srs_reg$coef[2,1]
        
      synth_samples$synth_v_beta[[i]]  <- (lm_srs_reg$coef[2,2]) ^ 2
        
      
      # ##PPD
      # stan_data <- list(
      #   n = n_synth,
      #   y = unlist(subset(orig_sample[sample_synth,], select = "Y")),
      #   x = unlist(subset(orig_sample[sample_synth,], select = "X"))
      # )
      # 
      # sim_est <- sampling(
      #   object = lm_posterior,
      #   data = stan_data,
      #   chains = 4,
      #   iter = 4000, warmup = 2000, thin = 5,
      #   # control = list(adapt_delta = 0.8, max_treedepth = 10),
      #   verbose = FALSE, refresh=0
      # )
      
      # PPD_samples <- rstan::extract(sim_est, permuted = TRUE)

      
      # plot(PPD_samples$beta2)
      #
      # length(apply(PPD_samples$syn_x, 2, mean)) #samples
      
     
      
      #rep_seq = dim(PPD_samples$syn_x)[1] -200/R * (R:1)
      #model synthetic dataset
      for (k in seq(1:R)) {
    
        # synth_samples$rep_ppd_y[i, k] <- mean(PPD_samples$syn_y[rep_seq[k],])
        # synth_samples$rep_ppd_v_y[i, k] <- var(PPD_samples$syn_y[rep_seq[k],])/n_synth  #(mean(apply(PPD_samples$syn_y,1, var)))/n_synth
        # 
        # synth_samples$rep_ppd_x[i, k] <- mean(PPD_samples$syn_x[rep_seq[k],])
        # synth_samples$rep_ppd_v_x[i, k] <- mean(PPD_samples$syn_x[rep_seq[k],])*(1-mean(PPD_samples$syn_x[rep_seq[k],]))/n_synth
        # 
        # synth_samples$rep_ppd_beta[i, k]  = summary(lm(PPD_samples$syn_y[rep_seq[k],] ~ PPD_samples$syn_x[rep_seq[k],]))$coef[2,1]
        # synth_samples$rep_ppd_v_beta[i, k]  = (summary(lm(PPD_samples$syn_y[rep_seq[k],] ~ PPD_samples$syn_x[rep_seq[k],]))$coef[2,2])^2

        ## Plug in without PPD
        
        rep_x = rbinom(n_synth, 1, q_l_x)
        
        rep_y = rnorm(n_synth) * lm_srs_reg$sigma + lm_srs_reg$coef[1,1] + lm_srs_reg$coef[2,1] * rep_x
        
        synth_samples$rep_y[i, k] <- mean(rep_y)
        synth_samples$rep_v_y[i, k] <- (1-n_synth/(c*n*L))*var(rep_y)/n_synth #(mean(apply(PPD_samples$syn_y,1, var)))/n_synth 
        
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
    # runs$wtreg_rep[j,] = apply(synth_samples$wtreg_boot, 2, mean)
    # runs$wtregvar_rep[j,] = apply(synth_samples$wtregvar_boot, 2, mean)
    
    
    runs$Qbar[j,] <- Qbar <- apply(Q,2,mean)
    runs$qbar[j,] <- qbar <- apply(q, 2, mean)
    runs$vbar[j,] <- vbar <- apply(v, 2, mean)
    
    runs$qbar_rep[j,] <- qbar_rep <- c(mean(synth_samples$rep_x), mean(synth_samples$rep_y), mean(synth_samples$rep_beta))
    runs$vbar_rep[j,] <- vbar_rep <- c(mean(synth_samples$rep_v_x), mean(synth_samples$rep_v_y), mean(synth_samples$rep_v_beta))
    
    runs$qbar_ppd_rep[j,] <- qbar_ppd_rep <- c(mean(synth_samples$rep_ppd_x), mean(synth_samples$rep_ppd_y), mean(synth_samples$rep_ppd_beta))
   runs$vbar_ppd_rep[j,] <- vbar_ppd_rep <- c(mean(synth_samples$rep_ppd_v_x), mean(synth_samples$rep_ppd_v_y), mean(synth_samples$rep_ppd_v_beta))
    
    ql_bar <- cbind(apply(synth_samples$rep_x,1,mean), apply(synth_samples$rep_y,1,mean), apply(synth_samples$rep_beta,1,mean))
    #wbar <- cbind(apply(synth_samples$rep_x,1,var), apply(synth_samples$rep_y,1,var), apply(synth_samples$rep_beta,1,var))
    wbar <- cbind(mean(apply(synth_samples$rep_x,1,var)), mean(apply(synth_samples$rep_y,1,var)), mean(apply(synth_samples$rep_beta,1,var)))
    ql_ppd_bar <- cbind(apply(synth_samples$rep_ppd_x,1,mean), apply(synth_samples$rep_ppd_y,1,mean), apply(synth_samples$rep_ppd_beta,1,mean))
    wbar_ppd <- cbind(mean(apply(synth_samples$rep_ppd_x,1,var)), mean(apply(synth_samples$rep_ppd_y,1,var)), mean(apply(synth_samples$rep_ppd_beta,1,var)))
    

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
    runs$vbar[j,tt] <- vbar_rep[tt] 
    runs$Bm_rep[j,tt] <- Bm_rep <- bm_rep - vbar_rep[tt] - wbar[tt]/R
  
    runs$Tm_rep[j,tt] <- tm_rep <- if_else(((1/M) + 1)*bm_rep - vbar_rep[tt] - wbar[tt]/R <= 0, 
                                          (1 + (2/M))*vbar_rep[tt] + (1/(R*M))*wbar[tt], 
                                          ((1/M) + 1)*bm_rep - vbar_rep[tt] - wbar[tt]/R)
    
    runs$Tm_rep_0[j,tt] <-  ((1/M) + 1)*bm_rep - vbar_rep[tt] - wbar[tt]/R
   
    runs$Tm_rep_02[j,tt] <-  ((1/M) + 1)*bm_rep - (1 + (1/R)) * vbar_rep[tt]
    
    #R=1
    runs$bm_rep_1[j,tt] <- bm_rep_1 <- var(q_rep_1[,tt])
    runs$Bm_rep_1[j,tt] <- Bm_rep_1 <- bm_rep_1 - 2*vbar_rep_1[tt]
    runs$Tm_rep_1_0[j,tt] <- ((1/M) + 1)*bm_rep_1 - 2*vbar_rep_1[tt]
    runs$Tm_rep_1[j,tt] <- tm_rep_1 <- if_else(((1/M) + 1)*bm_rep_1 - 2*vbar_rep_1[tt] <= 0, 
                                              (1 + (3/M))*vbar_rep_1[tt], 
                                              ((1/M) + 1)*bm_rep_1 - 2*vbar_rep_1[tt])
   ##ppd
    
    runs$bm_ppd_rep[j,tt] <- bm_ppd_rep <- sum((ql_ppd_bar[,tt] - qbar_ppd_rep[tt])^2)/(M-1)
    
    runs$wbar_ppd[j,tt] <- wbar_ppd[tt]
    
    runs$Bm_ppd_rep[j,tt] <- Bm_ppd_rep <- bm_ppd_rep - vbar_ppd_rep[tt] - wbar_ppd[tt]/R
    
    runs$Tm_ppd_rep[j,tt] <- tm_ppd_rep <- if_else(((1/M) + 1)*bm_ppd_rep - vbar_ppd_rep[tt] - wbar_ppd[tt]/R <= 0, 
                                           (1 + (2/M))*vbar_ppd_rep[tt] + (1/(R*M))*wbar_ppd[tt], 
                                           ((1/M) + 1)*bm_ppd_rep - vbar_ppd_rep[tt] - wbar_ppd[tt]/R)
    
    runs$Tm_ppd_rep_0[j,tt] <-  ((1/M) + 1)*bm_ppd_rep - vbar_ppd_rep[tt] - wbar_ppd[tt]/R

    
  }
    # ci_pop <- c(qnorm(0.025, mean = Qbar, sd = sqrt(Qvar)), qnorm(0.975, mean = Qbar, sd = sqrt(Qvar)))
    # ci_synth <- c(qnorm(0.025, mean = qbar, sd = sqrt(tm)), qnorm(0.975, mean = qbar, sd = sqrt(tm)))
    # ci_rep <- c(qnorm(0.025, mean = qbar_rep, sd = sqrt(tm_rep)), qnorm(0.975, mean = qbar_rep, sd = sqrt(tm_rep)))
    # ci_rep_1 <- c(qnorm(0.025, mean = qbar_rep_1, sd = sqrt(tm_rep_1)), qnorm(0.975, mean = qbar_rep_1, sd = sqrt(tm_rep_1)))
    # ci_srs <- c(qnorm(0.025, mean = qboot, sd = sqrt(qbootvar)), qnorm(0.975, mean = qboot, sd = sqrt(qbootvar)))
    # ci_ht <- confint(htmean)
    # ci_srssyn <- c(qnorm(0.025, mean = qsyn, sd = sqrt(qsynvar)), qnorm(0.975, mean = qsyn, sd = sqrt(qsynvar)))
    
    ci_pop <- cbind(qt(0.025, M-1) *sqrt(runs$Qvar[j,]) + runs$Qbar[j,], qt(0.975, M-1) *sqrt(runs$Qvar[j,]) + runs$Qbar[j,])
    ci_synth <- cbind(qt(0.025, M-1) *sqrt(runs$Tm[j,]) + runs$qbar[j,], qt(0.975, M-1) *sqrt(runs$Tm[j,]) + runs$qbar[j,])
    ci_rep <- cbind(qt(0.025, M-1) *sqrt(runs$Tm_rep[j,]) + runs$qbar_rep[j,], qt(0.975, M-1) *sqrt(runs$Tm_rep[j,]) + runs$qbar_rep[j,]) 
    ci_ppd_rep <- cbind(qt(0.025, M-1) *sqrt(runs$Tm_ppd_rep[j,]) + runs$qbar_ppd_rep[j,], qt(0.975, M-1) *sqrt(runs$Tm_ppd_rep[j,]) + runs$qbar_ppd_rep[j,]) 
    ci_rep_1 <- cbind(qt(0.025, M-1) *sqrt(runs$Tm_rep_1[j,]) + runs$qbar_rep_1[j,], qt(0.975, M-1) *sqrt(runs$Tm_rep_1[j,]) + runs$qbar_rep_1[j,]) 
    ci_srs <- cbind(qt(0.025, M-1) *sqrt(runs$qbootvar[j,]) + runs$qboot[j,], qt(0.975, M-1) *sqrt(runs$qbootvar[j,]) + runs$qboot[j,]) 
    ci_ht <- cbind(qt(0.025, M-1) *sqrt(runs$htbootvar[j,]) + runs$htboot[j,], qt(0.975, M-1) *sqrt(runs$htbootvar[j,]) + runs$htboot[j,])
    ci_srssyn <- cbind(qt(0.025, M-1) *sqrt(runs$qsynvar[j,]) + runs$qsyn[j,], qt(0.975, M-1) *sqrt(runs$qsynvar[j,]) + runs$qsyn[j,]) 
    ci_wtreg <- cbind(qt(0.025, M-1) *sqrt(runs$wtregvar[j,]) + runs$wtreg[j,], qt(0.975, M-1) *sqrt(runs$wtregvar[j,]) + runs$wtreg[j,])
    ci_wtreg_rep <- cbind(qt(0.025, M-1) *sqrt(runs$Tm_wtreg_rep[j,]) + runs$wtreg_rep[j,], qt(0.975, M-1) *sqrt(runs$Tm_wtreg_rep[j,]) + runs$wtreg_rep[j,])
    

    
    for (tt in 1:3){
    runs$result_pop[j, tt] <- if_else(ci_pop[tt,1] <=  true_val[tt] & true_val[tt] <= ci_pop[tt,2], "T", "F")
    runs$result_synth[j, tt] <- if_else(ci_synth[tt,1] <= true_val[tt] & true_val[tt]  <= ci_synth[tt,2], "T", "F")
    runs$result_rep[j, tt] <- if_else(ci_rep[tt,1] <= true_val[tt] & true_val[tt]  <= ci_rep[tt,2], "T", "F")
    runs$result_ppd_rep[j, tt] <- if_else(ci_ppd_rep[tt,1] <= true_val[tt] & true_val[tt]  <= ci_ppd_rep[tt,2], "T", "F")
    runs$result_rep_1[j, tt] <- if_else(ci_rep_1[tt,1] <= true_val[tt] & true_val[tt]  <= ci_rep_1[tt,2], "T", "F")
    runs$result_srs[j, tt] <- if_else(ci_srs[tt,1] <= true_val[tt] & true_val[tt]  <= ci_srs[tt,2], "T", "F")
    runs$result_ht[j, tt] <- if_else(ci_ht[tt,1] <= true_val[tt] & true_val[tt]  <= ci_ht[tt,2], "T", "F")
    runs$result_srssyn[j, tt] <- if_else(ci_srssyn[tt,1] <= true_val[tt] & true_val[tt]  <= ci_srssyn[tt,2], "T", "F")
    runs$result_wtreg[j, tt] <- if_else(ci_wtreg[tt,1] <= true_val[tt] & true_val[tt]  <= ci_wtreg[tt,2], "T", "F")
    runs$result_wtreg_rep[j, tt] <- if_else(ci_wtreg_rep[tt,1] <= true_val[tt] & true_val[tt]  <= ci_wtreg_rep[tt,2], "T", "F")
    }
    
    print(j)
  }

  output = NULL


  for (tt in 1:3){

    results <- data.frame(Method = c("Direct","HT", "SRSsyn", "Pseudo-Pop", "Pseudo-SRS", "SynRep-R", "SynRep-R-PPD", "SynRep-1","Wtreg","Wtreg-BS"),
                          QEst_bias = c(mean(runs$qboot[,tt]), mean(runs$htboot[,tt]), mean(runs$qsyn[,tt]), mean(runs$Qbar[,tt]), mean(runs$qbar[,tt]), mean(runs$qbar_rep[,tt]), mean(runs$qbar_ppd_rep[,tt]), mean(runs$qbar_rep_1[,tt]), mean(runs$wtreg[,tt]),mean(runs$wtreg_rep[,tt])) - true_val[tt],
                          QEst_bias_relpct = (c(mean(runs$qboot[,tt]), mean(runs$htboot[,tt]), mean(runs$qsyn[,tt]), mean(runs$Qbar[,tt]), mean(runs$qbar[,tt]), mean(runs$qbar_rep[,tt]), mean(runs$qbar_ppd_rep[,tt]),mean(runs$qbar_rep_1[,tt]), mean(runs$wtreg[,tt]),mean(runs$wtreg_rep[,tt])) - true_val[tt])/true_val[tt] * 100,
                          VarEst = c(mean(runs$qbootvar[,tt]), mean(runs$htbootvar[,tt]), mean(runs$qsynvar[,tt]), mean(runs$Qvar[,tt]), mean(runs$Tm[,tt]), mean(runs$Tm_rep[,tt]), mean(runs$Tm_ppd_rep[,tt]),mean(runs$Tm_rep_1[,tt]), mean(runs$wtregvar[,tt]), mean(runs$Tm_wtreg_rep[,tt])),
                          VarComp = c(var(runs$qboot[,tt]), var(runs$htboot[,tt]), var(runs$qsyn[,tt]), var(runs$Qbar[,tt]), var(runs$qbar[,tt]), var(runs$qbar_rep[,tt]), var(runs$qbar_ppd_rep[,tt]),var(runs$qbar_rep_1[,tt]), var(runs$wtreg[,tt]),var(runs$wtreg_rep[,tt])),
                          Var_sd = c(sd(runs$qbootvar[,tt]), sd(runs$htbootvar[,tt]), sd(runs$qsynvar[,tt]), sd(runs$Qvar[,tt]), sd(runs$Tm[,tt]), sd(runs$Tm_rep[,tt]), sd(runs$Tm_ppd_rep[,tt]),sd(runs$Tm_rep_1[,tt]), sd(runs$wtregvar[,tt]), sd(runs$Tm_wtreg_rep[,tt]))
    )


    results = results %>% mutate(Var_r = VarEst/VarComp,
                                 Var_r_dir = VarComp /var(runs$qboot[,tt]),
                                 Var_r_ht = VarComp /var(runs$htboot[,tt]),
                                 Var_r_srs = VarComp /var(runs$qsyn[,tt]),
                                 VarEst_raw = c(mean(runs$qbootvar[,tt]), mean(runs$htbootvar[,tt]), mean(runs$qsynvar[,tt]), mean(runs$Qvar[,tt]), mean(runs$Tm_0[,tt]), mean(runs$Tm_rep_0[,tt]), mean(runs$Tm_ppd_rep_0[,tt]),mean(runs$Tm_rep_1_0[,tt]), mean(runs$wtregvar[,tt]), mean(runs$Tm_wtreg_rep[,tt])),
                                 Var_r_raw = VarEst_raw/VarComp)


    ci_results <- data.frame(Method = c("Direct","HT", "SRSsyn", "Pseudo-Pop", "Pseudo-SRS", "SynRep-R", "SynRep-R-PPD", "SynRep-1","Wtreg","Wtreg-BS"),
                             Contained = c(mean(runs$result_srs[, tt] == "T"),
                                           mean(runs$result_ht[, tt] == "T"),
                                           mean(runs$result_srssyn[, tt] == "T"),
                                           mean(runs$result_pop[, tt] == "T"),
                                           mean(runs$result_synth[, tt] == "T"),
                                           mean(runs$result_rep[, tt] == "T"),
                                           mean(runs$result_ppd_rep[, tt] == "T"),
                                           mean(runs$result_rep_1[, tt] == "T"),
                                           mean(runs$result_wtreg[, tt] == "T"),
                                           mean(runs$result_wtreg_rep[, tt] == "T")),
                             Contained_raw = c(mean(runs$result_srs[runs$qbootvar[, tt] > 0, tt]== "T"),
                                               mean(runs$result_ht[runs$htbootvar[, tt] > 0, tt]== "T"),
                                               mean(runs$result_srssyn[runs$qsynvar[, tt] > 0, tt]== "T"),
                                               mean(runs$result_pop[runs$Qvar[, tt] > 0, tt]== "T"),
                                               mean(runs$result_synth[runs$Tm_0[, tt] > 0, tt]== "T"),
                                               mean(runs$result_rep[runs$Tm_rep_0[, tt] > 0, tt]== "T"),
                                               mean(runs$result_ppd_rep[runs$Tm_ppd_rep_0[, tt] > 0, tt]== "T"),
                                               mean(runs$result_rep_1[runs$Tm_rep_1_0[, tt] > 0, tt]== "T"),
                                               mean(runs$result_wtreg[runs$wtregvar[, tt] > 0, tt] == "T"),
                                               mean(runs$result_wtreg_rep[runs$Tm_wtreg_rep[, tt] >0, tt] == "T")),
                             prop_neg =c(mean(runs$qbootvar[, tt]<0), mean(runs$htbootvar[, tt]<0), mean(runs$qsynvar[, tt]<0), mean(runs$Qvar[, tt]<0), mean(runs$Tm_0[, tt]<0), mean(runs$Tm_rep_0[, tt]<0), mean(runs$Tm_ppd_rep_0[, tt]<0),mean(runs$Tm_rep_1_0[, tt]<0), mean(runs$wtregvar[,tt]<0),mean(runs$Tm_wtreg_rep[,tt]<0))
    )

    output[[tt]] = list(tt, results,ci_results)


  }
  return(output)

}

