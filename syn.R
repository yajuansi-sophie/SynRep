## Simulation from PPS

acs_w = readRDS("ACS2021_w.RDS")

syn = function(M, R, L, n){
  
  set.seed(3)
  
  N = length(acs_w)
  population <- data.frame(index = c(1:N),
                           X = acs_w,
                           Y = rep(NA, N),
                           p = rep(NA, N),
                           sample = rep(NA, N))
  
  population$Y <- rnorm(N, 20 + 0.2*population$X, 100)
  
  population$p <- sampling::inclusionprobabilities(population$X, n)
  population$wts <- 1/sampling::inclusionprobabilities(population$X, n)

  
  true_mean <- mean(population$Y)

  c <- 50
  run_num <- 5000
  
  runs <- list(qboot = rep(NA, run_num), #naive estimate from sample not accounting for weights
               qbootvar = rep(NA, run_num), #variance estimate for naive estimate from sample not accounting for weights
               htboot = rep(NA, run_num), #design-based estimate from sample
               htbootvar = rep(NA, run_num), #design-based estimate from sample
               qsyn = rep(NA, run_num), #estimate from synthetic data generated from original sample
               qsynvar = rep(NA, run_num), #variance estimate for qsyn
               Qbar = rep(NA, run_num), #pop level Q est
               qbar = rep(NA, run_num), #synth level Q est
               qbar_rep = rep(NA, run_num), #rep level Q est
               qbar_rep_1 = rep(NA, run_num),
               Bm = rep(NA, run_num), #pop level Bm est
               bm = rep(NA, run_num), #synth level bm est
               Bm_rep = rep(NA, run_num), #rep level Bm est
               Bm_rep_1 = rep(NA, run_num),
               vbar = rep(NA, run_num), #synth level vbar
               vbar_rep = rep(NA, run_num), #rep level vbar
               vbar_rep_1 = rep(NA, run_num),
               bm_rep = rep(NA, run_num), #rep level bm (var between qlbar's and qbar)
               bm_rep_1 = rep(NA, run_num),
               wbar = rep(NA, run_num), #mean of wl's
               Qvar = rep(NA, run_num), #pop level Qvar est
               Tm = rep(NA, run_num),  #synth level Qvar est
               Tm_rep = rep(NA, run_num), #rep level Qvar est
               Tm_rep_1 = rep(NA, run_num), #rep level Qvar est
               Tm_0 = rep(NA, run_num),  #possible negative var synthe level Qvar est
               Tm_rep_0 = rep(NA, run_num), #rep level Qvar est
               Tm_rep_1_0 = rep(NA, run_num), #rep level Qvar est
               result_pop = rep(NA, run_num),
               result_synth = rep(NA, run_num),
               result_rep = rep(NA, run_num),
               result_rep_1 = rep(NA, run_num),
               results_srs = rep(NA, run_num),
               result_ht = rep(NA, run_num),
               result_srssyn = rep(NA, run_num)
  ) 
  
  for (j in seq(1:run_num)) {
    
    #take sample pps
    orig_sample <- population[pps::ppss(population$X, n),]
    
    #record sample non-design-based estimate
    runs$qboot[[j]] <- qboot <- mean(orig_sample$Y)
    runs$qbootvar[[j]] <- qbootvar <- var(orig_sample$Y)/n
    
    #record design-based estimate
    dpps<- svydesign(id=~1, probs=orig_sample$p, data=orig_sample)
    
    htmean <- svymean(~Y, dpps)
    runs$htboot[[j]] <- coef(htmean)
    runs$htbootvar[[j]] <- SE(htmean)^2
    
    #simulate from original sample and take mean
    simple_synth <- rnorm(n, mean(orig_sample$Y), sd(orig_sample$Y))
    runs$qsyn[[j]] <- qsyn <- mean(simple_synth)
    runs$qsynvar[[j]] <- qsynvar <- var(simple_synth)/n
    
    #store as matrix - might be easier for sums
    synth_samples <- list(Q = rep(NA, M),
                          synth_q = rep(NA, M),
                          synth_v = rep(NA,M),
                          rep_q = matrix(data = NA, nrow = M, ncol = R),
                          rep_v = matrix(data = NA, nrow = M, ncol = R)
    )
    
    ### Bootstrap
    orig_sample$wts <- N*orig_sample$wts/sum(orig_sample$wts)
    dsgn <- svydesign(ids=~1, data = orig_sample, weights =~wts)
    dsgn.rw <- as.svrepdesign(design=dsgn, type="subbootstrap",replicates=M)
    repwt <- as.matrix(dsgn.rw$repweights)
    ########
    
    for (i in seq(1:M)) {
      #create bootstrap sample for obtaining weights and then create population
      n_synth <- n
 
      wts <- repwt[,i]*orig_sample$wts
      wtind <- wts!=0
      wts <- wts[wtind]
      polyaY <- orig_sample$Y[wtind]
 
      wts <- (c*n)*(wts/sum(wts))
      bootsamp_size <- length(wts)

      synpop <- double()
      
      
      for (l in seq(1:L)) {
        pop_synth <- c(polyapost::wtpolyap(polyaY, wts, c*n - bootsamp_size))
        synpop <- c(synpop, pop_synth)
      }
      
      sample_synth <- sample(synpop, n_synth)
      #calculate relevant quantities
      synth_samples$Q[[i]] <- mean(synpop)
      synth_samples$synth_q[[i]] <- q_l <- mean(sample_synth)
      synth_samples$synth_v[[i]] <- v_l <- (1-n_synth/(c*n*L))*(sd(sample_synth))^2/n_synth
      
      
      #model synthetic dataset
      for (k in seq(1:R)) {
        rep <- rnorm(n_synth, q_l, sd(sample_synth))
        synth_samples$rep_q[i, k] <- mean(rep)
        synth_samples$rep_v[i, k] <- (1-n_synth/(c*n*L))*(sd(rep))^2/n_synth
      }
    }
    
    
    #store information
    Q <- synth_samples$Q
    q <- synth_samples$synth_q
    v <- synth_samples$synth_v
    q_rep <- synth_samples$rep_q
    v_rep <- synth_samples$rep_v
    #R=1  
    q_rep_1 <- synth_samples$rep_q[,1]
    v_rep_1 <- synth_samples$rep_v[,1]
    
    runs$Qbar[[j]] <- Qbar <- mean(Q)
    runs$qbar[[j]] <- qbar <- mean(q)
    runs$vbar[[j]] <- vbar <- mean(v)
    runs$qbar_rep[[j]] <- qbar_rep <- mean(q_rep)
    runs$vbar_rep[[j]] <- vbar_rep <- mean(v_rep)
    ql_bar <- rowMeans(q_rep)
    #R=1 
    runs$qbar_rep_1[[j]] <- qbar_rep_1 <- mean(q_rep_1)
    runs$vbar_rep_1[[j]] <- vbar_rep_1 <- mean(v_rep_1)
    
    #compute parameters for posterior from pop level estimates
    runs$Bm[[j]] <- Bm <- (1/(M-1))*sum((Qbar - Q)^2)
    runs$Qvar[[j]] <- Qvar <- (1 + (1/M))*Bm
    
    #compute parameters for approx posterior from sample level estimates
    runs$bm[[j]] <- bm <- sum((q-qbar)^2)/(M-1)
    runs$Tm[[j]] <- tm <- if_else((1 + (1/M))*bm - vbar <= 0, 
                                  (1 + (2/M))*vbar, 
                                  (1 + (1/M))*bm - vbar)
    runs$Tm_0[[j]]  = (1 + (1/M))*bm - vbar
    
    
    #compute parameters for approx posterior from replicate sample estimates
    runs$bm_rep[[j]] <- bm_rep <- sum((ql_bar - qbar_rep)^2)/(M-1)
    runs$bm_rep_1[[j]] <- bm_rep_1 <- var(q_rep_1)
    
    wl <- apply(q_rep, 1, var)
    runs$wbar[[j]] <- wbar <- mean(wl)
    
    runs$Bm_rep[[j]] <- Bm_rep <- bm_rep - vbar_rep - wbar/R
    tm_rep <- ((1/M) + 1)*bm_rep - vbar_rep - wbar/R
    runs$Tm_rep[[j]] <- tm_rep <- if_else(((1/M) + 1)*bm_rep - vbar_rep - wbar/R <= 0, 
                                          (1 + (2/M))*vbar_rep + (1/(R*M))*wbar, 
                                          ((1/M) + 1)*bm_rep - vbar_rep - wbar/R)
    
    runs$Tm_rep_0[[j]] <-  ((1/M) + 1)*bm_rep - vbar_rep - wbar/R
    
    runs$Bm_rep_1[[j]] <- Bm_rep_1 <- bm_rep_1 - 2*vbar_rep_1
    tm_rep_1 <- ((1/M) + 1)*bm_rep_1 - 2*vbar_rep_1
    runs$Tm_rep_1[[j]] <- tm_rep_1 <- if_else(((1/M) + 1)*bm_rep_1 - 2*vbar_rep_1 <= 0, 
                                              (1 + (3/M))*vbar_rep_1, 
                                              ((1/M) + 1)*bm_rep_1 - 2*vbar_rep_1)
    
    runs$Tm_rep_1_0[[j]] <- ((1/M) + 1)*bm_rep_1 - 2*vbar_rep_1
    
    # ci_pop <- c(qnorm(0.025, mean = Qbar, sd = sqrt(Qvar)), qnorm(0.975, mean = Qbar, sd = sqrt(Qvar)))
    # ci_synth <- c(qnorm(0.025, mean = qbar, sd = sqrt(tm)), qnorm(0.975, mean = qbar, sd = sqrt(tm)))
    # ci_rep <- c(qnorm(0.025, mean = qbar_rep, sd = sqrt(tm_rep)), qnorm(0.975, mean = qbar_rep, sd = sqrt(tm_rep)))
    # ci_rep_1 <- c(qnorm(0.025, mean = qbar_rep_1, sd = sqrt(tm_rep_1)), qnorm(0.975, mean = qbar_rep_1, sd = sqrt(tm_rep_1)))
    # ci_srs <- c(qnorm(0.025, mean = qboot, sd = sqrt(qbootvar)), qnorm(0.975, mean = qboot, sd = sqrt(qbootvar)))
    # ci_ht <- confint(htmean)
    # ci_srssyn <- c(qnorm(0.025, mean = qsyn, sd = sqrt(qsynvar)), qnorm(0.975, mean = qsyn, sd = sqrt(qsynvar)))
    
    ci_pop <- c(qt(0.025, M-1) *sqrt(Qvar) + Qbar, qt(0.975, M-1) *sqrt(Qvar) + Qbar)
    ci_synth <- c(qt(0.025, M-1) *sqrt(tm) + qbar, qt(0.975, M-1) *sqrt(tm) + qbar)
    ci_rep <- c(qt(0.025, M-1) *sqrt(tm_rep) + qbar_rep, qt(0.975, M-1) *sqrt(tm_rep) + qbar_rep) 
    ci_rep_1 <- c(qt(0.025, M-1) *sqrt(tm_rep_1) + qbar_rep_1, qt(0.975, M-1) *sqrt(tm_rep_1) + qbar_rep_1) 
    ci_srs <- c(qt(0.025, M-1) *sqrt(qbootvar) + qboot, qt(0.975, M-1) *sqrt(qbootvar) + qboot) 
    ci_ht <- c(qt(0.025, M-1) *SE(htmean) + coef(htmean), qt(0.975, M-1) *SE(htmean) + coef(htmean))
    ci_srssyn <- c(qt(0.025, M-1) *sqrt(qsynvar) + qsyn, qt(0.975, M-1) *sqrt(qsynvar) + qsyn) 
    
    runs$result_pop[[j]] <- if_else(ci_pop[[1]] <= true_mean & true_mean <= ci_pop[[2]], "T", "F")
    runs$result_synth[[j]] <- if_else(ci_synth[[1]] <= true_mean & true_mean <= ci_synth[[2]], "T", "F")
    runs$result_rep[[j]] <- if_else(ci_rep[[1]] <= true_mean & true_mean <= ci_rep[[2]], "T", "F")
    runs$result_rep_1[[j]] <- if_else(ci_rep_1[[1]] <= true_mean & true_mean <= ci_rep_1[[2]], "T", "F")
    runs$result_srs[[j]] <- if_else(ci_srs[[1]] <= true_mean & true_mean <= ci_srs[[2]], "T", "F")
    runs$result_ht[[j]] <- if_else(ci_ht[[1]] <= true_mean & true_mean <= ci_ht[[2]], "T", "F")
    runs$result_srssyn[[j]] <- if_else(ci_srssyn[[1]] <= true_mean & true_mean <= ci_srssyn[[2]], "T", "F")
    print(j)
  }
  
  
  
  results <- data.frame(Method = c("Sample Direct","Sample HT", "SRS Syn", "Population", "Synthetic Sample", "Synthetic Replicate Sample", "Synthetic Replicate Sample (R=1)"),
                        QEst_bias = c(mean(runs$qboot), mean(runs$htboot), mean(runs$qsyn), mean(runs$Qbar), mean(runs$qbar), mean(runs$qbar_rep), mean(runs$qbar_rep_1)) - true_mean,
                        QEst_bias_relpct = (c(mean(runs$qboot), mean(runs$htboot), mean(runs$qsyn), mean(runs$Qbar), mean(runs$qbar), mean(runs$qbar_rep), mean(runs$qbar_rep_1)) - true_mean)/true_mean * 100,
                        VarEst = c(mean(runs$qbootvar), mean(runs$htbootvar), mean(runs$qsynvar), mean(runs$Qvar), mean(runs$Tm), mean(runs$Tm_rep), mean(runs$Tm_rep_1)),
                        VarComp = c(var(runs$qboot), var(runs$htboot), var(runs$qsyn), var(runs$Qbar), var(runs$qbar), var(runs$qbar_rep), var(runs$qbar_rep_1)))
  
  results = results %>% mutate(Var_r = VarEst/VarComp)
  
  results = results %>% mutate(VarEst_raw=c(mean(runs$qbootvar), mean(runs$htbootvar), mean(runs$qsynvar), mean(runs$Qvar), mean(runs$Tm_0), mean(runs$Tm_rep_0), mean(runs$Tm_rep_1_0))
  )
  
  ci_results <- data.frame(Method = c("Sample Direct","Sample HT", "SRS Syn", "Population", "Synthetic Sample", "Synthetic Replicate Sample", "Synthetic Replicate Sample (R=1)"),
                           Contained = c(length(runs$result_srs[runs$result_srs == "T"]),
                                         length(runs$result_ht[runs$result_ht == "T"]),
                                         length(runs$result_srssyn[runs$result_srssyn == "T"]),
                                         length(runs$result_pop[runs$result_pop == "T"]),
                                         length(runs$result_synth[runs$result_synth == "T"]),
                                         length(runs$result_rep[runs$result_rep == "T"]),
                                         length(runs$result_rep_1[runs$result_rep_1 == "T"]))/run_num,
                           Contained_raw = c(mean(runs$result_srs[runs$qbootvar > 0]== "T"),
                                             mean(runs$result_ht[runs$htbootvar > 0]== "T"),
                                             mean(runs$result_srssyn[runs$qsynvar > 0]== "T"),
                                             mean(runs$result_pop[runs$Qvar > 0]== "T"),
                                             mean(runs$result_synth[runs$Tm_0 > 0]== "T"),
                                             mean(runs$result_rep[runs$Tm_rep_0 > 0]== "T"),
                                             mean(runs$result_rep_1[runs$Tm_rep_1_0 > 0]== "T")),
                           prop_neg =c(mean(runs$qbootvar<0), mean(runs$htbootvar<0), mean(runs$qsynvar<0), mean(runs$Qvar<0), mean(runs$Tm_0<0), mean(runs$Tm_rep_0<0), mean(runs$Tm_rep_1_0<0))
  )
  
  return(list(results,ci_results))
  
}

