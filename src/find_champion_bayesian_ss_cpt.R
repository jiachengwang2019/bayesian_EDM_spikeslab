#  #######################################################################
#       File-Name:      find_champion_bayesian_ss_cpt.R
#       Version:        R 3.6.2
#       Date:           Aug, 21 2020
#       Author:         Jiacheng Wang <Jiacheng.Wang@nbcuni.com>
#       Purpose:        find champion bayesian function with spike and slab prior and change points
#       Input Files:    NONE
#       Output Files:   NONE
#       Data Output:    NONE
#       Previous files: NONE
#       Dependencies:   NONE
#       Required by:    NONE
#       Status:         IN PROGRESS
#       Machine:        NBCU laptop
#  #######################################################################
find_champion_bayesian_ss_cpt <- function(data,
                                          show_variable = 'SC3_Impressions',
                                          agg_timescale = 'Date',
                                          log_transformation = 1,
                                          OOS_start = OOS_start,
                                          regressors,
                                          nr_samples,
                                          nr_burnin = round(nr_samples / 4, 0),
                                          max_changepoints = 0,
                                          changepoint_dates = NA,
                                          changepoints_minseglen = 1){
  #' @description finds best bayesian forecast model with spike and slab prior for the given data according to out-of-sample date 
  #' and conditions provided on regressor candidates
  #' 
  #' @param data R data-frame, typically an output from `create_aggregated_data()` function (see `data_prep_function.R`)
  #' @param show_variable denotes which variable to fit the model to, default is `SC3_Impressions`
  #' @param agg_timescale denotes the time-scale of the aggregated data - choices are `Week`, `Date`, `Hour` or `Half_Hr`
  #' @param log_transformation denotes whether to use log-transformation to the variable, choices are 0 (no) and 1 (yes), default is 1
  #' @param OOS_start the first date for the out-of-sample, format has to be `as.Date()`
  #' @param regressors vector of all regressors to be used, include `trend` in this vector if you want to use drift/trend in the model
  #' @param nr_samples number of iterations for MCMC sample
  #' @param nr_burnin number of iterations to discard as burn in period in MCMC algorithm, default is one quarter of total iterations
  #' @param max_changepoints no of changepoints to be used in piecewise linear trend, default is 0 (not to use piecewise linear trend)
  #' @param changepoint_dates dates where trend behavior changes, default is `NA` (changepoints to be found using statistical means)
  #' @param changepoint_minseglen integer giving the minimum segment length (no. of observations between changes) for the changepoint function, default is 1

  #' @return list containing the following:
  #' - bayesian optimal model parameter choice
  #' - new data adding prediction and error based on the optimal model
  #' - R dataframe with all details and model fit and forecasts for train/test data
  #' - changepoint dates for piecewise linear function
  
  if (sum(is.na(changepoint_dates))==0 & length(changepoint_dates)>max_changepoints){
    stop("no of changepoint-dates is more than allowed no of changepoints")
  }
  # check if "trend" in the regressors
  if (!"trend" %in% regressors & max_changepoints > 0){
    stop("'trend' has to be included in regressors to allow for piecewise linear trend")
  }
  # check if "intercept" in the regressors
  if (!"intercept" %in% regressors & max_changepoints > 0){
    warning("'intercept' should be included in regressors to get better prediction")
  }
  
  # prepare train and test data 
  if (agg_timescale == "Week"){
    train <- data %>% filter(Week < OOS_start)
    test <- data %>% filter(Week >= OOS_start)
    cp_timescale <- "Week"
    date_lag =  date_diff(data$Week[which(data$Week< OOS_start)])
    date_lag_test = predict_date_diff(data$Week[which(data$Week< OOS_start)],data$Week[which(data$Week>=OOS_start)])
  } else{
    train <- data %>% filter(Date < OOS_start)
    test <- data %>% filter(Date >= OOS_start)
    cp_timescale <- "Date"
    date_lag = date_diff(data$Date[which(data$Date< OOS_start)])
    date_lag_test = predict_date_diff(data$Date[which(data$Date< OOS_start)],data$Date[which(data$Date>=OOS_start)])
  }
  
  # find train/test length and data
  train_periods = nrow(train)
  test_periods = nrow(test)
  
  if (log_transformation == 1){
    train_series = log(as.numeric(unlist(train[,show_variable])) + 0.01)
    test_series = log(as.numeric(unlist(test[,show_variable])) + 0.01)
  } else{
    train_series = as.numeric(unlist(train[,show_variable]))
    test_series = as.numeric(unlist(test[,show_variable]))
  }
  
  cv_range = c(0.001, 0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.25, 0.5, 0.75,0.9,1,2,3,4,5,10,11,12,20)
  cv_no = length(cv_range)
  
  reg_train <- train
  reg_test <- test
  full_formula <- paste(show_variable,paste(regressors[! regressors %in% c("trend","intercept")],collapse = "+"),sep = "~")
  full_model <- lm(formula(full_formula),data = train)
  regressors_subset <- names(full_model$coefficients)[-1]
  reg_train <- train[,regressors_subset]
  reg_test <- test[,regressors_subset]
  
  if ("intercept" %in% regressors){
    intercept <- rep(1,train_periods)
    reg_train <- cbind(intercept,reg_train)
    intercept <- rep(1,test_periods)
    reg_test <- cbind(intercept,reg_test)
  }
  if ("trend" %in% regressors){
    x1 = c(1:train_periods)
    reg_train <- cbind(reg_train,x1)
    colnames(reg_train)[ncol(reg_train)] = "trend"
    x2 = seq(train_periods+1,train_periods+test_periods,1)
    reg_test <- cbind(reg_test,x2)
    colnames(reg_test)[ncol(reg_test)] = "trend"
  }
  
  # adjust xreg matrices as needed
  idx = which(colSums(reg_train)==0)
  if (length(idx)>0){
    reg_train <- reg_train[,-idx]
    reg_test <- reg_test[,-idx]
  }
  regset <- names(reg_train)
  
  # prepare piecewise linear regression if adding any change points
  cp_date_out = NA
  cpt_no = 0
  knotpoints = as.list(rep(NA,max_changepoints))
  cpt_add = ifelse(max_changepoints > 0,1,0)
  cpt_add2 = ifelse(max_changepoints > 0,max_changepoints,1)
  if (max_changepoints > 0){
    cpt_no = c(rep(0,max_changepoints))
    cp_date_out = as.list(rep(NA,max_changepoints))
    if (is.na(changepoint_dates)==T){
      # experiment with at most the assumed no of "max_changepoints"
      for (i in 1:max_changepoints){
        options(warn=-1)
        changepoints <- cpt.mean(as.numeric(unlist(train[,show_variable])),
                                 method = "BinSeg",
                                 pen.value = "2*log(n)",
                                 Q = i,
                                 minseglen = changepoints_minseglen)
        options(warn=0)
        knotpoints[[i]] = cpts(changepoints)
        cpt_no[i] = length(knotpoints[[i]])
      }
    }else{
      knotpoints = which(train[,cp_timescale] %in% changepoint_dates)
      cpt_no = length(knotpoints)
    }
    train_cpt = as.list(rep(NA,max_changepoints))
    for (i in 1:max_changepoints){
      if (cpt_no[i]>0){
        train_cpt[[i]] = as.list(rep(NA, cpt_no[i]+1))
        for (j in 1:(cpt_no[i]+1)){
          if(j !=1 & j!=(cpt_no[i]+1)){
            train_cpt[[i]][[j]] = (knotpoints[[i]][j-1]+1):knotpoints[[i]][j]
          }else{
            train_cpt[[i]][[1]] = 1:knotpoints[[i]][1]
            train_cpt[[i]][[cpt_no[i]+1]] = (knotpoints[[i]][cpt_no[i]]+1):train_periods
          }
        }
        changepoint_dates = as.data.frame(train[knotpoints[[i]],cp_timescale])
        cp_date_out[[i]] = paste(as.character(changepoint_dates[,cp_timescale]),collapse = ",")
      }
    }
  }
  MAPE_train = c(rep(NA,cv_no*(1-cpt_add) + cv_no*max_changepoints*cpt_add))
  SMAPE_train = c(rep(NA,cv_no*(1-cpt_add) + cv_no*max_changepoints*cpt_add))
  MAPE_test = c(rep(NA,cv_no*(1-cpt_add) +cv_no*max_changepoints*cpt_add))
  SMAPE_test = c(rep(NA,cv_no*(1-cpt_add) + cv_no*max_changepoints*cpt_add))
  
  
  a1 = 2; a2 = 1; theta = .5; a = 1;  b = 1;  s = 1
  lm_model <- lm(train_series ~ as.matrix(reg_train)-1)
  beta_old<- coef(lm_model)
  pi_old <- c(1,rep(0,p-2), 1)
  sigma2_old <- var(predict(lm_model) - train_series)
  tau2_old <- 1
  theta_old <- 0.5
  
  res_param <- list()
  res_train_fit <- list()
  res_test_fit <- list()
  p <- ncol(as.matrix(reg_train))
  for (l in 1:length(cv_range)){
    cat('total:',length(cv_range),'+++ Now:',l,'\n')
    alpha1 = cv_range[l]
    variance_matrix_train = diag(train_periods)
    variance_matrix_test = diag(test_periods)
    Gamma = variance_matrix_train%*% (exponential_decay(date_lag,alpha1) + diag(dim(as.matrix(reg_train))[1]))
    GI <- qr.solve(Gamma)
    predict_cov_train = variance_matrix_train%*%exponential_decay(date_lag,alpha1) %*% variance_matrix_train
    predict_date_cov = variance_matrix_test%*%exponential_decay(date_lag_test,alpha1) %*% variance_matrix_train
    
    if (max_changepoints > 0){
      for (j in 1:max_changepoints){
        train_series_cpt = as.list(rep(NA,cpt_no[j]+1))
        train_series_cpt_outside = as.list(rep(NA,cpt_no[j]+1))
        reg_train_cpt = as.list(rep(NA,cpt_no[j]+1))
        reg_train_cpt_outside = as.list(rep(NA,cpt_no[j]+1))
        variance_cpt = as.list(rep(NA,cpt_no[j]+1))
        variance_cpt_outside = as.list(rep(NA,cpt_no[j]+1))
        covariance_interaction = as.list(rep(NA,cpt_no[j]+1))
        beta_cpt_old = as.list(rep(NA,cpt_no[j]+1))
        beta_cpt_outside_old = as.list(rep(NA,cpt_no[j]+1))
        ed1_inv = as.list(rep(NA,cpt_no[j]+1))
        ed3_inv = as.list(rep(NA,cpt_no[j]+1))
        ed5_inv = as.list(rep(NA,cpt_no[j]+1))
        res <- matrix(NA, nrow = nr_samples, ncol = 2*p*(cpt_no[[j]]+1) + 1 + 1 + 1)
        colnames(res) <- c(
          paste0('pi', seq(p*(cpt_no[[j]]+1))),
          paste0('beta', seq(p*(cpt_no[[j]]+1))),
          'sigma2', 'tau2', 'theta'
        )
        res[1, ] <- c(rep(c(1,rep(0,p-2), 1),cpt_no[[j]]+1), rep(as.vector(beta_old),cpt_no[[j]]+1), sigma2_old, 1, .5)
        train_fit <- matrix(NA, nrow = nr_samples-1, ncol = train_periods)
        test_fit <- matrix(NA, nrow = nr_samples-1, ncol = test_periods)
        for (k in 1:(cpt_no[j]+1)){
          train_series_cpt[[k]] = train_series[train_cpt[[j]][[k]]]
          train_series_cpt_outside[[k]] = train_series[-train_cpt[[j]][[k]]]
          reg_train_cpt[[k]] = reg_train[train_cpt[[j]][[k]],]
          reg_train_cpt_outside[[k]] = reg_train[-train_cpt[[j]][[k]],]
          variance_cpt[[k]] = GI[c(train_cpt[[j]][[k]]),c(train_cpt[[j]][[k]])]
          variance_cpt_outside[[k]] = GI[-train_cpt[[j]][[k]],-train_cpt[[j]][[k]]]
          covariance_interaction[[k]] = GI[train_cpt[[j]][[k]],-train_cpt[[j]][[k]]]
          ed1_inv[[k]] = t(reg_train_cpt[[k]])%*% variance_cpt[[k]] %*% as.matrix(reg_train_cpt[[k]])
          ed3_inv[[k]] = t(reg_train_cpt[[k]])%*% variance_cpt[[k]] %*% train_series_cpt[[k]]
          ed5_inv[[k]] = t(reg_train_cpt[[k]])%*%covariance_interaction[[k]]
        }
        for (i in 2:nr_samples){
          cat('total:',nr_samples,'+++ Now:',i,'\n')
          pi_prev <- res[i-1, seq(p*(cpt_no[[j]]+1))]
          beta_prev <- res[i-1, seq(p*(cpt_no[[j]]+1) + 1, 2*p*(cpt_no[[j]]+1))]
          sigma2_prev <- res[i-1, ncol(res) - 2]
          tau2_prev <- res[i-1, ncol(res) - 1]
          theta_prev <- res[i-1, ncol(res)]
          
          ## Start sampling from the conditional posterior distributions
          ##############################################################
          # sample theta from a Beta
          theta_new <-  (a + sum(pi_prev))/(a + sum(pi_prev) + b + sum(1 - pi_prev))
          # sample sigma2 from an Inverse-Gamma
          resid_cpt = as.list(rep(NA,cpt_no[j]+1))
          for (m in 1:(cpt_no[j]+1)){
            beta_cpt_old[[m]] = res[1,(p*(cpt_no[[j]]+1)+1+(m-1)*p):(p*(cpt_no[[j]]+1)+m*p)]
            resid_cpt[[m]] = train_series_cpt[[m]]-as.matrix(reg_train_cpt[[m]])%*% beta_cpt_old[[m]]
          }
          sigma2_new <- as.numeric((a2 + t(unlist(resid_cpt)) %*% GI %*% unlist(resid_cpt) / 2 + sum(beta_prev^2)/(2*tau2_prev))/(a1 + (train_periods + sum(pi_prev))/2-1))
          # sample tau2 from an Inverse Gamma
          tau2_new <- as.numeric(( s^2/2 + sum(beta_prev^2)/(2*sigma2_new))/(1/2 + 1/2 * sum(pi_prev) - 1))
          
          # sample beta from multivariate Gaussian
          beta_cpt_new = as.list(rep(NA,cpt_no[j]+1))
          for (ll in 1:(cpt_no[j]+1)){
            beta_cpt_cov <- qr.solve(ed1_inv[[ll]]/sigma2_new + diag(1/(tau2_new*sigma2_new), p))
            beta_cpt_new[[ll]] <- as.vector(beta_cpt_cov%*%((ed3_inv[[ll]] + ed5_inv[[ll]]%*%unlist(resid_cpt[-ll]))/sigma2_new))
          }
          
          # sample each pi_j in random order
          for (mm in 1:(cpt_no[j]+1)){
            for (nn in sample(seq(p)[2:(p-1)])){
              pi0 <- pi_prev[(1+p*(mm-1)):(p*mm)]
              pi0[nn] <- 0
              bp0 <- as.vector(t(beta_cpt_new[[mm]] * pi0))
              # compute the z variables and the conditional variance
              xj <- as.matrix(reg_train_cpt[[mm]])[, nn]
              z <- train_series_cpt[[mm]] - as.matrix(reg_train_cpt[[mm]]) %*% bp0
              cond_var <- t(xj) %*% variance_cpt[[mm]] %*% xj + 1/(tau2_new) 
              
              # compute chance parameter of the conditional posterior of pi_j (Bernoulli)
              l0 <- log(1 - theta_new)
              #cat('l0:',l0, '\n')
              l1 <- (
                log(theta_new) - .5 * log(tau2_new*sigma2_new) +
                  (t(z) %*% variance_cpt[[mm]] %*% xj + t(xj)%*%covariance_interaction[[mm]]%*% unlist(resid_cpt[-mm]))^2 / (2/sigma2_new*cond_var) + .5 * log(sigma2_new / cond_var)
              )
              #cat('l1:',l1, '\n')
              pi_prev[p*(mm-1)+nn] <- rbinom(1, 1, exp(l1) / (exp(l1) + exp(l0)))
            }
            beta_cpt_new[[mm]] = beta_cpt_new[[mm]]*pi_prev[(1+p*(mm-1)):(p*mm)]
          }
          pi_new <- pi_prev
          beta_new <- unlist(beta_cpt_new)
          res[i, ] <- c(pi_new, beta_new, sigma2_new, tau2_new, theta_new)
          # Estimate for the train & test
          # train
          resid_cpt_new = list()
          for (xx in 1:(cpt_no[j]+1)){
            resid_cpt_new[[xx]] = train_series_cpt[[xx]] - as.matrix(reg_train_cpt[[xx]])%*%as.vector(beta_cpt_new[[xx]])
          }
          resid_train_error <- sigma2_new*predict_cov_train%*%GI%*%as.vector(unlist(resid_cpt_new))
          evaluate_train_tmp = as.list(rep(NA,cpt_no[j]+1))
          for(yy in 1:(cpt_no[j]+1)){
            evaluate_train_tmp[[yy]] = as.matrix(reg_train_cpt[[yy]])%*%beta_cpt_new[[yy]]
          }
          evaluate_train = unlist(evaluate_train_tmp) + resid_train_error
          train_fit[i-1,] <- evaluate_train
          
          # test
          resid_error = sigma2_new*predict_date_cov%*%GI%*%unlist(resid_cpt)
          evaluate_test <- as.matrix(reg_test)%*%beta_cpt_new[[cpt_no[j]+1]]+ resid_error
          test_fit[i-1,] <- evaluate_test
        }
        if (nr_burnin > 0){
          parameter = res[-seq(nr_burnin), ]
          train_fit = train_fit[-seq(nr_burnin-1), ]
          test_fit = test_fit[-seq(nr_burnin-1), ]
        }else{
          parameter = res
          train_fit = train_fit
          test_fit = test_fit
        }
        train_predict <- apply(train_fit,2,mean)
        test_predict <- apply(test_fit,2,mean)
        if (log_transformation == 0){
          MAPE_train[max_changepoints*(l-1)+j] <- mape(train_series,train_predict)
          SMAPE_train[max_changepoints*(l-1)+j] <-smape(train_series,train_predict)
          MAPE_test[max_changepoints*(l-1)+j] <-mape(test_series,test_predict)
          SMAPE_test[max_changepoints*(l-1)+j] <-smape(test_series,test_predict)
        }else{
          MAPE_train[max_changepoints*(l-1)+j] <-mape(exp(train_series),exp(train_predict))
          SMAPE_train[max_changepoints*(l-1)+j] <-smape(exp(train_series),exp(train_predict))
          MAPE_test[max_changepoints*(l-1)+j] <-mape(exp(test_series),exp(test_predict))
          SMAPE_test[max_changepoints*(l-1)+j] <-smape(exp(test_series),exp(test_predict))
        }
        res_param[[max_changepoints*(l-1)+j]] <- parameter
        res_train_fit[[max_changepoints*(l-1)+j]] <- train_fit
        res_test_fit[[max_changepoints*(l-1)+j]] <- test_fit
      }
    }else{
      p <- ncol(as.matrix(reg_train))
      res <- matrix(NA, nrow = nr_samples, ncol = 2*p + 1 + 1 + 1)
      train_fit <- matrix(NA, nrow = nr_samples-1, ncol = train_periods)
      test_fit <- matrix(NA, nrow = nr_samples-1, ncol = test_periods)
      colnames(res) <- c(
        paste0('pi', seq(p)),
        paste0('beta', seq(p)),
        'sigma2', 'tau2', 'theta'
      )
      lm_model <- lm(train_series ~ as.matrix(reg_train)-1)
      res[1, ] <- c(1,rep(0,p-2), 1, coef(lm_model), var(predict(lm_model) - train_series), 1, .5)
      GI <- qr.solve(Gamma)
      XtGIX <- t(as.matrix(reg_train)) %*% GI %*% as.matrix(reg_train)
      XtGIy <- t(as.matrix(reg_train)) %*% GI %*% train_series
      # we start running the Gibbs sampler
      for (i in seq(2, nr_samples)) {
        # cat('total:',nr_samples,'+++Now:',i,'\n')
        # first, get all the values of the previous time point
        pi_prev <- res[i-1, seq(p)]
        beta_prev <- res[i-1, seq(p + 1, 2*p)]
        sigma2_prev <- res[i-1, ncol(res) - 2]
        tau2_prev <- res[i-1, ncol(res) - 1]
        theta_prev <- res[i-1, ncol(res)]
        
        ## Start sampling from the conditional posterior distributions
        ##############################################################
        
        # sample theta from a Beta
        #!! same as the previous general case
        # theta_new <- rbeta(1,a + sum(pi_prev), b + sum(1 - pi_prev) )
        theta_new <-  (a + sum(pi_prev))/(a + sum(pi_prev) + b + sum(1 - pi_prev))
        
        # sample sigma2 from an Inverse-Gamma
        #!! different from previous general case
        err <- train_series - as.matrix(reg_train) %*% beta_prev
        #sigma2_new <- rinvgamma(1, a1 + n_train/2, a2 + t(err) %*% GI %*% err / 2) 
        sigma2_new <- as.numeric((a2 + t(err) %*% GI %*% err / 2)/(a1 + train_periods/2-1))
        
        # sample tau2 from an Inverse Gamma
        # tau2_new <- 1 / rgamma(
        #   1, 1/2 + 1/2 * sum(pi_prev),
        #   s^2/2 + t(beta_prev) %*% beta_prev / (2*sigma2_new)
        # )
        tau2_new <- as.numeric(( s^2/2 + t(beta_prev) %*% beta_prev / (2*sigma2_new))/(1/2 + 1/2 * sum(pi_prev) - 1))
        
        # sample beta from multivariate Gaussian
        # !! different from previous general case
        beta_cov <- qr.solve(as.matrix.data.frame(XtGIX)/sigma2_new   + diag(1/(tau2_new*sigma2_new), p))
        beta_mean <- beta_cov %*% XtGIy * (1/sigma2_new)
        #beta_new <- mvtnorm::rmvnorm(1, beta_mean, beta_cov)
        beta_new <- as.vector(beta_mean)
        
        # sample each pi_j in random order
        for (j in sample(seq(p)[2:(p-1)])) {
          
          # get the betas for which beta_j is zero
          pi0 <- pi_prev
          pi0[j] <- 0
          bp0 <- as.vector(t(beta_new * pi0))
          
          # compute the z variables and the conditional variance
          xj <- as.matrix(reg_train)[, j]
          z <- train_series - as.matrix(reg_train) %*% bp0
          cond_var <- t(xj) %*% GI %*% xj + 1/(tau2_new) 
          
          # compute chance parameter of the conditional posterior of pi_j (Bernoulli)
          l0 <- log(1 - theta_new)
          l1 <- (
            log(theta_new) - .5 * log(tau2_new*sigma2_new) +
              (t(z) %*% GI %*% xj)^2 / (2/sigma2_new*cond_var) + .5 * log(sigma2_new / cond_var)
          )
          # sample pi_j from a Bernoulli
          pi_prev[j] <- rbinom(1, 1, exp(l1) / (exp(l1) + exp(l0)))
        }
        
        pi_new <- pi_prev
        # add new samples
        res[i, ] <- c(pi_new, beta_new*pi_new, sigma2_new, tau2_new, theta_new)
        # Estimate for the train & test
        beta_new <- beta_new*pi_new
        # train
        resid <- train_series - ts(as.matrix(reg_train)%*%as.vector(beta_new))
        resid_train_error <- sigma2_new*predict_cov_train%*%GI%*%resid
        evaluate_train <- as.matrix(reg_train)%*%as.vector(beta_new) + resid_train_error
        train_fit[i-1,] <- evaluate_train
        # test
        resid_error <- sigma2_new*predict_date_cov%*%GI%*%resid
        evaluate_test <- as.matrix(reg_test)%*%as.vector(beta_new) + resid_error
        test_fit[i-1,] <- evaluate_test
      }
      ########## Finished!
      if (nr_burnin > 0){
        parameter = res[-seq(nr_burnin), ]
        train_fit = train_fit[-seq(nr_burnin-1), ]
        test_fit = test_fit[-seq(nr_burnin-1), ]
      }else{
        parameter = res
        train_fit = train_fit
        test_fit = test_fit
      }
      train_predict <- apply(train_fit,2,mean)
      test_predict <- apply(test_fit,2,mean)
      if (log_transformation == 0){
        MAPE_train[l] <-mape(train_series,train_predict)
        SMAPE_train[l] <-smape(train_series,train_predict)
        MAPE_test[l] <-mape(test_series,test_predict)
        SMAPE_test[l] <-smape(test_series,test_predict)
      }else{
        MAPE_train[l] <-mape(exp(train_series),exp(train_predict))
        SMAPE_train[l] <-smape(exp(train_series),exp(train_predict))
        MAPE_test[l] <-mape(exp(test_series),exp(test_predict))
        SMAPE_test[l] <-smape(exp(test_series),exp(test_predict))
      }
    }
    res_param[[l]] <- parameter
    res_train_fit[[l]] <- train_fit
    res_test_fit[[l]] <- test_fit
  }
  MAPE_total = MAPE_train + MAPE_test
  SMAPE_total = SMAPE_train + SMAPE_test
  
  # find the optimal model!!
  if (max_changepoints >0){
    alpha_optimal = cv_range[ceiling(which.min(MAPE_test)/cpt_add2)]
    all_model_detail = c(rep(row.names(all_model),rep(cpt_add2,nrow(all_model))))
  }else{
    alpha_optimal <- cv_range[which.min(MAPE_test)]
    ## All model details
    all_model_detail = expand.grid(cv_range)
  }
  train_predict <- apply(res_train_fit[[which.min(MAPE_test)]],2,mean)
  test_predict <- apply(res_test_fit[[which.min(MAPE_test)]],2,mean)
  if (log_transformation == 0){
    train_predict_champion <- train_predict
    test_predict_champion <- test_predict
  }else{
    train_predict_champion <- exp(train_predict)
    test_predict_champion <- exp(test_predict)
  }
  ## get the champion_model function
  optimal_model = list()
  optimal_model$alpha <- alpha_optimal
  post_means <- apply(res_param[[which.min(MAPE_test)]], 2, mean)
  res_table_optimal <- cbind(
    post_means[grepl('beta', names(post_means))],
    post_means[grepl('pi', names(post_means))]
  )
  rownames(res_table_optimal) <- rep(regset, max_changepoints + 1)
  colnames(res_table_optimal) <- c('Post. Mean', 'Post. Inclusion')
  optimal_model$regressors <- res_table_optimal
  optimal_model$sigma2 <- as.numeric(post_means[length(post_means)-2])
  optimal_model$tau2 <- as.numeric(post_means[length(post_means)-1])
  optimal_model$theta <- as.numeric(post_means[length(post_means)])
  
  ## get the full_data_champion result
  combinedfit_champion <- as.data.frame(append(train_predict_champion , test_predict_champion))
  names(combinedfit_champion) = "Predict"
  full_data_champion <- as.data.frame(cbind(as.data.frame(data),combinedfit_champion)) 
  full_data_champion$Stream = show_variable
  full_data_champion$Error = ((full_data_champion[,show_variable] - full_data_champion$Predict)/full_data_champion[,show_variable])*100
  full_data_champion$APE = abs(full_data_champion$Error)
  detail_entry = data.frame(all_model_detail,
                            MAPE_train,
                            MAPE_test,
                            MAPE_total,
                            SMAPE_train,
                            SMAPE_test,
                            SMAPE_total,
                            log_transformation,
                            paste(regset,collapse = ","),
                            stringsAsFactors = T)
  names(detail_entry) = c('alpha',
                          "train_MAPE","test_MAPE","total_MAPE",
                          "train_SMAPE","test_SMAPE","total_SMAPE",
                          "log_transformation",
                          "regressors")
  detail_entry$champion = 0
  minidx_mape = which.min(detail_entry$test_MAPE)
  minidx_smape = which.min(detail_entry$test_SMAPE)
  if (minidx_mape != minidx_smape){
    detail_entry$champion[minidx_mape] = "MAPE_CHAMPION"
    detail_entry$champion[minidx_smape] = "SMAPE_CHAMPION"
  }else{
    detail_entry$champion[minidx_mape] = "CHAMPION"
  }
  detail_entry$ranking = 'Average'
  detail_entry$ranking[detail_entry$test_MAPE < 15] = 'Good'
  detail_entry$ranking[detail_entry$test_MAPE > 2*15] = 'Bad'
  return(list(champion_model = optimal_model,
              champion_result = full_data_champion,
              all_model_details = detail_entry,
              changepoint_dates = changepoint_dates))
}

#::::::::::::::::::::::::::: dependent functions
############## Dependent function
exponential_decay = function(date_diff, alpha){
  nrows = dim(date_diff)[1]
  ncols = dim(date_diff)[2]
  a = matrix(rep(0,nrows*ncols),nrow = nrows)
  a = exp(-alpha*date_diff)
  return(a)
}

date_diff = function(Date){
  nrow = length(Date)
  a = c(rep(0,nrow))
  b = matrix(rep(0,nrow*nrow),nrow=nrow)
  for (i in c(1:nrow)){
    a[i] = as.numeric(Date[i] - Date[1])
  }
  b[1,] = a
  for (j in c(2:nrow)){
    b[j,] = abs(a-a[j])
  }
  return(b)
}

predict_date_diff = function(Date1,Date2){
  nrow = length(Date2)
  ncol = length(Date1)
  a = matrix(rep(NA,nrow*ncol),nrow = nrow)
  for (i in c(1:nrow)){
    a[i,] = Date2[i]-Date1
  }
  return(a)
}

