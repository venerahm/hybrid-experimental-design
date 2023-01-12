#############################################################
# Functions modeled from Jamie Yap and Tianchen Qian's code #
#############################################################

# A utility function to find the change location of a vector.
# This is used to handle different length of observations among individuals.
### UPDATE: Allows for multiple probability values varying by group (responders vs nonresponders), for differing days within a study without invoking availability and includes output for graphing

find_change_location <- function(v){
  n <- length(v)
  if (n <= 1) {
    stop("The vector need to have length > 1.")
  }
  return(c(1, 1 + which(v[1:(n-1)] != v[2:n])))
}

binary_outcome_moderated_effect_resp <- function(
  dta,
  control_var,
  moderator,
  id_var,
  day_var = "studyday.0start",
  trt_var = "intervention",
  outcome_var = "outcome",
  avail_var = NULL,
  prob_treatment = "prob",
  significance_level = 0.05)
{
  ############## description ###############
  ##
  ## This function estimates the moderated treatment effect for binary outcome,
  ## and provides variance estimate, test statistics (t-test and F-test), and p-values.
  ##
  ## It incorporates two methods for small sample correction:
  ## 1) the usage of "Hat" matrix in the variance estimate (as in Mancl & DeRouen 2001)
  ## 2) the usage of t-distribution or F-distribution critical value with corrected degrees of freedom
  ##    (as in Liao et al. 2015)
  ##
  ## It is used as an internal function for SARA analysis wrapper functions.
  ##
  ## Note: dta needs to be sorted by id_var then day_var
  ##       (currently this is handled in each of the wrapper functions)
  
  ############## arguments ###############
  ##
  ## dta.............the data set in long format
  ## control_var...........vector of variable names used to reduce noise (Z in the model),
  ##                       could be NULL (no control covariates)
  ## moderator.............vector of variable names as effect modifiers (X in the model),
  ##                       could be NULL (no effect modifier)
  ## id_var................variable name for subject id (to distinguish between subjects in dta)
  ## day_var...............variable name for day in study
  ## trt_var...............variable name for treatment indicator
  ## outcome_var...........variable name for outcome variable
  ## avail_var.............variable name for availability variable
  ##                       NULL (default) means always-available
  ## prob_treatment........probability of treatment 
  ##### ALLOWS MULTIPLE VALUES #####
  ## significance_level....significance level for the hypothesis testing (default to 0.05)
  
  
  ############## return value ###############
  ##
  ## This function returns a list of the following components:
  ##
  ## beta_hat..............estimated beta
  ## alpha_hat.............estimated alpha
  ## beta_se...............standard error for beta_hat
  ## alpha_se..............standard error for alpha_hat
  ## beta_se_ssa...........standard error for beta_hat, with small sample correction (hat matrix)
  ## alpha_se..............standard error for alpha_hat, with small sample correction (hat matrix)
  ## test_result_t.........(two-sided) t-test result for each entry in beta_hat, which is a list consisting of test_stat, critical_value, p_value
  ## test_result_f.........F-test result for beta = 0, which is a list consisting of test_stat, critical_value, p_value
  ## varcov................estimated variance-covariance matrix for (beta_hat, alpha_hat)
  ## varcov_ssa............estimated variance-covariance matrix for (beta_hat, alpha_hat), with small sample correction (hat matrix)
  ## dims..................a list of p (dim of moderator) and q (dim of control_var), which includes the intercepts if added
  
  
  ############## part 1 :: preparation ###############
  
  # location of each individual's data in dta, when looping over individuals
  person_first_index <- find_change_location(dta[, id_var])
  if (length(person_first_index) != length(unique(dta[, id_var]))) {
    stop("The length of person_first_index doesn't equal the number of unique id_var's.")
  }
  person_data_location <- c(person_first_index, nrow(dta) + 1)
  
  # gather variables
  if (is.null(avail_var)) {
    avail <- rep(1, nrow(dta))
  } else {
    avail <- dta[, avail_var]
  }
  
  A <- dta[, trt_var]
  if (any(is.na(A[avail == 1]))) {
    stop("Treatment indicator is NA where availability = 1.")
  }
  A[avail == 0] <- 0
  
  cA <- A - dta[,prob_treatment] # centered A
  Xdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, moderator] ) )   # X design matrix, intercept added
  Zdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, control_var] ) ) # Z design matrix, intercept added
  Y <- dta[, outcome_var]
  
  n <- sample_size <- length(unique(dta[, id_var]))
  
  p <- ncol(Xdm) # dimension of beta
  q <- ncol(Zdm) # dimension of alpha
  
  Xnames <- c("Intercept", moderator)
  Znames <- c("Intercept", control_var)
  
  
  ############## part 2 :: estimate beta and alpha ###############
  
  estimating_function <- function(theta) {
    # function to be used in solver
    # theta is a vector of length: length(control_var) + 2
    beta <- as.matrix(theta[1:p])
    alpha <- as.matrix(theta[(p+1):(p+q)])
    
    exp_Zdm_alpha <- exp(Zdm %*% alpha)
    exp_negAXdm_beta <- exp(- A * (Xdm %*% beta))
    residual <- Y * exp_negAXdm_beta - exp_Zdm_alpha
    
    ef <- rep(NA, length(theta)) # value of estimating function
    for (i in 1:p) {
      ef[i] <- sum( residual * avail * cA * Xdm[, i])
    }
    for (i in 1:q) {
      ef[p + i] <- sum( residual * avail * exp_Zdm_alpha * Zdm[, i])
    }
    
    ef <- ef / sample_size
    return(ef)
  }
  
  # if come across solver error, try change the initial values here
  beta_initial <- rep(0, p)
  alpha_initial <- rep(0, q)
  root <- multiroot(estimating_function, c(beta_initial, alpha_initial), useFortran = FALSE)
  if (p == 1) {
    beta_root <- root$root[1]
  } else {
    beta_root <- as.matrix(root$root[1:p])
  }
  if (q == 1) {
    alpha_root <- root$root[(p+1):(p+q)]
  } else {
    alpha_root <- as.matrix(root$root[(p+1):(p+q)])
  }
  
  
  ############## part 3 :: estimate standard error by asymptotics ###############
  
  # Compute M_n matrix (M_n is the empirical expectation of the derivative of the estimating function)
  
  Mn_summand <- array(NA, dim = c(nrow(dta), p+q, p+q))
  for (it in 1:nrow(dta)) {
    
    # this is to make R code consistent whether X_it, Z_it contains more entries or is just 1.        
    if (p == 1) {
      Xbeta <- Xdm[it, ] * beta_root
    } else {
      Xbeta <- Xdm[it, ] %*% beta_root
    }
    if (q == 1) {
      Zalpha <- Zdm[it, ] * alpha_root
    } else {
      Zalpha <- Zdm[it, ] %*% alpha_root
    }
    
    exp_Zalpha <- as.vector(exp(Zalpha))
    exp_negAXbeta <- as.vector(exp(- A[it] * Xbeta))
    exp_2Zalpha <- as.vector(exp(2 * Zalpha))
    Mn_summand[it, 1:p, 1:p] <-
      - as.numeric(Y[it] * exp_negAXbeta * A[it] * cA[it]) * (Xdm[it, ] %o% Xdm[it, ]) * avail[it]
    Mn_summand[it, 1:p, (p+1):(p+q)] <-
      - as.numeric(cA[it] * exp_Zalpha) * (Xdm[it, ] %o% Zdm[it, ]) * avail[it]
    Mn_summand[it, (p+1):(p+q), 1:p] <-
      - as.numeric(Y[it] * exp_negAXbeta * A[it] * exp_Zalpha) * (Zdm[it, ] %o% Xdm[it, ]) * avail[it]
    Mn_summand[it, (p+1):(p+q), (p+1):(p+q)] <-
      as.numeric(Y[it] * exp_negAXbeta * exp_Zalpha - 2 * exp_2Zalpha) * (Zdm[it, ] %o% Zdm[it, ]) * avail[it]
  }
  Mn <- apply(Mn_summand, c(2,3), sum) / sample_size
  Mn_inv <- solve(Mn)
  
  # Compute \Sigma_n matrix (\Sigma_n is the empirical variance of the estimating function)
  
  exp_Zdm_alpha <- exp(Zdm %*% alpha_root)
  exp_negAXdm_beta <- exp(- A * (Xdm %*% beta_root))
  residual <- Y * exp_negAXdm_beta - exp_Zdm_alpha
  Sigman_summand <- matrix(NA, nrow = nrow(dta), ncol = p+q)
  for (i in 1:p) {
    Sigman_summand[, i] <- residual * avail * cA * Xdm[, i]
  }
  for (i in 1:q) {
    Sigman_summand[, p+i] <- residual * avail * exp_Zdm_alpha * Zdm[, i]
  }
  
  Sigman <- matrix(0, nrow = p+q, ncol = p+q)
  for (user in 1:sample_size) {
    rows_for_this_user <- person_data_location[user]:(person_data_location[user+1] - 1)
    # changed to allow for differing days within a study without invoking availability (forces vectors to be matrices)
    summand <- matrix(Sigman_summand[rows_for_this_user, ], nrow = length(rows_for_this_user))
    summand <- colSums(summand)
    summand <- summand %o% summand
    Sigman <- Sigman + summand
  }
  Sigman <- Sigman / sample_size
  
  # Compute the asymptotic variance matrix ( this is on the scale of \sqrt{n}(\hat{\beta} - \beta) )
  
  asymp_varcov <- Mn_inv %*% Sigman %*% t(Mn_inv)
  asymp_var <- diag(asymp_varcov)
  
  # get the standard error for beta and alpha from the asymptotic variance matrix
  
  beta_se <- sqrt(asymp_var[1:p] / sample_size)
  alpha_se <- sqrt(asymp_var[(p+1):(p+q)] / sample_size)
  
  
  ############## part 4 :: estimate standard error with small sample correction ###############
  
  # construct the new "meat" in the sandwich estimator (to replace \Sigma_n)
  
  meat <- 0
  for (user in 1:sample_size) {
    # preparation
    # the following 6 lines have been updated to allow for differing days within a study without invoking availability (forces vectors to be matrices)
    rows_for_this_user <- person_data_location[user]:(person_data_location[user+1] - 1)
    A_user <- matrix(A[rows_for_this_user], nrow = length(rows_for_this_user)) # A_i1 to A_iT
    cA_user <- matrix(cA[rows_for_this_user], nrow = length(rows_for_this_user)) # centered A
    Xdm_user <- matrix(Xdm[rows_for_this_user, ], nrow = length(rows_for_this_user)) # X_i1 to X_iT
    Zdm_user <- matrix(Zdm[rows_for_this_user, ], nrow = length(rows_for_this_user)) # Z_i1 to Z_iT
    Y_user <- matrix(Y[rows_for_this_user], nrow = length(rows_for_this_user)) # Y_i1 to Y_iT
    avail_user <- avail[rows_for_this_user]
    
    # this is to make R code consistent whether X_it, Z_it contains more entries or is just 1.    
    if (p == 1) {
      Xbeta <- Xdm_user * beta_root
    } else {
      Xbeta <- Xdm_user %*% beta_root
    }
    if (q == 1) {
      Zalpha <- Zdm_user * alpha_root
    } else {
      Zalpha <- Zdm_user %*% alpha_root
    }
    
    exp_Zalpha_plus_AXbeta <- as.vector(exp(Zalpha + A_user * Xbeta))
    exp_Zalpha_minus_AXbeta <- as.vector(exp(Zalpha - A_user * Xbeta))
    exp_negAXbeta <- as.vector(exp(- A_user * Xbeta))
    
    # compute r_i
    r_i <- matrix((Y_user - exp_Zalpha_plus_AXbeta) * avail_user)
    
    # compute D_i
    # updated to allow for differing days within a study without invoking availability (forces vectors to be matrices)
    D_i <- cbind( as.vector(exp_negAXbeta * cA_user) * Xdm_user, exp_Zalpha_minus_AXbeta * Zdm_user )
    # exp_Zalpha_minus_AXbeta * Xdm_user: multiply the vector exp_Abeta_plus_Xalpha to each column of Xdm_user
    
    # compute partial_e_i/partial_theta
    # updated to allow for differing days within a study without invoking availability (forces vectors to be matrices)
    de_i <- cbind( - as.vector(exp_Zalpha_plus_AXbeta * A_user) * Xdm_user * avail_user,
                   - exp_Zalpha_plus_AXbeta * Zdm_user * avail_user )
    
    # compute H_i
    H_i <- de_i %*% Mn_inv %*% t(D_i) / sample_size
    
    Ii_minus_Hi_inv <- solve(diag(length(rows_for_this_user)) - H_i)
    meat <- meat + t(D_i) %*% Ii_minus_Hi_inv %*% r_i %*% t(r_i) %*% t(Ii_minus_Hi_inv) %*% D_i
  }
  meat <- meat / sample_size
  
  # calculate asymptotic variance with small sample adjustment
  # "ssa" stands for small sample adjustment
  
  asymp_varcov_ssa <- Mn_inv %*% meat %*% t(Mn_inv)
  asymp_var_ssa <- diag(asymp_varcov_ssa)
  
  beta_se_ssa <- sqrt(asymp_var_ssa[1:p] / sample_size)
  alpha_se_ssa <- sqrt(asymp_var_ssa[(p+1):(p+q)] / sample_size)
  
  
  ############## part 5 :: p-value with small sample correction ###############
  
  # t test (two-sided -- note the use of significance_level/2)
  
  test_stat <- beta_root / beta_se_ssa
  critical_value <- qt(1 - significance_level/2, df = n - p - q) # two-sided
  p_val <- 2 * pt(abs(test_stat), df = n - p - q, lower.tail = FALSE) # two-sided
  names(test_stat) <- names(p_val) <- Xnames
  test_result_t <- list(test_stat = test_stat,
                        critical_value = critical_value,
                        p_value = p_val)
  
  # F test (two sided, by the nature of F-test)
  
  test_stat <- as.numeric( t(beta_root) %*% solve(asymp_varcov_ssa[1:p, 1:p] / sample_size) %*% beta_root )
  n <- sample_size
  critical_value <- qf((n-q-p) * (1-significance_level) / (p * (n-q-1)), df1 = p, df2 = n-q-p)
  p_val <- pf(test_stat, df1 = p, df2 = n-q-p, lower.tail = FALSE)
  test_result_f <- list(test_stat = test_stat,
                        critical_value = critical_value,
                        p_value = p_val)
  
  ############## part 6 :: generate output values with variable names ###############
  
  beta_hat <- as.vector(beta_root)
  names(beta_hat) <- names(beta_se) <- names(beta_se_ssa) <- Xnames
  alpha_hat <- as.vector(alpha_root)
  names(alpha_hat) <- names(alpha_se) <- names(alpha_se_ssa) <- Znames
  
  output = list(beta_hat = beta_hat, alpha_hat = alpha_hat,
                beta_se = beta_se, alpha_se = alpha_se,
                beta_se_ssa = beta_se_ssa, alpha_se_ssa = alpha_se_ssa,
                test_result_t = test_result_t,
                test_result_f = test_result_f,
                varcov = asymp_varcov / sample_size,
                varcov_ssa = asymp_varcov_ssa / sample_size,
                dims = list(p = p, q = q),
                sample_size = sample_size)
  
  # return(output)
  
  
  ############## part 7 :: generate table values ###############
  
  
  beta_hat = as.numeric(output$beta_hat)
  beta_se = as.numeric(output$beta_se)
  beta_se_ssa = as.numeric(output$beta_se_ssa)
  ttest_val = t(data.frame(output$test_result_t[1]))
  ttest_pval = t(data.frame(output$test_result_t[3]))
  beta_ci_lower = output$beta_hat - 1.96 * output$beta_se
  beta_ci_upper = output$beta_hat + 1.96 * output$beta_se
  beta_ci_lower_ssa = output$beta_hat - 1.96 * output$beta_se_ssa
  beta_ci_upper_ssa = output$beta_hat + 1.96 * output$beta_se_ssa
  risk_ratio = exp(output$beta_hat)
  
  beta_table = data.frame(rbind(beta_hat, beta_se, beta_se_ssa, risk_ratio, ttest_val, ttest_pval, beta_ci_lower, beta_ci_upper, beta_ci_lower_ssa, beta_ci_upper_ssa))
  
  names(beta_table) = Xnames
  
  beta_var_cov = output$varcov[1:p, 1:p]
  beta_var_cov_ssa = output$varcov_ssa[1:p, 1:p]
  studyday = unique(dta[, day_var])
  
  ############## part 8 : generate graphs ###############
  
  intercept = c(1)
  mod_mat = data.frame(crossing(intercept, unique(dta[moderator])))
  
  estimates = data.frame(as.matrix(mod_mat) %*% as.matrix(beta_hat))
  names(estimates) = c("est")
  
  var_mat = as.matrix(mod_mat) %*% beta_var_cov_ssa %*% t(as.matrix(mod_mat))
  
  sd = sqrt(diag(var_mat))
  
  cilow = estimates - 1.96 * sd
  names(cilow) = c("cilow")
  
  cihigh = estimates + 1.96 * sd
  names(cihigh) = c("cihigh")
  
  plotdf = cbind(estimates, cilow, cihigh, mod_mat)
  
  ############## part 9 : return output ###############
  
  beta_table_t = data.frame(t(beta_table))
  beta_table_t = beta_table_t[,c(1,3,9,10,5,6)]
  colnames(beta_table_t) = c("Estimate", "Robust SE", "95% CI Low", "95% CI High", "T", "Pr > |T|")
  
  listreturn = list(beta_table_t, plotdf, beta_var_cov_ssa)
  
  return(listreturn)
  
}



