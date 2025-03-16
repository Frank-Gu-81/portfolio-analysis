library(tidyverse)
library(tseries)
library(ShrinkCovMat)
library(quadprog)
options(digits = 4)


setwd("/Users/frank/Documents/portfolio-analysis")


# Question 2 ----

# excess returns of the 4 mutual funds is stored in `excess_returns`

# get risk free results from the tresury bill data
RF_Data <- read.csv("data/TB3MS.csv")
rff <- RF_Data$TB3MS
rfree <- (1 + rff/100)^(1/12)-1

# calculate the excess returns of the 4 mutual funds
excess_returns <- list()
get_excess_returns <- function(companyList) {
  for (ticker in companyList) {
    companyData <- get.hist.quote(
      instrument = ticker,
      start = "2019-12-01",
      end = "2024-12-31",
      quote = "AdjClose",
      compression = "m"
    )
    adjPrices <- as.vector(companyData)
    returns <- (adjPrices[-1] - adjPrices[-length(adjPrices)]) / adjPrices[-length(adjPrices)]
    excess_returns[[ticker]] <- returns - rfree
  }
  return(as.data.frame(excess_returns))
}

companyList <- c("PRRAX", "FSPHX", "VUIAX", "FSRBX")

excess_returns <- get_excess_returns(companyList)

excess_returns


# Question 3 ----

# excess returns of SP500 is stored in `sp500`

sp500_temp <- get.hist.quote(
  instrument = "^GSPC",
  start = "2019-12-01",
  end = "2024-12-31",
  quote = "AdjClose",
  compression = "m"
)
sp500_adj <- as.vector(sp500_temp)
sp500_returns <- (sp500_adj[-1] - sp500_adj[-length(sp500_adj)]) / sp500_adj[-length(sp500_adj)]
sp500 <- sp500_returns - rfree
sp500

# check for correctness, mean should be 0.009361 and standard deviation should be 0.05246
mean(sp500)
sd(sp500)



# Question 4 ----
sample_mu <- apply(excess_returns, 2, mean)
sample_sd <- apply(excess_returns, 2, sd)


sample_mu_df <- as.data.frame(sample_mu)
sample_mu_df <- sample_mu_df |> 
  rownames_to_column(var = "Fund")  |> 
  rename("mean" = "sample_mu")
sample_mu_df |> 
  knitr::kable()

sample_sd_df <- as.data.frame(sample_sd)
sample_sd_df <- sample_sd_df |> 
  rownames_to_column(var = "Fund")  |> 
  rename("sd" = "sample_sd")
sample_sd_df |> 
  knitr::kable()


# Question 5 ----
sample_cov <- cov(excess_returns)
sample_cov
sample_cov_df <- as.data.frame(sample_cov)
sample_cov_df |> 
  knitr::kable()



# Question 6 ----
cov_shrink <- shrinkcovmat.equal(t(excess_returns))$Sigmahat
cov_shrink

cov_shrink_df <- as.data.frame(cov_shrink)
cov_shrink_df |> 
  knitr::kable()


# Question 7 ----
get_alphas <- function(y) {
  summary(lm(y ~ sp500))$coefficients[1, 1]
}
get_betas <- function(y) {
  summary(lm(y ~ sp500))$coefficients[2, 1]
}
get_sd <- function(y) {
  summary(lm(y ~ sp500))$sigma
}
stks.alpha <- apply(excess_returns, 2, get_alphas)
stks.beta <- apply(excess_returns, 2, get_betas)
stks.s <- apply(excess_returns, 2, get_sd)
stks.Sig <- var(c(sp500)) * (stks.beta%*%t(stks.beta)) + diag(stks.s^2)
stks.Sig
stks.Sig.df <- as.data.frame(stks.Sig)
stks.Sig.df
sample_cov

sd_comp <- data.frame(diag(stks.Sig)^0.5, diag(sample_cov)^0.5) |> 
  rename("Sig Under SIM" = "diag.stks.Sig..0.5", "Sample Sig" = "diag.sample_cov..0.5")
sd_comp



# Question 8 ----
lambda.1 <- 5
ra.1 <- solve.QP(
  Dmat = lambda.1 * sample_cov,
  dvec = sample_mu,
  Amat = matrix(rep(1,4), 4, 1),
  bvec = 1,
  meq = 1
)
w_opt_ra.1 <- ra.1$solution
w_opt_ra_mean.1 <- sum(w_opt_ra.1 * sample_mu)
w_opt_ra_sd.1 <- (w_opt_ra.1 %*% sample_cov %*% w_opt_ra.1)^0.5
w_opt_ra_mean.1
w_opt_ra_sd.1


# Question 9 ----
lambda.2 <- 10
ra.2 <- solve.QP(
  Dmat = lambda.2 * sample_cov,
  dvec = sample_mu,
  Amat = matrix(rep(1,4), 4, 1),
  bvec = 1,
  meq = 1
)
w_opt_ra.2 <- ra.2$solution
w_opt_ra_mean.2 <- sum(w_opt_ra.2 * sample_mu)
w_opt_ra_sd.2 <- (w_opt_ra.2 %*% sample_cov %*% w_opt_ra.2)^0.5
w_opt_ra_mean.2
w_opt_ra_sd.2


# Question 10 ----
lambda.3 <- 20
ra.3 <- solve.QP(
  Dmat = lambda.3 * sample_cov,
  dvec = sample_mu,
  Amat = matrix(rep(1,4), 4, 1),
  bvec = 1,
  meq = 1
)
w_opt_ra.3 <- ra.3$solution
w_opt_ra_mean.3 <- sum(w_opt_ra.3 * sample_mu)
w_opt_ra_sd.3 <- (w_opt_ra.3 %*% sample_cov %*% w_opt_ra.3)^0.5
w_opt_ra_mean.3
w_opt_ra_sd.3


# Question 12 ----
w_tangency <- solve(sample_cov, sample_mu) / sum(solve(sample_cov, sample_mu))
w_tangency

w_tangency_mean <- sum(w_tangency * sample_mu)
w_tangency_sd <- (w_tangency %*% sample_cov %*% w_tangency)^0.5
w_tangency_mean
w_tangency_sd


# Question 13 Example 5.11 ----
wT_nn <- solve.QP(
  Dmat = 2 * sample_cov,
  dvec = rep(0, 4),
  Amat = cbind(sample_mu, diag(4)),
  bvec = c(1, 0, 0, 0, 0),
  meq = 1
)$solution

wT_nn <- wT_nn / sum(wT_nn)

wT_nn_mean <- sum(wT_nn * sample_mu)
wT_nn_sd <- (wT_nn %*% sample_cov %*% wT_nn)^0.5
wT_nn_mean
wT_nn_sd


# Question 14 ----
wT_shrink <- solve(cov_shrink, sample_mu) / sum(solve(cov_shrink, sample_mu))
wT_shrink

wT_shrink_mean <- sum(wT_shrink * sample_mu)
wT_shrink_sd <- (wT_shrink %*% cov_shrink %*% wT_shrink)^0.5
wT_shrink_mean
wT_shrink_sd


# Question 15 ----
# Compute Sharpe Ratios
sharpe_ra_5 <- (w_opt_ra_mean.1) / w_opt_ra_sd.1  # From (8)
sharpe_ra_10 <- (w_opt_ra_mean.2) / w_opt_ra_sd.2  # From (9)
sharpe_ra_20 <- (w_opt_ra_mean.3) / w_opt_ra_sd.3  # From (10)
sharpe_tangency <- (w_tangency_mean) / w_tangency_sd  # From (12)
sharpe_tangency_nn <- (wT_nn_mean) / wT_nn_sd  # From (13)
sharpe_tangency_shrink <- (wT_shrink_mean) / wT_shrink_sd  # From (14)

# Store results in a data frame
sharpe_ratios <- data.frame(
  Portfolio = c("Risk-Averse (λ=5)", "Risk-Averse (λ=10)", "Risk-Averse (λ=20)", 
                "Tangency Portfolio", "Non-Negative Tangency", "Shrinkage Tangency"),
  Sharpe_Ratio = c(sharpe_ra_5, sharpe_ra_10, sharpe_ra_20, 
                   sharpe_tangency, sharpe_tangency_nn, sharpe_tangency_shrink)
)
sharpe_ratios |> 
  knitr::kable()


# Question 16 ----
p_mu <- cbind(excess_returns, sp500) |> 
  apply(2, mean)
p_sd <- cbind(excess_returns, sp500) |> 
  apply(2, sd)
p_cov <- cov(cbind(excess_returns, sp500))

wT_p <- solve(p_cov, p_mu) / sum(solve(p_cov, p_mu))
wT_p_df <- as.data.frame(wT_p) |> 
  rownames_to_column(var = "Asset") |> 
  rename("Weight" = "wT_p")

wT_p_df |> 
  knitr::kable(
    caption = "Weight Allocation of Tangency Portfolio with 4 Funds and SP 500"
  ) |> 
  kableExtra::kable_styling(latex_options = "hold_position")


ra_16 <- solve.QP(
  Dmat=2*p_cov, 
  dvec=rep(0,5), 
  Amat=cbind(p_mu, diag(5)), 
  bvec=c(1,rep(0,5)), 
  meq=1
)
w_16 <- ra_16$solution/sum(ra_16$solution)
weight_16 <- formatC(w_16, format="g", digits=4)
weight_16_df <- as.data.frame(weight_16) |> 
  rownames_to_column(var = "Asset") |> 
  rename("Weight" = "weight_16")
weight_16_df$Asset <- c(colnames(excess_returns), "sp500")
weight_16_df |> 
  knitr::kable(
    caption = "Weight Allocation of Risk-Averse Portfolio with 4 Funds and SP 500"
  ) |> 
  kableExtra::kable_styling(latex_options = "hold_position")
  

sharpe_p <- (p_mu) / (p_sd) 


# Question 17 ----
p_betas <- apply(excess_returns, 2, get_betas)
p_betas_df <- as.data.frame(p_betas)
p_alphas <- apply(excess_returns, 2, get_alphas)
p_alphas_df <- as.data.frame(p_alphas)
p_params <- cbind(p_betas, p_alphas)

p_params |> 
  knitr::kable()


# Question 18 ----
f_pvalue <- function(y) {
  summary(lm(y ~ sp500))$coefficients[1, 4]
}
p_pvalue <- apply(excess_returns, 2, f_pvalue)
p_pvalue

p_pv_fdr <- p.adjust(p_pvalue, method = "fdr")
p_pv_fdr

min(p_pv_fdr)


# Question 19 ----
f_get_residual_se <- function(y) {
  summary(lm(y ~ sp500))$sigma
}
p_non_market_se <- apply(excess_returns, 2, f_get_residual_se)

p_market_se <- p_betas * sd(sp500)

p_risk_components <- cbind(p_non_market_se, p_market_se) |> 
  as.data.frame() |> 
  rename("Non-Market Risk" = "p_non_market_se", "Market Risk" = "p_market_se")
p_risk_components |> 
  knitr::kable()


# Question 20 ----
f_get_rsq <- function(y) {
  summary(lm(y ~ sp500))$r.squared
}
p_rsq <- apply(excess_returns, 2, f_get_rsq)
p_rsq_df <- as.data.frame(p_rsq) |> 
  rename("R-Squared" = "p_rsq")
p_rsq_df |>
  knitr::kable()

# Question 21 ----
p_eqw <- rowMeans(excess_returns)
alpha_p_eqw <- get_alphas(p_eqw)
beta_p_eqw <- get_betas(p_eqw)

mr_p_eqw <- beta_p_eqw * sd(sp500)
nmr_p_eqw <- f_get_residual_se(p_eqw)

ret_sd_p_eqw <- sd(p_eqw)

p_eqw_params <- cbind(beta_p_eqw, alpha_p_eqw, mr_p_eqw, nmr_p_eqw, ret_sd_p_eqw) |> 
  as.data.frame() |> 
  rename("Beta" = "beta_p_eqw", "Alpha" = "alpha_p_eqw", "Market Risk" = "mr_p_eqw", "Non-Market Risk" = "nmr_p_eqw", "Return SD" = "ret_sd_p_eqw")

p_eqw_params |>
  knitr::kable()


# Question 22 ----
# Add constraint that the portfolio is uncorrelated with S&P 500
w_mv_uncorrelated <- solve.QP(
  Dmat = 2 * sample_cov, 
  dvec = rep(0, 4), 
  Amat = cbind(p_betas, rep(1, 4)),  # Constraint matrix 
  bvec = c(0, 1),  # Enforce beta = 0 and sum(weights) = 1
  meq = 2
)$solution

# Compute expected return and standard deviation of the portfolio
uncorrelated_mean <- sum(w_mv_uncorrelated * sample_mu)
uncorrelated_sd <- (w_mv_uncorrelated %*% sample_cov %*% w_mv_uncorrelated)^0.5


# Question 23 ----
sharpe_uncorrelated <- (uncorrelated_mean) / (uncorrelated_sd)
sharpe_sp500 <- mean(sp500) / sd(sp500)
sharpe_23 <- cbind(sharpe_uncorrelated, sharpe_sp500) |> 
  as.data.frame() |> 
  rename(
    "Sharpe Ratio for SP500" = "sharpe_sp500",
    "Sharpe Ratio for Uncorrelated Portfolio" = "V1"
  )
sharpe_23 |>
  knitr::kable()


# Question 24 ----
# Compute portfolio returns using weights from part (22)
p_mv_uncorrelated <- as.matrix(excess_returns) %*% w_mv_uncorrelated  # Compute portfolio returns

alpha_p_mv <- get_alphas(p_mv_uncorrelated)  # Alpha
beta_p_mv <- get_betas(p_mv_uncorrelated)   # Beta
non_market_risk_mv <- f_get_residual_se(p_mv_uncorrelated)  # Residual standard deviation

# Compute Market Risk (Beta * SD of Market)
market_risk_mv <- beta_p_mv * sd(sp500)

# Compute R-squared: proportion of variance explained by the market
r_squared <- (beta_p_mv^2 * var(sp500)) / var(p_mv_uncorrelated)

# Display results in a data frame
market_model_params_24 <- data.frame(
  Alpha = round(alpha_p_mv, 6),
  Beta = round(beta_p_mv, 6),
  Market_Risk = round(market_risk_mv, 6),
  Non_Market_Risk = round(non_market_risk_mv, 6),
  R_Squared = round(r_squared, 6)
)

knitr::kable(market_model_params_24, caption = "Market Model Parameters for Minimum Variance Portfolio (Uncorrelated)")


# Question 25 ----
# Function to calculate Sharpe, Treynor, and Appraisal Ratios
calculate_ratios <- function(fund_returns, market_returns) {
  model <- lm(fund_returns ~ market_returns)
  alpha <- coef(model)[1]  # Alpha
  beta <- coef(model)[2]   # Beta
  total_sd <- sd(fund_returns)  # Total risk (std dev)
  market_sd <- beta * sd(market_returns)  # Market risk
  residual_sd <- summary(model)$sigma  # Non-market risk
  
  # Sharpe Ratio = Excess Returns / Total Risk
  sharpe_ratio <- (mean(fund_returns)) / total_sd
  
  # Treynor Ratio = Excess Returns / Beta
  treynor_ratio <- ifelse(beta != 0, (mean(fund_returns)) / beta, NA)
  
  # Appraisal Ratio = Alpha / Non-Market Risk
  appraisal_ratio <- ifelse(residual_sd != 0, alpha / residual_sd, NA)
  
  return(c(Sharpe_Ratio = sharpe_ratio, Treynor_Ratio = treynor_ratio, Appraisal_Ratio = appraisal_ratio))
}

# Compute ratios for each fund
ratios_list <- lapply(excess_returns, calculate_ratios, sp500)

# Convert results to a data frame
ratios_df <- do.call(rbind, ratios_list)
rownames(ratios_df) <- colnames(excess_returns)  # Fund names
ratios_df <- as.data.frame(ratios_df) |> 
  rename(
    "Sharpe Ratio" = "Sharpe_Ratio",
    "Treynor Ratio" = "Treynor_Ratio.market_returns",
    "Appraisal Ratio" = "Appraisal_Ratio"
  )

knitr::kable(ratios_df, caption = "Sharpe, Treynor, and Appraisal Ratios for Each Fund")







