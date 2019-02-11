library(susieR)
fitted <- list()
posterior <- list()
for (r in 1:ncol(data$Y)) {
  if ('Z' %in% names(data)) {
      data$Y[,r] = residuals(lm(data$Y[,r]~Z, na.action=na.exclude))
  }
  if (prior_var == 'auto') {
      fitted[[r]] <- susie_auto(data$X,data$Y[,r],L_max=100,tol=1e-3)
  } else if (prior_var == 0) {
      fitted[[r]] <- susie(data$X,data$Y[,r],
                               L=maxL,
                               max_iter=maxI,
                               estimate_residual_variance=TRUE,
                               estimate_prior_variance=TRUE,
                               null_weight=null_weight,
                               coverage=0.95,min_abs_corr=0.5,
                               optimV_method = V_method,
                               tol=1e-3)
  } else {
      fitted[[r]] <- susie(data$X,data$Y[,r],
                               L=maxL,
                               max_iter=maxI,
                               estimate_residual_variance=TRUE,
                               scaled_prior_variance=prior_var,
                               null_weight=null_weight,
                               coverage=0.95,min_abs_corr=0.5,
                               tol=1e-3)
  }
  posterior[[r]] <- summary(fitted[[r]])
}