#' @title Compare SuSiE fits to truth
#' @param sets a list of susie CS info from susie fit
#' @param pip probability for p variables
#' @param true_coef true regression coefficients
#' @return total the number of total CS
#' @return valid the number of CS that captures a true signal
#' @return size an array of size of CS
#' @return purity an array of purity of CS
#' @return top the number of CS whose highest PIP is the true causal
susie_scores = function(sets, pip, true_coef) {
  if (is.null(dim(true_coef))) beta_idx = which(true_coef!=0) 
  else beta_idx = which(apply(true_coef, 1, sum) != 0)
  cs = sets$cs
  if (is.null(cs)) {
    size = 0
    total = 0
    purity = 0
  } else {
    size = sapply(cs,length)
    purity = as.vector(sets$purity[,1])
    total = length(cs)
  }
  valid = 0
  top_hit = 0
  if (total > 0) {
    for (i in 1:total){
      if (any(cs[[i]]%in%beta_idx)) valid=valid+1
      set.idx = cs[[i]]
      highest.idx = which.max(pip[set.idx])
      if (set.idx[highest.idx]%in%beta_idx) top_hit=top_hit+1
    }
  }
  return(list(total=total, valid=valid, size=size, purity=purity, top=top_hit))
}

susie_scores_multiple = function(res, truth) {
  total = valid = size = purity = top = 0
  objective = vector()
  for (r in 1:length(res)) {
    out = susie_scores(res[[r]]$sets, res[[r]]$pip, truth[,r])
    total = total + out$total
    valid = valid + out$valid
    size = size + out$size
    purity = purity + out$purity
    top = top + out$top
    objective[r] = susieR::susie_get_objective(res[[r]])
  }
  return(list(total=total, valid=valid, size=size, purity=purity, top=top, objective=objective))
}