get_center <- function(k,n) {
  ## For given number k, get the range k surrounding n/2
  ## but have to make sure it does not go over the bounds
  if (is.null(k)) {
      return(1:n)
  }
  start = floor(n/2 - k/2)
  end = ceiling(n/2 + k/2)
  if (start<1) start = 1
  if (end>n) end = n
  return(start:end)
}
