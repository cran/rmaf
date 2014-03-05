qn <-
function(x, data = NULL, na.action = na.fail)  
{ 
  if (!is.null(data))
    x <- as.vector(na.action(data$x))
  else 
    x <- as.vector(na.action(x))
  if (is.matrix(x))
    stop("only univariate series are allowed")
  n <- as.integer(length(x))
  if (is.na(n)) 
    stop("invalid length(x)")
  xn <- (1:n)/n 
  lm.fit <- lm(x ~ 1 + xn + I(xn^2) + I(xn^3))
  denominator <- 4*(lm.fit[[1]][3])^2 + 12*lm.fit[[1]][3]*
    lm.fit[[1]][4] + 12*(lm.fit[[1]][4])^2
  q.n <- as.integer(floor((n)^(4/5)*(9/2)^(1/5)*
                            ((var(resid(lm.fit))/denominator)^(1/5))))
  qn <- ifelse(q.n < n, as.integer(min(q.n, n - q.n)), floor(n^(4/5)/2))
}
