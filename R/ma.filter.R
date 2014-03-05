ma.filter <-
function(xt, data = NULL, seasonal = FALSE, d = NULL,
                      na.action = na.fail, plot = FALSE)
{ 
  if (!is.null(data))
    xt <- na.action(as.ts(data$xt))
  else 
  xt <- na.action(as.ts(xt))
  if (is.matrix(xt))
    stop("only univariate series are allowed")
  n <- as.integer(length(xt))
  if (is.na(n)) 
    stop("invalid length(x)")
  
  trendest <- function(xt) 
  {
    xt <- as.vector(xt)
    n <- as.integer(length(xt))
    mhat <- NULL
    q.n <- ifelse(qn(xt) < n, as.integer(min(qn(xt), n - qn(xt))), 
                  floor(n^(4/5)/2))
    for (j in 1:n)  
    {  
      if (j >= q.n + 1 && j <= n - q.n)  
      {
        mhat[j] <- sum(xt[(j - q.n):(q.n + j)])/(2*q.n + 1)
      }  
      
      else if (j < q.n + 1)   
      {
        mhat[j] <- ((4*q.n^2 - 4*q.n*j + 6*q.n + 4*j^2-6*j + 2)*sum(xt[1:(q.n + j)]) - 
                    6*(q.n - j + 1)*sum(c((-j + 1):q.n)*xt[1:(q.n + j)]))/
                   ((q.n + j)*((q.n + j)^2-1))
      } 
      else if (j > n - q.n)
      {
        mhat[j] <- ((4*(n - j)^2 + 4*q.n*(q.n + j - n) + 2*(n + q.n - j))*
                    sum(xt[(j - q.n):n]) + 6*(q.n - n + j)*sum(c(-q.n:(n - j))*
                    xt[(j - q.n):n]))/((n + q.n -j + 2)*(n + q.n - j + 1)*(n + q.n - j))
      }
    }
    return(mhat)
  }
  
  seasonest <- function(st)
  {
    if (is.null(d))
      stop("please input seasonal period d")
    if (d %% 1 != 0)
      stop("seasonal period d must be an integer")
    if (d < 2L && n%%d != 0)
      stop("seasonal period must be at least 2 and proportional to length(xt)")
    
    s <- matrix(0,n/d,d)
    for (i in 1:d)
    {
      s[,i] <- trendest(st[seq(i,n,by = d)])
    }
    list(data = xt, trend = ts(trendest(xt), start = date.start,
         frequency = d), season = ts(as.vector(t(s)),start = date.start, 
         frequency = d), remainder = ts(xt - trendest(xt) - as.vector(t(s)),
         start = date.start, frequency = d),
         SSE = sum((xt - trendest(xt) - as.vector(t(s)))^2))
  }
  
  date.start <- as.vector(time(xt))[1]

  if (seasonal == FALSE)
    Estimate <- list(data = xt, trend = ts(trendest(xt), start = date.start), 
                     remainder = ts(xt - trendest(xt), start = date.start), 
                     SSE = sum((xt - trendest(xt))^2))
  else  
    Estimate <- seasonest(xt - trendest(xt))
  
  if (plot == TRUE)
  {
    m <- length(names(Estimate)) - 1
    op <- par(mfrow = c(m,1), mar = c(2,6,1,1) + 0.1)
    for (i in 1:m)
      plot(as.vector(time(Estimate[[1]])),Estimate[[i]], ylab = names(Estimate)[i],
           xlab = "", type = "l", col = i)
    par(op)
  }

  return(Estimate)
}
