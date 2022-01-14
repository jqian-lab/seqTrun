###functions to apply
fun_test <- function(a,l,x,r,w){
  fit <- fn_est(a,l,x,r,w)
  return(fit[[1]])
}
fun_covar_est <- function(a,var2,l,x,r){
  fit <- fun_covar(a,var2,l,x,r)
  return((fit[[1]]))
}
fun_variance_est <- function(s,var2,l,x,r,w){
  fit <- fun_variance(s,var2,l,x,r,w)
  return((fit[[1]]))
}
###function to estimate CDF of L
fn_est <- function(a,l,x,r,w){
  ## to estimate P(L<u)
  ## L<X<R,
  ## observable: (L,X,R)|X<R
  n<-length(x)
  fit <-  survival::survfit(survival::Surv(x, r, event=rep(1,n)) ~ 1,weights = w,
                            timefix = FALSE)
  s.r.x <- summary(fit,time=x)$surv
  s.r.x <- s.r.x[match(x,sort(x))]
  F.l <- 1/sum(w/s.r.x)*sum((l<=a)*w/s.r.x)
  #  bet <- sum(w/s.r.x)*n
  bet <- 1-n/sum(w/s.r.x) #truncation probability
  return(list(F.l,bet))
}
###function to estimate variance
fun_var <- function(s,l,x,r){
  n<-length(x)
  vec.temp <- cbind(x,l,r)[order(x),]
  x.o <- vec.temp[,1]
  l.o <- vec.temp[,2]
  r.o <- vec.temp[,3]
  fit <-  survival::survfit(survival::Surv(x, r, event=rep(1,n)) ~ 1, timefix = FALSE)
  #fit.o <- survival::survfit(survival::Surv(x.o, r.o, event=rep(1,n)) ~ 1,weights = w)
  fit.x <- survival::survfit(survival::Surv(-r, -x, event=rep(1,n)) ~ 1, timefix = FALSE)
  s.r.x <- summary(fit,time=x)$surv
  s.r.x <- s.r.x[match(x,sort(x))]
  s.r.x.o <- sort(s.r.x,decreasing = TRUE)
  F.x <- summary(fit.x)$surv
  F.x.o <- sort(F.x,decreasing = FALSE)
  F.x.o <- c(F.x.o[-1],1)
  H.x <- -log(s.r.x.o)
  h.x <- H.x-c(0,H.x[-length(H.x)])
  beta <- n*(1/sum(1/s.r.x))
  #if(length(F.x.o)!=length(s.r.x.o)){
  #  break
  #}
  temp <- h.x*beta/(s.r.x.o*F.x.o)
  diff <- x.o-c(0,x.o[-n])
  temp.diff <- temp*diff
  line.1.2 <- 1/sum(1/s.r.x)*sum((l<=s)/s.r.x)-ifelse(l.o<=s,1,0) #length=n
  line.7 <- (1-2*(1/sum(1/s.r.x)*sum((l<=s)/s.r.x)))
  line.9 <- (1/sum(1/s.r.x)*sum((l<=s)/s.r.x))^2
  line.3 = line.4  = line.8 = line.10 = part.1 <- rep(0,n)
  cum_sum <- cumsum(temp.diff)
  line.5 <- matrix(cum_sum, ncol = n, nrow = n)
  t_res <-t(line.5)
  line.5[lower.tri(line.5)] <- t_res[t(upper.tri(line.5))]
  #line.5.6.sum <- sum(line.5)

  for (i in 1:n) {
    line.4[i] <- 1/sum(1/s.r.x)*sum((ifelse(l==l.o[i]&x==x.o[i],1,0))/s.r.x.o[i])
    line.8[i] <- (ifelse(l.o[i]<=s,1,0)/s.r.x.o[i])*1/sum(1/s.r.x)*sum((ifelse(l==l.o[i]&x==x.o[i],1,0))/s.r.x.o[i])
    line.10[i] <- 1/sum(1/s.r.x)*sum((ifelse(l==l[i]&x==x[i],1,0))/s.r.x[i])# i-th component of delta_est
  }

  for (j in 1:n) {
    line.3[j] <- 1/sum(1/s.r.x)*sum((ifelse(l==l.o[j]&x==x.o[j],1,0))/s.r.x.o[j])
    #    part.1[j] <- sum(line.1.2[j]*line.1.2*line.3[j]*line.4*line.5.6.sum)
    part.1[j] <- sum(line.1.2[j]*line.1.2*line.3[j]*line.4*line.5[,j]) #key part
  }

  delta <- sum(1/(s.r.x)*line.10)
  var.est <- (sum(part.1)+beta*line.7*sum(line.8)+beta*delta*line.9)/n
  var.est <- ifelse(var.est>0,var.est,0)
  return(list(var.est))
}
###function to estimate covariance
fun_covar <- function(s,t,l,x,r){
  n<-length(x)
  vec.temp <- cbind(x,l,r)[order(x),]
  x.o <- vec.temp[,1]
  l.o <- vec.temp[,2]
  r.o <- vec.temp[,3]
  fit <-  survival::survfit(survival::Surv(x, r, event=rep(1,n)) ~ 1, timefix = FALSE)
  #fit.o <- survival::survfit(survival::Surv(x.o, r.o, event=rep(1,n)) ~ 1,weights = w)
  fit.x <- survival::survfit(survival::Surv(-r, -x, event=rep(1,n)) ~ 1, timefix = FALSE)
  s.r.x <- summary(fit,time=x)$surv
  s.r.x <- s.r.x[match(x,sort(x))]
  s.r.x.o <- sort(s.r.x,decreasing = TRUE)
  F.x <- summary(fit.x)$surv
  F.x.o <- sort(F.x,decreasing = FALSE)
  F.x.o <- c(F.x.o[-1],1)
  H.x <- -log(s.r.x.o)
  h.x <- H.x-c(0,H.x[-length(H.x)])
  beta <- n*(1/sum(1/s.r.x))
  #if(length(F.x.o)!=length(s.r.x.o)){
  #  break
  #}
  temp <- h.x*beta/(s.r.x.o*F.x.o)
  diff <- x.o-c(0,x.o[-n])
  temp.diff <- temp*diff
  #pos <- which(sort(u,decreasing = FALSE) <= y[i])
  line.1 <- 1/sum(1/s.r.x)*sum((l<=t)/s.r.x)-ifelse(l.o<=t,1,0)
  line.2 <- 1/sum(1/s.r.x)*sum((l<=s)/s.r.x)-ifelse(l.o<=s,1,0)
  line.6 <- (1-(1/sum(1/s.r.x)*sum((l<=t)/s.r.x)))
  line.9 <- (1/sum(1/s.r.x)*sum((l<=s)/s.r.x)) #revised after fixing a typo in p.10 of the supplementary
  line.8 <- ((1/sum(1/s.r.x))^2*sum((l<=s)/s.r.x)*sum((l<=t)/s.r.x))
  line.3 = line.4  = line.7 = line.10 = line.11 = part.1 <- rep(0,n)

  cum_sum <- cumsum(temp.diff)
  line.5 <- matrix(cum_sum, ncol = n, nrow = n)
  t_res <-t(line.5)
  line.5[lower.tri(line.5)] <- t_res[t(upper.tri(line.5))]
  #  line.5.sum <- sum(line.5)

  for (i in 1:n) {
    line.4[i] <- 1/sum(1/s.r.x)*sum((ifelse(l==l.o[i]&x==x.o[i],1,0))/s.r.x.o[i])
    line.7[i] <- (ifelse(l.o[i]<=s,1,0)/s.r.x.o[i])*1/sum(1/s.r.x)*sum((ifelse(l==l.o[i]&x==x.o[i],1,0))/s.r.x.o[i])
    line.10[i] <- (ifelse(l.o[i]<=t,1,0)/s.r.x.o[i])*1/sum(1/s.r.x)*sum((ifelse(l==l.o[i]&x==x.o[i],1,0))/s.r.x.o[i])
    line.11[i] <- 1/sum(1/s.r.x)*sum((ifelse(l==l[i]&x==x[i],1,0))/s.r.x[i])
  }

  for (j in 1:n) {
    line.3[j] <- 1/sum(1/s.r.x)*sum((ifelse(l==l.o[j]&x==x.o[j],1,0))/s.r.x.o[j])
    part.1[j] <- sum(line.1[j]*line.2*line.3[j]*line.4*line.5[,j])
  }
  delta <- sum(1/(s.r.x)*line.11)
  covar.est <- (sum(part.1)+beta*line.6*sum(line.7)+beta*delta*line.8-beta*line.9*sum(line.10))/n
  #covar.est <- ifelse(covar.est>0,covar.est,0)
  return(list(covar.est))
}
###function to estimate analysitical variance
fun_variance <- function(s,t,l,x,r,w){
  #s <- quantile(a,0.25)
  #t <- quantile(a,0.5)
  mu_s <- fn_est(s,l,x,r,w)[[1]]
  mu_t <- fn_est(t,l,x,r,w)[[1]]
  var_s <- fun_var(s,l,x,r)[[1]]
  var_t <- fun_var(t,l,x,r)[[1]]
  cov_st <- fun_covar(s,t,l,x,r)[[1]]
  res <- (mu_s^2/mu_t^2)*(var_s/mu_s^2-2*cov_st/(mu_s*mu_t)+var_t/mu_t^2)
  return(list(res))
}

#' Nonparametric Maximum Likelihood Estimation Under Sequential Truncation
#'
#' \code{seqTrunNPMLE} provides nonparametric maximum likelihood estimate for the distribution
#'      of \eqn{L} under the sequential truncation model \eqn{(L,X)} is independent of \eqn{R}
#'      given \eqn{L<X<R}, and \eqn{P(L<X)=1}. Variance estimation is also evaluated by a nonparametric
#'      bootstrap procedure and an analytical method based on asymptotic properties of the estimator.
#'
#' @param data dataframe containing \eqn{(L,X,R)} on columns 1,2,3 respectively.
#' @param a the point on which the cumulative distribution function (CDF) of \eqn{L} is evaluated.
#' @param t.a the upper limit of \eqn{L} to be evaluated.
#' @param B the number of boostrap replications.
#' @param V the indicator of whether analytical variance is computed.
#'
#' @return The output contains the following components:
#' \describe{
#'   \item{Evaluation point}{the points on which we estimate the CDF of \eqn{L}.}
#'   \item{Point estimate}{the estimate of the CDF of \eqn{L}.}
#'   \item{Bootstrap SE}{the standard error estimated by a bayesian bootstrap procedure.}
#'   \item{Analytical SE}{the standard error estimated by the proposed method.}
#'   \item{Truncation probability}{the probability of \eqn{X>R}.}
#' }
#' @references Betensky R. A., Qian, J. and Hou, J. (2022), Nonparametric and semiparametric
#'     estimation with sequentially truncated survival data. (Techincal Report)
#'
#' @importFrom stats rbeta sd
#' @export
#' @example inst/examples/ex_seqTrunNPMLE.R
###comprehesive estimation function
#data: data frame containing data points and interval bounds.
#a: The point on which the CDF is evaluated.
#t.a: The conditional upper bound of the evaluated points.
#B: The number of times of bootstrap procedure.
#V: The indicator of whether analytical variance is obtained.
seqTrunNPMLE <- function(data,a,t.a,B,V){
  ##vector to be used for estimation
  temp.dat <- data
  n <- dim(temp.dat)[1]
  l <- temp.dat[,1]
  x <- temp.dat[,2]
  r <- temp.dat[,3]
  w <- rep(1,length(x))
  npmle.temp <- matrix(0,B,length(a))
  beta <- fn_est(1,l,x,r,w)[[2]]
  ##point estimation
  point_est <- sapply(sort(a),fun_test,l=l,x=x,r=r,w=w)

  ##analytical variance
  #The analytical variance will be applied if V is TRUE.
  if(V==TRUE){
    F.l.cov.temp <- sapply(sort(a), fun_covar_est, var2=t.a,l=l,x=x,r=r)
    F.l.variance.temp <- sapply(sort(a), fun_variance_est, var2=t.a,l=l,x=x,r=r,w=w)
    F.l.variance.temp <- ifelse(F.l.variance.temp>0,F.l.variance.temp,0)
    F.l.se.temp.new<-sqrt(F.l.variance.temp)
  }

  ##bootstrap
  #The bootstrap procedure will be applied B times if B is larger than 1.
  if(B>1){
    for (j in 1:B) {
      set.seed(j)
      w <- 4*rbeta(length(x), 0.5, 1.5) #using beta dist.
      w <- w/sum(w)
      F.npmle.est.bp<-sapply(sort(a),fun_test,l=l,x=x,r=r,w=w)
      F.npmle.est.bp <- F.npmle.est.bp/F.npmle.est.bp[length(a)]
      npmle.temp[j,] <- F.npmle.est.bp
    }
    F.np.se.temp <- apply(npmle.temp, 2, sd)
  }
  out <- list(a,point_est,F.np.se.temp,F.l.se.temp.new,beta)
  names(out) <- c("Evaluation point","Point estimate","Bootstrap SE","Analytical SE","Truncation probability")
  return(out)
}
##
#Evaluation point: The points on which we estimate the CDF of L.
#Point estimate: The estimate of the CDF of L.
#Bootstrap SE: The standard error estimated by a bayesian bootstrap procedure.
#Analytical SE: The standard error estimated by the proposed method.
#Truncation probability: The estimated truncation probability P(X>R)
