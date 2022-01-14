###functions to apply
fun_est_semi <- function(a,l,x,r,w,v){
  fit <- fn_est_semi(a,l,x,r,w,v)
  return(fit[[1]])
}
fun_covar_est <- function(a,var2,l,x,r,v){
  fit <- fun_covar(a,var2,l,x,r,v)
  return((fit[[1]]))
}
fun_variance_est <- function(s,var2,l,x,r,w,v){
  fit <- fun_variance(s,var2,l,x,r,w,v)
  return((fit[[1]]))
}
###function to estimate theta
fun_theta <- function(x,r,v,u,w){
  ##the length of theta should be smaller than the length of v
  J <- length(v)-1
  theta <- rep(0,length(v)-1)
  ##estimate the theta based on the derivation
  j.sum <- NULL
  for (j in 1:J) {
    if(j!=J){
      j.sum <- rbind(j.sum,(ifelse(v[j+1]<=r&r<v[j+2],1,0)-ifelse(v[j+1]<=x&x<v[j+2],1,0)))
    }
    else{
      j.sum <- rbind(j.sum,(ifelse(v[j+1]<=r,1,0)-ifelse(v[j+1]<=x,1,0)))
    }
  }
  for (j in 1:J) {
    if(j!=J){
      theta[j] <- sum(w*ifelse(v[j]<=r&r<v[j+1],1,0))/(sum(w*ifelse(v[j]<=r&r<v[j+1],1,0)*(r-v[j]))-
                                                         sum(w*ifelse(v[j]<=x&x<v[j+1],1,0)*(x-v[j]))+
                                                         sum(w*apply(j.sum[j:J,], 2, sum))*(v[j+1]-v[j]))
    }
    else{
      theta[j] <- sum(w*ifelse(v[j]<=r&r<v[j+1],1,0))/(sum(w*ifelse(v[j]<=r&r<v[j+1],1,0)*(r-v[j]))-
                                                         sum(w*ifelse(v[j]<=x&x<v[j+1],1,0)*(x-v[j]))+
                                                         sum(w*j.sum[j,])*(v[j+1]-v[j]))
    }
  }
  ##calculate the cumulative hazard function based on the value of r
  h.theta <- matrix(0,length(u))
  for (i in 1:length(u)) {
    #i=48
    if(u[i]<min(v)){
      h.theta[i] <- 0
    }
    else{
      pos.temp <- max(which(v<=u[i]))
      if(pos.temp==1){
        h.theta[i] <- theta[pos.temp]*(u[i]-v[pos.temp])
      }
      else{
        h.theta[i] <- sum((v[1:pos.temp]-c(0,v[1:(pos.temp-1)]))*c(0,theta[1:(pos.temp-1)]))+
          theta[pos.temp]*(u[i]-v[pos.temp])
      }
    }
    #print(i)
  }
  ##calculate the survival function based on its property
  s.theta <- sort(exp(-h.theta),decreasing = TRUE)
  return(list(theta,h.theta,s.theta))
}
###function to estimate semiparametric CDF of L
fn_est_semi <- function(a,l,x,r,w,v){
  ## to estimate P(L<u)
  ## L<X<R,
  ## observable: (L,X,R)|X<R
  n<-length(x)
  theta <- rep(0,length(v)-1)
  fit <-  fun_theta(x,r,v,x,w)
  s.r.x <- fit[[3]]
  s.r.x <- s.r.x[match(x,sort(x))]
  #s.r.x[which(s.r.x==0)] <- min(s.r.x[which(s.r.x!=0)])
  F.l <- 1/sum(w/s.r.x)*sum((l<=a)*w/s.r.x)
  bet <- sum(w/s.r.x)*n
  bet <- 1-n/sum(w/s.r.x) #truncation probability
  return(list(F.l,bet))
}
###function to estimate variance
fun_var <- function(s,l,x,r,v){
  n<-length(x)
  #r <- rep(100000,n)
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
  #options(max.print = 10000)
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
fun_covar <- function(s,t,l,x,r,v){
  n<-length(x)
  vec.temp <- cbind(x,l,r)[order(x),]
  x.o <- vec.temp[,1]
  l.o <- vec.temp[,2]
  r.o <- vec.temp[,3]
  fit <-  survival::survfit(survival::Surv(x, r, event=rep(1,n)) ~ 1, timefix = FALSE)
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
fun_variance <- function(s,t,l,x,r,w,v){
  mu_s <- fn_est(s,l,x,r,w)[[1]]
  mu_t <- fn_est(t,l,x,r,w)[[1]]
  var_s <- fun_var(s,l,x,r,v)[[1]]
  var_t <- fun_var(t,l,x,r,v)[[1]]
  cov_st <- fun_covar(s,t,l,x,r,v)[[1]]
  res <- (mu_s^2/mu_t^2)*(var_s/mu_s^2-2*cov_st/(mu_s*mu_t)+var_t/mu_t^2)
  return(list(res))
}

#' Semiparametric Maximum Likelihood Estimation Under Sequential Truncation
#'
#' \code{seqTrunSPMLE} provides semiparametric maximum likelihood estimate for the distribution
#'      of \eqn{L} under the sequential truncation model \eqn{(L,X)} is independent of \eqn{R}
#'      given \eqn{L<X<R}, and \eqn{P(L<X)=1}. Variance estimation is also evaluated by a nonparametric
#'      bootstrap procedure.
#'
#' @param data dataframe containing \eqn{(L,X,R)} on columns 1,2,3 respectively.
#' @param a the point on which the cumulative distribution function (CDF) of \eqn{L} is evaluated.
#' @param v the intervals of the piecewise exponential distribution.
#' @param t.a the upper limit of \eqn{L} to be evaluated.
#' @param B the number of boostrap replications.
#'
#' @return The output contains the following components:
#' \describe{
#'   \item{Evaluation point}{the points on which we estimate the CDF of \eqn{L}.}
#'   \item{Point estimate}{the estimate of the CDF of \eqn{L}.}
#'   \item{Bootstrap SE}{the standard error estimated by a bayesian bootstrap procedure.}
#'   \item{Truncation probability}{the probability of \eqn{X>R}.}
#' }
#' @references Betensky R. A., Qian, J. and Hou, J. (2022), Nonparametric and semiparametric
#'     estimation with sequentially truncated survival data. (Techincal Report)
#'
#' @importFrom stats rbeta sd
#' @export
#' @example inst/examples/ex_seqTrunSPMLE.R
###comprehesive estimation function
#data: data frame containing data points and interval bounds.
#a: The point on which the CDF is evaluated.
#v: The intervals of the piecewise exponential distribution.
#t.a: The conditional upper bound of the evaluated points.
#B: The number of times of bootstrap procedure.
seqTrunSPMLE <- function(data,a,v,t.a,B){
  ##vector to be used for estimation
  temp.dat <- data
  n <- dim(temp.dat)[1]
  l <- temp.dat[,1]
  x <- temp.dat[,2]
  r <- temp.dat[,3]
  w <- rep(1,length(x))
  spmle.temp <- matrix(0,B,length(a))
  beta <- fn_est_semi(1,l,x,r,w,v)[[2]]
  ##point estimation
  point_est <- sapply(sort(a),fun_est_semi,l=l,x=x,r=r,w=w,v=v)

  ##bootstrap
  #The bootstrap procedure will be applied B times if B is larger than 1.
  if(B>1){
    for (j in 1:B) {
      set.seed(j)
      w <- 4*rbeta(length(x), 0.5, 1.5) #using beta dist.
      w <- w/sum(w)
      F.spmle.est.bp<-sapply(sort(a),fun_est_semi,l=l,x=x,r=r,w=w,v=v)
      F.spmle.est.bp <- F.spmle.est.bp/F.spmle.est.bp[length(a)]
      spmle.temp[j,] <- F.spmle.est.bp
    }
    F.sp.se.temp <- apply(spmle.temp, 2, sd)
  }
  out <- list(a,point_est,F.sp.se.temp,beta)
  names(out) <- c("Evaluation point","Point estimate","Bootstrap SE","Truncation probability")
  return(out)
}
##
#Evaluation point: The points on which we estimate the CDF of L.
#Point estimate: The estimate of the CDF of L.
#Bootstrap SE: The standard error estimated by a bayesian bootstrap procedure.
#Truncation probability: The estimated truncation probability P(X>R)

