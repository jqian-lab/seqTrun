fn_data <- function(l.shape,l.rate,u.start,u.end,r.shape,r.rate,n){
  temp.dat <- matrix(0,n,3)
  i<-1
  while (i <=n) {
    l <- rgamma(1,l.shape,l.rate)
    u <- runif(1,u.start,u.end)
    r <- rgamma(1,r.shape,r.rate)
    x <- round(l + u,5)
    if(x<=r&(x%in%temp.dat[,2]==FALSE)){
      temp.dat[i,] <- c(l,x,r)
      i<-i+1
    }
  }
  return(as.matrix(temp.dat))
}

l.shape <- 26
l.rate <- 0.5
u.start <- 0
u.end <- 3
n=200
r.shape <- 20
r.rate <- 0.5
dat <- fn_data(l.shape,l.rate,u.start,u.end,r.shape,r.rate,n)
p0=pgamma(58.09,l.shape,l.rate)
a=qgamma(seq(0.05,1,0.05)*p0,l.shape,l.rate)
t.a=qgamma(p0,l.shape,l.rate)
B=5
v=quantile(dat[,3],probs = seq(0,1,0.1))

seqTrunSPMLE(dat,a,v,t.a,B)


