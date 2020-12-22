##################################################################
# R code for Armillotta, Luati and Lupparelli 
##################################################################

##################################################################
# Estimation and plots from Covid-19 example
##################################################################

# function for loglikelihood and gradient of log-linear Poisson AR

logl.p <- function(delta){   
  
 g <- 0
  
 beta <- delta[1]   # parameters
 phi1 <- delta[2]
 theta1 <-  delta[3]
  
 eta <- vector()
  
 eta[1] <- 0      # starting value
  
 #log-lik.
 
 for (i in (r+1):n){
   
  eta[i] <- beta + phi1*log(y[i-1]+1) + theta1*eta[i-1]  
  g <- g + y[i]*eta[i]-exp(eta[i])-lfactorial(y[i])       
    
 }
  
 g <- -g/n
 return(g)
  
}

gradd.p <- function(delta){  
  
 gr <- rep(0, k)
 out <- matrix(0, k, k)
  
 beta <- delta[1]
 phi1 <- delta[2]
 theta1 <-  delta[3]
  
 eta <- vector()
  
 eta[1] <- 0
  
 for (i in (r+1):n){
  
  eta[i] <- beta + phi1*log(y[i-1]+1) + theta1*eta[i-1]
    
  # gradient
    
  gr[1] <- y[i]-exp(eta[i])
  gr[2] <- (y[i]-exp(eta[i]))*log(1+y[i-1])
  gr[3] <- (y[i]-exp(eta[i]))*eta[i-1]
    
  # outer product
    
  out <- out + gr%*%t(gr)
    
    
 }
  
 #gr <- -gr/n  # print gradient 
 out <- out/n   # print outer product
  
 return(out)
  
}


#-----------------------------------------------------------------

# function for QMLE estimation

ml <- function(k, y, n, d0, logl, gradd, mod, met){
  
 r <- max(p,q)
  
 # select different optimization algorithm
  
 if(met=="nm")

  num_mle <- optim(par=d0, fn=logl, method="Nelder-Mead",
                   hessian=TRUE,
                   control=list(trace=1, maxit=10000))
  
 if(met=="cg")
   
  num_mle <- optim(par=d0, fn=logl, method="CG",
                   hessian=TRUE,
                   control=list(trace=1, maxit=10000))
  
 if(met=="l-bfgs-b")
   
  num_mle <- optim(par=d0, fn=logl, method="L-BFGS-B",
                   lower=c(0, 0, 0, 0),
                   upper=c(Inf, Inf, Inf, Inf),
                   hessian=TRUE,
                   control=list(trace=1, maxit=100000))
  
 if(met=="bfgs")
   
  num_mle <- optim(par=d0, fn=logl, method="BFGS",
                   hessian=TRUE,
                   control=list(trace=1, maxit=10000))
  
 if(met=="constr"){
    
  ui <- rbind(c(0,-1,-1,0),c(0,-1,0,-1),c(0,1,-1,0))
  ci <- c(-1,-1,-1)
  num_mle <- constrOptim(theta=d0, f=logl,
                         ui=ui, ci=ci,
                         method="Nelder-Mead",
                         control=list(trace=1, maxit=10000))
  
  num_mle$hessian <- optimHess(num_mle$par, logl)
    
 }
 
 if(met=="constr2"){
    
  ui <- c(0, -1, -1)
  ci <- -0.999
  num_mle <- constrOptim(theta=d0, f=logl, 
                         ui=ui, ci=ci, method="Nelder-Mead",
                        control=list(trace=1, maxit=10000))
  
  num_mle$hessian <- hessian(logl, x=num_mle$par)
    
 }
  
 # store result of estimation
 # value of maximum (quasi) log-likelihood
 # compute IC
 
 num_par <- num_mle$par        
 num_l <- -num_mle$value       
 num_aic <- 2*k - 2*num_l      
 num_bic <- log(n)*k - 2*num_l
 
 #compute inverse Hessian (matrix J^(-1)) 
 #compute outer product (Matrix I)
 # sandiwhich estimator i.e. matrix J^(-1)IJ^(-1)
 # estimated variance
 # estimated SE
 # t test
  
 num_J <- solve(num_mle$hessian)                               
  
 num_I <- gradd(num_par)          
  
 num_var <- num_J%*%num_I%*%num_J/n 
 num_V <- diag(num_var)             
 num_V <- t(num_V)
 num_SE <- sqrt(num_V)              
 num_t <- num_par/num_SE
 num_p <- pt(abs(num_t), n-k, lower.tail=FALSE) 
  
 # compute fitted values
  
 x <- vector()
 x[1] <- 0
  
 if(mod=="log.p" || mod=="log.nb")
    
  for (i in (r+1):n){
      
   x[i] <- num_par[1] + num_par[2]*log(y[i-1]+1) +
           num_par[3]*x[i-1]
      
  }
  
 if(mod=="garma.p" || mod=="garma.nb"){
    
  y_ <- vector()
  gy <- vector()
  c <- 0.1
    
  y_[1] <- max(y[1], c)
  gy[1] <- log(y_[1])
    
  for (i in (r+1):n){
      
   y_[i] <- max(y[i], c)
   gy[i] <- log(y_[i])
   x[i] <- num_par[1] + num_par[2]*gy[i-1] +
           num_par[3]*(gy[i-1] - x[i-1])
      
  }
    
 }
  
 if(mod=="glarma.p"){
    
  lambda <- vector()
  lambda[1] <- 1
    
  for (i in (r+1):n){
      
   x[i] <- num_par[1] + num_par[2]*x[i-1] +
           num_par[3]*
           (y[i-1] - lambda[i-1])/sqrt(lambda[i-1])
      
   lambda[i] <- exp(x[i])
      
  }
    
 }
  
 if(mod=="glarma.nb"){
    
  lambda <- vector()
  var <- vector()
    
  lambda[1] <- 1
  var[1] <- 1+1/v
    
  for (i in (r+1):n){
      
   x[i] <- num_par[1] + num_par[2]*x[i-1] +
           num_par[3]*
           (y[i-1] - lambda[i-1])/sqrt(var[i-1])
      
   lambda[i] <- exp(x[i])
   var[i] <- lambda[i]*(1+lambda[i]/v)
      
  }
    
 }
  
  
 # compute residuals and print results
  
  
 mu <- exp(x)
  
 if(mod=="log.p" || mod=="garma.p" ||
    mod=="glarma.p" ||
    mod=="gen.p" ||
    mod=="gen1.p") 
   
  e <- (y-mu)/sqrt(mu)
 
 else e <- (y-mu)/sqrt(mu*(1+mu/v))
  
  
 return(list(est=num_par, std=num_SE, t.test=num_t,
             p.value=num_p, mean=mu, pred=x,
             res=e, max=num_l,
             aic=num_aic, bic=num_bic))
  
}


#---------------------------------------------------------------

# Replication of COVID-19 application

library(numDeriv)

jhu_url <- 
  paste("https://raw.githubusercontent.com/CSSEGISandData/", 
                 "COVID-19/master/csse_covid_19_data/",
                 "csse_covid_19_time_series/", 
                 "time_series_covid19_deaths_global.csv",
                 sep = "")   

# load data

jhu <- read.csv(jhu_url)   

# select time series for Italy

hub <- jhu[152,]           
hub <- hub[,-c(1:4)]        
hub <- as.vector(as.numeric(hub))

# remove the days before the pandemic outbreak

hub <- hub[-c(1:30)]        

timeh<-seq(1:length(hub))


# compute daily noncumulant deaths

hub1 <- hub[-length(hub)]
hub1 <- c(0, hub1)
co <- hub-hub1

# fix wrong data

co[177] <- 4       
co[125] <- 31

timeh<-seq(1:length(co))

n<-length(co) # length of the series

y <- co
time <- timeh
              

#select lag one models

p<-1
q<-1
r <- max(p,q)
k <- p+q+1 

# set initial values of the parameters

d0 <- rep(0.01, k)  


# estimate Poisson models

log.pois11 <- ml(k, y, n, d0, logl.p,
                 gradd.p, mod="log.p",
                 me="constr2")

garma.pois11 <- ml(k, y, n, d0, garmal.p,
                   garmagr.p, mod="garma.p",
                   me="constr2")

glarma.pois11 <- ml(k, y, n, d0, glarmal.p,
                    glarmagr.p, mod="glarma.p",
                    me="constr2")


# calibrate dispersion parameter

log.mean <- vector()
garma.mean <- vector()
glarma.mean <- vector()
gen.mean <- vector()

log.v1 <- 0
garma.v1 <- 0
glarma.v1 <- 0
gen.v1 <- 0

for(j in 1:n){
  
 log.mean[j] <- log.pois11$mean[j]
 garma.mean[j] <- garma.pois11$mean[j]
 glarma.mean[j] <- glarma.pois11$mean[j]
  
 log.v1 <- log.v1 + 
           ((y[j]- log.mean[j])^2 -
              log.mean[j])/( log.mean[j]^2)
 
 garma.v1 <- garma.v1 + 
           ((y[j]- garma.mean[j])^2 -
              garma.mean[j])/( garma.mean[j]^2)
 
 glarma.v1 <- glarma.v1 +
           ((y[j]- glarma.mean[j])^2-
              glarma.mean[j])/( glarma.mean[j]^2)
  
}

log.v1 <- log.v1/n
log.v1 <- 1/log.v1
garma.v1 <- garma.v1/n
garma.v1 <- 1/garma.v1
glarma.v1 <- glarma.v1/n
glarma.v1 <- 1/glarma.v1


# estimate NB models

v <- log.v1
log.nb11 <- ml(k, y, n, d0, logl.nb,
               gradd.nb, mod="log.nb",
               me="constr2")

v <- garma.v1
garma.nb11 <- ml(k, y, n, d0, garmal.nb,
                 garmagr.nb, mod="garma.nb",
                 me="constr2")

v <- glarma.v1
glarma.nb11 <- ml(k, y, n, d0, glarmal.nb,
                  glarmagr.nb, mod="glarma.nb",
                  me="constr2")


# Probabilistic calibration via PIT

log.pois.Px <- vector()
log.pois.Px1 <- vector()
log.nb.Px <- vector()
log.nb.Px1 <- vector()
log.pois.px <- vector()
log.nb.px <- vector()

garma.pois.Px <- vector()
garma.pois.Px1 <- vector()
garma.nb.Px <- vector()
garma.nb.Px1 <- vector()
garma.pois.px <- vector()
garma.nb.px <- vector()

glarma.pois.Px <- vector()
glarma.pois.Px1 <- vector()
glarma.nb.Px <- vector()
glarma.nb.Px1 <- vector()
glarma.pois.px <- vector()
glarma.nb.px <- vector()


for(i in 1:n){
  
 log.pois.Px[i] <- ppois(y[i], log.pois11$mean[i])
 
 log.pois.Px1[i] <- ppois(y[i]-1, log.pois11$mean[i])
 
 log.nb.Px[i] <- pnbinom(y[i], size=log.v,
                         prob=log.v/
                           (log.v+
                            log.nb11$mean[i]))
 
 log.nb.Px1[i] <- pnbinom(y[i]-1, size=log.v,
                          prob=log.v/
                            (log.v+
                            log.nb11$mean[i]))
 
 log.pois.px[i] <- dpois(y[i], log.pois11$mean[i])
 
 log.nb.px[i] <- dnbinom(y[i], size=log.v,
                          prob=log.v/
                           (log.v+
                            log.nb11$mean[i]))
  
 garma.pois.Px[i] <- ppois(y[i], garma.pois11$mean[i])
 
 garma.pois.Px1[i] <- ppois(y[i]-1,
                            garma.pois11$mean[i])
 
 garma.nb.Px[i] <- pnbinom(y[i], size=garma.v,
                           prob=garma.v/
                             (garma.v+
                              garma.nb11$mean[i]))
 
 garma.nb.Px1[i] <- pnbinom(y[i]-1, size=garma.v,
                            prob=garma.v/
                              (garma.v+
                               garma.nb11$mean[i]))
 
 garma.pois.px[i] <- dpois(y[i], garma.pois11$mean[i])
 
 garma.nb.px[i] <- dnbinom(y[i], size=garma.v,
                           prob=garma.v/
                             (garma.v+
                              garma.nb11$mean[i]))
  
 glarma.pois.Px[i] <- ppois(y[i], glarma.pois11$mean[i])
 
 glarma.pois.Px1[i] <- ppois(y[i]-1, glarma.pois11$mean[i])
 
 glarma.nb.Px[i] <- pnbinom(y[i], size=glarma.v,
                            prob=glarma.v/
                              (glarma.v+
                               glarma.nb11$mean[i]))
 
 glarma.nb.Px1[i] <- pnbinom(y[i]-1, size=glarma.v,
                             prob=glarma.v/
                               (glarma.v+
                                glarma.nb11$mean[i]))
 
 glarma.pois.px[i] <- dpois(y[i], glarma.pois11$mean[i])
 
 glarma.nb.px[i] <- dnbinom(y[i], size=glarma.v,
                            prob=glarma.v/
                              (glarma.v+
                               glarma.nb11$mean[i]))
  
}


# Sharpness via scoring rules

kk <- 1000                           

log.pois.logs <- - log(log.pois.px)
garma.pois.logs <- - log(garma.pois.px)
glarma.pois.logs <- - log(glarma.pois.px)

log.nb.logs <- - log(log.nb.px)
garma.nb.logs <- - log(garma.nb.px)
glarma.nb.logs <- - log(glarma.nb.px)

log.pois.norm <- rep(0, n)
garma.pois.norm <- rep(0, n)
glarma.pois.norm <- rep(0, n)

log.nb.norm <- rep(0, n)
garma.nb.norm <- rep(0, n)
glarma.nb.norm <- rep(0, n)

for(i in 1:n){
  
 for(j in 1:kk){
    
  log.pois.norm[i] <- log.pois.norm[i] +
                      dpois(j, log.pois11$mean[i])^2
    
  log.nb.norm[i] <- log.nb.norm[i] +
                    dnbinom(j, size=log.v, 
                            prob=log.v/
                              (log.v+
                               log.nb11$mean[i]))^2
    
  garma.pois.norm[i] <- garma.pois.norm[i] + 
                        dpois(j, garma.pois11$mean[i])^2
    
  garma.nb.norm[i] <- garma.nb.norm[i] +
                      dnbinom(j, size=garma.v,
                              prob=garma.v/
                                (garma.v+
                                garma.nb11$mean[i]))^2
    
  glarma.pois.norm[i] <- glarma.pois.norm[i] +
                         dpois(j, glarma.pois11$mean[i])^2
  
  glarma.nb.norm[i] <- glarma.nb.norm[i] +
                       dnbinom(j, size=glarma.v,
                               prob=glarma.v/
                                 (glarma.v+
                                  glarma.nb11$mean[i]))^2
    
    
 }
  
}

log.pois.qs <- - 2*log.pois.px +
                     log.pois.norm
garma.pois.qs <- - 2*garma.pois.px +
                     garma.pois.norm
glarma.pois.qs <- - 2*glarma.pois.px +
                      glarma.pois.norm

log.nb.qs <- - 2*log.nb.px + log.nb.norm
garma.nb.qs <- - 2*garma.nb.px + garma.nb.norm
glarma.nb.qs <- - 2*glarma.nb.px + glarma.nb.norm

log.pois.sphs <- - log.pois.px / 
                     sqrt(log.pois.norm)
garma.pois.sphs <- - garma.pois.px /
                      sqrt(garma.pois.norm)
glarma.pois.sphs <- - glarma.pois.px /
                        sqrt(glarma.pois.norm)

log.nb.sphs <- - log.nb.px /
                   sqrt(log.nb.norm)
garma.nb.sphs <- - garma.nb.px /
                    sqrt(garma.nb.norm)
glarma.nb.sphs <- - glarma.nb.px /
                     sqrt(glarma.nb.norm)

log.diff.p <- rep(0, n)
log.diff.nb <- rep(0, n)
diffx <- rep(0, n)
log.rps.p <- rep(0, n)
log.rps.nb <- rep(0, n)

garma.diff.p <- rep(0, n)
garma.diff.nb <- rep(0, n)
garma.rps.p <- rep(0, n)
garma.rps.nb <- rep(0, n)

glarma.diff.p <- rep(0, n)
glarma.diff.nb <- rep(0, n)
glarma.rps.p <- rep(0, n)
glarma.rps.nb <- rep(0, n)

for(i in 1:n){

 for(j in 1:kk){
    
  log.diff.p[i] <- ppois(j, log.pois11$mean[i])
  log.diff.nb[i] <- pnbinom(j, size=log.v,
                            prob=log.v/
                              (log.v+
                               log.nb11$mean[i]))
  diffx[i] <- ifelse(y[i]<=j,1,0)
  
  log.rps.p[i]  <- log.rps.p[i] +
                   (log.diff.p[i] - diffx[i])^2
  
  log.rps.nb[i]  <- log.rps.nb[i] + 
                    (log.diff.nb[i] - diffx[i])^2
    
  garma.diff.p[i] <- ppois(j, garma.pois11$mean[i])
  garma.diff.nb[i] <- pnbinom(j, size=garma.v,
                              prob=garma.v/
                                (garma.v+
                                 garma.nb11$mean[i]))
  
  garma.rps.p[i]  <- garma.rps.p[i] +
                     (garma.diff.p[i] - diffx[i])^2
  
  garma.rps.nb[i]  <- garma.rps.nb[i] +
                      (garma.diff.nb[i] - diffx[i])^2
    
  glarma.diff.p[i] <- ppois(j, glarma.pois11$mean[i])
  
  glarma.diff.nb[i] <- pnbinom(j, size=glarma.v,
                               prob=glarma.v/
                                 (glarma.v+
                                    glarma.nb11$mean[i]))
  glarma.rps.p[i]  <- glarma.rps.p[i] +
                      (glarma.diff.p[i] - diffx[i])^2
  glarma.rps.nb[i]  <- glarma.rps.nb[i] +
                      (glarma.diff.nb[i] - diffx[i])^2
    
 }
  
}

log.pois.rps <- log.rps.p
garma.pois.rps <- garma.rps.p
glarma.pois.rps <- glarma.rps.p

log.nb.rps <- log.rps.nb
garma.nb.rps <- garma.rps.nb
glarma.nb.rps <- glarma.rps.nb

log.pois.dss <- (y-log.pois11$mean)^2/log.pois11$mean +
                 2*log(log.pois11$mean)
log.nb.dss <- (y-log.nb11$mean)^2/
              (log.nb11$mean*
                 (1+log.nb11$mean/log.v)) +
                  2*log(log.nb11$mean*
                          (1+log.nb11$mean/log.v))

log.pois.ses <- (y-log.pois11$mean)^2
log.nb.ses <- (y-log.nb11$mean)^2

garma.pois.dss <- (y-garma.pois11$mean)^2/
                  garma.pois11$mean +
                  2*log(garma.pois11$mean)

garma.nb.dss <- (y-garma.nb11$mean)^2/
                (garma.nb11$mean*
                   (1+garma.nb11$mean/garma.v)) +
                    2*log(garma.nb11$mean*
                            (1+garma.nb11$mean/garma.v))

garma.pois.ses <- (y-garma.pois11$mean)^2
garma.nb.ses <- (y-garma.nb11$mean)^2

glarma.pois.dss <- (y-glarma.pois11$mean)^2/
                    glarma.pois11$mean +
                    2*log(glarma.pois11$mean)

glarma.nb.dss <- (y-glarma.nb11$mean)^2/
                 (glarma.nb11$mean*
                    (1+glarma.nb11$mean/glarma.v)) +
                     2*log(glarma.nb11$mean*
                             (1+glarma.nb11$mean/
                                glarma.v))

glarma.pois.ses <- (y-glarma.pois11$mean)^2
glarma.nb.ses <- (y-glarma.nb11$mean)^2


# PRINT RESULTs

#mle results

round(c(log.pois11$est, 0 , log.pois11$aic,
        log.pois11$bic), 3)

round(log.pois11$std,3)

round(c(garma.pois11$est, 0, garma.pois11$aic,
        garma.pois11$bic), 3)

round(garma.pois11$std,3)

round(c(glarma.pois11$est, 0, glarma.pois11$aic,
        glarma.pois11$bic), 3)

round(glarma.pois11$std,3)

round(c(log.nb11$est, log.v1, log.nb11$aic,
        log.nb11$bic), 3)

round(log.nb11$std,3)

round(c(garma.nb11$est, garma.v1, garma.nb11$aic,
        garma.nb11$bic), 3)

round(garma.nb11$std,3)

round(c(glarma.nb11$est, glarma.v1, glarma.nb11$aic,
        glarma.nb11$bic), 3)

round(glarma.nb11$std,3)


#result of scoring rules

# logarithmic score  # quadratic score 
# spherical score  # ranked probability score
# Dawid-Sebastiani score


round(c(mean(log.pois.logs), mean(log.pois.qs),
        mean(log.pois.sphs), mean(log.pois.rps),
        mean(log.pois.dss)), 4)

round(c(mean(log.nb.logs), mean(log.nb.qs),
        mean(log.nb.sphs), mean(log.nb.rps),
        mean(log.nb.dss)), 4)

round(c(mean(garma.pois.logs), mean(garma.pois.qs),
        mean(garma.pois.sphs), mean(garma.pois.rps),
        mean(garma.pois.dss)), 4)

round(c(mean(garma.nb.logs), mean(garma.nb.qs),
        mean(garma.nb.sphs), mean(garma.nb.rps),
        mean(garma.nb.dss)), 4)

round(c(mean(glarma.pois.logs), mean(glarma.pois.qs),
        mean(glarma.pois.sphs), mean(glarma.pois.rps),
        mean(glarma.pois.dss)), 4)
round(c(mean(glarma.nb.logs), mean(glarma.nb.qs),
        mean(glarma.nb.sphs), mean(glarma.nb.rps),
        mean(glarma.nb.dss)), 4)


# plot time series and acf

par(mfrow=c(2,2))
plot(y ~ time, type = "l", xlab="Time", ylab="counts",    
     main="Daily COVID-19 deaths in Italy")
acf<- acf(y, xlab="Time", ylab="ACF",
          main="ACF COVID-19 deaths in Italy")
acf<- acf(log.pois11$res, xlab="Time", ylab="ACF",
          main="ACF standardized residuals Pois Log-AR")
acf<- acf(log.nb11$res, xlab="Time", ylab="ACF",
          main="ACF standardized residuals NB Log-AR")

# plot PIT

par(mfrow=c(2,3))
pit(y, log.pois.Px, log.pois.Px1, n.bins=10,
    y.max=3, my.title="PIT Poisson log-AR")
pit(y, garma.pois.Px, garma.pois.Px1, n.bins=10,
    y.max=3, my.title="PIT Poisson GARMA")
pit(y, glarma.pois.Px, glarma.pois.Px1, n.bins=10, 
    y.max=3, my.title="PIT Poisson GLARMA")
pit(y, log.nb.Px, log.nb.Px1, n.bins=10,
    y.max=3, my.title="PIT NB log-AR")
pit(y, garma.nb.Px, garma.nb.Px1, n.bins=10,
    y.max=3, my.title="PIT NB GARMA")
pit(y, glarma.nb.Px, glarma.nb.Px1, n.bins=10,
    y.max=3, my.title="PIT NB GLARMA")

#-----------------------------------------------------------
