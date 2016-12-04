library(reshape)
library(ggplot2)
library("yuima")
dyn.load("gaborlib")
H0 <- 100 #Definitely pushes the edges out. But creates wide swings --====-- 75 was good!
a <- (1/7)#70 works well for almost all 
b <- (1/25)
I <- complex(real = 0, imaginary = 1)
czero <- complex(real = 0, imaginary =0)
mu <- 0.5
sigma <- 0.2
theta <- 0.1
#Obs <- c(500, 5000, 50000)
Obs <- c(10000)
mcoeff <- 0.4 #Less than the exact modulus of continuity.
br <- .075 #Fraction of endpoints to discard
t0 <- 0
T <- 1
K0 <- ceiling((T - t0)/a)+1
#N <- 100 #Number of paths to simulate
N <- 1
coeffmat <- matrix(complex(real = 0, imaginary =0), ncol=2*K0+1, nrow=(2*H0)+1)
coeffmatc <- matrix(complex(real = 0, imaginary =0), ncol=2*K0+1, nrow=(2*H0)+1)
ise <- matrix(0, ncol =1)
getvols <- function(dx2, ti) {
                nti <- length(ti)
                cmat <- matrix(coeffmatc, ncol =1)
                sigmas <- matrix(czero, ncol =1, nrow=nti)
                res <- .C("getsigmas",
                    h0 = as.integer(H0),
                    k0 = as.integer(K0),
                    xsq = as.double(dx2),
                    reac = as.double(Re(cmat)), #no longer needed
                    imac = as.double(Im(cmat)), #no longer beeded
                    t = as.double(ti),
                    realsig = as.double(Re(sigmas)),
                    imagsig = as.double(Im(sigmas)),
                    nt = as.integer(length(ti)))
                    return (res$realsig)
                }

adata <- c()
gdata <- c()
cdata <- c()
odata <- c()
getmise <- function(sisqf, simf, nobs, outdata){
    actvol <- c()
    volhat <- c()
    mise <- c()
    biassq <- c()
    varsq <- c()
    ti <- c()
    j = 1
   
    for (obs in nobs) {
        #generate observation times
        ti <- seq(from = t0, to =T, by = (T - t0)/obs)[-(obs+1)] 
        nti <- length(ti)
        et <- ceiling(nti*br)
        actvol <- matrix(0, ncol = nti - 2*et + 1, nrow = N)
        volhat <- matrix(0, ncol = nti - 2*et + 1, nrow = N)
        for (i in 1:N){
            #Simulate  data. sigma is constant here
            data <- simf(obs) 
            #get deltaXsquared
            dx2 <- matrix((diff(data))^2, ncol=1)
            #get threshold
            u_n <- obs^{-(mcoeff)}
            dx2 <- matrix(ifelse(dx2 > u_n^2,0, dx2), ncol=1)            
            if(sisqf == "ou")
                actvol[i,] <- matrix(sigma**2, ncol= nti -2*et+1,nrow=1)
            if(sisqf == "gbm")
                actvol[i,] <- (((sigma*data))**2)[et:(nti-et)]
            if(sisqf == "cir")
                actvol[i,] <- ((sigma**2)*data)[et:(nti-et)]
            if(sisqf == "wbm")
                actvol[i,] <- ((sigma*(cos(data)/(1+data**2)))**2)[et:(nti-et)]

            volhat[i,] <- getvols(dx2,ti)[et:(nti-et)]
        }   
        biassq[j] <- mean((apply(volhat - actvol,c(2),mean))**2)
        varsq[j] <- mean(apply(volhat,c(2),
                               function(x) { n <- length(x); 
                               return (((n-1)/n)*var(x))}
                               ))
        mise[j] <- mean(apply((volhat - actvol)**2, c(1), mean))
        j<- j+1 #j indexes nobs
    }
   
    outdata = cbind(Times=ti[et:(nti-et)], Estimate= volhat[N,],Actual= actvol[N,])
    res <- cbind (mise, biassq,varsq)
    colnames(res)<-c(paste(sisqf, "-m"), paste(sisqf, "-b"), paste(sisqf, "-v"))
    return (list(res, outdata))
}

oufunc <- function (obs) {
    mod <- setModel(
                   drift=c("theta*(mu -  x)"), 
                   diffusion="sigma", 
                   jump.coeff="1", 
                   measure=list(intensity="5", 
                                df=list("dnorm(z,0,0.2)")), 
                   measure.type="CP", 
                   solve.variable="x", 
                   xinit="1")
    X<-(simulate(mod,sampling=setSampling(n=obs), true.p=list(theta=theta,mu = mu, sigma=sigma)))
    b<- (get.zoo.data(X))
    return (ts(b[[1]]))
}
gbmfunc <- function (obs) {
    mod <- setModel(
                   drift=c("mu * x"), 
                   diffusion="sigma * x", 
                   jump.coeff="1", 
                   measure=list(intensity="2", 
                                df=list("dnorm(z,0,0.15)")), 
                   measure.type="CP", 
                   solve.variable="x", 
                   xinit="1")
    X<-(simulate(mod,sampling=setSampling(n=obs), true.p=list(mu=mu,sigma=sigma)))
    b<- (get.zoo.data(X))
    return (ts(b[[1]]))
}
wbmfunc <- function (obs) {
    mod <- setModel(
                   drift=c("theta*(mu - sin(x))"), 
                   diffusion="sigma*(cos(x)/(1+x^2))", 
                   jump.coeff="1", 
                   measure=list(intensity="1", 
                                df=list("dnorm(z,0,0.2)")), 
                   measure.type="CP", 
                   solve.variable="x", 
                   xinit="1")
    X<-(simulate(mod,sampling=setSampling(n=obs), true.p=list(mu=mu,sigma=sigma)))
    b<- (get.zoo.data(X))
    return (ts(b[[1]]))
}
 gbmfunc <- function (obs) {
    mod <- setModel(
                   drift=c("mu * x"), 
                   diffusion="sigma * x", 
                   jump.coeff="1", 
                   measure=list(intensity="1", 
                                df=list("dnorm(z,0,0.2)")), 
                   measure.type="CP", 
                   solve.variable="x", 
                   xinit="1")
    X<-(simulate(mod,sampling=setSampling(n=obs), true.p=list(mu=mu,sigma=sigma)))
    b<- (get.zoo.data(X))
    return (ts(b[[1]]))
} 
cirfunc <- function (obs) {
    mod <- setModel(
                   drift=c("theta * (mu - x)"), 
                   diffusion="sigma * x^(0.5)", 
                   jump.coeff="1", 
                   measure=list(intensity="5", 
                                df=list("dnorm(z,0,0.05)")), 
                   measure.type="CP", 
                   solve.variable="x", 
                   xinit="1")
    X<-(simulate(mod,sampling=setSampling(n=obs), true.p=list(mu=mu,sigma=sigma, theta=theta)))
    b<- (get.zoo.data(X))
    return (ts(b[[1]]))
} 

#set.seed(123)
#wmise <- getmise("wbm", wbmfunc, Obs, adata)
#set.seed(123)
#omise <- getmise("ou", oufunc, Obs,odata)
set.seed(123)
gmise <- getmise("gbm", gbmfunc, Obs,gdata)
#set.seed(312)
#cmise <- getmise("cir", cirfunc, Obs,cdata)
#saveRDS(matrix(wmise,ncol=1),"WBMmise.rds")
#saveRDS(matrix(omise,ncol=1),"OUmise.rds")
#saveRDS(matrix(gmise,ncol=1),"GBMmise.rds")
#saveRDS(matrix(cmise,ncol=1),"CIRmise.rds")
X<-data.frame(gmise[[2]])
dl <- melt(X, id = "Times")
pa <- ggplot(data=dl, aes(x=Times, y=value, colour=variable)) + geom_line()
print(pa)

