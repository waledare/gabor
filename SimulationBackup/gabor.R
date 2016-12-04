library("Sim.DiffProc")
library("sde")
#library("reshape2")
#library("ggplot2")
dyn.load("gaborlib")
H0 <- 250 #Definitely pushes the edges out. But creates wide swings 
a <- (1/3) #The code is designed to take this into account implicitly
b <- (1/3)
I <- complex(real = 0, imaginary = 1)
czero <- complex(real = 0, imaginary =0)
mu <- 0.5
sigma <- 0.2
rev <- 0.1
init <- 0.8
#Obs <- c(500, 5000, 50000)
Obs <- c(200000)
br <- .05 #Number of data in each path
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
                        
            if(sisqf == "abm" | sisqf == "ou")
            actvol[i,] <- matrix(sigma**2, ncol= nti -2*et+1,nrow=1)
            if(sisqf == "gbm")
                actvol[i,] <- (((sigma*data))**2)[et:(nti-et)]
            if(sisqf == "cir")
                actvol[i,] <- ((sigma**2)*data)[et:(nti-et)]

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

abmfunc <- function (obs) ABM(N=obs, t0 = t0, T=T, x0 = init, theta=mu, sigma=sigma)
oufunc <- function (obs) OU(N=obs, t0 = t0, T=T, x0 = init, r=mu, sigma=sigma)
gbmfunc <- function (obs) Sim.DiffProc::GBM(N=obs, t0 = t0, T=T, x0 = init, theta=mu, sigma=sigma)
cirfunc <- function (obs) {
    return (as.matrix(sde.sim(X0=init,N=obs,t0=t0,T=T, theta=c(rev, mu, sigma), model ="CIR"))) 
}
amise <- getmise("abm", abmfunc, Obs, adata)
omise <- getmise("ou", oufunc, Obs,odata)
gmise <- getmise("gbm", gbmfunc, Obs,gdata)
cmise <- getmise("cir", cirfunc, Obs,cdata)
#saveRDS(matrix(amise,ncol=1),"ABMmise.rds")
#saveRDS(matrix(omise,ncol=1),"OUmise.rds")
#saveRDS(matrix(gmise,ncol=1),"GBMmise.rds")
#saveRDS(matrix(cmise,ncol=1),"CIRmise.rds")
#nt <- length(ti)
#print(plot(ti[br:(nt-br)],sig[1][br:(nt-br)]))

#
#mises <- cbind(amise,omise,gmise,cmise)
#miseout <- (apply(mises, c(1),num))


