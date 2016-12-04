library("Sim.DiffProc")
library("sde")
dyn.load("gaborlib")
H0 <- 50 #Definitely pushes the edges out. But creates wide swings 
a <- (1/3) #The code is designed to take this into account implicitly
b <- (1/5)
I <- complex(real = 0, imaginary = 1)
mu <- 0.5
sigma <- 0.2
rev <- 0.1
init <- 0.8
Obs <- c(500, 5000, 50000)
br <- .05 #Number of data in each path
t0 <- 0
T <- 1
K0 <- ceiling((T - t0)/a)+1
N <- 100 #Number of paths to simulate
coeffmat <- matrix(complex(real = 0, imaginary =0), ncol=2*K0+1, nrow=(2*H0)+1)
coeffmatc <- matrix(complex(real = 0, imaginary =0), ncol=2*K0+1, nrow=(2*H0)+1)
ise <- matrix(0, ncol =1)

getmise <- function(sisqf, simf, nobs){
    mise <- c()
    j = 1
    for (obs in nobs) {
        for (i in 1:N){
            #Simulate  data. sigma is constant here
            data <- simf(obs) 
            #get deltaXsquared
            dx2 <- matrix((diff(data))^2, ncol=1)
            #generate observation times
            ti <- seq(from = t0, to =T, by = (T - t0)/obs)[-(obs+1)] 
            nti <- length(ti)
            et <- ceiling(nti*br)
            sigmas <- matrix(complex(real = 0, imaginary =0), ncol =1, nrow=nti)
            getvols <- function() {
                cmat <- matrix(coeffmatc, ncol =1)
                res <- .C("getsigmas",
                    h0 = as.integer(H0),
                    k0 = as.integer(K0),
                    xsq = as.double(dx2),
                    reac = as.double(Re(cmat)),
                    imac = as.double(Im(cmat)),
                    t = as.double(ti),
                    realsig = as.double(Re(sigmas)),
                    imagsig = as.double(Im(sigmas)),
                    nt = as.integer(length(ti)))
                    return (res$realsig)
                }   

            if(sisqf == "abm" | sisqf == "ou")
                sisq <- sigma**2
            if(sisqf == "gbm")
                sisq <- ((sigma*data)[-length(data)])**2
            if(sisqf == "cir")
                sisq <- (sigma**2)*data[-length(data)]

            ise[i] <- mean(((getvols() - sisq)**2)[et:(nti-et)])
        }   
        mise[j] <- mean(ise)
        j<- j+1
    }
    return (mise)
}

abmfunc <- function (obs) ABM(N=obs, t0 = t0, T=T, x0 = init, theta=mu, sigma=sigma)
oufunc <- function (obs) OU(N=obs, t0 = t0, T=T, x0 = init, r=mu, sigma=sigma)
gbmfunc <- function (obs) Sim.DiffProc::GBM(N=obs, t0 = t0, T=T, x0 = init, theta=mu, sigma=sigma)
cirfunc <- function (obs) {
    return (as.matrix(sde.sim(X0=init,N=obs,t0=t0,T=T, theta=c(rev, mu, sigma), model ="CIR"))) 
}
amise <- getmise("abm", abmfunc, Obs)
omise <- getmise("ou", oufunc, Obs)
gmise <- getmise("gbm", gbmfunc, Obs)
cmise <- getmise("cir", cirfunc, Obs)
#saveRDS(matrix(amise,ncol=1),"ABMmise.rds")
#saveRDS(matrix(omise,ncol=1),"OUmise.rds")
#saveRDS(matrix(gmise,ncol=1),"GBMmise.rds")
#saveRDS(matrix(cmise,ncol=1),"CIRmise.rds")
#nt <- length(ti)
#print(plot(ti[br:(nt-br)],sig[1][br:(nt-br)]))
print(amise)
print(omise)
print(gmise)
print(cmise)

num <- function(x) return (paste("\num[round-precision=2,round-mode=figures]{",x, "} &"))
mises <- cbind(amise,omise,gmise,cmise)
miseout <- (apply(mises, c(1),num))
#cat(miseout,"miseout.txt")
