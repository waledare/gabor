library("Sim.DiffProc")
H0 <- 25
a <- 1 #The code is designed to take this into account implicitly
b <- (1/3)
I <- complex(real = 0, imaginary = 1)
gbmmu <- 0.5
gbmsigma <- 0.2
gbminit <- 0.8
simN <- 40000
gbmt0 <- 0
gbmT <- 5
K0 <- ceiling((gbmT - gbmt0)/a)+1
#Gabor frame generator
g <- function (x){
    r <- 0
    if (x <1 & x >= 0)
        r <- x
    if (x < 2 & x >= 1)
       r <- (2 -x)
    return (r)
}
#Dual Gabor frame generator
gt <- function (x){
    r <- 0
    if (x <0 & x >= -1)
        r <- (2/3)*(x+1)
    if (x < 2 & x >= 0)
       r <- (1/3)*(2 -x)
    return (r)
}
#Gabor frame elements
ghk <- function (h, k , x){
    return (exp(2*pi*I*h*b*x)*g(x - k))
}
#Dual Gabor frame elements
gthk <- function(h,k,x){
    return (exp(2*pi*I*h*b*x)*gt(x - k))
}
#Simulate GBM data. sigma is constant here
gbmdata <- GBM(N=simN, t0 = gbmt0, T=gbmT, x0 = gbminit, theta=gbmmu, sigma=gbmsigma) 
#get deltaXsquared
gbmdx2 <- matrix((diff(gbmdata))^2, ncol=1)
#generate observation times
ti <- seq(from = gbmt0, to =gbmT, by = (gbmT - gbmt0)/simN)[-(simN+1)] 

fnhkti <- function (f,h,k,ti){
    sapply(ti, function(x) {f(h,k,x)})
}
coeffmat <- matrix(0, ncol=2*K0+1, nrow=(2*H0)+1)


getcoefmat <- function(){
    mat <- matrix(0, ncol=2*K0+1, nrow=(2*H0)+1)
    i <- 0
    for (h in (-H0):H0){
        i<- i+1
        j <- 0
        for (k in (-K0):K0){
            j<- j + 1
            res <- ((Conj(fnhkti(gthk,h,k,ti))%*%gbmdx2));
            mat[i,j] <-  res
        }
    }
    return (mat)
}
coeffmat <- getcoefmat()

shat <- function(t) {
    buck <- matrix(0, ncol=2*K0+1, nrow=(2*H0)+1)
    i <- 0
    for (h in (-H0):H0){
        i<- i+1
        j <- 0
        for (k in (-K0):K0){
            j<- j + 1
            buck[i,j] <- coeffmat[i,j]*ghk(h,k,t)
        }
    }
    return (sum(buck))
}
print(plot(ti, x<-sapply(ti,shat)))
print((gbmT-1)*Re(mean((x-(gbmsigma*gbmdata[-(simN+1)])^2)[1000:49000])))
