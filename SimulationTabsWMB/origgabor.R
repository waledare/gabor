library("Sim.DiffProc")
dyn.load("gaborlib")
H0 <- 300 #Definitely pushes the edges out. But creates wide swings 
a <- (1/3) #The code is designed to take this into account implicitly
b <- (1/5)
I <- complex(real = 0, imaginary = 1)
br <- 500
gbmmu <- 0.5
gbmsigma <- 0.2
gbminit <- 0.8
simN <- 1500
gbmt0 <- 0
gbmT <- 1
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
#gbmdata <- GBM(N=simN, t0 = gbmt0, T=gbmT, x0 = gbminit, theta=gbmmu, sigma=gbmsigma) 
gbmdata <- ABM(N=simN, t0 = gbmt0, T=gbmT, x0 = gbminit, theta=gbmmu, sigma=gbmsigma) 
#get deltaXsquared
gbmdx2 <- matrix((diff(gbmdata))^2, ncol=1)
#generate observation times
ti <- seq(from = gbmt0, to =gbmT, by = (gbmT - gbmt0)/simN)[-(simN+1)] 


fnhkti <- function (f,h,k,ti){
    sapply(ti, function(x) {f(h,k,x)})
}
coeffmat <- matrix(complex(real = 0, imaginary =0), ncol=2*K0+1, nrow=(2*H0)+1)
coeffmatc <- matrix(complex(real = 0, imaginary =0), ncol=2*K0+1, nrow=(2*H0)+1)
sigmas <- matrix(complex(real = 0, imaginary =0), ncol =1, nrow=length(ti))

getvols <- function() {
    i <- complex(real = 0, imaginary =1)
    cmat <- matrix(coeffmatc, ncol =1)
    res <- .C("getsigmas",
       h0 = as.integer(H0),
       k0 = as.integer(K0),
       xsq = as.double(gbmdx2),
       reac = as.double(Re(cmat)),
       imac = as.double(Im(cmat)),
       t = as.double(ti),
       realsig = as.double(Re(sigmas)),
       imagsig = as.double(Im(sigmas)),
       nt = as.integer(length(ti)))
    #cmat <- res$reac +  res$imac * i
    #cmat <- matrix(mat, ncol=2*K0+1, nrow=(2*H0)+1)
    return (res)
}
res <- getvols();


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
#coeffmat <- getcoefmat()
#coeffmatc <- getcoefmatc()

shat <- function(t) {
    buck <- matrix(0, ncol=2*K0+1, nrow=(2*H0)+1)
    i <- 0
    for (h in (-H0):H0){
        i<- i+1
        j <- 0
        for (k in (-K0):K0){
            j<- j + 1
            buck[i,j] <- coeffmatc[i,j]*ghk(h,k,t)
        }
    }
    return (sum(buck))
}
#print(plot(ti, x<-sapply(ti,shat)))
x<-res$realsig
nt <- length(ti)
print(plot(ti[br:(nt-br)],x[br:(nt-br)]))

#print((gbmT-1)*Re(mean((x-(gbmsigma*gbmdata[-(simN+1)])^2)[1000:49000])))
#print(coeffmat)
#print(coeffmatc)
