taamMPCI <-
function(x,LSL,USL,Target,alpha){

if(missing(Target)) { # Estimating Target 
     Target<-LSL+(USL-LSL)/2
     }
 
if(missing(alpha)) { # Setting alpha 
     alpha<-0.0027
     }

p<-ncol(x)
n<-nrow(x)
means<-matrix(0,nrow=p,ncol<-1)

for(i in 1:p){
     means[i,1]<-mean(x[,i])
     }
s<-cov(x)
UPL<-matrix(0,nrow=p,ncol=1)
LPL<-matrix(0,nrow=p,ncol=1)

VTR<-2*(prod((USL-LSL)/2))*pi^(p/2)/(p*gamma(p/2)) # Vol. Tolerance Region   
VE<-det(s)^0.5*((pi*qchisq(1-alpha,p))^(p/2))/gamma(p/2+1) # Vol. Estimated 99.73% Process Region
Cp<-VTR/VE
D<-(1+n/(n-1)*(t(Target-means)%*%solve(s)%*%(Target-means)))^0.5
MCpm<-Cp/D

outList = list ("Taam et al. (1993) Multivariate Capability Index (MCpm)","MCpm"=MCpm)
print(outList)
}

