shahMPCI <-
function(x,LSL,USL,Target,alpha){

if(missing(Target)) { # Estimating Target 
     Target<-LSL+(USL-LSL)/2
     }
 
if(missing(alpha)) { # Setting alpha 
     alpha<-0.0027
     } 

###First Component
p<-ncol(x)
m<-nrow(x)
means<-matrix(0,nrow=p,ncol<-1)

for(i in 1:p){
     means[i,1]<-mean(x[,i])
     }
s<-cov(x)
UPL<-matrix(0,nrow=p,ncol=1)
LPL<-matrix(0,nrow=p,ncol=1)

s1<-solve(s)

for(i in 1:p){
     s2<-solve(s)

     for (k in 1:nrow(s1)){
         for (j in 1:nrow(s1)){

             if( k==i | j==i){s2[k,j]<-NaN}
         }
     }
     s3<-matrix((s2[sapply(s2, function(s2) !any(is.na(s2)))]),nrow=p-1)

     LPL[i,1]<-means[i,1]-sqrt((det(s3)*(qchisq((1-alpha),df=p)))/det(solve(s)))
     UPL[i,1]<-means[i,1]+sqrt((det(s3)*(qchisq((1-alpha),df=p)))/det(solve(s)))
  }

spec<-cbind(USL,LSL)
dif.spec<-USL-LSL
dif.proc<-UPL-LPL

d1<-1
d2<-1
for(i in 1:p){
d1<-dif.spec[i]*d1
d2<-dif.proc[i]*d2
}
CpM<-(d1/d2)^(1/p) 

###Second Component
PV<-1-pf((t(Target-means)%*%solve(s)%*%(Target-means))*((m*(m-p))/(p*(m-1))),p,m-p)

###Third Component
a<-0
proc<-cbind(LPL,UPL)

for( i in 1:nrow(spec)){
     if((spec[i,1]>proc[i,1]) | (spec[i,2]<proc[i,2]))(a<-a+1)
     }
ifelse (a==0,LI<-1,LI<-0)

###Vector
outList = list ("Shahriari(1995) Multivariate Capability Vector","CpM"=CpM,"PV"=PV,"LI"=LI)
print(outList)

}

