mpci<- function(index,x,LSL,USL,Target,npc,alpha,Method,perc,graphic,...){   
require(ellipse)

if (missing(index)){
       stop("Attention: Specify the index to use (e.g.: index=shah, index=taam, etc. )")}
     
 if (index=="shah") {
	 
 
     if(missing(Target)) { # Estimating Target 
         Target<-LSL+(USL-LSL)/2
       }
 
     if(missing(alpha)) { # Setting alpha 
         alpha<-0.0027
         } 

     ###First Component
     p<-ncol(x)
     m<-nrow(x)
     Xmv<-means<-matrix(0,nrow=p,ncol<-1)

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

     spec<-cbind(LSL,USL)
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
     a <- 0
     proc <- cbind(LPL,UPL)

     for( i in 1:nrow(spec)){
         if((spec[i,1] > proc[i,1]) | (spec[i,2] < proc[i,2]))(a<-a + 1)
         }
     ifelse (a == 0,LI <- 1,LI <- 0)

     ###Vector
     outList = list ("Shahriari et al. (1995) Multivariate Capability Vector","CpM"=CpM,"PV"=PV,"LI"=LI)
     print(outList)
  

	} 

	 	 

 if (index=="taam") {	 
	 
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
	 
	 
 if (index=="wang") {	 
	 
     SL<-cbind(LSL,USL) # matrix of specifications
     TOL<-USL-LSL # difference in specifications
     m<-length(LSL) #dimensions

     SS<-cov(x) # covariance matrix 
     n<-nrow(x)
     Xmv<-colMeans(x) #mean vector
     
     if(missing(perc)){
	   perc<-0.8
	     }
	 
     if(missing(Target)) { # Estimating Target 
          Target<-LSL+TOL/2
          }
 
     Ue<-eigen(SS, symmetric=TRUE, EISPACK = TRUE)$vectors    # eigenvectors
     DDe<-eigen(SS, symmetric=TRUE,EISPACK = TRUE)$values     # eigenvalues

     if(!missing(npc)) { 
          if(npc<=0 | npc>m | !is.numeric(npc) | npc != as.integer(npc) | length(npc) > 1){
              stop("Attention: the number of principal components (npc) must be a integer between 1 
          and the number of quality characteristics")
              }
          }


     if(missing(npc)) {    #number of principal components
          #Modified Algorithm of Rencher,A.C.(2002) Methods of Multivariate Analysis. John Wiley and Sons.
          #12.6 DECIDING HOW MANY COMPONENTS TO RETAIN
  
 
         if(missing(alpha)) { # Setting alpha 
              alpha<-0.05
              }
 
         if(!missing(alpha)) {  
                 if(alpha<0 | alpha>1){
                     stop("Attention: the significance level (alpha) must be between 0 and 1")
                     }
                 } 

         if(missing(Method)) { # Method to select the number of principal components
             Method<-"Percentage"   
             }
     
         if(!missing(Method)) { 
                 if((Method<=0 | Method>5) & (Method!="Percentage"& Method!="Average"& Method!="Scree"& Method!="Bartlett.test"& Method!="Anderson.test")){
                     stop("Attention: the Method must be a integer between 1 and 5 or one of the followings:
                     Percentage, Average, Scree, Bartlett.test or Anderson.test")
                     }
                 }
     
      ########################### 
         if (Method=="Percentage"|Method==1) {
             DDeV<-matrix(0,3,ncol=length(DDe)) # matrix of the Proportion of explained variance by the PCs
         DDeV[1,]<-seq(1:length(DDe))
         for (i in 1:length(DDe)){
             DDeV[2,i]<-DDe[i]/sum(DDe)
             ifelse(i==1,DDeV[3,i]<-DDeV[2,1],DDeV[3,i]<- DDeV[2,i]+ DDeV[3,i-1])
             }
         npc<-1   #number of Principal Component(npc)
         j<-1
         while (DDeV[3,j]<perc) {
             npc<-npc+1
             j<-j+1
             }  
         } 

     ###########################
         if (Method=="Average"|Method==2){
             mv<-matrix(0,nrow=1,ncol=length(DDe)) #mean vector
             for(i in 1:length(DDe)) {
                 ifelse (DDe[i]>mean(DDe),mv[1,i]<-1,mv[1,i]<-0)
                 }
             npc<-sum(mv)
             }
 
     ##########################
         if (Method=="Scree"|Method==3){ 
             plot(DDe,main="Scree graph for eigenvalue",xlab="Eigenvalue number",ylab="Eigenvalue size",type="o",pch=1,lty=1)
             cat("\n","Enter the number of principal component(npc) according to the scree graph:: ","\n")
             npc<-scan(n=1) 
         
             if(npc<=0 | npc>m | !is.numeric(npc) | npc != as.integer(npc) | length(npc) > 1){
                 stop("Attention: the number of principal components (npc) must be a integer between 1 
                 and the number of quality characteristics")
                 }
  
             }

     ##########################
 
     if (Method=="Bartlett.test"|Method==4){
         chi.p<-matrix(0,1,(length(DDe))) # (practical) chi-squared value
         chi.t<-matrix(0,1,(length(DDe))) # (theoretical) chi-squared value
         npc.matrix<-matrix(0,1,(length(DDe)))

         for (i in 1:length(DDe)){
             DD<-DDe[i:m]
             chi.p[i]<-(n-(2*m+11)/6)*((m-i+1)*log(mean(DD))-sum(log(DD)))
             chi.t[i]<- qchisq(1-alpha,((length(DD)-1)*(length(DD)+2)/2)) 
             }

         for (i in 1:length(chi.p)){
             ifelse(chi.p[i]>chi.t[i],npc.matrix[i]<-1,npc.matrix[i]<-0)
             }

         npc<-sum(npc.matrix) # number of principal components
         if (npc==0){
     stop("There are no difference between principal components according to Bartlett's Test.
     Please use another method (1,2,3 or 5)")}
  
     }

     ##########################
     if (Method=="Anderson.test"|Method==5){
     
         chi.p<-matrix(0,1,(length(DDe))) # (practical) chi-squared value
         chi.t<-matrix(0,1,(length(DDe))) # (theoretical) chi-squared value
         npcmatrix<-matrix(0,1,(length(DDe)))

         for (i in 1:length(DDe)){
             DD<-DDe[i:m]
             chi.p[i]<-(n-1)*length(DD)*log(sum(DD)/length(DD))-(n-1)*sum(log(DD))
             chi.t[i]<- qchisq(1-alpha,((length(DD)-1)*(length(DD)+2)/2)) 
     }

         for (i in 1:length(chi.p)){
             ifelse(chi.p[i]>chi.t[i],npcmatrix[i]<-1,npcmatrix[i]<-0)
             }

         npc<-sum(npcmatrix) # number of principal components
         if (npc==0){
     stop("There are no difference between principal components according to Anderson's Test.
     Please use another method (1,2,3 or 4)")}
     }

     }

     SLpce<-t(Ue)%*%SL  #matrix of spec. limits for each PC (rows:LSLPC, USLPC)

     ##vector of sample means for each PC
     Xmvpce<-t(Ue)%*%Xmv

     ##vector of target values for each PC
     Targetpce<-t(Ue)%*%Target 

     ##COMPUTING MATRIX OF Cp 
     numCpe<-abs(SLpce[,2]-SLpce[,1]) 

     ##this is matrix of the Cppc
     MatCpe<-numCpe/(6*sqrt(DDe))

     ##COMPUTING MATRIX OF Cpm 
     MatCpme<-numCpe/(6*sqrt(DDe+(Xmvpce-Targetpce)^2))

     ##COMPUTING MATRIX OF Cpk

     ##this is matrix of the UPPER
     Cppcue<-abs((SLpce[,2]-Xmvpce))/(3*sqrt(DDe)) 
     ##this is matrix of the LOWER
     Cppcle<-abs((Xmvpce-SLpce[,1]))/(3*sqrt(DDe)) 
     appo1e<-cbind(Cppcle,Cppcue)
     MatCpke<-matrix(NA, nrow=m, ncol=1)

     ##COMPUTING MATRIX OF Cpmk 

     ##this is matrix of the UPPER
     Cpmpcue<-abs((SLpce[,2]-Xmvpce))/(3*sqrt(DDe+(Xmvpce-Targetpce)^2)) 
     #this is matrix of the LOWER
     Cpmpcle<-abs((Xmvpce-SLpce[,1]))/(3*sqrt(DDe+(Xmvpce-Targetpce)^2)) 
     appo2e<-cbind(Cpmpcle,Cpmpcue)
     MatCpmke<-matrix(NA, nrow=m, ncol=1)

     ##  MatCpke MatCpmke
     for  (j in 1:m) {
                  MatCpke[j,]<-min(appo1e[j,]);
                  MatCpmke[j,]<-min(appo2e[j,]);
                 }

     ##WANG-CHEN
     MCpe<-(prod(MatCpe[1:npc]))^(1/npc)
     MCpke<-(prod(MatCpke[1:npc]))^(1/npc)
     MCpme<-(prod(MatCpme[1:npc]))^(1/npc)
     MCpmke<-(prod(MatCpmke[1:npc]))^(1/npc)

     outList = list ("Wang and Chen (1998) Multivariate Process Capability Indices(PCI) based on PCA",
     "number of principal components"=npc,"MCp"=MCpe,"MCpk"=MCpke,"MCpm"=MCpme,"MCpmk"=MCpmke)
     print(outList)
      
     }

     
	 
if (index=="xeke") {	 

     SL<-cbind(LSL,USL) # matrix of specifications
TOL<-USL-LSL # difference in specifications
m<-length(LSL) #dimensions

SS<-cov(x) # covariance matrix
n<-nrow(x)
Xmv<-colMeans(x) #mean vector

     if(missing(perc)){
	   perc<-0.8
	     }
		 
if(missing(Target)) { # Estimating Target 
     Target<-LSL+TOL/2
     }
 
Ue<-eigen(SS, symmetric=TRUE, EISPACK = TRUE)$vectors    # eigenvectors
DDe<-eigen(SS, symmetric=TRUE,EISPACK = TRUE)$values     # eigenvalues

if(!missing(npc)) { 
     if(npc<=0 | npc>m | !is.numeric(npc) | npc != as.integer(npc) | length(npc) > 1){
         stop("Attention: the number of principal components (npc) must be a integer between 1 
     and the number of quality characteristics")
         }
     }


if(missing(npc)) {    #number of principal components
     #Modified Algorithm of Rencher,A.C.(2002) Methods of Multivariate Analysis. John Wiley and Sons.
     #12.6 DECIDING HOW MANY COMPONENTS TO RETAIN
  
 
     if(missing(alpha)) { # Setting alpha 
         alpha<-0.05
         }
 
 if(!missing(alpha)) {  
         if(alpha<0 | alpha>1){
             stop("Attention: the significance level (alpha) must be between 0 and 1")
             }
         } 

     if(missing(Method)) { # Method to select the number of principal components
         Method<-"Percentage"   
         }
     
 if(!missing(Method)) { 
              if((Method<=0 | Method>5) & (Method!="Percentage"& Method!="Average"& Method!="Scree"&
  Method!="Bartlett.test"& Method!="Anderson.test")){
                 stop("Attention: the Method must be a integer between 1 and 5 or one of the followings:
         Percentage, Average, Scree, Bartlett.test or Anderson.test")
             }
         }
     
      ########################### 
     if (Method=="Percentage"|Method==1) {
         DDeV<-matrix(0,3,ncol=length(DDe)) # matrix of the Proportion of explained variance by the PCs
         DDeV[1,]<-seq(1:length(DDe))
         for (i in 1:length(DDe)){
             DDeV[2,i]<-DDe[i]/sum(DDe)
             ifelse(i==1,DDeV[3,i]<-DDeV[2,1],DDeV[3,i]<- DDeV[2,i]+ DDeV[3,i-1])
             }
         npc<-1   #number of Principal Component(npc)
         j<-1
         while (DDeV[3,j]<perc) {
             npc<-npc+1
             j<-j+1
             }  
         } 

     ###########################
     if (Method=="Average"|Method==2){
         mv<-matrix(0,nrow=1,ncol=length(DDe)) #mean vector
         for(i in 1:length(DDe)) {
             ifelse (DDe[i]>mean(DDe),mv[1,i]<-1,mv[1,i]<-0)
         }
         npc<-sum(mv)
         }
 
     ##########################
     if (Method=="Scree"|Method==3){
         plot(DDe,main="Scree graph for eigenvalue",xlab="Eigenvalue number",ylab="Eigenvalue size",type="o",pch=1,lty=1)
         cat("\n","Enter the number of principal component(npc) according to the scree graph:: ","\n")
         npc<-scan(n=1) 
         
 if(npc<=0 | npc>m | !is.numeric(npc) | npc != as.integer(npc) | length(npc) > 1){
             stop("Attention: the number of principal components (npc) must be a integer between 1 
         and the number of quality characteristics")
             }
  
 }

     ##########################
 
     if (Method=="Bartlett.test"|Method==4){
         chi.p<-matrix(0,1,(length(DDe))) # (practical) chi-squared value
         chi.t<-matrix(0,1,(length(DDe))) # (theoretical) chi-squared value
         npc.matrix<-matrix(0,1,(length(DDe)))

         for (i in 1:length(DDe)){
             DD<-DDe[i:m]
             chi.p[i]<-(n-(2*m+11)/6)*((m-i+1)*log(mean(DD))-sum(log(DD)))
             chi.t[i]<- qchisq(1-alpha,((length(DD)-1)*(length(DD)+2)/2)) 
             }

         for (i in 1:length(chi.p)){
             ifelse(chi.p[i]>chi.t[i],npc.matrix[i]<-1,npc.matrix[i]<-0)
             }

         npc<-sum(npc.matrix) # number of principal components
         if (npc==0){
     stop("There are no difference between principal components according to Bartlett's Test.
     Please use another method (1,2,3 or 5)")}
  
 }

     ##########################
     if (Method=="Anderson.test"|Method==5){
     
         chi.p<-matrix(0,1,(length(DDe))) # (practical) chi-squared value
         chi.t<-matrix(0,1,(length(DDe))) # (theoretical) chi-squared value
         npcmatrix<-matrix(0,1,(length(DDe)))

         for (i in 1:length(DDe)){
             DD<-DDe[i:m]
             chi.p[i]<-(n-1)*length(DD)*log(sum(DD)/length(DD))-(n-1)*sum(log(DD))
             chi.t[i]<- qchisq(1-alpha,((length(DD)-1)*(length(DD)+2)/2)) 
     }

         for (i in 1:length(chi.p)){
             ifelse(chi.p[i]>chi.t[i],npcmatrix[i]<-1,npcmatrix[i]<-0)
             }

         npc<-sum(npcmatrix) # number of principal components
         if (npc==0){
     stop("There are no difference between principal components according to Anderson's Test.
     Please use another method (1,2,3 or 4)")}
     }

    }

SLpce<-t(Ue)%*%SL  #matrix of spec. limits for each PC (rows:LSLPC, USLPC)

##vector of sample means for each PC
Xmvpce<-t(Ue)%*%Xmv

##vector of target values for each PC
Targetpce<-t(Ue)%*%Target 

##COMPUTING MATRIX OF Cp 
numCpe<-abs(SLpce[,2]-SLpce[,1]) 


##this is matrix of the Cppc
MatCpe<-numCpe/(6*sqrt(DDe))

##COMPUTING MATRIX OF Cpm 
MatCpme<-numCpe/(6*sqrt(DDe+(Xmvpce-Targetpce)^2))

##COMPUTING MATRIX OF Cpk 

##this is matrix of the UPPER
Cppcue<-abs((SLpce[,2]-Xmvpce))/(3*sqrt(DDe)) 

##this is matrix of the LOWER
Cppcle<-abs((Xmvpce-SLpce[,1]))/(3*sqrt(DDe)) 

appo1e<-cbind(Cppcle,Cppcue)
MatCpke<-matrix(NA, nrow=m, ncol=1)

##COMPUTING MATRIX OF Cpmk 

##this is matrix of the UPPER
Cpmpcue<-abs((SLpce[,2]-Xmvpce))/(3*sqrt(DDe+(Xmvpce-Targetpce)^2)) 

#this is matrix of the LOWER
Cpmpcle<-abs((Xmvpce-SLpce[,1]))/(3*sqrt(DDe+(Xmvpce-Targetpce)^2)) 

appo2e<-cbind(Cpmpcle,Cpmpcue)
MatCpmke<-matrix(NA, nrow=m, ncol=1)

##  MatCpke MatCpmke
for  (j in 1:m) {
                  MatCpke[j,]<-min(appo1e[j,]);
                  MatCpmke[j,]<-min(appo2e[j,]);
                 }
##SUM OF WEIGHTS
su_we<-sum(DDe[1:npc])

##Xekalaki-Perakis
vMXCpe<-MatCpe*DDe
MXCpe<-(sum(vMXCpe[1:npc]))/su_we

vMXCpke<-MatCpke*DDe
MXCpke<-(sum(vMXCpke[1:npc]))/su_we

vMXCpme<-MatCpme*DDe
MXCpme<-(sum(vMXCpme[1:npc]))/su_we

vMXCpmke<-MatCpmke*DDe
MXCpmke<-(sum(vMXCpmke[1:npc]))/su_we

outList = list ("Xekalaki and Perakis (2002) Multivariate Process Capability Indices(PCI) based on PCA",
"number of principal components"=npc,
"MCp"=MXCpe,"MCpk"=MXCpke,"MCpm"=MXCpme,"MCpmk"=MXCpmke)
print(outList)
     }
	 
	

if (index=="wangw") {		

  SL<-cbind(LSL,USL) # matrix of specifications
TOL<-USL-LSL # difference in specifications
m<-length(LSL) #dimensions

SS<-cov(x) # covariance matrix
n<-nrow(x)
Xmv<-colMeans(x) #mean vector

     if(missing(perc)){
	   perc<-0.8
	     }

if(missing(Target)) { # Estimating Target 
     Target<-LSL+TOL/2
     }
 
Ue<-eigen(SS, symmetric=TRUE, EISPACK = TRUE)$vectors    # eigenvectors
DDe<-eigen(SS, symmetric=TRUE,EISPACK = TRUE)$values     # eigenvalues

if(!missing(npc)) { 
     if(npc<=0 | npc>m | !is.numeric(npc) | npc != as.integer(npc) | length(npc) > 1){
         stop("Attention: the number of principal components (npc) must be a integer between 1 
     and the number of quality characteristics")
         }
     }


if(missing(npc)) {    #number of principal components
     #Modified Algorithm of Rencher,A.C.(2002) Methods of Multivariate Analysis. John Wiley and Sons.
     #12.6 DECIDING HOW MANY COMPONENTS TO RETAIN
  
 
     if(missing(alpha)) { # Setting alpha 
         alpha<-0.05
         }
 
 if(!missing(alpha)) {  
         if(alpha<0 | alpha>1){
             stop("Attention: the significance level (alpha) must be between 0 and 1")
             }
         } 

     if(missing(Method)) { # Method to select the number of principal components
         Method<-"Percentage"   
         }
     
 if(!missing(Method)) { 
              if((Method<=0 | Method>5) & (Method!="Percentage"& Method!="Average"& Method!="Scree"&
  Method!="Bartlett.test"& Method!="Anderson.test")){
                 stop("Attention: the Method must be a integer between 1 and 5 or one of the followings:
         Percentage, Average, Scree, Bartlett.test or Anderson.test")
             }
         }
     
      ########################### 
     if (Method=="Percentage"|Method==1) {
         DDeV<-matrix(0,3,ncol=length(DDe)) # matrix of the Proportion of explained variance by the PCs
         DDeV[1,]<-seq(1:length(DDe))
         for (i in 1:length(DDe)){
             DDeV[2,i]<-DDe[i]/sum(DDe)
             ifelse(i==1,DDeV[3,i]<-DDeV[2,1],DDeV[3,i]<- DDeV[2,i]+ DDeV[3,i-1])
             }
         npc<-1   #number of Principal Component(npc)
         j<-1
         while (DDeV[3,j]<perc) {
             npc<-npc+1
             j<-j+1
             }  
         } 

     ###########################
     if (Method=="Average"|Method==2){
         mv<-matrix(0,nrow=1,ncol=length(DDe)) #mean vector
         for(i in 1:length(DDe)) {
             ifelse (DDe[i]>mean(DDe),mv[1,i]<-1,mv[1,i]<-0)
         }
         npc<-sum(mv)
         }
 
     ##########################
     if (Method=="Scree"|Method==3){
         plot(DDe,main="Scree graph for eigenvalue",xlab="Eigenvalue number",ylab="Eigenvalue size",type="o",pch=1,lty=1)
         cat("\n","Enter the number of principal component(npc) according to the scree graph:: ","\n")
         npc<-scan(n=1) 
         
 if(npc<=0 | npc>m | !is.numeric(npc) | npc != as.integer(npc) | length(npc) > 1){
             stop("Attention: the number of principal components (npc) must be a integer between 1 
         and the number of quality characteristics")
             }
  
 }

     ##########################
 
     if (Method=="Bartlett.test"|Method==4){
         chi.p<-matrix(0,1,(length(DDe))) # (practical) chi-squared value
         chi.t<-matrix(0,1,(length(DDe))) # (theoretical) chi-squared value
         npc.matrix<-matrix(0,1,(length(DDe)))

         for (i in 1:length(DDe)){
             DD<-DDe[i:m]
             chi.p[i]<-(n-(2*m+11)/6)*((m-i+1)*log(mean(DD))-sum(log(DD)))
             chi.t[i]<- qchisq(1-alpha,((length(DD)-1)*(length(DD)+2)/2)) 
             }

         for (i in 1:length(chi.p)){
             ifelse(chi.p[i]>chi.t[i],npc.matrix[i]<-1,npc.matrix[i]<-0)
             }

         npc<-sum(npc.matrix) # number of principal components
         if (npc==0){
     stop("There are no difference between principal components according to Bartlett's Test.
     Please use another method (1,2,3 or 5)")}
  
 }

     ##########################
     if (Method=="Anderson.test"|Method==5){
     
         chi.p<-matrix(0,1,(length(DDe))) # (practical) chi-squared value
         chi.t<-matrix(0,1,(length(DDe))) # (theoretical) chi-squared value
         npcmatrix<-matrix(0,1,(length(DDe)))

         for (i in 1:length(DDe)){
             DD<-DDe[i:m]
             chi.p[i]<-(n-1)*length(DD)*log(sum(DD)/length(DD))-(n-1)*sum(log(DD))
             chi.t[i]<- qchisq(1-alpha,((length(DD)-1)*(length(DD)+2)/2)) 
     }

         for (i in 1:length(chi.p)){
             ifelse(chi.p[i]>chi.t[i],npcmatrix[i]<-1,npcmatrix[i]<-0)
             }

         npc<-sum(npcmatrix) # number of principal components
         if (npc==0){
     stop("There are no difference between principal components according to Anderson's Test.
     Please use another method (1,2,3 or 4)")}
 }

     }

SLpce<-t(Ue)%*%SL  #matrix of spec. limits for each PC (rows:LSLPC, USLPC)

##vector of sample means for each PC
Xmvpce<-t(Ue)%*%Xmv

##vector of target values for each PC
Targetpce<-t(Ue)%*%Target 

##COMPUTING MATRIX OF Cp 
numCpe<-abs(SLpce[,2]-SLpce[,1])  

##this is matrix of the Cppc
MatCpe<-numCpe/(6*sqrt(DDe))

##COMPUTING MATRIX OF Cpm 
MatCpme<-numCpe/(6*sqrt(DDe+(Xmvpce-Targetpce)^2))

##COMPUTING MATRIX OF Cpk 

##this is matrix of the UPPER
Cppcue<-abs((SLpce[,2]-Xmvpce))/(3*sqrt(DDe)) 
##this is matrix of the LOWER
Cppcle<-abs((Xmvpce-SLpce[,1]))/(3*sqrt(DDe)) 
appo1e<-cbind(Cppcle,Cppcue)
MatCpke<-matrix(NA, nrow=m, ncol=1)

##COMPUTING MATRIX OF Cpmk 

##this is matrix of the UPPER
Cpmpcue<-abs((SLpce[,2]-Xmvpce))/(3*sqrt(DDe+(Xmvpce-Targetpce)^2)) 
#this is matrix of the LOWER
Cpmpcle<-abs((Xmvpce-SLpce[,1]))/(3*sqrt(DDe+(Xmvpce-Targetpce)^2)) 
appo2e<-cbind(Cpmpcle,Cpmpcue)
MatCpmke<-matrix(NA, nrow=m, ncol=1)

##  MatCpke MatCpmke
for  (j in 1:m) {
                  MatCpke[j,]<-min(appo1e[j,]);
                  MatCpmke[j,]<-min(appo2e[j,]);
                 }

##SUM OF WEIGHTS
su_we<-sum(DDe[1:npc])

##CH Wang
vMWCpe<-(MatCpe)^(DDe)
MWCpe<-(prod(vMWCpe[1:npc]))^(1/su_we)

vMWCpke<-(MatCpke)^(DDe)
MWCpke<-(prod(vMWCpke[1:npc]))^(1/su_we)

vMWCpme<-(MatCpme)^(DDe)
MWCpme<-(prod(vMWCpme[1:npc]))^(1/su_we)

vMWCpmke<-(MatCpmke)^(DDe)
MWCpmke<-(prod(vMWCpmke[1:npc]))^(1/su_we)

outList = list ("Wang(2005) Multivariate Process Capability Indices(PCI) based on PCA",
"number of principal components"=npc,
"MCp"=MWCpe,"MCpk"=MWCpke,"MCpm"=MWCpme,"MCpmk"=MWCpmke)
print(outList)
     }
###############
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

###############Ploting

if(!missing(graphic)) { 
Xmv<-colMeans(x)
 if(length(Xmv)==2 & graphic==TRUE ){
SS<-cov(x) # covariance matrix
 
 ####
hi<-(USL[1]-LSL[1])/2
lo<-(USL[2]-LSL[2])/2
Xm<-colMeans(rbind(LSL,USL))
n=201

    d2 <- (hi-lo)*(hi+lo)                 
    ang <- 2 * pi * seq(0,1, len = 201)
    r <- lo*hi / sqrt(lo^2 + d2 * sin(ang)^2)
    pro <- r * cbind(cos(ang), sin(ang))
    al <- alpha * pi/180
    xx<-pro %*% rbind(c(cos(al), sin(al)), c(-sin(al), cos(al))) + cbind(rep(Xm[1],n),rep(Xm[2],n))
	plot(xx,xlab="Var 1",ylab="Var 2", lty=2, type="n") ; polygon(xx,border=4)
 #####
 
lines(ellipse(SS,centre=Xmv,level=0.9973),type="l")
points(x)
rect( xleft<-LSL[1],xright<-LSL[2],yleft<-USL[1],yright<-USL[2], border=2)
rect( xleft<-LPL[1],xright<-LPL[2],yleft<-UPL[1],yright<-UPL[2], border=3) 
legend(locator(1),cex=0.7,c("Process Region","Tolerance Region","Modified Process Region", "Modified Tolerance Region"),lty=c(1,1,1,1),col=c(1,2,3,4))

 print("CpM index of Shahriari et al. (1995) is the ratio of the Tolerance Region and the Modified Process Region")
 print("MCpm index of Taam et al. (1993) is the ratio of the ellipsoids: Modified Tolerance Region and the Process Region ")
 
 }
}	 
	 
}
