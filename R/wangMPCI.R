wangMPCI <-
function(x,LSL,USL,Target,npc,alpha,Method,...){

SL<-cbind(LSL,USL) # matrix of specifications
TOL<-USL-LSL # difference in specifications
m<-length(LSL) #dimensions

SS<-cov(x) # covariance matrix
n<-nrow(x)
Xmv<-colMeans(x) #mean vector

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
         while (DDeV[3,j]<0.8) {
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

outList = list ("Wang & Chen(1998) Multivariate Process Capability Indices(PCI) based on PCA",
"number of principal components"=npc,"MCp"=MCpe,"MCpk"=MCpke,"MCpm"=MCpme,"MCpmk"=MCpmke)
print(outList)

}

