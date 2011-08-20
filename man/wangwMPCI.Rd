\name{wangwMPCI}
\alias{wangwMPCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wang (2005) Multivariate Process Capability Indices.
}
\description{
Perform the Multivariate Process Capability Indices based on Principal Component Analysis (PCA)
according to Wang (2005).
}
\usage{
wangwMPCI(x,LSL,USL,Target,npc,alpha,Method,...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{  
  \item{x}{
matrix of the quality characteristics.
}
  \item{LSL}{
is the vector of the Lower Specification Limits
}
  \item{USL}{
is the vector of the Upper Specification Limits
}
  \item{Target}{
is the vector of the target of the process
}
  \item{npc}{
is the number of Principal Components
}
  \item{alpha}{
the significance level
}
  \item{Method}{
is the Method used to decide how many Principal Components
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
When the parameter Target is not specified, then is estimated of centered way as Target=LSL+(USL-LSL)/2.
A priori, the user can choose the number of Principal Components using: summary(princomp(x))

In order to work with the number of Principal Components (npc), the user can:

Specify the value of npc.

Else, choose between five methods described below (without introduce a value of npc).

If npc and Method is not specified the function works with Method 1 (default Method=1).

On the other hand, if npc and Method are both specified the function computes the MPCI using the 
number of principal components. 


Can be used one of the five following methods (modified of Rencher (2002)):

Method 1 or Percentage: guarantee at least the 80 percent of Cumulative Proportion of explained variance.

Method 2 or Average: works with the principal components whose eigenvalues are greater than the average of the
eigenvalues.

Method 3 or Scree: using the Scree Graph (plot Eigenvalue number vs. Eigenvalue size) user can choose the 
principal components to retain.

Method 4 or Bartlett.test: implement the Bartlett's test to ignore the principal components not significantly 
different to the rest.

Method 5 or Anderson.test: compute the Anderson's test to ignore the principal components not significantly 
different to the each other.   

}
\value{
The function return the followings indices: Mcp, MCpk, MCpm and MCpmk.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Anderson, T W (1963). "Asymptotic theory for principal component analysis". 
The Annals of Mathematical Statistics, 34, 122-148.

Jackson, J E (1991). "A user's guide to principal components".John Wiley and Sons.

Rencher,A C (2002) "Methods of Multivariate Analysis". Second edition. John Wiley and Sons.

Wang, C (2005) "Constructing multivariate process capability indices for short-run production." 
The International Journal of Advanced Manufacturing Technology, 26, 1306-1311.
}
\author{
Michele Scagliarini, Edgar Santos Fernandez
}
\note{
The results can be different according to the Methods to compute the npc.
The value of alpha used in Indices based on PCA (usually 0.05) is different of the value 
in shahMPCI and taamMPCI (conventionally = 0.0027).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
wangMPCI, xekeMPCI
For other alternatives to select npc, see nFactors package or Jolliffe,I.T. (2002).
Principal Component Analysis. Second edition. Springer Series in Statistics.
}
\examples{
######Example 1########

data(dataset1)
x<-dataset1
Target<-c(30, 70, 15, 12, 120) 
LSL<-c(24, 60, 10, 8, 100)
USL<-c(36, 80, 20, 16, 140)
alpha<-0.05

# Computing the Wang (2005) Multivariate Process Capability Indices 
# using the Anderson's test to select how many principal components.
wangwMPCI(x,LSL,USL,
         Target,
		 alpha=alpha,
		 Method=5)


#Computing the Wang (2005) using the Method 2 (Average)
wangMPCI(x,LSL,USL,
         Target,
		 alpha=alpha,
		 Method=2)

######Example 2########

Target<-c(2.2, 304.8, 304.8) 
LSL<-c(2.1, 304.5, 304.5)
USL<-c(2.3, 305.1, 305.1)
x<-matrix(c(2.196,2.184,2.135,2.140,2.119,2.163,2.145,2.209,2.227,2.277
,304.728,304.704,304.713,304.721,304.724,304.670,304.699,304.791,304.737
,304.859,304.798,304.746,304.680,304.719,304.767,304.792,304.753,304.816
,304.754,304.822),nrow=10)

#Computing the Wang (2005) Multivariate Process Capability Indices
wangwMPCI(x,LSL,USL,
          Target,
		  npc=2)

#Working with the MCpm index
a<-wangwMPCI(x,LSL,USL,
         Target,
		 npc=2)$MCpm; print(a)

#Computing the index using only x and the specification limits (LSL and USL).
rm(Target);rm(npc) #Removing the parameters Target and npc previously defined.
wangwMPCI(x,LSL,USL)
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

