\name{taamMPCI}
\alias{taamMPCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Taam et al. (1993) Multivariate Capability Index (MCpm)
}
\description{
Perform the Multivariate Process Capability Index according to Taam et al. (1993)
}
\usage{
taamMPCI(x,LSL,USL,Target,alpha)

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
 \item{alpha}{
is the proportion of nonconforming products (conventionally = 0.0027)
}
}
\details{
When the parameter Target is not specified, then is estimated of centered way as Target=LSL+(USL-LSL)/2. 
On the other hand, if alpha is missing the function assume alpha=0.0027
}
\value{
Return the value of MCpm. When the process mean vector equals the target vector, and the index has the value 1, then 99.73 percent of the process values lie within the modified tolerance region.

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Taam W, Subbaiah P, Liddy W (1993). "A note on multivariate capability indices." Journal of Applied Statistics, 20, 339-351.
}
\author{
Edgar Santos Fernandez, Michele Scagliarini.

}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
shahMPCI
}
\examples{
alpha<-0.0027
Target<-c(2.2, 304.8, 304.8) 
LSL<-c(2.1, 304.5, 304.5)
USL<-c(2.3, 305.1, 305.1)
x<-matrix(c(2.196,2.184,2.135,2.140,2.119,2.163,2.145,2.209,2.227,2.277
,304.728,304.704,304.713,304.721,304.724,304.670,304.699,304.791,304.737
,304.859,304.798,304.746,304.680,304.719,304.767,304.792,304.753,304.816
,304.754,304.822),nrow=10)

#Computing the Taam el al. (1993) Multivariate Capability Vector
a<-taamMPCI(x,LSL,USL,Target,alpha)

#Computing the index using only x and the specification limits (LSL and USL). 
#The function estimate the Target and assume alpha=0.0027
rm(Target);rm(alpha) #Removing the parameters Target and alpha 
                     #previously defined.
taamMPCI(x,LSL,USL)
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

