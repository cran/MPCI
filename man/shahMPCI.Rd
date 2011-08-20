\name{shahMPCI}
\alias{shahMPCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Shahriari et al. (1995) Multivariate Capability Vector.
}
\description{
Perform the Multivariate Process Capability Vector according to Shahriari et al. (1995)
}
\usage{
shahMPCI(x,LSL,USL,Target,alpha)

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
The function return a vector of three componenets (CpM; PV; LI). 
The First (CpM) is a ratio of areas or volumes. A value of CpM higher than one indicate that the modified process region is smaller than the engineering tolerance region. 
The second component of the vector is defined as the significance level of a Hotelling T2 statistic computed under the assumption that the center of the engineering specifications is considered to be the true underlying mean of the process.
Values of PV close to zero indicate that the center of the process is far from the engineering target value.
The third component of the vector summarizes a comparison of the location of the modified process region and the tolerance region. It indicates whether any part of the modified process region falls outside the engineering specifications.

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Shahriari H, Hubele N, Lawrence F (1995). "A multivariate process capability vector." Proceedings of the 4th Industrial Engineering Research Conference, 1, 304-309.
}
\author{
Edgar Santos Fernandez, Michele Scagliarini.

}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
taamMPCI
}
\examples{
alpha <- 0.0027
Target <- c(2.2, 304.8, 304.8) 
LSL <- c(2.1, 304.5, 304.5)
USL <- c(2.3, 305.1, 305.1)
x <- matrix(c(2.196,2.184,2.135,2.140,2.119,2.163,2.145,2.209,2.227,2.277
,304.728,304.704,304.713,304.721,304.724,304.670,304.699,304.791,304.737
,304.859,304.798,304.746,304.680,304.719,304.767,304.792,304.753,304.816
,304.754,304.822),nrow = 10) # Matrix with three columns (quality 
# characteristics) and ten rows (number observations) 

#Computing the Shahriari et al. (1995) Multivariate Capability Vector
shahMPCI(x,LSL,USL,Target,alpha)

#Working with the CpM Index
q <- shahMPCI(x,LSL,USL,Target,alpha)$CpM; print(q)

#Computing the index using only x and the specification limits (LSL and USL). 
#The function estimate the Target and assume alpha = 0.0027
rm(Target);rm(alpha) #Removing the parameters Target and alpha 
                     #previously defined.
shahMPCI(x,LSL,USL)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

