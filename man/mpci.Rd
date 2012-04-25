\name{mpci}
\alias{mpci}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate Process Capability Index
}
\description{
Performs capability indices for multivariate process in which the variables are correlated.

Present the most accepted indices such as: 

Shahriari et al. (1995) Multivariate Capability Vector.

Taam et al. (1993) Multivariate Capability Index (MCpm).

and the followings based on Principal Component Analysis (PCA):

Wang and Chen (1998) Multivariate Process Capability Indices (MPCI),

Xekalaki and Perakis (2002) MPCI.

Wang (2005) MPCI

}
\usage{

 mpci(index,x,LSL,USL,Target,npc,alpha,Method,perc,graphic,...)

}

\arguments{
  \item{index}{
is the name of the type of index to use. For instance: index="shah" or index="taam"
}
\item{x}{
matrix of the quality characteristics.
}
  \item{LSL}{
is the vector of the Lower Specification Limits.
}
  \item{USL}{
is the vector of the Upper Specification Limits.
}
  \item{Target}{
is the vector of the target of the process.
}
  \item{npc}{
is the number of principal component to use. Only for indices based on PCA.
}
  \item{alpha}{
for the Shahriari et al. (1995) and Taam et al. (1993) indices, alpha is the proportion of nonconforming products (conventionally = 0.0027)

In the case of the indices based on Principal Component Analysis (PCA) (Wang and Chen (1998), Xekalaki and Perakis (2002), Wang (2005)), alpha is the the significance level.
}
  \item{Method}{
is the Method used to decide how many Principal Components (works only with PCA indices).
}
  \item{perc}{
is the percent of Cumulative Proportion of explained variance specified by the user in Method 1 (only for PCA indices).
}
  \item{graphic}{
allows in bivariate case the graphical representation. 
}
  \item{...}{
others parameters
}
}
\details{
When the parameter Target is not specified, then is estimated of centered way as Target=LSL+(USL-LSL)/2. 

On the other hand, if alpha is missing the function assume alpha=0.0027 for the first two indices 

and 0.05 for the others.

A priori, the user can choose the number of Principal Components using: summary(princomp(x))


In order to work with the number of Principal Components (npc), the user can:

Specify the value of npc.

Else, choose between five methods described below (without introduce a value of npc).

If npc and Method is not specified the function works with Method 1 (default Method=1).

In this case, the user could specify other desired percentage different of the 80%(e.g.: perc=0.90).

On the other hand, if npc and Method are both specified the function computes the MPCI using the 
number of principal components. 


Can be used one of the five following methods (modified of Rencher (2002)):

Method 1 or Percentage: guarantees at least the 80 percent of Cumulative Proportion of explained variance.
(or other percent specified by the user with the argument perc) 

Method 2 or Average: works with the principal components whose eigenvalues are greater than the average of the
eigenvalues.

Method 3 or Scree: using the Scree Graph (plot Eigenvalue number vs. Eigenvalue size) user can choose the 
principal components to retain.

Method 4 or Bartlett.test: implements the Bartlett's test to ignore the principal components not significantly 
different to the rest.

Method 5 or Anderson.test: computes the Anderson's test to ignore the principal components not significantly 
different to the each other.   


}
\value{
Using the index="shah" the function return a vector of three componenets (CpM; PV; LI). 
The First (CpM) is a ratio of areas or volumes. A value of CpM higher than one indicate that the modified process region is smaller than the engineering tolerance region. 
The second component of the vector is defined as the significance level of a Hotelling T2 statistic computed under the assumption that the center of the engineering specifications is considered to be the true underlying mean of the process.
Values of PV close to zero indicate that the center of the process is far from the engineering target value.
The third component of the vector summarizes a comparison of the location of the modified process region and the tolerance region. It indicates whether any part of the modified process region falls outside the engineering specifications.

When index="taam" the function return the value of MCpm. When the process mean vector equals the target vector, and the index has the value 1, then 99.73 percent of the process values lie within the modified tolerance region.

For values of index="wang", "xeke" or "wangw" function return the followings indices: Mcp, MCpk, MCpm and MCpmk.
When all indices are greater than 1, the process is capable.

}
\references{
Anderson, T W (1963). "Asymptotic theory for principal component analysis". 
The Annals of Mathematical Statistics, 34, 122-148.

Chen H (1994). "A multivariate process Capability index over a rectangular solid zone." Statistica Sinica, 4, 749-758.

Hubele N, Shahriari H, Cheng C (1991). "A bivariate process Capability vector". In JB. Keats and DC.Montgomery (Eds.) "Statistical Process Control in Manufacturing".Marcel Dekker,New York.

Jackson, J E (1991). "A user's guide to principal components". John Wiley & Sons. New York.

Rencher,A C (2002) "Methods of Multivariate Analysis". 2nd Ed. John Wiley & Sons. New York.

Santos-Fernandez E, Scagliarini M (2012). "MPCI: An R Package for Computing Multivariate Process Capability Indices". Journal of Statistical Software, 47(7), 1-15, URL http://www.jstatsoft.org/v47/i07/.

Scagliarini M (2011). "Multivariate process capability using principal component analysis in the presence of measurement errors." AStA Advances in Statistical Analysis, 95, 346-357.

Shahriari H, Hubele N, Lawrence F (1995). "A multivariate process capability vector." Proceedings of the 4th Industrial Engineering Research Conference, 1, 304-309.

Taam W, Subbaiah P, Liddy W (1993). "A note on multivariate capability indices." Journal of Applied Statistics, 20, 339-351.

Wang, C (2005) "Constructing multivariate process capability indices for short-run production." The International Journal of Advanced Manufacturing Technology, 26, 1306-1311.

Wang F, Chen J (1998). "Capability index using principal components analysis." Quality Engineering, 11, 21-27.

Wang F, Hubele N, Lawrence F, Miskulin J, Shahriari H (2000). "Comparison of three multivariate process capability indices." Journal of Quality Technology, 32, 263-275.

Xekalaki E, Perakis M (2002). "The Use of principal component analysis in the assessment of process capability indices." 
Proceedings of the Joint Statistical Meetings of the American Statistical Association, The Institute of Mathematical Statistics, The Canadian Statistical Society. New York.

}
\author{
Edgar Santos-Fernandez, Michele Scagliarini.
}
\note{
The results can be different according to the Methods to compute the npc.
The value of alpha used in Indices based on PCA (usually 0.05) is different of the value 
in index="shah" or "taam" (conventionally = 0.0027).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For other alternatives to select npc, see nFactors package or Jolliffe,I.T. (2002).
Principal Component Analysis.Second edition. Springer Series in Statistics.
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

mpci(index="shah",x,LSL,USL,Target,alpha)

#Computing the Shahriari et al. (1995) index using only x and the
#specification limits (LSL and USL). 
#The function estimate the Target and assume alpha = 0.0027
rm(Target);rm(alpha) #Removing the parameters Target and alpha 
                     #previously defined.
mpci(index="shah",x,LSL,USL)

#Computing the Taam el al. (1993) Multivariate Capability Vector
alpha <- 0.0027
Target <- c(2.2, 304.8, 304.8) 
a <- mpci(index="taam",x,LSL,USL,Target,alpha)



# Using the indices based on based on PCA

data("dataset1")
x<-dataset1
Target <- c(30, 70, 15, 12, 120) 
LSL <- c(24, 60, 10, 8, 100)
USL <- c(36, 80, 20, 16, 140)
alpha <- 0.05

# Computing the Wang and Chen (1998) Multivariate Process Capability Indices (MPCI)
# using the Bartlett's test to select how many principal components.
mpci(index="wang",x,LSL,USL,
         Target,
		 alpha = alpha,
		 Method = 4)

# Wang and Chen (1998) index using the default method (Method 1 or
# Percentage) and value of alpha
mpci(index="wang",x,LSL,USL)

# Xekalaki and Perakis (2002) Multivariate Process Capability Indices (MPCI)
# using the Bartlett's test to select how many principal components.
mpci(index="xeke",x,LSL,USL,
         Target,
		 alpha = alpha,
		 Method = 4)

# Xekalaki and Perakis (2002) using the default method (Method 1 or
# Percentage) and value of alpha
mpci(index="xeke",x,LSL,USL,Target)

# Wang (2005) Multivariate Process Capability Indices 
# using the Anderson's test to select how many principal components.
mpci(index="wangw",x,LSL,USL,
         Target,
		 alpha = alpha,
		 Method = 5)

#Computing the Wang (2005) using the Method 2 (Average)
mpci(index="wangw",x,LSL,USL,
         Target,
		 alpha = alpha,
		 Method = 2)


}

\keyword{mpci}

