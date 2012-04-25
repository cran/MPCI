\name{MPCI-package}
\alias{MPCI-package}
\alias{MPCI}
\docType{package}
\title{
Multivariate Process Capability Indices (MPCI)
}
\description{
The aim of the MPCI package is the approach to the proposed Multivariate Process Capability Indices (MPCI).
Perform the followings indices:

Shahriari et al. (1995) Multivariate Capability Vector.

Taam et al. (1993) Multivariate Capability Index (MCpm).

and the followings based on Principal Component Analysis (PCA):

Wang and Chen (1998) Multivariate Process Capability Indices (MPCI),

Xekalaki and Perakis (2002) MPCI.

Wang (2005) MPCI

}
\details{
\tabular{ll}{
Package: \tab MPCI\cr
Type: \tab Package\cr
Version: \tab 1.0.4\cr
Date: \tab 2011-04-25\cr
License: GPL(>=2)\cr
LazyLoad: yes
}
}
\author{
Edgar Santos-Fernandez, Michele Scagliarini.

Maintainer: Edgar Santos-Fernandez <edgar.santos@etecsa.cu>
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
\keyword{ package }
\seealso{ MEWMA and MCUSUM package}
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

#Performing the Taam et al. (1993) Multivariate Capability Index (MCpm)
mpci(index="taam",x,LSL,USL,Target,alpha)

#Computing the Wang and Chen (1998) Multivariate Process Capability Indices(MPCI)
alpha = 0.05
mpci(index="wang", x,LSL,USL,Target)

#Performing the Xekalaki and Perakis (2002) MPCI.
mpci(index="xeke",x,LSL,USL,
          Target,
		  Method = 2)

#Computing the Wang (2005) MPCI
mpci(index="wangw",x,LSL,USL,
          Target,
		  Method = 1,
		  alpha = 0.05)

}
