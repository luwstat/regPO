\name{po_fit}
\alias{po_fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Regression analysis of left-truncated arbitrarily censored data or arbitrarily censored data under the proportional odds model
}
\description{
Fits the semi-parametric maximum likelihood estimator for left-truncated arbitrarily censored data or just arbitrarily censored data under the proportional odds modelusing EM algorithm.
}
\usage{
po_fit(L,R,truncation = F,C,x,order = 3,equal_space = T,nknot, myknots,conv_cri = 1e-7)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{L}{Lower bound of the observed intervals, a length n vector.}
  \item{R}{Upper bound of the observed intervals, a length n vector.}
  \item{truncation}{Logical, if TRUE the left trunction times for each observation should be provide in argument C. The default if FALSE.}
  \item{C}{Left truncation times of the observations, a length n vector.}
  \item{x}{Covariate matrix, an n x p matrix.}
  \item{order}{The order of I-spline basis function.}
  \item{equal_space}{Logical, if TRUE knots are placed evenly across the observed intervals based on the input data set; if FALSE knots should be specified by the user. see myknots. The default is TRUE.}
  \item{nknot}{The total number of knots.}
  \item{myknots}{Knots specified by the user.}
  \item{conv_cri}{The convergence criterion of the EM algorithm.}
}
\details{
L and R are left and right ends of the censoring interval. This allows for left censored (L == 0), right censored (R == inf), uncensored (L == R) along with interval censored observations (0 < L < R < inf).

The EM algorithm converges when the sum of square of the difference in the parameter estimates is less than conv_cri.
}
\value{
\item{beta}{Estimates of regression coefficients.}
\item{beta_se}{The standard errors of regression coefficients estimate.}
\item{CI}{Confidence intervals for the regression coefficients.}
\item{spline_coef}{Estimates of spline coefficients.}
\item{knots}{The knots for I-spline.}
\item{AIC}{The Akaike information criterion. }
\item{BIC}{The Bayesian information criterion.}
\item{Baseline_Surv}{Plot of the estimated baseline survival curve.}
\item{Baseline_hazard}{Plot of the estimated baseline hazard rate curve.}
\item{Baseline_odds}{Plot of the estimated baseline odds curve.}
}
\references{
Lu Wang & Lianming Wang. Regression analysis of  arbitrarily censored and left-truncated data under the proportional odds model
}
\author{
Wang,L. and Wang,L.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(icenReg)
data("IR_diabetes")
gender <- as.matrix(as.numeric(IR_diabetes[,'gender'] == 'male'))
po_fit(L = IR_diabetes[,'left'],R = IR_diabetes[,'right'],x = gender,order = 3,nknot = 10)
  }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
