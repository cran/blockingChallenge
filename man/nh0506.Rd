\name{nh0506}
\alias{nh0506}
\docType{data}
\title{
Smoking and homocysteine levels in NHANES 2005-2006.
}
\description{
This is an example data set for \code{makeblocks}. For analysis of this 
data see Pimentel et al. (2016), Rosenbaum (2018).
}
\usage{data("nh0506")}
\format{
  A data frame with 2475 observations on the following 19 variables.
  \describe{
    \item{\code{SEQN}}{NHANES id number}
    \item{\code{female}}{1 if female, 0 if male}
    \item{\code{age}}{age in years, >=20}
    \item{\code{black}}{1 if black, 0 otherwise}
    \item{\code{hispanic}}{1 if hispanic, 0 otherwise}
    \item{\code{education}}{Education}
    \item{\code{povertyr}}{Ratio of family income to the poverty level, capped at 5x}
    \item{\code{creactiveprotein}}{creactive protein}
    \item{\code{homocysteine}}{homocysteine}
    \item{\code{cotinine}}{cotinine in blood}
    \item{\code{cadmium}}{cadmium in blood}
    \item{\code{lead}}{lead in blood}
    \item{\code{bmi}}{Body mass index}
    \item{\code{cigs100life}}{1 if smoked more than 100 cigarettes in lifetime, 0 otherwise}
    \item{\code{smokenow}}{1 if smokes now, 0 otherwise}
    \item{\code{cigsdays30}}{Days smoked in last 30 days: 0 if never smoker, 30 if daily smoker}
    \item{\code{cigsperday30}}{Cigarettes smoked per day in last 30 days}
    \item{\code{tobacco5days}}{1 = used tobacco in the last 30 days, 0 otherwise}
    \item{\code{dailysmoker}}{1 = daily smoker, 0 = never smoker}
  }
}
\details{
Data from NHANES 2005-2006 concerning homocysteine levels in daily smokers and never
smokers, aged 20 and older. Daily smokers smoked every day for the last 30 days, smoking
an average of at least 10 cigarettes per day. Never smokers smoked fewer than 100 cigarettes in
their lives, do not smoke now, and had no tobacco use in the previous 5 days.}
\source{
NHANES, the US National Health and Nutrition Examination Survey, 2005-2006.
}
\references{
Pimentel, S. D., Small, D. S. and Rosenbaum, P. R. (2016). Constructed second control groups and
attenuation of unmeasured biases. Journal of the American Statistical Association, 111, 1157-1167.

Rosenbaum, P. R. (2018). Sensitivity analysis for stratified comparisons in an 
observational study of the effect of smoking on homocysteine levels. Annals of Applied Statistics,
to appear.
}
\examples{
data(nh0506)
}
\keyword{datasets}
