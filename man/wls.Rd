\name{wls}
\alias{wls}
\docType{data}
\title{
WLS data on Catholic schooling, Catholic religion, and access to a Catholic school.
}
\description{
Does Catholic schooling causally effect outcomes, e.g. wages, later in the life?
This data is extracted from a longitudinal study of high school graduates
in 1957 in Wisconsin.}
\usage{data("wls")}
\format{
  A data frame with 4469 observations on the following 24 variables.
  \describe{
    \item{\code{idpub}}{Public id from wls.}
    \item{\code{gwiiq_bm}}{Student's IQ}
    \item{\code{edfa57q}}{Father's education}
    \item{\code{edmo57q}}{Mother's education}
    \item{\code{bmpin1}}{Parent's income}
    \item{\code{ocsf57}}{occupation score}
    \item{\code{ocpf57}}{occupational prestige score}
    \item{\code{incg400}}{Parent's income in the top 1\%}
    \item{\code{incg250}}{Parent's income in the top 5\%}
    \item{\code{relfml}}{Catholic religion.  TRUE = Catholic; FALSE = Non-Catholic}
    \item{\code{res57}}{A numeric code of area of residence of the graduate.}
    \item{\code{hsdm57}}{TRUE = Catholic schooling; FALSE = Public schooling.}
    \item{\code{edfa57q.miss}}{logical vector indicating whether edfa57q is missing.}
	\item{\code{edmo57q.miss}}{logical vectors indicating whether edmo57q is missing.}
	\item{\code{bmpin1.miss}}{logical vectors indicating whether bmpin1 is missing.}
	\item{\code{ocsf57.miss}}{logical vectors indicating whether ocsf57 is missing.}
	\item{\code{ocpf57.miss}}{ logical vectors indicating whether ocpf57 is missing.}
	\item{\code{res57Dic}}{Urban or Rural. TRUE = Urban; FALSE = Rural.}
	\item{\code{edfa57q.NoNA}}{Father's education with NAs replaced by the average.}
    \item{\code{edmo57q.NoNA}}{Mother's education with NAs replaced by the average.}
    \item{\code{bmpin1.NoNA}}{Parent's income with NAs replaced by the average.}
    \item{\code{ocsf57.NoNA}}{Occupation score with NAs replaced by the average.}
    \item{\code{ocpf57.NoNA}}{Occupational prestige score with NAs replaced by the average.}
    \item{\code{yrwg74}}{1974 Graduate's wages and salaries in \$100's last year.}
  }
}
\details{
This data set is constructed out of publicly available variables of Wisconsin
Longitudinal Study.  Only male graduates are considered. For more detailed 
description of the variables see https:// www.ssc.wisc.edu/ wlsresearch/ documentation/.}
\source{
Wisconsin Longitudinal Study, https://www.ssc.wisc.edu/wlsresearch/
}
\references{
Herd, Pamela, Deborah Carr, and Carol Roan. 2014. "Cohort Profile: Wisconsin 
Longitudinal Study (WLS)." International Journal of Epidemiology 43:34--41 
PMCID: PMC3937969}
\examples{
data(wls)
}
\keyword{datasets}
