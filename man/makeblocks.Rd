\name{makeblocks}
\alias{makeblocks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates strata of size \code{bsize} aiming to minimize within strata distance.
}
\description{
From a distance matrix of \code{n} points (or from a data set, in 
which case the distances are calculated using the \code{smahal} function) and
block size, \code{bsize}, this function outputs a stratification of the points
of blocks of size \code{bsize}.  
}
\usage{
makeblocks(distmat, bsize, Itr = 30, .data = NULL, vars, maxItr = 200, 
					verbose = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distmat}{
A numeric distance matrix of size \code{n} by \code{n}. Can be created using 
\code{smahal} function.
}
  \item{bsize}{
The size of each block. 
}
  \item{Itr}{
Number of iterations of the randomization algorithm to be run, (initially).
See 'Usage Tricks' section below.
}
  \item{.data}{
An optional data frame. Recommended.
}
  \item{vars}{
A character vector of column names/indices to be used to create the distance 
matrix. Ignored, if \code{distmat} is provided. If not provided all variables of .data are used.
}
  \item{maxItr}{
Maximum number of iterations of the randomization algorithm to be run, (initially).
See 'Usage Tricks' section below.
}
  \item{verbose}{
0, 1 or 2 for the amount to detail to be printed during the run of the function.
}
  \item{\dots}{
Additional arguments. May be ignored.
}
}
\details{
At the minimum the function requires \code{distmat} and \code{bsize}. In 
this case it will return a blocking structure as a \code{floor(n/bsize)}
by \code{bsize} matrix.  Each entry of the matrix will be a unit assigned 
to that block.

The argument .data of the data indices is not required but recommended for 
convenience. If .data is provided but not 'distmat', the distance matrix will be 
created by a call to \code{smahal} function to \code{.data[,vars]}.  If both
\code{distmat} and \code{.data} are provided then \code{distmat} is used for
blocking and the strata index, from 1 to \code{floor(n/bsize)}, for the data
points are added as an extra column to the data set. This makes the 
diagnostics process simpler, see examples below.

IMPORTANT NOTE: In order to perform matching, \code{makeblocks} requires the 
user to load the optmatch (>= 0.9-1) package separately. The manual loading is 
required due to software license issues. If the package is not loaded the 
\code{makeblocks} command will fail with an error saying the 'optmatch' package 
is not present. Reference to 'optmatch' is given below.
}

\section{Usage Tricks}{
The function depends on an iterative algorithm. The iterative algorithm 
is usually pretty smart and can figure out the number of iterations required 
to get close to the optimal solution with high probability. 

To make the function make all decisions
on its own, set \code{maxItr} to a  very high value, or Inf, and do not specify \code{Itr}.
If the problem is easy, it should automatically stop much earlier than \code{maxItr}.

For a more interactive use of the function one can start with a relatively 
small value of \code{Itr}, e.g., 20 or 30, and then it should prompt how many more 
iterations would be needed to maximize get close to the optimal solution. The
user can choose how many more iterations to run. Some run log will be printed 
to further help the user make this decision, for more detailed logs use 
\code{verbose=1} or \code{verbose=2}.

To preset the number of iterations, say 200, provide \code{Itr=200} and 
\code{maxItr=200}.

Since the algorithm is a randomization algorithm, the outcomes will not be 
reproducible unless the seed is set, see \code{set.seed}. 
}

\section{Challenge}{
The pdf version of this manual shows the math equations clearly; the 
html version may be difficult to read.

The problem of finding an optimal stratification is NP-hard, in general. This 
algorithm combines multiple heuristics to solve the following problem.

Let d is a \code{n} by \code{n} matrix of distances, i.e. \code{d[ij] >= 0}, 
\code{d[ij]=d[ji]} and \code{d[ij] <= d[ik]+d[kj]}, for all \code{i, j, k}. Fix
integer \code{b}.  Find \code{B[1], \ldots, B[floor(n/b)]}, \code{floor(n/b)}
exclusive subsets of \code{\{1,\ldots,n\}} of size \code{b} such that it solves 
the following optimization problem. 
\deqn{\min_{B[1], \ldots, B[\lfloor(n/b)\rfloor]; B[l]\cap B[l']=\emptyset; |B[l]|=b; \forall l,l'}  
\sum_{l=1}^{\lfloor(n/b)\rfloor} \sum_{i,j\in B[l]} d[ij]. }{
min_{B[1], \ldots, B[\lfloor(n/b)\rfloor]; B[l]\cap B[l']=\emptyset; |B[l]|=b; for all l,l'}  
\sum_{l=1}^{floor(n/b)} \sum_{i,j\in B[l]} d[ij].
}

The challenge has two versions.

\subsection{(1)}{ Construct another algorithm that beats this algorithm across multiple problems,
decided by \code{n}, \code{bsize}.  Pathological examples can be insightful in 
this challenge.
}
\subsection{(2)}{ Provide theoretical results. If you set \code{maxItr=}\eqn{{n \choose \lfloor(n/b)\rfloor}}
the solution will be a 2-approximation (with a small and fixable caveat). The 
heuristic algorithm believes that it can do not too bad with much less 
computation. How much worse is this heuristic algorithm? What would be a decent 
compromise, the number of iterations \code{Itr}, to balance computational burden 
and approximation error? Answers to these questions are useful both universally 
and in average cases.
}
}
\value{
A list consisting of the following elements
\item{strata}{ A blocking structure as a \code{floor(n/bsize)}
by \code{bsize} matrix.  Each entry being the units in that block. The unit names
are based on the row names of \code{.data} or the row and column names of 
\code{distmat}, or 1:nrow(distmat).  The rownames of \code{strata} are the 
pivot elements of the blocks.}
\item{.data}{If \code{.data} is provided, then the same data set is returned 
with an extra column 'strata'. \code{.data$strata} are the strata index 
for the corresponding data point.  The data points not assigned to any strata 
is indexed with a negative number.
}
}

\author{
Bikram Karmakar
}

\references{
Hansen, B. B. and Klopfer, S. O. (2006) Optimal full matching and related designs via network flows,
JCGS 15 609--627.

Karmakar, B., Small, D. S. and Rosenbaum, P. R. (2018). 
Reinforced designs in observational studies of treatment effects: Multiple instruments plus controls as evidence factors.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## IN THESE EXAMPLES YOU NEED TO LOAD THE R PACKAGE 
#  'optmatch'.  OTHERWISE, THEY WILL RUN IN A CODE CHECK MODE.
# library(optmatch)
###### 		Usage 1		##########
data(nh0506)
\dontshow{
	nh0506 <- nh0506[1:13,]
}
res <- makeblocks(nh0506, bsize=4, maxItr=15, 
         vars=c("female","age","black","hispanic","education","povertyr"))

		
\dontrun{
###### Some diagnostics
v = "age"
## Boxplots
temp = res$.data
temp[[paste0(v,'.instrata')]] <- temp[[v]] 
meanbystrat <- aggregate(temp[[v]], by =list(strata=temp$strata), mean)
for(b in unique(temp$strata)){
     who <- which( temp$strata == b )
     temp[[paste0(v,'.instrata')]][who] <- temp[[paste0(v,'.instrata')]][who] 
                        - meanbystrat[meanbystrat[,1]==b,2]
}

temp[[paste0(v,'.instrata')]] = temp[[paste0(v,'.instrata')]] + mean(temp[[v]])

boxplot(v~x,
data = data.frame(v=c(temp[[v]],temp[[paste0(v,'.instrata')]][temp$strata>0]), 
               x=c(rep('prestratification', nrow(temp)), 
                        rep('poststratification', sum(temp$strata>0)))))

## Anova
summary(aov(age~factor(strata), data=res$.data[res$.data$strata>0,]))
}


###### 		Usage 2		##########
distmat <- smahal(nh0506[,c("female","age","black","hispanic",
								"education","povertyr")])
res <- makeblocks(distmat, .data=nh0506, bsize=4, maxItr=15)


\dontrun{
# Other usages
## Internally calls 'smahal' to create the distances.
res <- makeblocks(nh0506[,c("female","age","black","hispanic",
                                           "education","povertyr")], 4)

## Returns the blocking structure in a matrix form, if data is not provided.
distmat <- smahal(nh0506[,c("female","age","black",
                                  "hispanic","education","povertyr")])
res <- makeblocks(distmat, bsize=4, maxItr=20)
}


###### 		Usage 3		##########
data(wls)

\dontrun{
data(wls)
library(optmatch)
wls4match <- wls

## This code replicates the blocking algorithm used in the paper
##	Karmakar, Small, and Rosenbaum (2018).

## Create the distance matrix

distmat1 <- smahal(wls4match[,"gwiiq_bm"]) 				## IQ

## Father's and mother's edu and parent's income
distmat2 <- smahal(wls4match[,c("edfa57q.NoNA", "edmo57q.NoNA", "bmpin1.NoNA",	
              ## Father's and mother's edu and parent's income
              "incg400", "incg250")])				
              ## Indicators for income in the top 5% and 1%			
			  
## occupation score
distmat2.2 <- smahal(wls4match[,c("ocsf57.NoNA", "ocpf57.NoNA")])		

## missing indicators
distmat3 <- smahal(wls4match[,c("edfa57q.miss", "edmo57q.miss", 	
              "bmpin1.miss", "ocsf57.miss", "ocpf57.miss")])

## The IQ = gwiiq_bm is given more weight.
##  parents' education and parent's income 
distmat = distmat1*10+6*distmat2+3*distmat2.2+2*distmat3

## creating the blocks.  This can take about 30min to run.
##	May take more time depending of the computation power of the system.
set.seed(0841)
res.20.2 = makeblocks(distmat, bsize=25, Itr=250, maxItr=250, .data=wls4match)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Challenge }% use one of  RShowDoc("KEYWORDS")
\keyword{ Stratification algorithm }% __ONLY ONE__ keyword per line
