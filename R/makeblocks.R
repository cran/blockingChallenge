makeblocks <-
function(distmat, bsize, Itr=30, .data=NULL, vars, maxItr = 200, verbose=0, ...){
#if (!requireNamespace("optmatch", quietly = TRUE)) {
#stop("Error: package optmatch (>= 0.9-1) not loaded. You must install optmatch first and agree to the terms of its license.")
#}

stopifnot(length(bsize)==1)
stopifnot(bsize==round(bsize))
stopifnot(bsize>=3)

if(missing(distmat)){
stopifnot(!is.null(.data))
stopifnot(is.data.frame(.data)|is.matrix(.data))
stopifnot((dim(.data)[1])>(2*bsize))
if(is.null(rownames(.data)))
rownames(.data) = 1:nrow(.data)
if(!missing(vars)){
stopifnot(all(vars %in% colnames(.data)))

distmat = smahal(.data[,vars])
rownames(distmat) = rownames(.data)
colnames(distmat) = rownames(.data)
} 
if(missing(vars)){
distmat = smahal(.data)
rownames(distmat) = rownames(.data)
colnames(distmat) = rownames(.data)
}
}

if(!missing(distmat)){
stopifnot(is.data.frame(distmat)|is.matrix(distmat))
if(dim(distmat)[1]!=dim(distmat)[2]){
if(!is.null(.data)){
.data = .data
} else { .data = distmat }

if(is.null(rownames(.data)))
rownames(.data) = 1:nrow(.data)
if(!missing(vars)){
stopifnot(all(vars %in% colnames(.data)))

distmat = smahal(.data[,vars])
rownames(distmat) = rownames(.data)
colnames(distmat) = rownames(.data)
}
if(missing(vars)){
distmat = smahal(.data)
rownames(distmat) = rownames(.data)
colnames(distmat) = rownames(.data)
}
}
if(is.null(rownames(distmat))) rownames(distmat) = 1:nrow(distmat)
stopifnot((dim(distmat)[1])>(2*bsize))
}

nobs<-dim(distmat)[1] #number of observations
nblocks<-floor(nobs/bsize) #number of blocks
units <- rownames(distmat)

if(!missing(maxItr) & !missing(Itr)) {
if(maxItr==Itr){
 Itr4est = -1
} else if(Itr > 30){
 Itr4est = 30
} else Itr4est = Itr

runallItr = TRUE
} else { 
if(missing(maxItr)){
runallItr = FALSE

if(missing(Itr)){
Itr = 30
Itr4est <- Itr
 } else if(Itr > 30){
Itr4est = 30
 } else Itr4est = Itr
} else {
runallItr = TRUE

Itr = min(maxItr, 30)
Itr4est <- Itr
}
}

prev.strat = NA
prev.strat.dist = 10*bsize*nblocks*max(distmat)

strat.dist <- c()


itr <- 1
Itr_est <- -1
while(1){
#for(itr in 1:Itr){
cat("Iteration no.", itr,":: ")
# Randomly select the first unit in each block
if(verbose>0) cat("\n\tSelecting random starts: ")
 
first <- sample(units, nblocks, replace=FALSE)

current.strat <- matrix(NA, nblocks, bsize)
rownames(current.strat) <- first

distmat.temp1 <- distmat[first, setdiff(units, first)]

if(eval(parse(text="requireNamespace('optmatch', quietly = TRUE)"))){
om <- eval(parse(text="suppressWarnings(pairmatch(distmat.temp1, controls=bsize-1))"))
} else {
if(itr==1) { 
cat("\n")
cat("Error: package optmatch (>= 0.9-1) not loaded. You must install optmatch first and agree to the terms of its license.")
cat("\nRunning code check mode.")
cat("\n")
}

om <- paste0("1.", c(sample(1:length(first)), sample(rep(1:length(first), bsize-1))))
om <- c(om, rep(NA, (length(units)%%bsize)))
om = factor(om)
names(om) = c(first, colnames(distmat.temp1))
}

invisible(gc())
stopifnot(all(units %in% names(om)))
stopifnot(length(levels(om)) == nblocks)

rm(distmat.temp1)
for(b in rownames(current.strat))
current.strat[b,] = na.omit(names(om)[om==om[b]])

sinks <- names(om)[is.na(om)]

current.strat.dist = sum(sapply(first, function(b) sum(distmat[b, current.strat[b,]], na.rm=TRUE)))

if(verbose>0) cat("\n\t First match created.  Average within strata distance = ", current.strat.dist/nblocks,"\n")
if(verbose>0) cat("\tTrying local optimization.")

ntry = 0
nrogue = 1

while(1){
##identify the furthest block unit in each block
first.nxt.itr = first
rogue.unit = c()
rogue.unitid = list()
rogue.dist = 0

 ntry = ntry+1
if(verbose>1) cat("\n\t\t ",ntry, ": ( nrogue=",nrogue,")")
for(b in rownames(current.strat)){
unitid = order(distmat[b, current.strat[b,]], decreasing=TRUE)[1:nrogue]#which.max(distmat[b, current.strat[b,]])[1]
names(unitid) = current.strat[b,unitid]
if(distmat[b, tail(names(unitid),1)] == 0){
first.nxt.itr = setdiff(first.nxt.itr, b)
} else {
rogue.unit = c(rogue.unit, names(unitid))
rogue.dist = rogue.dist + sum(distmat[b, names(unitid)])
rogue.unitid[[b]] = unitid
}
}

if(length(first.nxt.itr)==0) break;#########

rogue.unit = c(rogue.unit, sinks)# add the unmatched units in the last try to the rogue ones 

#distmat.temp2 <- distmat[first.nxt.itr, rogue.unit]
distmat.temp2 <- distmat[first.nxt.itr, rogue.unit]
for(b in first.nxt.itr)
distmat.temp2[b,rogue.unit] = colMeans(distmat[current.strat[b,-rogue.unitid[[b]]], rogue.unit])


if(eval(parse(text="requireNamespace('optmatch', quietly = TRUE)"))){
om2 <- eval(parse(text="suppressWarnings(pairmatch(distmat.temp2, controls=nrogue))"))

} else { 
om2 <- paste0("1.", c(sample(1:length(first.nxt.itr)), sample(rep(1:length(first.nxt.itr), nrogue))))
om2 <- c(om2, rep(NA, (sum(dim(distmat.temp2))-length(first.nxt.itr)*(1+nrogue)) ))
om2 = factor(om2)
names(om2) = c(first.nxt.itr, colnames(distmat.temp2))
}

invisible(gc())
stopifnot(names(om2)[1:length(first.nxt.itr)] == first.nxt.itr)
stopifnot(length(levels(om2)) == length(first.nxt.itr))
rm(distmat.temp2)

new.dist = 0
matched.unit = list()
for(b in first.nxt.itr){
matched.unit[[b]] = names(na.omit(om2[om2==om2[b]]))[-1]
new.dist = new.dist + sum(distmat[b, matched.unit[[b]]])
}

if(verbose>0) cat(" ", new.dist/nblocks, rogue.dist/nblocks)
if(new.dist < rogue.dist){
cat(" . ")
for(b in first.nxt.itr)
current.strat[b, rogue.unitid[[b]]] = matched.unit[[b]]
sinks = names(om2[is.na(om2)])
current.strat.dist = current.strat.dist - rogue.dist + new.dist
}

if(new.dist >= rogue.dist){
if(nrogue >= 2){ 
nrogue = 1
break;
} else{ nrogue = nrogue+1 }
}
}

if(verbose>0) cat("\n\tFinal  average within block distance of iteration", itr,"::",current.strat.dist/nblocks,".\n")
strat.dist <- c(strat.dist,current.strat.dist)

#print(prev.strat.dist)
#print(current.strat.dist)
if(current.strat.dist < prev.strat.dist){
prev.strat.dist = current.strat.dist
#print(prev.strat.dist)
prev.strat = current.strat
}
cat("best average within block distance::", prev.strat.dist/nblocks,"\n") 


if(itr == Itr4est){
C_K <- 1 + (max(strat.dist)/min(strat.dist)-1)/(2*(max(strat.dist)/min(strat.dist)))^2
pC <- mean(strat.dist <= C_K*min(strat.dist))
pC <- pC/2
eps <- 0.1
Itr_est <- ceiling( (eps^(-1)-1)*(pC^(-1)-1) )
#print(Itr_est )
}
itr <- itr + 1

if(itr > Itr & Itr_est > Itr){
itext <- paste("\nWe suggest at least", (ceiling(Itr_est)-Itr), "more iteration(s) to get potentially better result. 
Enter (y) to Continue, (n) to Stop, or enter the number of more iterations: ")

if(!runallItr){
input <- ifelse(interactive(),readline(itext),
{cat(itext); readLines("stdin",n=1)})
} else input <- min(maxItr, Itr_est) - Itr

if(grepl("n", input)){
 break; } else if(grepl("y", input)) {
Itr <- Itr_est
cat("\nRunning", Itr-(itr-1), "more iteration(s).\n\n")
if(Itr4est<30){ 
Itr4est = ifelse(Itr>30, 30, Itr)
} else if(Itr4est>=30) Itr_est = Itr-1

} else{
if(Itr4est<30) Itr4est = 30
Itr <- Itr+round(as.numeric(input))
if(Itr<itr){ 
break; 
} else {
if(!runallItr) cat("\nRunning", Itr-(itr-1), "more iteration(s)\n\n")
if(Itr4est<30){ Itr4est = ifelse(Itr>30, 30, Itr)
} else if(Itr4est>=30) Itr_est = Itr-1
}
}

} else if(itr > Itr) break;

}

cat("\n")
## summary F statistics
bk = rep(NA, nobs)
names(bk) = units
for(b in 1:nrow(prev.strat))
bk[prev.strat[b,]] = b

unblock = sum(is.na(bk))
  if (unblock>0) bk[is.na(bk)]<-(-(1:unblock))

who<-outer(bk,bk,"==")
inside<-sum(distmat[who])
outside<-sum(distmat[!who])

ndist<-bsize*(bsize-1)*nblocks

#print(c(mean.inside=inside/ndist, mean.outside = outside/(nobs*nobs-(nobs+ndist))))
#return(strat.dist)

if(!is.null(.data)){
.data = cbind(.data, strata=NA)
.data[names(bk),'strata'] = bk

cat("A new column, strata, is added to the data with strata ids.\n")
return(list(strata = prev.strat, .data = .data)) 
##mean.inside=inside/ndist, mean.outside = outside/(nobs*nobs-(nobs+ndist)) ))

}

if(is.null(.data))
return(list(strata = prev.strat))##, mean.inside=inside/ndist, mean.outside = outside/(nobs*nobs-(nobs+ndist))))

}
