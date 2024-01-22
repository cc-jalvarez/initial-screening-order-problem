############
# 
# Initial Screening Order simulations
#
############

### general settings

# cleaning environment
rm(list = ls())

# load packages (and install them, if needed)
llibrary = function(packages) {
  new.packages = packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  for( p in packages) {
    require(p, character.only = TRUE)
  }  
}

# required packages
llibrary( c("TruncatedNormal", "matrixStats", "copula", "doParallel", "latex2exp") ) 
registerDoParallel(cores=detectCores()/4*3) # using 3/4 of cores

# repeatability
set.seed(0)

# set current and output directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
outdir = "../output/" # output directory

# settings in plots
par(mar=c(4,4,1.2,1))
redcol = "#e41a1c"
bluecol = "#377eb8"
greencol = "#4daf4a"
greycol = 8

### simulation global settings

# digit precision in experiments
seqs_psi = seq(.3, .8, .01)
seqs_q = seq(0, .6, .05)

### random generation of applicant scores

# rscorenormv0: truncated normal distribution
vectd0 = Vectorize(function(x) dtmvnorm(x, .5, .02, lb=0, ub=1))
rscorenormv0 = function(n) rtmvnorm(n, .5, .02, lb=0, ub=1)
rscorenormv0(10) # example

# rscorenormv1: a variant with less qualified than unqualified applicants
vectd1 = Vectorize(function(x) dtmvnorm(x, .8, .05, lb=0, ub=1))
rscorenormv1 = function(n) rtmvnorm(n, .8, .05, lb=0, ub=1)
rscorenormv1(10) # example

# rscorenormv2: another variant with mostly qualified applicants
vectd2 = Vectorize(function(x) dtmvnorm(x, 1, .05, lb=0, ub=1))
rscorenormv2 = function(n) rtmvnorm(n, 1, .05, lb=0, ub=1)
rscorenormv2(10) # example

# plot densities
par(mar=c(3,4,2,1))
par(xpd=TRUE)
curve(vectd0(x), xlim=c(0, 1), ylim=c(0, 4), xlab="", ylab="", las=1, col=redcol, lwd=2, lty=1)
mtext("x", side = 1, line = 1.7, font=2, cex=1.2)
mtext("density", side = 2, line = 2.3, font=2, cex=1.2)
box(lwd=2)
curve(vectd1(x), xlim=c(0, 1), ylim=c(0, 3), ylab="density", col=bluecol, add=T, lwd=2, lty=2)
curve(vectd2(x), xlim=c(0, 1), ylim=c(0, 4), ylab="density", col=greencol, add=T, lwd=2, lty=4)
lab_fig1 = c('tN(0.5, 0.02)', 'tN(0.8, 0.05)', 'tN(1, 0.05)')
legend(0.1, 4.58, pch=16, lab_fig1, ncol=3, 
       lty=c(1,2,4), bty='n', col=c(redcol, bluecol, greencol), cex=1.1, pt.cex=1, lwd=2,
       text.width=strwidth("longlonglonglong"))
dev.copy2pdf(file=paste(outdir,"fig1.pdf", sep=''), width = 6, height = 4)

### random generation of rankings, consisting of candidate order, scores and isprotected

#' @param n number of candidates
#' @param s random score generation
#' @param prob fraction of protected
#' @param corr spearman correlation between order and scores (NA = independent, -1 = scores descending)
rscore = function(n, s=rscorenormv0, prob=0.2, corr=NA) {
  applicants = 1:n
  nprot = round(n*prob)
  isprotected = rbinom(n, 1, prob) #sample(c(rep(1, nprot), rep(0,n-nprot)))
  scores = s(n)
  # rearrange scores to achieve corr
  if( !is.na(corr)) {
    sscores = sort(scores)
    norm.cop = normalCopula(corr)
    v = rCopula(n, norm.cop)
    o1 = order(v[,1])
    o2 = order(v[,2])
    a1 = rep(0, length(scores))
    a2 = rep(0, length(scores))
    a1[o1] = sscores
    a2[o2] = sscores
    scores = a1[order(a2)]
  }
  list(applicants=applicants, isprotected=isprotected, scores=scores)
}
# example: independent
rscore(5)
# example: corr=-1 (ordered by descending score)
res = rscore(5, corr=-1)
cat('spearman_corr =', cor(res$applicants, res$scores, method="spearman"))
res

### utility notions, modeled as list of 3 elements
###    f(s) = utility function
###    cmp(u1, u2, k) = comparison between utility u1 and u2
###    isVector = T if utility function returns a vector

# utility as sum and cmp as ratio
uSumRatio = list(f = function(s) if(length(s)==0) NA else sum(s$scores), 
                 cmp = function(u1, u2, k) u2/u1, 
                 isVector=FALSE)
# utility as set and cmp as Jaccard distance
uJSim = list(f = function(s) if(length(s)==0) NULL else s$applicants, 
             cmp = function(u1, u2, k) {
               if(!(is.vector(u1) & is.vector(u2)))
                 return(NA)
               intersection = length(intersect(u1, u2)) 
               union = length(u1) + length(u2) - intersection
               return (intersection/union) 
             }, 
             isVector=TRUE)
# utility as proportion of protected
uProtProp = list(f = function(s) if(length(s)==0) NA else mean(s$isprotected), 
                 cmp = function(u1, u2, k) u2, 
                 isVector=FALSE)

### fatigue functions

curve(dnorm(x, 0, .005), xlim=c(-.6, .6), ylim=c(0, 8))
curve(dnorm(x, 0, .05), add=T)
curve(dnorm(x, 0, .5), add=T)
eps1_fatigue = function(v, t) pmin(1, pmax(0, v + rnorm(length(v), 0, .005*t))) # use t%%120 to model a break every 2 hours
eps1_fatigue(0.3, 1); eps1_fatigue(0.3, 20)
eps1_fatigue(c(0.3, 0.4), c(1, 2))

eps2_fatigue = function(v, t) pmin(1, pmax(0, v + rnorm(length(v), -.005*t, .001*t))) # use t%%120 to model a break every 2 hours
eps2_fatigue(0.3, 1); eps2_fatigue(0.3, 20)
eps2_fatigue(c(0.3, 0.4), c(1, 2))

### algorithms

# computing bestk through ExaminationSearch
#' @param r ranking data
#' @param k number of desired candidates
#' @param q quota of protected (NA = no quota)
bestk = function(r, k, q=0, fatigue=NULL) {
  o = order(if(is.null(fatigue)) r$scores else fatigue(r$scores, 1:length(r$scores)-1), 
          decreasing=TRUE)
  if(q==0) {
    o = o[1:k]
  } else {
    minprot = max(1,round(k*q))
    if(minprot>sum(r$isprotected))
      return(list()) # not enough protected
    nother = k-minprot
    posprot = which(r$isprotected[o]==1)[1:minprot]
    posother = (1:length(o))[-posprot][1:nother]
    o = o[sort(c(posprot, posother))]
  } 
  list(applicants=r$applicants[o], isprotected=r$isprotected[o], scores=r$scores[o])
}
# example:
res = rscore(500)
bestk(res, 6, q = 0)
bk = bestk(res, 6, q = 0.1); bk
mean(bk$isprotected)
bestk(res, 6, q = 0.5, fatigue=eps1_fatigue)
bestk(res, 300, q = 1) # no solution

# computing goodk through CascadeSearch
#' @param r ranking data
#' @param k number of selected candidates
#' @param q quota of protected (NA = no quota)
goodk = function(r, k, psi, q=0, fatigue=NULL) {
  if(is.null(fatigue)) {
    # more efficient version if no fatigue
    o = which(r$scores >= psi)
    if(length(o)<k)
      return(list()) # not enough candidates
    if(q==0) {
      o = o[1:k]
    } else {
      minprot = max(1, round(k*q))
      nother = k-minprot
      posprot = which(r$isprotected[o]==1)
      if(length(posprot)<minprot)
        return(list()) # not enough protected
      posprot = posprot[1:minprot]
      posother = (1:length(o))[-posprot][1:nother]
      o = o[sort(c(posprot, posother))]
    }
    return(list(applicants=r$applicants[o], isprotected=r$isprotected[o], scores=r$scores[o]))
  }
  # fatigue is not NULL
  minprot = ifelse(is.na(q), 0, max(1,round(k*q)))
  maxother = k-minprot
  n = length(r$scores)
  t = 1
  o = c()
  scores = c()
  for(i in 1:n) {
    if(maxother==0 & r$isprotected[i]==0)
      next
    score = fatigue(r$scores[i], t-1)
    t = t + 1
    if(score>= psi) {
      k = k - 1
      if(r$isprotected[i] & minprot>0)
          minprot = minprot-1
      else
          maxother = maxother-1
      o = append(o, i)
      scores = append(scores, r$scores[i]) # append ground-truth scores, not fatigued ones
      if(k==0)
        break
    }
  }
  if(k>0)
    return(list()) # not enough candidates
  list(applicants=r$applicants[o], isprotected=r$isprotected[o], scores=scores)
}
# example:
goodk(res, 6, psi=0.7, q = 0)
goodk(res, 6, psi=0.7, q = 0.5)
goodk(res, 6, psi=0.7, q = 0.5, fatigue=eps1_fatigue)
goodk(res, 6, psi=0.99, q = 0.5) # no solution

### simulations

# single run
#' @param n number of candidates
#' @param prob fraction of protected
#' @param corr spearman correlation between order and scores (NA = independent, -1 = scores descending)
#' @param genscore random generation function
#' @param k number of selected candidates
#' @param q quota of protected (NA = no quota)
#' @param util utility notion
simRun = function(n, prob, corr, genscore, k, q, util, fatigue) {
  # generate scores
  r = rscore(n, s=genscore, prob=prob, corr=corr)
  # best-k algorithm
  bestk_res = bestk(r, k, q=q)
  # utility of best-k solution
  ubestk_res = util$f(bestk_res)
  # utility of good-k solution av variation of psi
  ugoodk_res_s = sapply(seqs_psi, function(psi) util$f(goodk(r, k=k, psi=psi, q=q, fatigue=fatigue)))
  # difference in utility at distance from s_kth
  if (util$isVector) {
    if(is.matrix(ugoodk_res_s)) 
      res = apply(ugoodk_res_s, 2, function(ugoodk_res) util$cmp(ubestk_res, ugoodk_res, k))
    else
      res = unlist(lapply(ugoodk_res_s, function(ugoodk_res) util$cmp(ubestk_res, ugoodk_res, k)))
  }
  else
    res = util$cmp(ubestk_res, ugoodk_res_s, k)
  res
}

simRunBest = function(n, prob, corr, genscore, k, util, fatigue) {
  # generate scores
  r = rscore(n, s=genscore, prob=prob, corr=corr)
  # utility of best-k without fatigue at variation of q
  ubestk_res_s = sapply(seqs_q, function(q) util$f(bestk(r, k=k, q=q)))
  # utility of best-k with fatigue at variation of q
  ugoodk_res_s = sapply(seqs_q, function(q) util$f(bestk(r, k=k, q=q, fatigue=fatigue)))
  # difference in utility at distance from s_kth
  res = util$cmp(ubestk_res_s, ugoodk_res_s, k)
  res
}

# main simulation algorithms
#' @param n number of canidates
#' @param prob fraction of protected
#' @param corr spearman correlation between order and scores (NA = independent, -1 = scores descending)
#' @param genscore random generation function
#' @param k number of selected candidates
#' @param q quota of protected (NA = no quota)
#' @param util utility notion
sim = function(n, iter, prob, corr, genscore, k, q, util, fatigue) {
  res = foreach(id = 1:iter, 
                .combine = cbind, .export=c("simRun", "rscore", "bestk", "goodk", "seqs_psi"), #ls(globalenv()),
                .packages = c("TruncatedNormal", "matrixStats", "copula")) %dopar% {
                  simRun(n, prob=prob, corr=corr, genscore=genscore, k=k, q=q, util=util, fatigue=fatigue)
                }
  list(mean=rowMeans(res, na.rm=TRUE), sd=rowSds(res, na.rm=TRUE))
}

simBest = function(n, iter, prob, corr, genscore, k, util, fatigue) {
  res = foreach(id = 1:iter, 
                .combine = cbind, .export=c("simRunBest", "rscore", "bestk", "goodk", "seqs_q"), #ls(globalenv()),
                .packages = c("TruncatedNormal", "matrixStats", "copula")) %dopar% {
                  simRunBest(n, prob=prob, corr=corr, genscore=genscore, k=k, util=util, fatigue=fatigue)
                }
  list(mean=rowMeans(res, na.rm=TRUE), sd=rowSds(res, na.rm=TRUE))
}

# sequential version (only for debug)
simSeq = function(n, iter, prob, corr, genscore, k, q, util, fatigue) {
  res = replicate(iter, simRun(n, prob=prob, corr=corr, genscore=genscore, k=k, q=q, util=util, fatigue=fatigue))
  list(mean=rowMeans(res, na.rm=TRUE), sd=rowSds(res, na.rm=TRUE))
}

# experimental settings
expbase = function(ns, ks, iter=1000, corr=c(NA,NA,NA), prob=0.2, genscore=rscorenormv0, yline=2.7,
                   q=c(0,0,0), util=uSumRatio, fatigue=NULL, ylab, leglab=NA, legdis=1.06, ylim=c(), output=NA) {
  res1 = sim(n=ns[1], k=ks[1], iter=iter, corr=corr[1], prob=prob, genscore=genscore[[1]], q=q[1], util=util, fatigue=fatigue)
  res2 = sim(n=ns[2], k=ks[2], iter=iter, corr=corr[2], prob=prob, genscore=genscore[[2]], q=q[2], util=util, fatigue=fatigue)
  res3 = sim(n=ns[3], k=ks[3], iter=iter, corr=corr[3], prob=prob, genscore=genscore[[3]], q=q[3], util=util, fatigue=fatigue)
  ylim = if(length(ylim)==0) range(res1$mean, res2$mean, res3$mean, na.rm=TRUE) else ylim
  plot(seqs_psi, res1$mean, col=redcol, las=1, type="l", lty=1, 
       ylim=ylim, xlab="", ylab="", lwd=2)
  mtext(expression(bold(psi)), side = 1, line = 1.7, font=2, cex=1.2)
  mtext(ylab, side = 2, line = yline, font=2, cex=1.2)
  box(lwd=2)
  points(seqs_psi, res2$mean, pch=16, cex=0.3, col=bluecol, lwd=2, type="l", lty=2)
  points(seqs_psi, res3$mean, pch=16, cex=0.3, col=greencol, lwd=2, type="l", lty=4)
  legend(seqs_psi[1]*1.08, ylim[2]*legdis, pch=16, leglab, ncol=3, 
         lty=c(1,2,4), bty='n', col=c(redcol, bluecol, greencol), cex=1.1, pt.cex=1, lwd=2,
         text.width=strwidth("longlonglonglonglong"))
  if(!is.na(output))
    dev.copy2pdf(file=paste(outdir, output, sep=''), width = 6, height = 4)
  list(res1=res1, res2=res2, res3=res3)
}


# experimental settings
expbaseBest = function(ns, ks, iter=1000, corr=c(NA,NA,NA), prob=0.2, genscore=rscorenormv0, yline=2.7,
                   util=uSumRatio, fatigue=NULL, ylab, leglab=NA, legdis=1.06, ylim=c(), output=NA) {
  res1 = simBest(n=ns[1], k=ks[1], iter=iter, corr=corr[1], prob=prob, genscore=genscore[[1]], util=util, fatigue=fatigue)
  res2 = simBest(n=ns[2], k=ks[2], iter=iter, corr=corr[2], prob=prob, genscore=genscore[[2]], util=util, fatigue=fatigue)
  res3 = simBest(n=ns[3], k=ks[3], iter=iter, corr=corr[3], prob=prob, genscore=genscore[[3]], util=util, fatigue=fatigue)
  #boxplot.matrix(res1, use.cols = FALSE)
  ylim = if(length(ylim)==0) range(res1$mean, res2$mean, res3$mean, na.rm=TRUE) else ylim
  plot(seqs_q, res1$mean, col=redcol, las=1, type="l", lty=1, 
       ylim=ylim, xlab="", ylab="", lwd=2)
  mtext("q", side = 1, line = 1.7, font=2, cex=1.2)
  mtext(ylab, side = 2, line = yline, font=2, cex=1.2)
  box(lwd=2)
  points(seqs_q, res2$mean, pch=16, cex=0.3, col=bluecol, lwd=2, type="l", lty=2)
  points(seqs_q, res3$mean, pch=16, cex=0.3, col=greencol, lwd=2, type="l", lty=4)
  legend(seqs_q[1]*1.08, ylim[2]*legdis, pch=16, leglab, ncol=3, 
         lty=c(1,2,4), bty='n', col=c(redcol, bluecol, greencol), cex=1.1, pt.cex=1, lwd=2,
         text.width=strwidth("longlonglonglonglong"))
  if(!is.na(output))
    dev.copy2pdf(file=paste(outdir, output, sep=''), width = 6, height = 4)
  list(res1=res1, res2=res2, res3=res3)
}

exp1 = function(iter=1000, corr=NA, q=.5, util=uSumRatio, fatigue=NULL, ylim=c(), 
                ylab='ratio to baseline', legdis=1.06, output=NA) {
  leglab = lab_fig1
  ns = rep(120,3)
  ks = rep(6,3)
  genscore = c(rscorenormv0, rscorenormv1, rscorenormv2)
  expbase(ns=ns, ks=ks, iter=iter, corr=rep(corr,3), genscore=genscore, q=rep(q,3), util=util, legdis=legdis,
          fatigue=fatigue, ylab=ylab, leglab=leglab, ylim=ylim, output=output)
}

exp1Best = function(iter=1000, corr=NA, util=uSumRatio, fatigue=NULL, ylim=c(), 
                ylab='ratio to baseline', legdis=1.06, output=NA) {
  leglab = lab_fig1
  ns = rep(120,3)
  ks = rep(6,3)
  genscore = c(rscorenormv0, rscorenormv1, rscorenormv2)
  expbaseBest(ns=ns, ks=ks, iter=iter, corr=rep(corr,3), genscore=genscore, util=util, legdis=legdis,
          fatigue=fatigue, ylab=ylab, leglab=leglab, ylim=ylim, output=output)
}

res=exp1(iter=10000, ylim=c(.65, 1), q=.5, output="rq05v012.pdf")
res=exp1(iter=10000, ylim=c(0, 1), q=.5, legdis=1.16, util=uJSim, ylab = 'Jaccard similarity', output="jq05v012.pdf")
#res=exp1(iter=10000, ylim=c(0, 1), q=.5, util=uProtProp, ylab = 'proportion of protected', output="pq05v012.pdf")

res=exp1(iter=1000, fatigue=eps1_fatigue, ylim=c(.65, 1), q=.5)

res=exp1(iter=10000, fatigue=eps1_fatigue, ylim=c(.65, 1), q=.5, output="rq05v012f1.pdf")
res=exp1(iter=10000, fatigue=eps2_fatigue, ylim=c(.65, 1), q=.5, output="rq05v012f2.pdf")
res=exp1(iter=10000, fatigue=eps1_fatigue, ylim=c(0, 1), q=.5, legdis=1.16, util=uJSim, ylab = 'Jaccard similarity', output="jq05v012f1.pdf")

res=exp1Best(iter=1000, fatigue=eps1_fatigue, ylim=c(.65, 1), output="rq05v012f1bk.pdf")
res=exp1Best(iter=10000, fatigue=eps2_fatigue, ylim=c(.65, 1), output="rq05v012f2bk.pdf")

exp2 = function(iter=1000, corr=NA, q=0, util=uSumRatio, fatigue=NULL, ylim=c(), 
                genscore=rscorenormv0, ylab='ratio to baseline', legdis=1.06, output=NA) {
  leglab = c('n=120, k=6', 'n=400, k=20', 'n=30, k=6')
  ns = c(120, 400, 30)
  ks = c(6, 20, 6)
  genscore = c(genscore, genscore, genscore)
  expbase(ns=ns, ks=ks, iter=iter, corr=rep(corr,3), genscore=genscore, q=rep(q,3), util=util, legdis=legdis,
          fatigue=fatigue, ylab=ylab, leglab=leglab, ylim=ylim, output=output)
}

res=exp2(iter=10000, ylim=c(.65, 1), q=.5, output="rq05r3v0.pdf")
#res=exp2(iter=10000, ylim=c(.65, 1), q=.5, genscore=rscorenormv1)
#res=exp2(iter=10000, ylim=c(.65, 1), q=.5, genscore=rscorenormv2)
#res=exp2(iter=10000, ylim=c(0, 1), q=.5, legdis=1.16, util=uJSim, ylab = 'Jaccard similarity')


exp3 = function(iter=1000, q=0, util=uSumRatio, fatigue=NULL, ylim=c(),  
                genscore=rscorenormv0, ylab='ratio to baseline', legdis=1.052, output=NA) {
  leglab = c(TeX(r'($\rho=-1)'), TeX(r'($\rho=-0.8)'), TeX(r'($\rho=-0.5)'))
  ns = rep(120,3)
  ks = rep(6,3)
  genscore = c(genscore, genscore, genscore)
  corr = c(-1, -0.8, -0.5)
  expbase(ns=ns, ks=ks, iter=iter, corr=corr, genscore=genscore, q=rep(q,3), util=util, legdis=legdis,
          fatigue=fatigue, ylab=ylab, leglab=leglab, ylim=ylim, output=output)
}

res=exp3(iter=10000, ylim=c(.65, 1), q=.5, output="rq05c3v0.pdf")
#res=exp3(iter=10000, ylim=c(.75, 1), q=.5, genscore=rscorenormv1)
res=exp3(iter=10000, ylim=c(.65, 1), q=.5, genscore=rscorenormv2, output="rq05c3v2.pdf")
#res=exp3(iter=10000, ylim=c(0, 1), q=.5, util=uJSim, legdis=1.15, ylab = 'Jaccard similarity')
#res=exp3(iter=10000, ylim=c(0, 1), q=.5, genscore=rscorenormv1, util=uJSim, legdis=1.15, ylab = 'Jaccard similarity', output="jq05c3v1.pdf")
#res=exp3(iter=10000, ylim=c(0, 1), q=.5, genscore=rscorenormv1, util=uJSim, legdis=1.15, ylab = 'Jaccard similarity')

res=exp3(iter=10000, fatigue=eps1_fatigue, ylim=c(.65, 1), q=.5, output="rq05c3v0f1.pdf")
res=exp3(iter=10000, fatigue=eps1_fatigue, ylim=c(.65, 1), q=.5, genscore=rscorenormv2, output="rq05c3v2f1.pdf")
res=exp3(iter=10000, fatigue=eps2_fatigue, ylim=c(.65, 1), q=.5, output="rq05c3v0f2.pdf")

exp4 = function(iter=1000, corr=NA, util=uSumRatio, fatigue=NULL, ylim=c(), yline=2.7,
                genscore=rscorenormv0, ylab='ratio to baseline', legdis=1.052, output=NA) {
  leglab = c('q=0', 'q=0.25', 'q=0.5')
  ns = rep(400,3)
  ks = rep(20,3)
  genscore = c(genscore, genscore, genscore)
  q = c(0, 0.25, 0.5)
  expbase(ns=ns, ks=ks, iter=iter, corr=corr, genscore=genscore, q=q, util=util, legdis=legdis,
          fatigue=fatigue, ylab=ylab, leglab=leglab, ylim=ylim, yline=yline, output=output)
}

res=exp4(iter=10000, ylim=c(0, 1), util=uProtProp, legdis=1.15, yline=2.3,
         ylab = expression(bold('f(S'^k*''['good']*')')), output="pq3v0.pdf")
res=exp4(iter=10000, ylim=c(.65, 1), output="rq3v0.pdf")
res=exp4(iter=10000, ylim=c(.65, 1), genscore=rscorenormv2, output="rq3v2.pdf")
#res=exp4(iter=10000, ylim=c(0, 1), genscore=rscorenormv2, util=uJSim, legdis=1.15, ylab = 'Jaccard similarity')

#res=exp4(iter=10000, ylim=c(0, 1), fatigue=eps1_fatigue, util=uProtProp, legdis=1.15, yline=2.3,
#         ylab = expression(bold('f(S'^k*''['good']*')')))
#res=exp4(iter=10000, fatigue=eps1_fatigue, ylim=c(.65, 1))

