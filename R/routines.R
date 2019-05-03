# given two data frames of somatic and germline locations, annotate
# the somatic data frame with the position of the nearest germline entry.
find.nearest.germline <- function(som, germ, chrs=c(1:22,'X')) {
    som$nearest.het <- NA
    for (chr in chrs) {
        gpos <- germ$pos[germ$chr == chr]
        spos <- som$pos[som$chr == chr]
        gidx <- findInterval(spos, gpos)
        nearest.idx <- ifelse(abs(gpos[gidx] - spos) <= abs(gpos[gidx+1] - spos), gidx, gidx+1)
        som$nearest.het[som$chr == chr] <- gpos[nearest.idx]
    }
    som
}

muttype.map <- c(
    'A>C'='T>G',
    'A>G'='T>C',
    'A>T'='T>A',
    'C>A'='C>A',
    'C>G'='C>G',
    'C>T'='C>T',
    'G>A'='C>T',
    'G>C'='C>G',
    'G>T'='C>A',
    'T>A'='T>A',
    'T>C'='T>C',
    'T>G'='T>G'
)

get.3mer <- function(df) {
    require(BSgenome)
    require(BSgenome.Hsapiens.UCSC.hg19)

    comp <- c('A', 'C', 'G', 'T')
    names(comp) <- c('T', 'G', 'C', 'A')

    x <- df

    if (!('muttype' %in% colnames(x))) {
        cat("adding mutation types..\n")
        x$muttype <- muttype.map[paste(x$refnt, x$altnt, sep = ">")]
    }

    x$ctx <- getSeq(BSgenome.Hsapiens.UCSC.hg19,
                    names=paste("chr", x$chr, sep=''),
                    start=x$pos-1, end=x$pos+1, as.character=TRUE)
    x$ctx.rc <- sapply(strsplit(x$ctx, ""),
                    function(s) paste0(comp[s[c(3,2,1)]], collapse=''))

    x$type.and.ctx <- ifelse(x$refnt == 'C' | x$refnt == 'T',
                       paste0(x$ctx, ":", x$muttype),
                       paste0(x$ctx.rc, ":", x$muttype))
    x
}

# somatic and hsnps must have 'af' and 'dp' columns
get.fdr.tuning.parameters <- function(somatic, hsnps, bins=20, random.seed=0)
{
    cat(sprintf("estimating bounds on somatic mutation rate (seed=%d)..\n",
        random.seed))
    # fcontrol -> estimate.somatic.burden relies on simulations to
    # estimate the artifact:mutation ratio.
    set.seed(random.seed)

    max.dp <- as.integer(quantile(hsnps$dp, prob=0.8))
    fcs <- lapply(0:max.dp, function(dp)
        fcontrol(germ.df=hsnps[hsnps$dp == dp,],
                som.df=somatic[somatic$dp == dp,],
                bins=bins)
    )
    fc.max <- fcontrol(germ.df=hsnps[hsnps$dp > max.dp,],
                som.df=somatic[somatic$dp > max.dp,],
                bins=bins)
    fcs <- c(fcs, list(fc.max))
    cat(sprintf("        profiled hSNP and somatic VAFs at depths %d .. %d\n",
        0, max.dp))

    burden <- as.integer(
        c(sum(sapply(fcs, function(fc) fc$est.somatic.burden[1])),  # min est
          sum(sapply(fcs, function(fc) fc$est.somatic.burden[2])))  # max est
    )

    cat(sprintf("        estimated callable somatic mutation burden range (%d, %d)\n",
        burden[1], burden[2]))
    cat("          -> using MAXIMUM burden\n")
    list(bins=bins, burden=burden, fcs=fcs, max.dp=max.dp)
}


apply.fdr.tuning.parameters <- function(somatic, fdr.tuning) {
    somatic$popbin <- ceiling(somatic$af * fdr.tuning$bins)
    somatic$popbin[somatic$dp == 0 | somatic$popbin == 0] <- 1

    nt.na <- mapply(function(dp, popbin) {
        idx = min(dp, fdr.tuning$max.dp+1) + 1
        if (is.null(fdr.tuning$fcs[[idx]]$pops))
            c(0.1, 0.1)
        else
            fdr.tuning$fcs[[idx]]$pops$max[popbin,]
    }, somatic$dp, somatic$popbin)

    nt.na
}

genotype.somatic <- function(gatk, gatk.lowmq, sc.idx, bulk.idx,
    sites.with.ab, sc.cigars, bulk.cigars, cigar.training, cigar.emp.score,
    fdr.tuning, spikein=FALSE,
    cap.alpha=TRUE, cg.id.q=0.05, cg.hs.q=0.05, random.seed=0, target.fdr=0.1,
    bulkref=bulk.idx+1, bulkalt=bulk.idx+2, scref=sc.idx+1, scalt=sc.idx+2,
    min.sc.alt=0, min.sc.dp=0, min.bulk.dp=0)
{
    call.fingerprint <- as.list(environment())
    # almost all of this information is saved in the results
    dont.save <- c('gatk', 'gatk.lowmq', 'sites.with.ab',
        'sc.cigars', 'bulk.cigars')
    call.fingerprint <- call.fingerprint[!(names(call.fingerprint) %in% dont.save)]

    # N.B.: there used to be some randomness in this function, but I
    # believe all of it now resides elsewhere. To be safe, I've left
    # this set.seed so that results can be reproduced.
    set.seed(random.seed)

    cat("step 1: preparing data\n")
    gatk$muttype <- muttype.map[paste(gatk$refnt, gatk$altnt, sep=">")]
    gatk$dp <- gatk[,scalt] + gatk[,scref]
    gatk$af <- gatk[,scalt] / gatk$dp
    gatk$bulk.dp <- gatk[,bulkalt] + gatk[,bulkref]

    # sites only has columns 'chr','pos','refnt','altnt', which match gatk.
    # so this merge call is really just subsetting gatk.
    somatic <- merge(gatk, sites.with.ab, all.y=TRUE)
    cat(sprintf("        %d somatic SNV candidates\n", nrow(somatic)))
    # choose the AB nearest to the AF of each candidate
    # af can be NA if the site has 0 depth
    somatic$gp.mu <- ifelse(!is.na(somatic$af) & somatic$af < 1/2,
        -abs(somatic$gp.mu), abs(somatic$gp.mu))
    somatic$ab <- 1/(1+exp(-somatic$gp.mu))


    cat("step 2: computing p-values for filters\n")
    cat("        allele balance consistency\n")
    somatic$abc.pv <-
        mapply(abc2, altreads=somatic[,scalt], gp.mu=somatic$gp.mu,
            gp.sd=somatic$gp.sd, factor=1, dp=somatic$dp)

    cat("        lysis artifacts\n")
    somatic$lysis.pv <-
        mapply(test2, altreads=somatic[,scalt], gp.mu=somatic$gp.mu,
            gp.sd=somatic$gp.sd, dp=somatic$dp, div=2)

    cat("        MDA artifacts\n")
    somatic$mda.pv <-
        mapply(test2, altreads=somatic[,scalt], gp.mu=somatic$gp.mu,
            gp.sd=somatic$gp.sd, dp=somatic$dp, div=4)


    cat(sprintf("step 3: tuning FDR = %0.3f\n", target.fdr))
    if (cap.alpha)
        cat(sprintf("        cap.alpha=TRUE: alpha <= %0.3f enforced despite artifact prevalence\n", target.fdr))

    nt.na <- apply.fdr.tuning.parameters(somatic, fdr.tuning)

    cat("        lysis artifact FDR\n")
    somatic <- cbind(somatic,
                lysis=t(mapply(match.fdr2,
                            gp.mu=somatic$gp.mu, gp.sd=somatic$gp.sd,
                            dp=somatic$dp,
                            nt=nt.na[1,], na=nt.na[2,],
                            target.fdr=target.fdr, div=2,
                            cap.alpha=cap.alpha)))

    cat("        MDA artifact FDR\n")
    somatic <- cbind(somatic,
                mda=t(mapply(match.fdr2,
                            gp.mu=somatic$gp.mu, gp.sd=somatic$gp.sd,
                            dp=somatic$dp,
                            nt=nt.na[1,], na=nt.na[2,],
                            target.fdr=target.fdr, div=4,
                            cap.alpha=cap.alpha)))


    cat("step 5: applying optional alignment filters\n")
    if (missing(gatk.lowmq)) {
        cat("        WARNING: skipping low MQ filters. will increase FP rate\n")
        somatic$lowmq.test <- TRUE
    } else {
        cat(sprintf("        attaching low MQ data..\n"))
        lmq <- gatk.lowmq[,c('chr', 'pos', 'refnt', 'altnt',
                             colnames(gatk.lowmq)[c(scref,scalt,bulkref,bulkalt)])]
        somatic <- merge(somatic, lmq, by=c('chr', 'pos', 'refnt', 'altnt'),
            all.x=TRUE, suffixes=c('', '.lowmq'))
        # At the moment only removing sites that have bulk support when
        # the MQ threshold is lowered.
        cn <- paste0(colnames(gatk.lowmq)[bulkalt], '.lowmq')
        # In spikein mode, known hSNPs are being treated like somatics,
        # so it is necessary to short circuit this test.
        somatic$lowmq.test <- is.na(somatic[,cn]) | somatic[,cn] == 0 | spikein
    }

    somatic <- merge(somatic, sc.cigars, by=c('chr', 'pos'), all.x=T)
    somatic <- merge(somatic, bulk.cigars, by=c('chr', 'pos'), all.x=T,
        suffixes=c('', '.bulk'))
    somatic$id.score.y <- somatic$ID.cigars / somatic$dp.cigars
    somatic$id.score.x <- somatic$ID.cigars.bulk / somatic$dp.cigars.bulk
    somatic$id.score <- cigar.emp.score(training=cigar.training, test=somatic, which='id')
    somatic$hs.score.y <- somatic$HS.cigars / somatic$dp.cigars
    somatic$hs.score.x <- somatic$HS.cigars.bulk / somatic$dp.cigars.bulk
    somatic$hs.score <- cigar.emp.score(training=cigar.training, test=somatic, which='hs')

    cat("        Excessive indel CIGAR ops\n")
    somatic$cigar.id.test <-
        somatic$id.score > quantile(cigar.training$id.score, prob=cg.id.q, na.rm=T)
    cat("        Excessive clipped read CIGAR ops\n")
    somatic$cigar.hs.test <-
        somatic$hs.score > quantile(cigar.training$hs.score, prob=cg.hs.q, na.rm=T)
    

    somatic$dp.test <- somatic$dp >= min.sc.dp & somatic$bulk.dp >= min.bulk.dp

    cat("step 6: calling somatic SNVs\n")
    somatic$pass <-
        somatic$abc.pv > 0.05 &
        somatic$lysis.pv <= somatic$lysis.alpha &
        somatic$mda.pv <= somatic$mda.alpha &
        somatic$cigar.id.test & somatic$cigar.hs.test &
        somatic$lowmq.test & somatic$dp.test & somatic[,scalt] >= min.sc.alt
    cat(sprintf("        %d passing somatic SNVs\n", sum(somatic$pass)))
    cat(sprintf("        %d filtered somatic SNVs\n", sum(!somatic$pass)))

    # return non-data parameters and the output
    # the final RDA files are already large enough to be painful to work with
    return(c(call.fingerprint, list(somatic=somatic)))
}

# probability distribution of seeing y variant reads at a site with
# depth d with the estimate (gp.mu, gp.sd) of the local AB. model is
#     Y|p ~ Bin(d, 1/(1+exp(-b))),    b ~ N(gp.mu, gp.sd)
# the marginal distribution of Y (# variant reads) requires integrating
# over b, which has no closed form solution. we approximate it using
# gauss-hermite quadrature. increasing the number of nodes in ghd
# increases accuracy.
# 'factor' allows changing the relationship of ab -> af
# NOTE: for many calls, is best for the caller to compute ghd once
# and supply it to successive calls.
dreads <- function(ys, d, gp.mu, gp.sd, factor=1, ghd=gaussHermiteData(128)) {
    require(fastGHQuad)
    
    sapply(ys, function(y)
        ghQuad(function(x) {
                # ghQuad requires the weight function on x to be exp(-x^2)
                # (which is NOT a standard normal)
                b <- sqrt(2)*gp.sd*x + gp.mu
                exp(dbinom(y, size=d, prob=1/(factor*(1+exp(-b))), log=TRUE) - log(pi)/2)
            }, ghd
        )
    )
}

# determine the beta (power) as a function of alpha (FP rate)
# note that this is subject to heavy integer effects, specifically
# when D is low or when af is close to 0 or 1
# div controls the error model
estimate.alphabeta2 <- function(gp.mu, gp.sd, dp=30, alphas=10^seq(-4,0,0.05), div=2) {
    td <- data.frame(dp=0:dp,
        mut=dreads(0:dp, d=dp, gp.mu=gp.mu, gp.sd=gp.sd),
        err1=dreads(0:dp, d=dp, gp.mu=gp.mu, gp.sd=gp.sd, factor=div),
        err2=dreads(0:dp, d=dp, gp.mu=-gp.mu, gp.sd=gp.sd, factor=div))

    # (td[,4]+td[,5])/2 corresponds to an even mixture of errors on the p
    # and 1-p alleles.
    # XXX: i really don't know how I feel about this mixture. but use it
    # for now, maybe later come up with a better prior for the two error
    # modes.
    td$err <- (td$err1 + td$err2)/2
    td <- td[order(td$err),]
    td$cumerr <- cumsum(td$err)

    # FP probability assuming mutation occurs on this allele (wrt mu)
    fps <- sapply(alphas, function(alpha) {
        reject <- td$cumerr <= alpha
        # because there are two potential H0 models, a single alpha
        # does not make sense. Instead, use the maximum alpha to stay
        # conservative.
        max.alpha <- max(c(0, td$cumerr[td$cumerr <= alpha]))
        c(alpha=max.alpha, beta=sum(c(0, td$mut[reject == TRUE])))
    })

    list(td=td, alphas=fps[1,], betas=fps[2,], input.alphas=alphas)
}

# NOT equivalent to test(div=1)
abc2 <- function(altreads, gp.mu, gp.sd, dp, factor=1) {
    pmf <- dreads(0:dp, d=dp, gp.mu=gp.mu, gp.sd=gp.sd, factor=factor)
    p.value <- sum(pmf[pmf <= pmf[altreads + 1]])
    c(p.value=p.value)  # just names the return
}

test2 <- function(altreads, gp.mu, gp.sd, dp, div) {
    # equal mixture of errors in cis and trans
    err <- (dreads(0:dp, d=dp, gp.mu=gp.mu, gp.sd=gp.sd, factor=div) +
            dreads(0:dp, d=dp, gp.mu=-gp.mu, gp.sd=gp.sd, factor=div))/2

    p.value <- sum(err[err <= err[altreads + 1]])
    c(p.value=p.value)  # just names the return
}


bin.afs <- function(afs, bins=20) {
    sq <- seq(0, 1, 1/bins)
    x <- findInterval(afs, sq, left.open=TRUE)
    # af=0 will be assigned to bin 0 because intervals are (.,.]
    x[x==0] <- 1
    tx <- tabulate(x, nbins=bins)
    names(tx) <- apply(cbind(sq[-length(sq)], sq[-1]), 1, mean)
    tx
}

# the "rough" interval is not a confidence interval. it is just a
# heuristic that demonstrates the range of  reasonably consistent
# somatic mutation burdens
# WARNING! this is the *callable somatic burden*. if you want an
# estimate of the total somatic burden genome-wide, adjust this
# number by the fraction of genome represented in the somatic
# candidate set.
estimate.somatic.burden <- function(fc, min.s=1, max.s=5000, n.subpops=10, display=FALSE, rough.interval=0.99) {
    sim <- function(n.muts, g, s, n.samples=1000, diagnose=FALSE) {
        samples <- rmultinom(n=n.samples, size=n.muts, prob=g)
        if (diagnose) {
            boxplot(t(samples))
            lines(s, lwd=2, col=2)
        }
        mean(apply(samples, 2, function(col) all(col <= s)))
    }

    # determine how often a sample of N somatic mutations from the
    # proper het distribution (germlines) "fits" inside the somatic
    # candidate distribution.
    # always try nmut=1-100
    srange <- c(1:100, seq(101, max.s, length.out=n.subpops))
    srange <- srange[srange < max.s]
    fraction.embedded <- sapply(srange, sim, g=fc$g, s=fc$s)
    # max(c(0,... in some cases, there will be no absolutely no overlap
    # between hSNP VAFs and somatic VAFs, leading to fraction.embedded=0
    # for all N.  In this case, return 0 and let the caller decide what
    # to do.
    min.burden <- max(c(0, srange[fraction.embedded >= 1 - (1 -
        rough.interval)/2]))
    max.burden <- max(c(0, srange[fraction.embedded >= (1 - rough.interval)/2]))
    c(min=min.burden, max=max.burden)
}

# {germ,som}.df need only have columns named dp and af
# estimate the population component of FDR
# ignore.100 - ignore variants with VAF=100
fcontrol <- function(germ.df, som.df, bins=20, rough.interval=0.99) {
    germ.afs <- germ.df$af[!is.na(germ.df$af) & germ.df$af > 0]
    som.afs <- som.df$af[!is.na(som.df$af)]
    g <- bin.afs(germ.afs, bins=bins)  # counts, not probabilities
    s <- bin.afs(som.afs, bins=bins)

    # when fcontrol is used on small candidate sets (i.e., when controlling
    # for depth), there may be several 0 bins in s.
#    g <- g*1*(s > 0)
    # XXX: i don't like this. it doesn't quite match the principle, but
    # rather addresses a limitation of the heuristic method.

    if (length(s) == 0 | all(g == 0))
        return(list(est.somatic.burden=c(0, 0),
             binmids=as.numeric(names(g)),
             g=g, s=s, pops=NULL))

    # returns (lower, upper) bound estimates
    approx.ns <- estimate.somatic.burden(fc=list(g=g, s=s),
        min.s=1, max.s=sum(s), n.subpops=min(sum(s), 100),
        rough.interval=rough.interval)

cat(sprintf("fcontrol: dp=%d, max.s=%d (%d), n.subpops=%d, min=%d, max=%d\n",
germ.df$dp[1], nrow(som.df), sum(s), min(nrow(som.df),100),
as.integer(approx.ns[1]), as.integer(approx.ns[2])))
    pops <- lapply(approx.ns, function(n) {
#        nt <- pmax(n*(g/sum(g))*1*(s > 0), 0.1)
nt <- pmax(n*(g/sum(g)), 0.1)
        # ensure na > 0, since FDR would be 0 for any alpha for na=0
        # XXX: the value 0.1 is totally arbitrary and might need to be
        # more carefully thought out.
        na <- pmax(s - nt, 0.1)
        cbind(nt=nt, na=na)
    })

    return(list(est.somatic.burden=approx.ns,
         binmids=as.numeric(names(g)),
         g=g, s=s, pops=pops))
}

# SUMMARY
# given a target FDR, try to match the FDR as best as we can without
# exceeding it.
#
# DETAILS
# for a given (af, dp), find the relation (alpha, beta) to determine
# what sensitivities and specificities are achievable.  using the
# estimated numbers of true mutations nt and artifacts na in the
# population, calculate the corresponding FDRs (alpha,beta,FDR).
# return (alpha, beta, FDR) such that
# FDR = arg max { FDR : FDR <= target.fdr }
# to break ties (when more than one (alpha,beta) produce the same
# FDR), select the triplet with minimal alpha and maximal beta.
# cap.alpha - because there is often such a strong spike at low VAF
#             amongst somatic candidate sites, the FDR control
#             procedure sometimes results in prior error rates of <1%
#             amongst high VAF candidate sSNVs. In this case, it will
#             often accept variants that are very clearly more consistent
#             with the artifact model(s) than with a true mutation. The
#             cap.alpha parameter helps to address this issue by requiring
#             alpha <= target FDR.
#             To be more clear: when nt >> na, FDR matching can pass
#             variants that have very high alpha (e.g., cases with
#             alpha ~ 0.4, beta ~ 0.7). This puts significant pressure
#             on the accuracy of the nt and na estimates.
match.fdr2 <- function(gp.mu, gp.sd, dp, nt, na, target.fdr=0.1, div=2, cap.alpha=TRUE) {
    alphabeta <- estimate.alphabeta2(gp.mu=gp.mu, gp.sd=gp.sd, dp=dp, div=div)
    x <- data.frame(alpha=alphabeta$alphas, beta=alphabeta$betas,
               fdr=ifelse(alphabeta$alphas*na + alphabeta$betas*nt > 0,
                          alphabeta$alphas*na / (alphabeta$alphas*na + alphabeta$betas*nt), 0))
    x <- rbind(c(0, 0, 0), x)  # when no parameter meets target fdr
    # use (a,) from the highest FDR <= target.fdr
    x <- x[x$fdr <= target.fdr,] 
    if (cap.alpha)
        x <- x[x$alpha <= target.fdr,]
    x <- x[x$fdr == max(x$fdr),]  # may be more than one row
    # breaking ties by minimum alpha, maximum beta
    x <- x[x$alpha == min(x$alpha),]
    x <- x[x$beta == max(x$beta),]
    unlist(x[1,])
}




###############################################################################
# various plotting routines
###############################################################################


plot.ab <- function(ab) {
    layout(matrix(1:4,nrow=2,byrow=T))
    td <- ab$td[order(ab$td$dp),]
    plot(td$dp, td$mut, type='l', ylim=range(td[,c('mut', 'err1', 'err2')]),
        xlab="Depth", ylab="Model probabiblity")
    lines(td$dp, pmax(td$err1,td$err2), lty='dotted', col=2)
    plot(x=ab$input.alphas, y=ab$alphas, xlab="Requested alpha", ylab="Estimated alpha", log="xy")
    plot(ab$alphas, ab$betas, xlab='FP rate', ylab='Power')
    abline(h=0, lty='dotted')
    plot(ab$alphas, ab$betas, log='x', xlab='log(FP rate)', ylab='Power')
    abline(h=0, lty='dotted')
}


# for 96 dimensional mut sigs
mutsig.cols <- rep(c('deepskyblue', 'black', 'firebrick2', 'grey', 'chartreuse3', 'pink2'), each=16)

plot.3mer <- function(x, no.legend=FALSE, ...) {
    bases <- c("A", 'C', 'G', 'T')

    # need to make a table of all possible contexts because they may not
    # be observed after filtering.
    t <- rep(0, 96)
    names(t) <- paste0(rep(bases, each=4),
                      rep(c('C', 'T'), each=48),
                      rep(bases, times=4),
                      ":",
                      rep(c('C', 'T'), each=48),
                      ">",
                      c(rep(c('A', 'G', 'T'), each=16),
                        rep(c('A', 'C', 'G'), each=16)))
    print(table(x$ctx))
    t2 <- table(x$type.and.ctx)
    t[names(t2)] <- t2
    tn <- do.call(rbind, strsplit(names(t), ":"))
    t <- t[order(tn[,2])]
    print(t)
    p <- barplot(t, las=3, col=mutsig.cols, names.arg=tn[order(tn[,2]), 1], space=0.5, border=NA, ...)
    abline(v=(p[seq(4,length(p)-1,4)] + p[seq(5,length(p),4)])/2, col='grey')
    if (!no.legend)
        legend('topright', ncol=2, legend=sort(unique(tn[,2])),
            fill=mutsig.cols[seq(1, length(mutsig.cols), 16)])
}


plot.fcontrol <- function(fc) {
    layout(matrix(1:(1+length(fc$pops)), nrow=1))
    plot(fc$binmids, fc$g/sum(fc$g),
        ylim=range(c(fc$g/sum(fc$g), fc$s/sum(fc$s))),
        type='l', lwd=2)
    lines(fc$binmids, fc$s/sum(fc$s), col=2, lwd=2)

    for (i in 1:length(fc$pops)) {
        pop <- fc$pops[[i]]
        barplot(names.arg=fc$binmids, t(pop), col=1:2,
            main=sprintf('Assumption: ~%d true sSNVs', sum(round(pop[,1],0))), las=2)
        legend('topright', fill=1:2, legend=c('Ntrue', 'Nartifact'))
    }
}


# using the (alpha, beta) relationships and fcontrol population
# estimations, determine average sensitivity per AF with a
# (theoretically) controlled FDR
plot.fdr <- function(fc, dps=c(10,20,30,60,100,200), target.fdr=0.1, div=2) {
    afs <- fc$binmids
    layout(matrix(1:(3*length(fc$pops)), nrow=3))
    for (i in 1:length(fc$pops)) {
        # from fcontrol: pops[[i]] has rows corresponding to afs
        pop <- fc$pops[[i]]
        l <- lapply(dps, function(dp) {
            sapply(1:length(afs), function(i)
                match.fdr(afs[i], dp, nt=pop[i,1], na=pop[i,2],
                    target.fdr=target.fdr, div=div)
            )
        })
        matplot(x=afs, sapply(l, function(ll) ll[3,]), type='l', lty=1,
            main=sprintf("Assuming %d true sSNVs", sum(pop[,1])),
            xlab="AF (binned)", ylab="FDR", ylim=c(0, 1.1*target.fdr))
        abline(h=target.fdr, lty='dotted')
        matplot(x=afs, sapply(l, function(ll) ll[1,]), type='l', lty=1,
            xlab="AF (binned)", ylab="log(alpha)", log='y', ylim=c(1e-5,1))
        abline(h=10^-(1:5), lty=2)
        matplot(x=afs, sapply(l, function(ll) ll[2,]), type='l', lty=1,
            xlab="AF (binned)", ylab="Power", ylim=0:1)
    }
}


plot.ssnv.region <- function(chr, pos, alt, ref, fits, fit.data, upstream=5e4, downstream=5e4, gp.extend=1e5, n.gp.points=100, blocks=50) {
    d <- fit.data[fit.data$chr==chr & fit.data$pos >= pos - upstream &
                  fit.data$pos <= pos + downstream,]

    cat("estimating AB in region..\n")
    # ensure that we estimate at exactly pos
    est.at <- c(seq(pos - upstream, pos-1, length.out=n.gp.points/2), pos,
                seq(pos+1, pos + downstream, length.out=n.gp.points/2))
    fit.chr <- fits[[chr]]
    cat("WARNING: if this function produces NA predictions, try reducing the value of 'blocks'\n")
    # infer.gp can fail when using chunk>1.
    # this does not affect the calling code, just this plot
    gp <- infer.gp(ssnvs=data.frame(pos=est.at),
        fit=fit.chr, hsnps=fit.data[fit.data$chr == chr,],
        chunk=blocks, flank=gp.extend, max.hsnp=500)
    gp$pos <- est.at
    plot.gp.confidence(df=gp, add=FALSE)
    points(d$pos, d$hap1/d$dp, pch=20, ylim=0:1)
    af <- alt/(alt+ref) # adjust to closest AB
    gp.at <- gp[gp$pos == pos,]$gp.mu
    ab <- 1/(1+exp(-gp.at))
    af <- ifelse(abs(af - ab) <= abs(af - (1-ab)), af, 1-af)
    points(pos, af, pch=20, cex=1.25, col=2)
}


# NOTE: the 95% probability interval is in the NORMAL space, not the
# fraction space [0,1]. After the logistic transform, the region in
# fraction space may not be an HDR.
plot.gp.confidence <- function(pos, gp.mu, gp.sd, df, sd.mult=2,
    logspace=FALSE, tube.col=rgb(0.9, 0.9, 0.9),
    line.col='black', tube.lty='solid', line.lty='solid', add=TRUE)
{
    if (!missing(df)) {
        pos <- df$pos
        gp.mu <- df$gp.mu
        gp.sd <- df$gp.sd
    }

    # maybe the analytical solution exists, but why derive it
    if (!logspace) {
        cat("transforming to AF space...\n")
        sd.upper <- 1 / (1 + exp(-(gp.mu + sd.mult*gp.sd)))
        sd.lower <- 1 / (1 + exp(-(gp.mu - sd.mult*gp.sd)))
        gp.mu <- 1 / (1 + exp(-gp.mu))
    } else {
        sd.upper <- gp.mu + sd.mult*gp.sd
        sd.lower <- gp.mu - sd.mult*gp.sd
    }

    if (!add) {
        plot(NA, NA, xlim=range(pos), ylim=0:1)
    }

    polygon(c(pos, rev(pos)), c(gp.mu, rev(sd.lower)),
        col=tube.col, border=line.col, lty=tube.lty)
    polygon(c(pos, rev(pos)), c(gp.mu, rev(sd.upper)),
        col=tube.col, border=line.col, lty=tube.lty)
    lines(pos, gp.mu, lwd=2, col=line.col, lty=line.lty)
}



#############################################################################
# routines for joint calling
#############################################################################

joint.caller <- function(dfs, hdfs) {
    jdf <- jhelper(dfs)
    jhdf <- jhelper(hdfs)
    cutoffs <- get.joint.cutoffs(jhdf, n=length(dfs))
    jdf$pass <- apply.joint.cutoffs(jdf, cutoffs)
    dfs <- lapply(dfs, function(df) {
        # af > 1/dp is a roundabout way to enforce alt>1, because
        # at this point in the pipeline the column corresponding to
        # alt counts for this sample has been lost.
        df$jpass <- as.logical(jdf$pass * (df$af > 0) &
            df$dp.test & df$lowmq.test & df$af > 1/df$dp)
        df$jabc <- jdf$jabc
        df
    })
    dfs
}



jhelper <- function(dfs, gp.sd.penalty=0.1, gp.sd.cutoff=1) {
    n=length(dfs)
    cps <- 4 # columns per sample
    jdf <- do.call(cbind, lapply(dfs,function(df) df[,c('af','dp','abc.pv','gp.sd')]))

    jdf$ncells <- rowSums(!apply(jdf[,seq(1,n*cps,cps)], 1:2, is.nan) & jdf[,seq(1,n*cps,cps)]>0)

    # column reordering makes it easy to recover af, dp, pv
    jdf$jabc <- apply(jdf[,c(seq(1,n*cps,cps), seq(2,n*cps,cps), seq(3,n*cps,cps), seq(4,n*cps,cps))], 1,
        function(row) {
            afs <- row[1:n]
            dps <- row[(n+1):(2*n)]
            pvs <- row[(2*n+1):(3*n)]
            sds <- row[(3*n+1):(4*n)]
            n.penalties <- sum(sds[dps>0 & afs>0] >= gp.sd.cutoff)
            prod(pvs[dps>0 & afs>0]) * prod(rep(gp.sd.penalty, n.penalties))
    })
    jdf
}

# n - number of samples
get.joint.cutoffs <- function(jdf, n, q=0.1) {
    counts <- sapply(2:n, function(i) sum(jdf$ncells == i))
    if (any(counts) == 0) {
        cat("counts for ncells:\n")
        for (i in 2:n)
            cat(sprintf("ncells=%2d, count=%d\n", i, counts[i-1]))
        stop("Found ncell count=0 (see above).  Try increasing the number of hSNP spikeins")
    }

    cutoffs <- sapply(2:n, function(i)
        quantile(jdf$jabc[jdf$ncells == i], prob=q))
    # ncells=0,1: don't pass any of these. they will be determined only by
    # single sample statistics
    data.frame(ncells=0:n, cutoff=c(Inf, Inf, cutoffs))
}


apply.joint.cutoffs <- function(jdf, cutoffs) {
    mapply(function(ncells, jabc)
        jabc >= cutoffs$cutoff[cutoffs$ncells == ncells],
        ncells=jdf$ncells, jabc=jdf$jabc)
}
