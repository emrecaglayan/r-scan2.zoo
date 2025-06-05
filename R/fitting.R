# n.tiles - if NULL, then downsampling is not performed (uses full set of
#           hSNP training sites).

abmodel.fit.one.chrom = function(path, sc.sample, chrom, genome.seqinfo, n.chunks = 1, 
    n.logp.samples.per.chunk = 20000, hsnp.tilesize = 100, n.tiles = 250, 
    refine.n.steps = 4, refine.top.n = 50, logp.samples.per.step = 20000, 
    alim = c(-7, 2), blim = c(2, 4), clim = c(-7, 2), dlim = c(2, 
        6)) 
{
    chrom <- as.character(chrom)
    if (!(chrom %in% seqnames(genome.seqinfo))) {
        stop(paste0("invalid chromosome name '", chrom, "'\n"))
    }
    region <- as(genome.seqinfo, "GRanges")[chrom, ]
    hsnps <- read.training.hsnps(path = path, sample.id = sc.sample, 
        region = region)
    hsnps <- abmodel.downsample.hsnps(hsnps, hsnp.tilesize = hsnp.tilesize, 
        n.tiles = n.tiles, verbose = TRUE)
        
    # Assign max(scalt, scref) to scalt
    maxes = pmax(hsnps$scalt, hsnps$scref)
    mins = pmin(hsnps$scalt, hsnps$scref)
    hsnps$scalt = maxes
    hsnps$scref = mins

    # Assign max(scalt, scref) to phased.hap1
    hsnps$phased.hap1 = maxes
    hsnps$phased.hap2 = mins


    refine.records <- list()
    for (i in 1:refine.n.steps) {
        cat("Chrom:", chrom, "refinement step", i, "n.chunks", 
            n.chunks, "\n")
        refine.record <- abmodel.refine.parameter.space(hsnps, 
            n.chunks = n.chunks, n.logp.samples.per.chunk = n.logp.samples.per.chunk, 
            top.n = refine.top.n, hsnp.tilesize = hsnp.tilesize, 
            alim = alim, blim = blim, clim = clim, dlim = dlim)
        refine.records[[i]] <- refine.record
        alim <- refine.record$new.param.space[, "a"]
        blim <- refine.record$new.param.space[, "b"]
        clim <- refine.record$new.param.space[, "c"]
        dlim <- refine.record$new.param.space[, "d"]
    }
    list(n.chunks = n.chunks, n.logp.samples.per.chunk = n.logp.samples.per.chunk, 
        n.tiles = n.tiles, hsnp.tilesize = hsnp.tilesize, n.steps = refine.n.steps, 
        top.n = refine.top.n, refine.records = refine.records, 
        fit = refine.records[[refine.n.steps]]$logp.samples[1, 
            , drop = FALSE])
}


abmodel.downsample.hsnps <- function(hsnps, hsnp.tilesize=100, n.tiles=250, verbose=TRUE) {
    if (is.null(n.tiles)) {
        if (verbose) cat('skipping downsampling; n.tiles=NULL\n')
        return(hsnps)
    }
    # Now subset by tiles. This is critical for compute efficiency (>10-fold
    # reduction in runtime with very similar results).
    hsnps$tile.id <- head(rep(1:nrow(hsnps), each=hsnp.tilesize), nrow(hsnps))
    max.tile <- max(hsnps$tile.id)
    tiles.to.use <- sort(sample(max.tile, min(max.tile, n.tiles)))
    if (verbose)
        cat(sprintf("downsampling %d hSNPs (%d tiles) to %d hSNPs (%d tiles)\n",
            nrow(hsnps), max(hsnps$tile.id),
            sum(hsnps$tile.id %in% tiles.to.use), length(tiles.to.use)))
    hsnps[hsnps$tile.id %in% tiles.to.use]
}


# Default values for a/b/c/dlim correspond to the pre-refinment
# parameter space usually searched.
abmodel.refine.parameter.space <-
    function(hsnps, n.chunks=1, n.logp.samples.per.chunk=20000, top.n=50, hsnp.tilesize=100,
    alim=c(-7, 2), blim=c(2, 4), clim=c(-7, 2), dlim=c(2, 6))
{
    cat(sprintf('    param space: a=(%0.3f,%0.3f), b=(%0.3f,%0.3f), c=(%0.3f,%0.3f), d=(%0.3f,%0.3f)\n',
        alim[1], alim[2], blim[1], blim[2], clim[1], clim[2], dlim[1], dlim[2]))

    progressr::with_progress({
        p <- progressr::progressor(along=1:n.chunks)
        p(amount=0, class='sticky')
        logp.samples <- do.call(rbind, future.apply::future_lapply(1:n.chunks, function(i) {
            p(amount=0)
            ctx <- abmodel.approx.ctx(x=hsnps$pos, y=hsnps$phased.hap1,
                d=hsnps$phased.hap1 + hsnps$phased.hap2,
                hsnp.chunksize=hsnp.tilesize)
            ret <- abmodel.sample(n=n.logp.samples.per.chunk,
                alim=alim, blim=blim, clim=clim, dlim=dlim,
                ctx=ctx, seed=NULL)
            p()
            ret
        }, future.seed=0))
        # Randomness here only influences what (a,b,c,d) points are chosen
        # in the parameter space window.
    })

    # Order the samples so that the best fit is always the first row
    logp.samples <- logp.samples[order(logp.samples$logp, decreasing=TRUE),]

    new.param.space <- abmodel.shrink.parameter.space(logp.samples, top.n=top.n)
    list(logp.samples=logp.samples, new.param.space=new.param.space)
}


abmodel.shrink.parameter.space <- function(logp.samples, top.n=50) {
    dn <- dimnames(logp.samples)
    logp.samples <- as.matrix(logp.samples)
    # swapping columns (1,2) and (3,4) to force b < d.
    # ifelse returns a value the same shape as the first argument
    logi.mat <- matrix(rep(logp.samples[,2] < logp.samples[,4], times=5), ncol=5)
    logp.samples <- as.data.frame(ifelse(logi.mat, logp.samples, logp.samples[,c(3,4,1,2,5)]))
    dimnames(logp.samples) <- dn
    logp.samples <- logp.samples[order(logp.samples[,5], decreasing=TRUE),]

    # Use the top 'top.n' logp values to build a new parameter range
    logp.samples[,2] <- log10(logp.samples[,2])   # b and d bounds are in log10 space
    logp.samples[,4] <- log10(logp.samples[,4])
    bounds <- apply(head(logp.samples[,-5], top.n), 2, range)
    colnames(bounds) <- colnames(logp.samples)[-5]

    bounds
}


# Allocate all necessary vectors for the C fitting code.
# Trying to avoid memory copying since these matrices can be large.
abmodel.approx.ctx <- function(x, y, d, hsnp.chunksize=100) {
    n <- hsnp.chunksize

    # numeric() vectors are initialized to 0 by default.
    list(hsnp.chunksize=as.integer(hsnp.chunksize),
         x=x,
         y=y,
         d=d,
         U=numeric(n),
         V=numeric(n),
         B=numeric(n),
         sqrtW=numeric(n),
         K=numeric(n*n),
         A=numeric(n*n))
}

# ctx is the set of working memory space allocated above
# IMPORTANT: the main return of this C function is the approximate
# logp. However, the buffers U, V, B, sqrtW, K and A all contain
# information about the LAST BLOCK of the Laplace approximation.
# This means that calling this function with hsnp.chunksize=length(x)
# (or y or d) returns, among other things, the mode of the Laplace
# approx.
# USE THIS WITH CAUTION!
abmodel.approx.logp <- function(a, b, c, d, ctx,
    max.it=as.integer(50), verbose=FALSE) {
    result <- .Call("laplace_approx_chunk_cpu",
        ctx$hsnp.chunksize, c(a, b, c, d),
        ctx$x, ctx$y, ctx$d, as.integer(length(ctx$x)),
        ctx$U, ctx$V, ctx$B, ctx$sqrtW,
        ctx$K, ctx$A,
        max.it, verbose,
        PACKAGE="scanZooTrial1")
    return(result)
}

# if seed=NULL, don't set it. this assumes that the caller is properly
# setting the seed to be both: (1) different between parallelized chunks and
# (2) reproducible
abmodel.sample <- function(n=1000, alim=c(-7,2), blim=c(2,4), clim=c(-7,2), dlim=c(2,6),
    ctx, seed=0, max.it=50, verbose=FALSE) {

    if (!is.null(seed)) {
        set.seed(seed)
    }
    max.it <- as.integer(max.it)

    params <- data.frame(a=runif(n, min=alim[1], max=alim[2]),
        b=10^runif(n, min=blim[1], max=blim[2]),
        c=runif(n, min=clim[1], max=clim[2]),
        d=10^runif(n, min=dlim[1], max=dlim[2]))

    logps <- mapply(abmodel.approx.logp, a=params$a, b=params$b, c=params$c, d=params$d,
        MoreArgs=list(ctx=ctx, max.it=max.it, verbose=verbose))

    return(cbind(params, logp=logps))
}


# must match the kernel function in src/laplace.c
K.func <- function(x, y, a, b, c, d) exp(a - (x - y)^2 / b^2) + exp(c - (x-y)^2 / d^2)

# From Rasmussen & Williams 2006.  Calculates the conditional distn of
# the latent GP given the observed points without inverting K.  Since
# the latent GP is MVN, "computing the distribution" only requires
# solving for the mean and covariance.
alg3.2.2 <- function(a, b, c, d, ctx, Xnew) {
    # (Using my notation): this conditional distribution is B|Y,X,D.
    # Approximated by the Laplace method, B|Y,X,D ~ MVN(mode, covmat).
    # mode is the maximizer of the nonapproximate distn and
    # covmat=(K + W^-1)^-1, where W is the Hessian of log p(Y|B).
    mode <- ctx$B[1:length(ctx$x)]
    # this is stored in ctx, but I'm not sure that the C code creates
    # a matrix in the format R expects.
    K <- outer(ctx$x, ctx$x, K.func, a=a, b=b, c=c, d=d)
    covK <- outer(ctx$x, Xnew, K.func, a=a, b=b, c=c, d=d)

    # Infer the GP at some new positions Xnew
    # ( Y - d...) is del log p(Y|B=mode)
    mean.new <- t(covK) %*% (ctx$y - ctx$d * exp(mode) / (1 + exp(mode)))

    # v satisfies: v^T v = k(X_sSNV, X)^T (K + W^-1)^-1 k(X_sSNV, X)
    W <- ctx$d * exp(mode) / (1 + exp(mode))^2
    sqrtW <- sqrt(W)
    L <- t(chol(diag(length(ctx$x)) + outer(sqrtW, sqrtW) * K))
    v <- forwardsolve(L, outer(sqrtW, rep(1, ncol(covK))) * covK)

    cov.new <- outer(Xnew, Xnew, K.func, a=a, b=b, c=c, d=d) - t(v) %*% v
    if (any(is.na(mean.new)) | any(is.na(cov.new)))
        stop("NA values predicted")

    list(mean=mean.new, cov=cov.new)
}

# form a large enough block around the set of variants "vars",
# infer the Laplace-approximate distribution of B|Y at the training
# sites within the block, then infer the same approximate distribution
# on B*|Y*, the balances at the candidate variant sites.
# returns the mean and variance of the GP at the candidate sites.
    # WARNING: THIS IS ONLY GUARANTEED TO WORK WHEN ssnvs IS A SINGLE ROW
infer.gp.block <- function(ssnvs, fit, hsnps, ctx, flank=1e5, max.hsnps=150, verbose=FALSE) {
    a <- fit$a
    b <- fit$b
    c <- fit$c
    dparam <- fit$d

    # make a window of [ssnv position - flank, ssnv position + flank]
    # then trim it down to a maximum of max.hsnps in each direction.
    # Remember that ctx has already been allocated assuming its window
    # will be no larger than 2*max hsnps.
    # WARNING: THIS IS ONLY GUARANTEED TO WORK WHEN ssnvs IS A SINGLE ROW
    # the reason is there is  no way to bound the 'middle' term here for
    # an arbitrary list of positions.
    middle <- 0
    right <- findInterval(range(ssnvs$pos)[2], hsnps$pos)
    down <- findInterval(range(ssnvs$pos)[2] + flank, hsnps$pos)
    middle <- down
    down <- min(down, right + max.hsnps)
    left <- findInterval(range(ssnvs$pos)[1], hsnps$pos)
    up <- findInterval(range(ssnvs$pos)[1] - flank, hsnps$pos)
    middle <- middle - up + 1
    up <- max(up, left - max.hsnps)
    window <- c(up, down)

    d <- hsnps[max(window[1], 1):min(window[2], nrow(hsnps)),]
    if (verbose) {
        print(middle)
        print(window)
        cat(sprintf("infer.gp.block: window=%d-%d, %d nearby hets\n",
            min(d$pos), max(d$pos), nrow(d)))
        cat(sprintf("positions:"))
        print(ssnvs$pos)
    }

    # approx. distn of B|Y at the training sites
    ctx$x <- d$pos
    ctx$y <- d$phased.hap1
    ctx$d <- d$phased.hap1 + d$phased.hap2
    abmodel.approx.logp(a=a, b=b, c=c, d=dparam, ctx=ctx)

    # insert the position of the variant to be tested
    z2 <- alg3.2.2(a=a, b=b, c=c, d=dparam, ctx=ctx, Xnew=ssnvs$pos)
    data.frame(gp.mu=z2$mean, gp.sd=sqrt(diag(z2$cov)))
}


# Update: 3/24/2020: now allows inference at a single hSNP by leaving the
# hSNP out of the fitting set.
# "chunks" here are NOT the 250 hSNP blocks used in parameter fitting.
# "ssnvs" are the candidate sSNVs. the data frame need only have a 'pos'
#        column, but should only contain candidates from one chromosome
# "hsnps" should be the phased hSNPs used for fitting, but again only
#        from one chromosome corresponding to ssnvs.
# spikein - if set to TRUE, then ssnvs is expected to be a subset of hsnps.
#           for each spikein snp, AB is estimated by temporarily leaving
#           the single hsnp out of the training set.
#           WARNING: spikein is only meant to be used with chunk=1!
infer.gp <- function(ssnvs, fit, hsnps, chunk=2500, flank=1e5, max.hsnps=150,
    verbose=FALSE, spikein=FALSE) {

    if (verbose) cat(sprintf("mode=%s\n", ifelse(spikein, 'spikein', 'somatic')))
    if (spikein) {
        if (chunk != 1)
            stop("infer.gp: can only run in spikein mode with chunk=1\n")
        ssnv.is.hsnp <- ssnvs$pos %in% hsnps$pos
        cat("infer.gp: building ssnvs <-> hsnps map\n")
        hsnp.map <- sapply(1:nrow(ssnvs), function(i) {
            if (!ssnv.is.hsnp[i]) NA
            else which(hsnps$pos == ssnvs$pos[i])
        })
        #if (length(hsnp.map) != nrow(ssnvs))
            #stop(sprintf("For spike-in mode, 'ssnvs' (%d rows) must be a subset of 'hsnps' (%d rows), but only %d hsnps are in ssnvs", nrow(ssnvs), nrow(hsnps), length(hsnp.map)))
        #ssnvs <- hsnps[ssnvs,,drop=FALSE]
        cat(sprintf("infer.gp: performing %d leave-1-out hSNP AB estimations\n", nrow(ssnvs)))  
    }

    nchunks <- ceiling(nrow(ssnvs)/chunk)

    ctx <- abmodel.approx.ctx(c(), c(), c(), hsnp.chunksize=2*max.hsnps + 10)
    do.call(rbind, lapply(1:nchunks, function(i) {
        if (i %% 100 == 0)
            cat(sprintf("infer.gp: progress: finished %d of %d sites (%0.1f%%)\n",
                i, nchunks, 100*i/nchunks))
        start <- 1 + (i-1)*chunk
        stop <- min(i*chunk, nrow(ssnvs))
        h <- hsnps
        if (spikein) {
            if (ssnv.is.hsnp[i])
                # ssnv i is hsnp hsnp.map[i], so remove it from the training set
                h <- hsnps[-hsnp.map[i],]
        }
        infer.gp.block(ssnvs[start:stop,,drop=FALSE],
            fit, h,
            ctx=ctx, flank=flank, max.hsnps=max.hsnps,
            verbose=verbose)
    }))
}



# special case of infer.gp with chunksize=1 and without any assumptions
# about ssnvs and hsnps overlapping
#
# There are several reasons to use chunksize=1.
# 1. Most importantly, when hSNPs are included in sSNVs (i.e., to use the
#    leave-one-out strategy), they must be handled one at a time.  It would
#    be incorrect to simultaneously score any somatic candidates with an
#    hSNP left out (because of lost information) and it would also be
#    incorrect to simultaneously score any other hSNPs.
# 2. Compute efficiency. infer.gp.block requires matrix inversion, which is
#    >O(n^2), so it is beneficial to make each matrix as small as possible
#    without sacrificing accuracy.
infer.gp1 <- function(ssnvs, fit, hsnps, flank=1e5, max.hsnps=150,
    verbose=FALSE)
{
    ssnv.is.hsnp <- ssnvs$pos %in% hsnps$pos
    if (verbose) {
        cat(sprintf("infer.gp: %d/%d loci scheduled for AB estimation are training hSNPs\n",
            sum(ssnv.is.hsnp), nrow(ssnvs)))
        cat(sprintf("infer.gp: using leave-1-out strategy for AB estimation at hSNPs\n"))
    }

    if (verbose) cat("infer.gp: building ssnvs <-> hsnps map\n")
    hsnp.map <- sapply(1:nrow(ssnvs), function(i) {
        if (!ssnv.is.hsnp[i]) NA
        else which(hsnps$pos == ssnvs$pos[i])
    })


    progressr::with_progress({
        p <- progressr::progressor(along=1:(nrow(ssnvs)/100))
        ctx <- abmodel.approx.ctx(c(), c(), c(), hsnp.chunksize=2*max.hsnps + 10)
        ret <- sapply(1:nrow(ssnvs), function(i) {
            h <- hsnps
            if (ssnv.is.hsnp[i])
                # ssnv i is hsnp hsnp.map[i], so remove it from the training set
                h <- hsnps[-hsnp.map[i],]
            # returns (gp.mu, gp.sd) when chunk=1
            ret <- infer.gp.block(ssnvs[i,,drop=FALSE], fit, h,
                    ctx=ctx, flank=flank, max.hsnps=max.hsnps, verbose=FALSE)
            if (i %% 100 == 1) p()
            c(gp.mu=ret$gp.mu, gp.sd=ret$gp.sd)
        })
    }, enable=verbose)
    t(ret)
}
