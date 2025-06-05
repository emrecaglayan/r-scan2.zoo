# Retrieve the "high quality" mutations used to build the true mutation
# signature.  Nothing special here, this function just makes it easier
# to have a consistent definition for high quality mutations across the
# package code.
get.high.quality.mutations <- function(gatk, muttype=c('snv', 'indel')) {
    mt <- match.arg(muttype)
    gatk[muttype == mt & pass == TRUE]
}


# SCAN2 mutation signature rescue is only applied to sites that would
# pass normally if target.fdr=0.5 (i.e., 50% false discovery rate, which
# is very high; recommended target.fdr is 0.01 (1%)).
reduce.table <- function(gatk, target.fdr) {
    compute.filter.reasons(gatk[static.filter & lysis.fdr <= 0.5 & mda.fdr <= 0.5], target.fdr=target.fdr)
}


# Modifies 'o' by reference. The returned o is not a particularly valid
# SCAN2 object.
#
# N.B. there are currently no reasons to separate SNVs and indels here
# because the various test columns already incoporate differences in
# calling parameters.
compute.filter.reasons <- function(gatk, target.fdr = o@call.mutations$target.fdr) 
{
    m <- gatk[, .(pass = !pass, abc.test, lysis.test = !is.na(lysis.fdr) & 
        lysis.fdr <= target.fdr, mda.test = !is.na(mda.fdr) & 
        mda.fdr <= target.fdr, cigar.id.test, cigar.hs.test, 
        lowmq.test, dp.test, min.sc.alt.test, csf.test)]
    m[is.na(m)] <- FALSE
    gatk[, `:=`(filter.reasons, apply(!m, 1, function(row) paste(colnames(m)[row], 
        collapse = "&")))]
}

mutsig.rescue.one <- function(object, artifact.sig, true.sig,
    target.fdr=object@call.mutations$target.fdr,
    rescue.target.fdr=0.01, muttype=c('snv', 'indel'))
{
    mt <- match.arg(muttype)

    # All work in this function will be done on a copy of the object with a much, much
    # smaller GATK table.  Results will be joined back at the end.
    tmpgatk <- copy(reduce.table(object@gatk, target.fdr=target.fdr))

    sigtype <- if (mt == 'snv') sbs96 else id83
    mutsigs <- sigtype(tmpgatk[muttype == mt & filter.reasons == 'lysis.test']$mutsig)

    sigscores <- get.sig.score(mutsigs=mutsigs,
        artifact.sig=artifact.sig, true.sig=true.sig)

    # it doesn't seem to be possible to use a column assigned by := for another
    # assignment in the same data.table statement.  i.e., to combine all of these
    # into a single statement.
    tmpgatk[muttype == mt & filter.reasons == 'lysis.test',
        rweight := as.numeric(10^-sigscores$postp[mutsig])]  # as.numeric: get rid of table class
    tmpgatk[muttype == mt & filter.reasons == 'lysis.test', rescue.fdr := 
        lysis.pv / (lysis.pv + lysis.beta * rweight * nt/na)]

    # rescue refers uniquely to rescued sites, even though regularly PASSed sites
    # would also meet these criteria.
    tmpgatk[muttype == mt & filter.reasons == 'lysis.test', rescue := 
        !pass & rescue.fdr <= rescue.target.fdr]
    data.table::setkey(tmpgatk, chr, pos, refnt, altnt)  # probably should already be this way

    # Now join the results back to the main (much larger) table.
    # This modifies object by reference, no need to return it.
    object@gatk[tmpgatk, on=.(chr, pos, refnt, altnt),
        c('rweight', 'rescue.fdr', 'rescue') := list(i.rweight, i.rescue.fdr, i.rescue)]


    # Compute the signature homogeneity test w.r.t. the true signature provided
    # to this function.
    sig.homogeneity.test <- sig.homogeneity.test(object, true.sig, muttype)

    list(rescue.target.fdr=rescue.target.fdr,
        sig.homogeneity.test=sig.homogeneity.test,
        postp = sigscores$postp,
        test.spectrum=as.spectrum(mutsigs),
        true.sig=true.sig,
        artifact.sig=artifact.sig,
        nmuts=length(mutsigs),
        weight.true=sigscores$weight.true,
        weight.artifact=sigscores$weight.artifact,
        relative.error=sigscores$rel.error)
}


# Produces a weight vector of the same length as the provided signatures
# (which must all be of the same length). High values (>0) indicate high
# likelihood of originating from the artifact signature; low values (<0)
# indicate high likelihood of originating from the true signature.
#
# mutsigs - mutation signature channel values for each mutation under
#    consideration for rescue.  DO NOT call as.spectrum() before calling
#    get.sig.score().
get.sig.score <- function(mutsigs, true.sig, artifact.sig, eps=0.001) {
    # pracma::lsqnonneg needs a vector; a 'table' doesn't cut it
    test.spectrum <- as.numeric(as.spectrum(mutsigs))

    sigs <- cbind(true.sig, artifact.sig)
    weights <- pracma::lsqnonneg(sigs, test.spectrum)$x
    rel.error <- norm(test.spectrum - as.vector(sigs %*% weights), type='2') / norm(test.spectrum, type='2')

    # Force weight > 0 so ratios always exist; weights can sometimes
    # be small so adding eps may have a large effect; but this is
    # somewhat controlled by using as.spectrum() on test.spectrum,
    # which converts the spectrum into a probability distn (i.e.,
    # sum(all channels) = 1).
    weights <- weights + eps
    # if we assume there are only 2 possible generating signatures,
    # then dividing weights by sum(weights) would convert the fit
    # weights into probabilities. but there's no numerical need for
    # that since we only ever consider the ratio of the weights.
    # a better future method might consider the remaining error in
    # the fit to encapsulate all other unknown processes.
    postp <- log10(artifact.sig*weights[2]) - log10(true.sig*weights[1])
    list(postp=postp, weight.true=weights[1], weight.artifact=weights[2], rel.error=rel.error)
}


sig.homogeneity.test <- function(object, true.sig, muttype=c('snv', 'indel')) {
    muttype <- match.arg(muttype)
    true.muts <- get.high.quality.mutations(object@gatk, muttype=muttype)
    if (muttype == 'snv')
        true.muts <- table(sbs96(true.muts$mutsig))
    if (muttype == 'indel')
        true.muts <- table(id83(true.muts$mutsig))
    
    sig.homogeneity.test.vs.sig(true.muts, true.sig)
}


sig.homogeneity.test.vs.sig <- function(true.muts, true.sig, n.samples=1e5, seed=10) {
    set.seed(10)   # for reproducibility
    logp.cell <- dmultinom(true.muts, size=sum(true.muts), prob=true.sig, log=TRUE)
    randoms <- stats::rmultinom(n.samples, sum(true.muts), prob=true.sig)
    logps <- apply(randoms, 2, dmultinom, size=sum(true.muts), prob=true.sig, log=TRUE)
    mean(logps < logp.cell)
}
