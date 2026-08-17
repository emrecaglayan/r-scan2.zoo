setClassUnion('null.or.df', c('NULL', 'data.frame'))
setClassUnion('null.or.dt.or.raw', c('NULL', 'data.table', 'raw'))
setClassUnion('null.or.Seqinfo', c('NULL', 'Seqinfo'))
setClassUnion('null.or.GRanges', c('NULL', 'GRanges'))
setClassUnion('null.or.list', c('NULL', 'list'))
setClassUnion('null.or.logical', c('NULL', 'logical'))
setClassUnion('null.or.character', c('NULL', 'character'))

# if adding new slots to the class, be sure to update concat(...)
# to combine the slots properly.
setClass("SCAN2", slots=c(
    region='null.or.GRanges',
    genome.string='character',
    genome.seqinfo='null.or.Seqinfo',
    single.cell='character',
    bulk='character',
    gatk="null.or.dt.or.raw",
    integrated.table.path='null.or.character',
    resampled.training.data="null.or.list",
    ab.fits='null.or.df',
    static.filter.params='null.or.list',
    ab.estimates='null.or.df',
    mut.models='null.or.df',
    cigar.data='null.or.df',
    excess.cigar.scores='null.or.list',
    fdr.prior.data='null.or.list',
    fdr='null.or.list',
    call.mutations='null.or.list',
    depth.profile='null.or.list',
    mutburden='null.or.list',
    mutsig.rescue='null.or.list',
    spatial.sensitivity='null.or.list'))


# "flip" the GP mu around 0 to best match the AF of each candidate mutation. 
match.ab <- function(af, gp.mu) {
    ifelse(af < 1/2, -abs(gp.mu), abs(gp.mu))
}

# currently will always fire when SCAN2 is built from chunked objects.
# firing isn't a problem, just causes a warning message.
# need to define what the "full genome" is (in terms of GRanges objects)
# for the recognized genomes.
check.chunked <- function(object, message) {
    if (!is.null(object@region))
        warning(message)
}

genome.string.to.seqinfo.object = function(genome = c("hs37d5", "hg38", "mm10", "mm39", "rheMac10", "panTro6", "calJac4", "musFur1")) 
{
    genome <- match.arg(genome)
    if (genome == "hs37d5") {
        return(GenomeInfoDb::Seqinfo(genome = "GRCh37.p13"))
    }
    else if (genome == "hg38") {
        return(GenomeInfoDb::Seqinfo(genome = "hg38"))
    }
    else if (genome == "mm10") {
        return(GenomeInfoDb::Seqinfo(genome = "mm10"))
    }
    else if (genome == "rheMac10") {
        return(GenomeInfoDb::Seqinfo(genome = "rheMac10"))
    }
    else if (genome == "mm39") {
        return(GenomeInfoDb::Seqinfo(genome = "mm39"))
    }
    else if (genome == "panTro6") {
        return(GenomeInfoDb::Seqinfo(genome = "panTro6"))
    }
    else if (genome == "calJac4") {
        return(GenomeInfoDb::Seqinfo(genome = "calJac4"))
    }
    else if (genome == "musFur1") {
        return(GenomeInfoDb::Seqinfo(genome = "musFur1"))
    }
    else {
        stop("unsupported genome string")
    }
}

genome.string.to.bsgenome.object = function(genome = c("hs37d5", "hg38", "mm10", 'rheMac10', 'mm39', 'panTro6', 'calJac4', 'musFur1'))                                 
{
    genome <- match.arg(genome)
    if (genome == "hs37d5") {
        require(BSgenome.Hsapiens.1000genomes.hs37d5)                          
        genome <- BSgenome.Hsapiens.1000genomes.hs37d5                         
    }
    else if (genome == "hg38") {
        require(BSgenome.Hsapiens.UCSC.hg38)                                   
        genome <- BSgenome.Hsapiens.UCSC.hg38                                  
    }
    else if (genome == "mm10") {
        require(BSgenome.Mmusculus.UCSC.mm10)                                  
        genome <- BSgenome.Mmusculus.UCSC.mm10                                 
    }
    else if (genome == "rheMac10") {
        require(BSgenome.Mmulatta.UCSC.rheMac10)                                  
        genome <- BSgenome.Mmulatta.UCSC.rheMac10                               
    }
    else if (genome == "mm39") {
        require(BSgenome.Mmusculus.UCSC.mm39)                                  
        genome <- BSgenome.Mmusculus.UCSC.mm39                               
    }
    else if (genome == "panTro6") {
        require(BSgenome.Ptroglodytes.UCSC.panTro6)                                  
        genome <- BSgenome.Ptroglodytes.UCSC.panTro6                                 
    }
    else if (genome == "calJac4") {
        require(BSgenome.Cjacchus.UCSC.calJac4)                                  
        genome <- BSgenome.Cjacchus.UCSC.calJac4                                 
    }
    else if (genome == "musFur1") {
        require(BSgenome.Mfuro.UCSC.musFur1)                                  
        genome <- BSgenome.Mfuro.UCSC.musFur1                                 
    }
    else {
        stop("unsupported genome string")                                      
    }
    genome
}


genome.string.to.tiling = function (genome, tilewidth = 1e+07, scaffolds = '') 
{
    sqi <- genome.string.to.seqinfo.object(genome)
    sqi = sqi[scaffolds]
    sqi_levels = seqlevels(sqi)
    chroms.to.tile = sqi_levels

    grs <- GenomicRanges::tileGenome(seqlengths = sqi[chroms.to.tile], 
        tilewidth = tilewidth, cut.last.tile.in.chrom = TRUE)
    grs
}



make.scan = function(single.cell, bulk, genome = c("hs37d5", "hg38", "mm10", 'rheMac10', 'mm39', 'panTro6', 'calJac4', 'musFur1'), region = NULL) 
{
    genome <- match.arg(genome, choices = c("hs37d5", "hg38", "mm10", 'rheMac10', 'mm39', 'panTro6', 'calJac4', 'musFur1'))
    new("SCAN2", single.cell = single.cell, bulk = bulk, genome.string = genome, 
        genome.seqinfo = genome.string.to.seqinfo.object(genome), 
        region = region, gatk = NULL, ab.fits = NULL, ab.estimates = NULL, 
        mut.models = NULL, cigar.data = NULL, excess.cigar.scores = NULL, 
        static.filter.params = NULL, fdr.prior.data = NULL, fdr = NULL, 
        call.mutations = NULL, mutburden = NULL, mutsig.rescue = NULL, 
        spatial.sensitivity = NULL)
}


setValidity("SCAN2", function(object) {
    if (length(object@genome.string) != 1)
        return("must provide exactly one genome name")

    if (!is.null(object@genome.seqinfo)) {
        if (class(object@genome.seqinfo) != 'Seqinfo')
            stop('@genome.seqinfo must be of class Seqinfo')
    }

    if (length(object@single.cell) != 1)
        return("must provide exactly one single cell sample name")

    if (length(object@bulk) != 1)
        return("must provide exactly one bulk sample name")

    if (!is.null(object@region)) {
        if (!is(object@region, 'GRanges')) 
            stop('"region" must be a GRanges object')
        # This isn't possible if we're to combine chunks for the final
        # pipeline steps. Need a check that allows either single region
        # GRanges or the expected full autosome set.
        #
        # Disabling for now, but this isn't safe in general.
        #if (length(object@region) != 1) {
           #stop('"region" must contain exactly 1 range')
        #}
    }

    if (!is.null(object@gatk)) {
        if (ncol(object@gatk) != 13) {
            return("@gatk must have 13 columns")
        }
        if (!all(colnames(object@gatk)[1:7] ==
                c('chr', 'pos', 'refnt', 'altnt', 'mq', 'mqrs'))) {
            return("@gatk is improperly formatted")
        }
    }

    if (!is.null(object@ab.fits)) {
        # XXX: should also check that relevant chromosomes are present
        if (!all(colnames(object@ab.fits == c('a', 'b', 'c', 'd', 'logp'))))
            return("@ab.fits is improperly formatted")
    }

    if (!is.null(object@ab.estimates)) {
        if (!all(c('ab', 'gp.mu', 'gp.sd') %in% colnames(object@gatk)))
            return("ab.estimates were not properly imported")
    }

    if (!is.null(object@mut.models)) {
        if (!all(c('abc.pv', 'lysis.pv', 'lysis.beta', 'mda.pv', 'mda.beta')
            %in% colnames(object@gatk)))
            return("mut.models were not properly calculated")
    }

    # XXX: add validator for rest of slots
    return(TRUE)
})


# various steps in the pipeline require certain slots to be filled.
# this function will ensure all slots in 'slots: (character)' are not
# NULL and will generate an error message otherwise.
check.slots <- function(object, slots, abort=TRUE) {
    error.occurred = FALSE
    for (s in slots) {
        if (is.null(slot(object, s))) {
            error.occurred = TRUE
            if (s == 'gatk')
                cat("must import GATK read counts first (see: read.integrated.table())\n")
            if (s == 'ab.fits')
                cat('must import AB model parameters (see: add.ab.fits())\n')
            if (s == 'ab.estimates')
                cat("must import allele balance model estimates first (see: add.ab.estimates())\n")
            if (s == 'cigar.data')
                cat("must import CIGAR data first (see: add.cigar.data())\n")
            if (s == 'excess.cigar.scores')
                cat("must compute excess CIGAR scores first (see: compute.excess.cigar.scores())\n")
            if (s == 'static.filter.params')
                cat("must apply static site filters first (see: add.static.filters())\n")
            if (s == 'fdr.prior.data')
                cat("must compute or import FDR priors first (see: compute.fdr.prior.data())\n")
        }
    }
    if (error.occurred & abort)
        stop("One or more required slots are missing. See above for details.")
    return(error.occurred)
}


# to use: e.g., x <- do.call(concat, xs)
# where xs is a list of SCAN2 chunks. the chunks must not overlap. one day this
# might be enforced programatically

setGeneric("concat", function(...) standardGeneric("concat"))
setMethod("concat", signature="SCAN2", function(...) {
	args <- list(...)
	if (length(args) == 0) stop("tried to concat() 0 objects")
	if (length(args) == 1) return(args[[1]])
	init <- args[[1]]
	ret <- make.scan(init@single.cell, init@bulk, genome = init@genome.string, region = reduce(do.call(c, lapply(args, function(a) a@region))))
	ret@gatk <- rbindlist(lapply(args, function(a) a@gatk))
	ensure.same <- function(l, slot.name, var.name, var.name2) {
	    if (missing(var.name) & !missing(var.name2)) 
		stop("var.name2 only valid if var.name also specified")
	    if (all(sapply(l, function(element) is.null(slot(element, slot.name))))) {
		return()
	    }
	    if (!missing(var.name)) {
		if (missing(var.name2)) {
		    if (!all(sapply(l, function(element) all(slot(l[[1]], slot.name)[[var.name]] == slot(element, slot.name)[[var.name]])))) 
		        stop(paste("list of SCAN2 chunks cannot be concatenated; slot", slot.name, "does not have consistent values for", var.name))
		}
		else {
		    if (!all(sapply(l, function(element) all(slot(l[[1]], slot.name)[[var.name]][[var.name2]] == slot(element, slot.name)[[var.name]][[var.name2]])))) 
		        stop(paste("list of SCAN2 chunks cannot be concatenated; slot", slot.name, "does not have consistent values for", var.name))
		}
	    }
	    else {
		if (!all(sapply(l, function(element) all(slot(l[[1]], slot.name) == slot(element, slot.name))))) 
		    stop(paste("list of SCAN2 chunks cannot be concatenated; slot", slot.name, "does not have consistent values"))
	    }
	}
	ensure.same(args, "genome.string")
	ensure.same(args, "single.cell")
	ensure.same(args, "bulk")
	for (mt in c("snv", "indel")) {
	    ensure.same(args, "static.filter.params", mt, "min.sc.alt")
	    ensure.same(args, "static.filter.params", mt, "min.sc.dp")
	    ensure.same(args, "static.filter.params", mt, "min.bulk.dp")
	    ensure.same(args, "static.filter.params", mt, "max.bulk.alt")
	    ensure.same(args, "static.filter.params", mt, "max.bulk.af")
	    ensure.same(args, "static.filter.params", mt, "cg.id.q")
	    ensure.same(args, "static.filter.params", mt, "cg.hs.q")
	    ensure.same(args, "static.filter.params", mt, "disable.hs.filter")
	}
	ret@static.filter.params <- init@static.filter.params
	ensure.same(args, "cigar.data", "sc.path")
	ensure.same(args, "cigar.data", "bulk.path")
	ret@cigar.data <- data.frame(sc.sites = sum(sapply(args, function(a) ifelse(is.null(a@cigar.data), 0, a@cigar.data$sc.sites))), bulk.sites = sum(sapply(args, function(a) ifelse(is.null(a@cigar.data), 0, a@cigar.data$bulk.sites))))
	for (mt in c("snv", "indel")) {
	    ensure.same(args, "excess.cigar.scores", mt, "legacy")
	    ensure.same(args, "excess.cigar.scores", mt, "null.sites")
	}
	ret@excess.cigar.scores <- init@excess.cigar.scores
	ret@excess.cigar.scores$snv$sites <- sum(sapply(args, function(a) ifelse(is.null(a@excess.cigar.scores$snv), 0, a@excess.cigar.scores$snv$sites)))
	ret@excess.cigar.scores$indel$sites <- sum(sapply(args, function(a) ifelse(is.null(a@excess.cigar.scores$indel), 0, a@excess.cigar.scores$indel$sites)))
	for (mt in c("snv", "indel")) {
	    ensure.same(args, "fdr.prior.data", mt, "bins")
	    ensure.same(args, "fdr.prior.data", mt, "max.dp")
	    ensure.same(args, "fdr.prior.data", mt, "candidates.used")
	    ensure.same(args, "fdr.prior.data", mt, "hsnps.used")
	    ensure.same(args, "fdr.prior.data", mt, "nt.tab")
	    ensure.same(args, "fdr.prior.data", mt, "na.tab")
	    ensure.same(args, "fdr.prior.data", mt, "mode")
	}
	ret@fdr.prior.data <- init@fdr.prior.data
	ret@fdr <- NULL
	if (any(sapply(args, function(a) !is.null(a@fdr)))) {
	    ret@fdr <- list()
	    for (mt in c("snv", "indel")) {
		ensure.same(args, "fdr", mt, "mode")
		ret@fdr[[mt]] <- data.frame(mode = init@fdr[[mt]]$mode, sites = sum(sapply(args, function(a) ifelse(is.null(a@fdr[[mt]]), 0, a@fdr[[mt]]$sites))))
	    }
	}
	if (any(sapply(args, function(a) any(as.matrix(init@ab.fits) != as.matrix(a@ab.fits))))) stop("@ab.fits must be identical for all concat() elements")
	ret@ab.fits <- init@ab.fits
	ret@ab.estimates <- data.frame(sites = sum(sapply(args, function(a) ifelse(is.null(a@ab.estimates), 0, a@ab.estimates$sites))))
	ret@mut.models <- data.frame(sites = sum(sapply(args, function(a) ifelse(is.null(a@mut.models), 0, a@mut.models$sites))))
	ret
})


read.integrated.table = function(object, path, quiet = FALSE) 
{
    object@gatk <- read.and.annotate.integrated.table(path = path, 
        sample.id = object@single.cell, region = object@region, 
        quiet = quiet)
    object@integrated.table.path <- path

    id <- object@gatk[, paste(chr, pos, refnt, altnt)]
    dupid <- duplicated(id)
    if (sum(dupid)) {
        stop(paste("found", sum(dupid), "duplicate (chr,pos,refnt,altnt) instances in integrated table. This is likely a bug, please report it."))
    }
    object
}


#setGeneric("read.integrated.table", function(object, path, quiet=FALSE)
#    standardGeneric("read.integrated.table"))
#setMethod("read.integrated.table", "SCAN2", function(object, path, quiet=FALSE) {
#    object@gatk <- read.and.annotate.integrated.table(path=path, sample.id=object@single.cell,
#        region=object@region, quiet=quiet)
#    object@integrated.table.path <- path
#
#    id <- object@gatk[,paste(chr, pos, refnt, altnt)]
#    dupid <- duplicated(id)
#    if (sum(dupid)) {
#        stop(paste('found', sum(dupid), 'duplicate (chr,pos,refnt,altnt) instances in integrated table. This is likely a bug, please report it.'))
#    }
#
#    object
#})



# Add AB Gaussian process parameter fits for each chromosome.
# These take a long time to compute, so the snakemake pipeline will almost
# always be necessary here.
setGeneric("add.ab.fits", function(object, path)
    standardGeneric("add.ab.fits"))
setMethod("add.ab.fits", "SCAN2", function(object, path) {
    object@ab.fits <- get(load(path))
    object
})



# Using tile subsampling, AB model fitting can feasibly be performed on
# a single multicore node.
# AB fits must be computed using all training sites on a chromosome.
# AB fits are produced by sampling 20,000 random parameter 4-tuples. The
# highest logP values are considered the best fits. The random sampling
# is iteratively refined in 'n.steps' steps by restricting random parameter
# sampling to a smaller subset of the space.
#   - samples.per.chunk, n.chunks: compute 'samples.per.chunk' random
#         parameter samplings in 'n.chunks' independent (possibly parallel)
#         threads. samples.per.chunk*n.chunks should be kept at 20,000.
#   - refine.n.steps: number of iterations in which the (a,b,c,d) parameter space
#         is refined and sampled.
#   - refine.top.n: use the top 'top.n' parameter values (by logP) to create the
#         next refined parameter space
#   - n.tiles and hsnp.tilesize: the log-likelihood of (a,b,c,d|training hSNPs)
#         is approximated by breaking the hSNPs into non-overlapping tiles of
#         size hsnp.tilesize. This is necessary because the approximation requires
#         inverting the (#hSNPs x #hSNPs) covariance matrix. Furthermore, the
#         additional information about (a,b,c,d) provided by each tile of hSNPs
#         becomes less and less as more tiles are added. We have found that ~200
#         tiles (=20,000 hSNPs) is a good trade-off between compute time and
#         accuracy.
#   - alim, blim, clim, dlim - starting bounds for the (a,b,c,d) parameter
#         space.
setGeneric("compute.ab.fits", function(object, path, chroms, #=1:22,
    n.cores=future::availableCores(),
    logp.samples.per.step=20000, refine.n.steps=4, refine.top.n=50,
    n.tiles=250, hsnp.tilesize=100,
    alim=c(-7, 2), blim=c(2, 4), clim=c(-7, 2), dlim=c(2, 6))
    standardGeneric("compute.ab.fits"))
setMethod("compute.ab.fits", "SCAN2", function(object, path, chroms, #=1:22,
    n.cores=future::availableCores(),
    logp.samples.per.step=20000, refine.n.steps=4, refine.top.n=50,
    n.tiles=250, hsnp.tilesize=100,
    alim=c(-7, 2), blim=c(2, 4), clim=c(-7, 2), dlim=c(2, 6)) 
{
    cat("using", n.cores, "cores\n")
    n.chunks <- 100  # 4*n.cores  # using a multiple of n.cores gives a smoother progress bar
    # using n.cores and ceiling can lead to (very small) differences in number
    # of samples taken, and thus slightly different results. just use a very large
    # number of chunks (like 100)
    if (logp.samples.per.step %% n.chunks != 0)
        stop(sprintf('logp.samples.per.step must be a multiple of n.chunks (%d)', n.chunks))
    samples.per.chunk <- ceiling(logp.samples.per.step / n.chunks)

    # Check all chroms up front so the loop doesn't die after a significant amount of work
    not.in <- chroms[!(chroms %in% seqnames(object@genome.seqinfo))]
    if (length(not.in) > 0) {
        cat("the following chromosomes are not recognized:\n")
        print(not.in)
        cat("valid chromosomes names for genome", object@genome.string, 'are:\n')
        print(seqnames(object@genome.seqinfo))
        stop('invalid chromosomes, see above for details')
    }

    chrom.refine.records <- setNames(lapply(chroms, abmodel.fit.one.chrom,
        path=path, sc.sample=object@single.cell,
        genome.seqinfo=object@genome.seqinfo,
        hsnp.tilesize=hsnp.tilesize, n.tiles=n.tiles,
        refine.n.steps=refine.n.steps, n.chunks=n.chunks,
        n.logp.samples.per.chunk=samples.per.chunk), chroms)

    chrom.refine.records
})



# Actually calculate the AB estimates.
# Currently, computes AB at ALL SITES rather than just somatic candidates.
# This will be useful in the future for mosaics and perhaps for using more
# germline hSNPs to model what somatic mutations should look like.
#
# IMPORTANT: AB estimation benefits greatly from chunking. However,
# AB estimation uses a window of 100kb up and downstream from every site at
# which AB is being estimated. For sites at the edge of each chunk, data
# either upstream or downstream will not be available. To solve this,
# training data must be read in AGAIN with the 100kb flanking regions added.
setGeneric("compute.ab.estimates", function(object, quiet=FALSE)
    standardGeneric("compute.ab.estimates"))
setMethod("compute.ab.estimates", "SCAN2", function(object, quiet=FALSE) {
    check.slots(object, c('gatk', 'ab.fits'))

    training.hsnps <- get.training.sites.for.abmodel(object=object,
        region=object@region, integrated.table.path=object@integrated.table.path, quiet=quiet)

    ab <- compute.ab.given.sites.and.training.data(
        sites=object@gatk[,.(chr,pos,refnt,altnt)],
        training.hsnps=training.hsnps,
        ab.fits=object@ab.fits, quiet=quiet)

    # Add an extra 'ab' column that transforms gp.mu (-Inf, Inf) -> [0, 1]
    object@gatk[, c('ab', 'gp.mu', 'gp.sd') := 
        list(1/(1+exp(-ab[,'gp.mu'])), ab[,'gp.mu'], ab[,'gp.sd'])]
    object@ab.estimates <- data.frame(sites=nrow(ab))
    object
})



setGeneric("compute.models", function(object, verbose=TRUE)
    standardGeneric("compute.models"))
setMethod("compute.models", "SCAN2", function(object, verbose=TRUE) {
    check.slots(object, c('gatk', 'ab.estimates'))

    if (nrow(object@gatk) > 0) {
        matched.gp.mu <- match.ab(af=object@gatk$af, gp.mu=object@gatk$gp.mu)
        pvb <- compute.pvs.and.betas(object@gatk$scalt, object@gatk$dp,
                                    matched.gp.mu, object@gatk$gp.sd, verbose=verbose)
    } else {
        pvb <- data.table(abc.pv=numeric(0), lysis.pv=numeric(0), lysis.beta=numeric(0), mda.pv=numeric(0), mda.beta=numeric(0))
    }

    object@gatk[, c('abc.pv', 'lysis.pv', 'lysis.beta', 'mda.pv', 'mda.beta') := pvb]
    object@mut.models <- data.frame(sites=nrow(pvb))
    object
})


setGeneric("compute.fdr.prior.data", function(object, mode=c('legacy', 'new'), quiet=FALSE)
    standardGeneric("compute.fdr.prior.data"))
setMethod("compute.fdr.prior.data", "SCAN2", function(object, mode=c('legacy', 'new'), quiet=FALSE)
{
	check.slots(object, c("gatk", "static.filter.params"))
	check.chunked(object, "compute.fdr.prior must be called on a SCAN2 object containing all sites, not a chunked parallelized object")
	mode <- match.arg(mode)
	if (!("static.filter" %in% colnames(object@gatk))) stop("must compute static filters before running compute.fdr.prior.data")
	muttypes <- c("snv", "indel")
	object@fdr.prior.data <- c(setNames(lapply(muttypes, function(mt) {
	    if (mode == "legacy") {
		sfp <- object@static.filter.params[["snv"]]
		cand <- object@gatk[muttype == mt & balt <= sfp$max.bulk.alt & bulk.gt == "0/0" & scalt >= sfp$min.sc.alt & dp >= sfp$min.sc.dp & bulk.dp >= sfp$min.bulk.dp & (is.na(balt.lowmq) | balt.lowmq <= sfp$max.bulk.alt)]
		hets = object@gatk[muttype == mt & training.site == TRUE & scalt >= sfp$min.sc.alt]
	    } else {
		sfp <- object@static.filter.params[[mt]]
		cand <- object@gatk[muttype == mt & somatic.candidate == TRUE & scalt >= sfp$min.sc.alt & dp >= sfp$min.sc.dp]
		hets <- object@gatk[muttype == mt & training.site == TRUE & scalt >= sfp$min.sc.alt & dp >= sfp$min.sc.dp & bulk.dp >= sfp$min.bulk.dp]
	    }
	    c(compute.fdr.prior.data.for.candidates(candidates = cand, hsnps = hets, random.seed = 0, quiet = quiet), mode = mode)
	}), muttypes), list(mode = mode))
	object
})


setGeneric("compute.fdr", function(object, path, mode=c('legacy', 'new'), quiet=FALSE)
    standardGeneric("compute.fdr"))
setMethod("compute.fdr", "SCAN2", function(object, path, mode=c('legacy', 'new'), quiet=FALSE)
{
    mode <- match.arg(mode)
    check.slots(object, c("gatk", "ab.estimates", "mut.models"))
    if (!missing(path) & !is.null(slot(object, "fdr.prior.data"))) 
        stop("fdr.prior.data is already loaded; path to new fdr.prior.data would overwrite")
    if (!missing(path)) 
        object@fdr.prior.data <- get(load(path))
    check.slots(object, "fdr.prior.data")
    muttypes <- c("snv", "indel")
    object@fdr <- c(setNames(lapply(muttypes, function(mt) {
        object@gatk[muttype == mt, `:=`(c("nt", "na"), estimate.fdr.priors(.SD, 
            object@fdr.prior.data[[mt]], use.ghet.loo = FALSE))]
        object@gatk[muttype == mt & training.site == TRUE & scalt >= 
            object@static.filter.params[[mt]]$min.sc.alt, `:=`(c("nt", 
            "na"), estimate.fdr.priors(.SD, object@fdr.prior.data[[mt]], 
            use.ghet.loo = TRUE))]
        if (mode == "legacy") {
            sfp <- object@static.filter.params[["snv"]]
            object@gatk[muttype == mt & (resampled.training.site | 
                (balt == 0 & bulk.gt == "0/0" &
                  scalt >= sfp$min.sc.alt & dp >= sfp$min.sc.dp)), 
                `:=`(c("lysis.fdr", "mda.fdr"), compute.fdr.legacy(altreads = scalt, 
                  dp = dp, gp.mu = match.ab(af = af, gp.mu = gp.mu), 
                  gp.sd = gp.sd, nt = nt, na = na, verbose = !quiet))]
            sites <- nrow(object@gatk[muttype == mt & balt == 
                0 & bulk.gt == "0/0" & scalt >= 
                sfp$min.sc.alt & dp >= sfp$min.sc.dp])
        } else if (mode == "new") {
            object@gatk[muttype == mt, `:=`(c("lysis.fdr", "mda.fdr"), 
                list(lysis.pv * na/(lysis.pv * na + lysis.beta * 
                  nt), mda.pv * na/(mda.pv * na + mda.beta * 
                  nt)))]
            sites <- nrow(object@gatk[muttype == mt])
        }
        list(sites = sites)
    }), muttypes), list(mode = mode))
    object
})


setGeneric("add.static.filter.params", function(object, config.path,
    muttype=c('snv', 'indel'), min.sc.alt=2, min.sc.dp=6,
    max.bulk.alt=0, min.bulk.dp=11, exclude.dbsnp=TRUE, cg.id.q=0.05, cg.hs.q=0.05)
        standardGeneric("add.static.filter.params"))
setMethod("add.static.filter.params", "SCAN2",
function(object, config.path, muttype=c('snv', 'indel'),
    min.sc.alt=2, min.sc.dp=6, max.bulk.alt=0, min.bulk.dp=11,
    exclude.dbsnp=TRUE, cg.id.q=0.05, cg.hs.q=0.05)
{
    if (!missing(muttype) & !missing(config.path))
        stop('exactly one of muttype or config.path must be specified')

    muttype <- match.arg(muttype)

    # if a config file is given, it is expected to have explicit values for all settings
    if (!missing(config.path)) {
        yaml <- yaml::read_yaml(config.path)
        for (mt in c('snv', 'indel')) {
            object@static.filter.params[[mt]] <- list(
                min.sc.alt=yaml[[paste0(mt, '_min_sc_alt')]],
                min.sc.dp=yaml[[paste0(mt, '_min_sc_dp')]],
                max.bulk.alt=yaml[[paste0(mt, '_max_bulk_alt')]],
                max.bulk.af=yaml[[paste0(mt, '_max_bulk_af')]],
                min.bulk.dp=yaml[[paste0(mt, '_min_bulk_dp')]],
                # exclude.dbsnp, cg.id.q and cg.hs.q are not user configurable
                # at the moment. the defaults in this function's signature are
                # always used.
                exclude.dbsnp=exclude.dbsnp,
                cg.id.q=cg.id.q,
                cg.hs.q=cg.hs.q
            )
        }
    } else {
        object@static.filter.params[[muttype]] <- list(
            min.sc.alt=min.sc.alt, min.sc.dp=min.sc.dp,
            max.bulk.alt=max.bulk.alt, max.bulk.af=1,  # ignored, cede control to max.bulk.alt
            min.bulk.dp=min.bulk.dp,
            exclude.dbsnp=exclude.dbsnp, cg.id.q=cg.id.q, cg.hs.q=cg.hs.q)
    }
    object
})


# legacy mode - id.score/hs.score both must be strictly > the id/hs.score.q cutoffs.
#               this was later found to be problematic for non-BWA MEM aligners that
#               generate no such CIGAR ops, causing all points to have the same scores.
#               in that case, allowing >= cutoff effectively disables the filter, which
#               is necessary.

compute.static.filters = function(object, mode=c('new', 'legacy'))
{
    check.slots(object, c("gatk", "cigar.data", "mut.models"))
    mode <- match.arg(mode)
    object@static.filter.params$mode <- mode
    for (mt in c("snv", "indel")) {
        sfp <- object@static.filter.params[[mt]]
        object@gatk[muttype == mt, `:=`(c("cigar.id.test", "cigar.hs.test", 
            "lowmq.test", "dp.test", "abc.test", "min.sc.alt.test", 
            "max.bulk.alt.test", "max.bulk.af.test", "csf.test"), list(id.score >= object@excess.cigar.scores[[mt]]$id.score.q, 
            hs.score >= object@excess.cigar.scores[[mt]]$hs.score.q, 
            is.na(balt.lowmq) | balt.lowmq <= sfp$max.bulk.alt, 
            dp >= sfp$min.sc.dp & bulk.dp >= sfp$min.bulk.dp, 
            abc.pv > 0.05, scalt >= sfp$min.sc.alt, balt <= sfp$max.bulk.alt, 
            bulk.af <= sfp$max.bulk.af, muttype == 'snv' | max.out <= 2))]
        if (mode == "legacy") {
            object@gatk[muttype == mt, `:=`(c("cigar.id.test", 
                "cigar.hs.test"), list(id.score > object@excess.cigar.scores[[mt]]$id.score.q, 
                hs.score > object@excess.cigar.scores[[mt]]$hs.score.q))]
        }
        object@gatk[muttype == mt, `:=`(static.filter, cigar.id.test & 
            cigar.hs.test & lowmq.test & dp.test & abc.test & 
            min.sc.alt.test & max.bulk.alt.test & max.bulk.af.test & csf.test)]
    }
    object

}

#setGeneric("compute.static.filters", function(object, mode=c('new', 'legacy'))
#        standardGeneric("compute.static.filters"))
#setMethod("compute.static.filters", "SCAN2", function(object, mode=c('new', 'legacy')) {
#    check.slots(object, c('gatk', 'cigar.data', 'mut.models'))
#    mode <- match.arg(mode)
#
#    object@static.filter.params$mode <- mode
#
#    for (mt in c('snv', 'indel')) {
#        sfp <- object@static.filter.params[[mt]]
#        object@gatk[muttype == mt, c('cigar.id.test', 'cigar.hs.test', 'lowmq.test',
#                'dp.test', 'abc.test', 'min.sc.alt.test',
#                'max.bulk.alt.test', 'max.bulk.af.test',
#                'dbsnp.test', 'csf.test') :=
#                    list(id.score >= object@excess.cigar.scores[[mt]]$id.score.q,
#                        hs.score >= object@excess.cigar.scores[[mt]]$hs.score.q,
#                        is.na(balt.lowmq) | balt.lowmq <= sfp$max.bulk.alt,
#                        dp >= sfp$min.sc.dp & bulk.dp >= sfp$min.bulk.dp,
#                        abc.pv > 0.05,
#                        scalt >= sfp$min.sc.alt,
#                        balt <= sfp$max.bulk.alt,
#                        bulk.af <= sfp$max.bulk.af,
#                        !sfp$exclude.dbsnp | dbsnp == '.',
#                        muttype == 'snv' | unique.donors <= 1 | max.out <= 2)]
#
#        if (mode == 'legacy') {
#            object@gatk[muttype == mt, c('cigar.id.test', 'cigar.hs.test') := 
#                    list(id.score > object@excess.cigar.scores[[mt]]$id.score.q,
#                        hs.score > object@excess.cigar.scores[[mt]]$hs.score.q)]
#        }
#
#        object@gatk[muttype == mt, static.filter :=
#            cigar.id.test & cigar.hs.test & lowmq.test & dp.test &
#            abc.test & min.sc.alt.test & max.bulk.alt.test & max.bulk.af.test &
#            dbsnp.test & csf.test]
#    }
#    object
#})


# Do not change any static filter parameter manually if any of fdr.priors, fdr or
# call.mutations have been run. The output of these routines depends on the
# initial set of candidates (the fraction of true positives within this set is
# estimated), so changing the candidate set (by changing the static filters) will
# change the result.
#
# This function is intended for those who wish to experiment with changing
# parameters; it is not feasible to completely rerun the SCAN2 tool, or even the
# final call_mutations step many 10s or 100s of times to explore the param space.
#
# N.B. if you want to update target.fdr, this is not the function you are looking
#      for. target.fdr can be updated without the cumbersome recalculation of FDR
#      priors/FDR rates. Use call.mutations().
#
# `new.params` - must be a 2 element list, named 'snv' and 'indel', each of which
#                are also lists. These lists may contain named members corresponding
#                to already existing static filter parameters.
#
# XXX.mode - override the mode used for stage XXX. If not specified, the original
#            mode (recorded in the SCAN2 object) will be used.
#   - `compute.static.filters.mode`
#   - `fdr.prior.mode`
#   - `fdr.mode`
setGeneric("update.static.filter.params", function(object, fdr.prior.mode, fdr.mode, compute.static.filters.mode, new.params=list(snv=list(), indel=list()), quiet=FALSE)
    standardGeneric("update.static.filter.params"))
setMethod("update.static.filter.params", "SCAN2", function(object, fdr.prior.mode, fdr.mode, compute.static.filters.mode, new.params=list(snv=list(), indel=list()), quiet=FALSE) {
    if (!all(names(new.params) %in% c('snv', 'indel')))
        stop('new.params must be a list containing at least one list named "snv" or "indel"')

    for (mt in c('snv', 'indel')) {
        sfp <- object@static.filter.params[[mt]]
        new <- new.params[[mt]]
        for (p in names(new)) {
            if (p %in% names(sfp)) {
                if (new[p] != sfp[[p]]) {
                    if (!quiet) cat(paste0("    (",mt,") updating ", p, ": ", sfp[[p]], " -> ", new[[p]], '\n'))
                    object@static.filter.params[[mt]][[p]] <- new[[p]]
                }
            } else {
                if (!quiet) cat(paste0("    (",mt,") ignoring unrecognized parameter ", p, '\n'))
            }
        }
    }

    # Update somatic candidate loci and recompute static filters
    sfp <- object@static.filter.params
    annotate.gatk.candidate.loci(object@gatk,
        snv.min.bulk.dp=sfp$snv$min.bulk.dp,
        snv.max.bulk.alt=sfp$snv$max.bulk.alt,
        snv.max.bulk.af=sfp$snv$max.bulk.af,
        indel.min.bulk.dp=sfp$indel$min.bulk.dp,
        indel.max.bulk.alt=sfp$indel$max.bulk.alt,
        indel.max.bulk.af=sfp$indel$max.bulk.af)
    object <- compute.static.filters(object,
        mode=ifelse(missing(compute.static.filters.mode), object@static.filter.params$mode, compute.static.filters.mode))

    if (!is.null(slot(object, 'fdr.prior.data'))) {
        if (!quiet) cat("    voiding N_T/N_A FDR prior calculations\n")
        old.mode <- object@fdr.prior.data$mode
        object@fdr.prior.data <- NULL
        object@gatk[, nt := NA]
        object@gatk[, na := NA]
        object <- compute.fdr.prior.data(object,
            mode=ifelse(missing(fdr.prior.mode), old.mode, fdr.prior.mode),
            quiet=quiet)
    }
    if (!is.null(slot(object, 'fdr'))) {
        if (!quiet) cat("    voiding per-locus FDR calculations\n")
        old.mode <- object@fdr$mode
        object@fdr <- NULL
        object@gatk[, lysis.fdr := NA]
        object@gatk[, mda.fdr := NA]
        object <- compute.fdr(object,
            mode=ifelse(missing(fdr.mode), old.mode, fdr.mode),
            quiet=quiet)
    }
    if (!is.null(slot(object, 'call.mutations'))) {
        if (!quiet) cat("    recalling all mutations\n")
        object@gatk[, pass := NA]
        object@gatk[, training.pass := NA]
        object <- call.mutations(object)
    }
    if (!is.null(slot(object, 'mutburden'))) {
        if (!quiet) cat("    recalling all mutations\n")
        object@mutburden <- NULL
        object <- compute.mutburden(object)
    }
    if (!is.null(slot(object, 'mutsig.rescue'))) {
        if (!quiet) cat("    voiding mutation signature-based rescue calls\n")
        object@mutsig.rescue <- NULL
        object@gatk[, rescue := NA]
    }
    if (!is.null(slot(object, 'spatial.sensitivity'))) {
        if (!quiet) cat("    voiding spatial sensitivity estimates\n")
        object@spatial.sensitivity <- NULL
    }

    object 
})

# When reading the integrated table: read all metadata columns
# and only the 3 genotype/count columns corresponding to `sample.id`. This can
# save significant memory overhead in large (100+ cell) projects.  This optimized
# reading is handled by read.integrated.table.1sample().
#
# Next, annotate the integrated table with information corresponding to
# `sample.id`. This includes training site definitions (which varies from cell
# to cell depending on whether there is no data at the site in that single cell;
# gt=./.) and assigning single cell ref/alt read counts to phased haplotypes.

read.and.annotate.integrated.table <- function(path, sample.id, region=NULL, quiet=FALSE, ncores = 10) {

	tr <- read.integrated.table.1sample(path, sample.id, region = region, quiet = quiet)
	
	if(nrow(tr) == 0){
		return(tr)
	}
	
	setindex(tr, resampled.training.site)
	tr[, `:=`(dp, scalt + scref)]
	tr[, `:=`(af, scalt/dp)]
	data.table::setkey(tr, chr, pos, refnt, altnt)
	setindex(tr, muttype)
	sc.gt <- tr[[sample.id]]

	# Modification: Instead of properly phased snps, hsnps are chosen based on their af and depth using binomial distribution

	# Remove homozygous
	tokeep1 = which(!(tr$balt == 0 | tr$bulk.gt == '1/1' | tr$bulk.gt == '0/0'))

	require(parallel)
	tokeep1_perc <- quantile(tokeep1, probs = seq(0, 1, by = 0.025))
	indexes <- as.numeric( sapply(tokeep1_perc, function(p) which.min(abs(tokeep1 - p))) )
	indexes = tokeep1[indexes]
	
	resL = mclapply(tokeep1, mc.cores = ncores, function(x){
			if(x %in% indexes){print(paste0('Choosing likely hSNPs without phasing: ', round(100 * (x / max(tokeep1))), '%'))}
			res = binom.confint(tr[x,bulk.dp] * 0.5, tr[x,bulk.dp], 0.99, method = 'bayes')
			interval_check = tr[x,bulk.af] > res$lower & tr[x,bulk.af] < res$upper
			nan_check = is.nan(res$lower) | is.nan(res$upper)
			!(nan_check) & interval_check
			})

	if(length(resL) == 0){
		return(NULL)
	} else{
	    tokeep2 = tokeep1[which(unlist(resL))]
	    tr[, `:=`(hetSNP, .I %in% tokeep2)]
	    tr[, `:=`(training.site, hetSNP == T & sc.gt != "./." & bulk.gt != 
		"./.")]
	    tr[, `:=`(hetSNP, NULL)]
	    return(tr)
	}

}


read.training.hsnps <- function(path, sample.id, region=NULL, quiet=FALSE) {
    read.and.annotate.integrated.table(path=path, sample.id=sample.id, region=region, quiet=quiet)[training.site == TRUE & muttype == 'snv']
}


# After all site-specific information is computed, somatic mutations can
# be called.
#
# Somatic mutation calling involves determining FDR prior distributions
# (which are defined by the set of candidate sites), calculating FDR
# heuristics and finally passing sites based on a set of filters.
#
# The most important thing to note about calling mutations is
# that calls depend on the set of candidate sites defined by static
# filter parameters. That is, if you change something like the minimum alt
# reads per single cell, you will in turn change the set of possible
# candidate sites and need to rerun this function.
#
# N.B. because call.mutations must not be called on chunked SCAN2 objects,
# mode=legacy is not computationally feasible.
setGeneric("call.mutations", function(object, target.fdr=0.01, quiet=FALSE)
        standardGeneric("call.mutations"))
setMethod("call.mutations", "SCAN2", function(object, target.fdr=0.01, quiet=FALSE) {
    check.slots(object, c('gatk', 'static.filter.params', 'mut.models', 'excess.cigar.scores', 'fdr.prior.data', 'fdr'))

    # Try to handle indel calling when the cross-sample filter does not meet
    # the usual requirements (>1 unique individual).
    #   1. If there is only one unique individual but several cells from that
    #      individual, then just discard indels shared across cells. With a
    #      few cells present this likely is a decent filter.
    #   2. If there is only one cell from one individual, indel calling should
    #      not be attempted.
    suppress.shared.indels <- max(object@gatk$unique.donors, na.rm=TRUE) == 1
    # <= 2: bulks are actually counted in unique.cells, so the case where
    # there is only one cell from one individual have cause max(unique.cells)=2
    suppress.all.indels <- max(object@gatk$unique.cells, na.rm=TRUE) <= 2

    # Pass somatic mutations
    object@gatk[, pass := static.filter == TRUE &
        (muttype == 'snv' | (!suppress.all.indels & (!suppress.shared.indels | unique.cells == 1))) &
        lysis.fdr <= target.fdr & mda.fdr <= target.fdr]

    # Pass germline heterozygous sites using L-O-O for sensitivity estimation
    object@gatk[training.site == TRUE,
        training.pass := 
            # same tests as static.filter EXCEPT lowmq.test, max.bulk.alt.test,
            # which mostly will fail at training het germline sites.
            (cigar.id.test & cigar.hs.test & dp.test & abc.test & min.sc.alt.test) &
            lysis.fdr <= target.fdr & mda.fdr <= target.fdr]

    # NAs are especially present in legacy output where not all sites
    # have CIGAR data, not all sites have FDR values and sites not in
    # the cross sample panel (which were originally removed by merge())
    # may have NAs (though they shouldn't; 0-fills should be used).
    object@call.mutations <- c(
        as.list(unlist(object@gatk[muttype == 'snv',
            .(snv.pass=as.integer(sum(pass, na.rm=TRUE)),
              snv.resampled.training.pass=as.integer(sum(resampled.training.site & training.pass, na.rm=TRUE)))])),
        as.list(unlist(object@gatk[muttype == 'indel',
            .(indel.pass=as.integer(sum(pass, na.rm=TRUE)),
              indel.resampled.training.pass=as.integer(sum(resampled.training.site & training.pass, na.rm=TRUE)))])),
        list(target.fdr=target.fdr,
            suppress.shared.indels=suppress.shared.indels,
            suppress.all.indels=suppress.all.indels))
    object
})


setGeneric("add.depth.profile", function(object, depth.path)
        standardGeneric("add.depth.profile"))
setMethod("add.depth.profile", "SCAN2", function(object, depth.path) {
    vars <- load(depth.path)  # should contain two objects: dptab and clamp.dp
    if (!all(c('dptab', 'clamp.dp') %in% vars))
        stop('depth.profile RDA file expected to contain dptab and clamp.dp objects')

    object@depth.profile <- list(
        dptab=dptab,
        clamp.dp=clamp.dp
    )
    object
})


setGeneric("is.compressed", function(object) standardGeneric("is.compressed"))
setMethod("is.compressed", "SCAN2", function(object) {
    is.raw(object@gatk)
})


setGeneric("compress", function(object) standardGeneric("compress"))
setMethod("compress", "SCAN2", function(object) {
    if (!is.compressed(object)) {
        object@gatk <- qs::qserialize(object@gatk)
    } else {
        warning('object already compressed, skipping compression')
    }
    object
})


setGeneric("decompress", function(object) standardGeneric("decompress"))
setMethod("decompress", "SCAN2", function(object) {
    if (is.compressed(object)) {
        object@gatk <- qs::qdeserialize(object@gatk)
    } else {
        warning('object already decompressed, skipping decompression')
    }
    object
})
