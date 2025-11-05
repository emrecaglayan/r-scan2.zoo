# Wrappers and expected formats for read.table.1sample
gatk.meta.cols <- c(
    chr='character',
    pos='integer',
    refnt='character',
    altnt='character'
)
read.gatk.1sample <- function(path, sample.id, region=NULL, quiet=FALSE) {
    read.table.1sample(path=path, sample.id=sample.id, region=region,
        meta.cols=gatk.meta.cols, quiet=quiet)
}

integrated.table.meta.cols <- c(
    chr='character',
    pos='integer',
    refnt='character',
    altnt='character',
    mq = 'numeric',
    mqrs = 'numeric',
    bulk.gt='character',
    bref='integer',
    balt='integer',
    bulk.dp='integer',
    bulk.af='numeric',
#    bulk.binom.prob='numeric',
    tref='integer',
    talt='integer',
#    tcells='integer',
    muttype='character',
    mutsig='character',
    balt.lowmq='integer',
    nalleles='integer',
    unique.donors='integer',
    unique.cells='integer',
    unique.bulks='integer',
    max.out='integer',
    sum.out='integer',
    sum.bulk='integer',
    somatic.candidate='logical',
#    training.site='logical',
    resampled.training.site='logical'
)

read.integrated.table.1sample = function(path, sample.id, region = NULL, quiet = FALSE) 
{
    read.table.1sample(path = path, sample.id = sample.id, region = region, 
        meta.cols = integrated.table.meta.cols, quiet = quiet)
}


# Some extra work to make sure we only read in the part of the
# table relevant to one sample. Otherwise, memory can become
# an issue for projects with 10s-100s of cells.
#
# meta.cols - columns that should be read in addition to columns
#             specific to `sample.id`.
#
# region can be a GRanges object with a single interval to read only
# a subset of the GATK table. The table is tabix indexed, so this can
# be done quickly.

read.table.1sample <- function(path, sample.id, meta.cols, region=NULL, quiet=FALSE){
    if (!quiet) cat("Importing GATK table..\n")

	if(quiet == F){
		print(meta.cols)
	}
	
    n.meta.cols <- length(meta.cols)

    # Step 1: just get the header and detect the columns corresponding to sample.id
    tf <- Rsamtools::TabixFile(path)
    open(tf)
    header <- read.tabix.header(tf)
    close(tf)
    col.strings <- strsplit(header, '\t')[[1]]
    tot.cols <- length(col.strings)
    sample.idx <- which(col.strings == sample.id)

    if (!quiet) {
        cat("Selecting columns:\n")
        for (i in 1:length(col.strings)) {
            if (i <= n.meta.cols) {
                cat(sprintf("    (%d)", i), col.strings[i], '\n')
            } else if (any(i == sample.idx + 0:2)) {
                cat(sprintf("    (%d)", i), col.strings[i], '[single cell]\n')
            }
        }
    }
    
    # remove colClasses option
    gatk <- read.tabix.data(path=path, region=region, header=header, quiet=quiet)

    new.sample.idx <- which(colnames(gatk) == sample.id)
    colnames(gatk)[new.sample.idx+1:2] <- c('scref', 'scalt')

    # Rearrange columns so that the single cell triplet is first, then bulk triplet
    cols.to.keep <- col.strings[1:n.meta.cols]
    cols.to.keep <- c(cols.to.keep, sample.id, c('scref', 'scalt'))

    gatk <- gatk[,..cols.to.keep]
    gatk$chr = as.character(gatk$chr)

    gatk
}


# annotate `gatk.meta` with 'bulk.gt', 'bulk.dp', 'bulk.af', 'bref' and 'balt' taken from `gatk`,
# which are the bulk genotype string assigned by HaplotypeCaller, bulk depth, bulk VAF, the number
# of ref supporting bulk reads and mutation supporting bulk reads.
#
# also add the total sum of alt/ref reads across all single cells in the batch.
#
# sc.samples - list of all sample IDs supplied to SCAN2 as --sc-bam arguments
#
# legacy - SCANSNV's legacy behavior was to sum alt-supporting reads from ALL BAMs
#          except the named bulk sample to determine somatic candidate sites.
#          This can cause problems if the user supplies additional single cells
#          not for primary analysis (like large sets of low-depth cells) or
#          additional bulks that are not the primary bulk comparison (i.e., bulks
#          from other tissues).
#
# only `gatk.meta` is modified (by reference).
annotate.gatk.counts <- function(gatk.meta, gatk, bulk.sample, sc.samples, legacy=FALSE, quiet=FALSE) {
    col.strings <- colnames(gatk)
    tot.cols <- length(col.strings)
    bulk.idx <- Inf
    if (!missing(bulk.sample))
        bulk.idx <- which(col.strings == bulk.sample)

    if (!quiet) {
        cat("Selecting columns:\n")
        for (i in 1:length(col.strings)) {
            if (any(i == bulk.idx + 0:2)) {
                cat(sprintf("    (%d)", i), col.strings[i], '[bulk]\n')
            }
        }
    }

    gatk.meta[, c('bulk.gt', 'bref', 'balt') :=
        list(gatk[[bulk.idx]], gatk[[bulk.idx+1]], gatk[[bulk.idx+2]])]
    gatk.meta[, c('bulk.dp', 'bulk.af') := list(bref+balt, balt/(bref+balt))]

    refs <- c()
    alts <- c()
    if (legacy) {
        # legacy SCANSNV behavior used every alt read count
        alts <- which(colnames(gatk) == 'alt')
        refs <- which(colnames(gatk) == 'ref')
    } else {
        #sample.col.idxs <- seq(8, ncol(gatk), 3)    # all sample column IDs
        # adding the bulk sample so that the `- bref` and `- balt` calculations below can
        # be applied whether legacy mode is chosen or not.
        sample.col.idxs <- which(colnames(gatk) %in% c(bulk.sample, sc.samples))

        refs <- sample.col.idxs + 1
        alts <- sample.col.idxs + 2
    }
    # data.table is very finicky about accessing variables in the calling scope.
    # especially bad when data.tables are nested in data.table formulae.
    tref.var <- rowSums(as.matrix(gatk[, ..refs])) - gatk.meta$bref
    talt.var <- rowSums(as.matrix(gatk[, ..alts])) - gatk.meta$balt

    # Annotate the total number of single cell alt and ref reads across all SC samples.
    gatk.meta[, c('tref', 'talt') := list(tref.var, talt.var)]
    gatk.meta
}


# `gatk` is a data.table, so all of these updates happen by reference.
# no need to return the result.
#
# This can be slow with add.mutsig, particularly for indels. Best used
# on chunked data.
annotate.gatk <- function(gatk, genome.string, add.mutsig=TRUE) {
    data.table::setkey(gatk, chr, pos, refnt, altnt)

    # Determine SNV/indel status and then annotate mutation signature channels
    gatk[, muttype := ifelse(nchar(refnt) == 1 & nchar(altnt) == 1, 'snv', 'indel')]
    # allow for fast selection of SNVs or indels
    setindex(gatk, muttype)

    if (add.mutsig) {
        # This is the only place in scan2 where BSgenome is necessary.
        # This is important because BSgenome has several dependencies (in
        # particular, rtracklayer) that uses ~350 MB of RAM over the other
        # scan2 dependencies. This mem usage is multiplied when using
        # library(future), plan(multicore) for parallelization so is quite
        # significant.
        #
        # N.B. using the integrated table workflow, we should only ever
        # have to annotate these values once.
        require(BSgenome)
        
        gatk[muttype == 'snv',
            mutsig := get.3mer(chr=chr, pos=pos, refnt=refnt, altnt=altnt,
                               genome=genome.string.to.bsgenome.object(genome.string))]
	       
        gatk_indel <- gatk[muttype == 'indel'] 	
        if (nrow(gatk_indel) != 0) {
          chs_list <- lapply(1:nrow(gatk_indel), function(i) {
  	    chr=gatk_indel$chr[i]
  	    pos=gatk_indel$pos[i]
  	    refnt=gatk_indel$refnt[i]
  	    altnt=gatk_indel$altnt[i]
            if(i%%100==0){print(paste0('Number of variants processed: ', i)); flush.console()}
   	    chs=attribution_of_indels(genome=genome.string.to.bsgenome.object(genome.string),
  			in_CHROM=chr, in_POS=pos, in_REF=refnt, in_ALT=altnt) })
          gatk[muttype == 'indel', mutsig := unlist(chs_list)] 
        }  
     }
}


# 'gatk' is a data.table to be modified by reference
annotate.gatk.lowmq <- function(gatk, path, bulk, region, quiet=FALSE) {
    lowmq <- read.gatk.1sample(path, sample.id=bulk, region=region, quiet=quiet)
    data.table::setkey(lowmq, chr, pos, refnt, altnt)

    gatk[lowmq, on=.(chr,pos,refnt,altnt), balt.lowmq := i.scalt]
}


# 'gatk' is a data.table to be modified by reference
# 'phasing.path' points to a phased VCF output by either SHAPEIT or EAGLE
# with sample column 'phasedgt'. The VCF file must have a header line
# beginning with #CHROM, as is usual for the VCF spec.
#
# IMPORTANT: this VCF contains no information about single cells, only bulk.
#
# to allow the very simple table join to work well, we assume phasing was
# restricted to biallelic sites for this individual. while not ideal, the
# performance of such phasing has been good enough.
#
# To annotate phasing status (phased.hap1, phased.hap2) the estimated
# phase needs to be joined to a single cell table with alt and ref read
# counts. Based on phase, (alt,ref) maps to either (hap1,hap2) or (hap2,hap1).
annotate.gatk.phasing <- function(gatk, phasing.path, region, quiet=FALSE) {
    # VCF column 1 should be named "CHROM"
    phase.data <- read.tabix.data(path=phasing.path, region=region, quiet=quiet,
        colClasses=list(character='CHROM'))

    if (ncol(phase.data) != 10)
        stop('expected single-sample VCF format with 10 columns')

    # phasing.path is a single sample standard VCF, so we specify the column
    # format explicitly.
    colnames(phase.data) <- c('chr', 'pos', 'dbsnp', 'refnt', 'altnt', 'qual', 'filter', 'info', 'format', 'phasedgt')

    # This assumes "GT" is the first element of the GT string format, which isn't
    # guaranteed but is the case for our data.

    unrecognized <- setdiff(unique(phase.data$phasedgt), c('1|0', '0|1', './.'))
    if (length(unrecognized) > 0)
        stop(paste('phasing genotypes expected to be either 0|1, 1|0 or ./., but found', unrecognized, collapse='\n'))

    # First join the phase genotype (string is either 0|1, 1|0 or ./. if no call)
    gatk[phase.data, on=.(chr,pos,refnt,altnt), phased.gt := i.phasedgt]
}


# if 'panel.path' is NULL, then columns with dummy counts of 0 (so as to
# not throw off filters) will be joined. this is easily discernable from
# real data since unique.donors can never be <1 for any non-ref site in
# this sample (the donor would include this one).
annotate.gatk.panel <- function(gatk, panel.path, region=NULL, quiet=FALSE) {
    if (!is.null(panel.path)) {
        panel <- read.tabix.data(path=panel.path, region=region, quiet=quiet,
            colClasses=list(character='chr'))  # force chromosome to be interpreted as a string

        # update by reference
        # N.B. only join by chromosome and position here. Because the panel often
        # involves 50+ cells and 10+ individuals, it is common for sites that
        # would be biallelic in this sample to become multiallelic when considering
        # all samples.
        # At multiallelic sites, the panel counts actually do not tell us which of
        # the alleles is supported in the other cells/subjects, only that there are
        # non-reference alleles supported in those other cells/subjects. One day
        # we will hopefully avoid this by decomposing alleles via vcfallelicprimitives
        # or some similar tool.
        gatk[panel, on=.(chr,pos),
            c('nalleles', 'unique.donors', 'unique.cells', 'unique.bulks', 'max.out', 'sum.out', 'sum.bulk') :=
                list(i.nalleles, i.unique.donors, i.unique.cells, i.unique.bulks, i.max.out, i.sum.out, i.sum.bulk)]
    } else {
        gatk[,
            c('nalleles', 'unique.donors', 'unique.cells', 'unique.bulks', 'max.out', 'sum.out', 'sum.bulk') :=
                list(0, 0, 0, 0, 0, 0, 0)]
    }
}


# Adds a columns `somatic.candidate` to `gatk` by reference. Loci with
# somatic.candidate=TRUE pass some of the basic requirements to be considered
# a mutation locus, meaning the site may be a candidate in at least one single
# cell.
#
# XXX: any filter not dependent on specific single cell info (like min.bulk.dp)
# should really be applied here.

annotate.gatk.candidate.loci <- function(gatk, snv.min.bulk.dp, snv.max.bulk.alt, snv.max.bulk.af, indel.min.bulk.dp, indel.max.bulk.alt, indel.max.bulk.af, mode = c("new", "legacy"))
{
    mode <- match.arg(mode)
    if (mode == "new") {
        snv.allow.bulk.gt <- snv.max.bulk.alt > 0 | snv.max.bulk.af > 
            0
        indel.allow.bulk.gt <- indel.max.bulk.alt > 0 | indel.max.bulk.af > 
            0
        gatk[, `:=`(somatic.candidate, ((muttype == "snv" & balt <= 
            snv.max.bulk.alt & bulk.dp >= snv.min.bulk.dp & !is.na(bulk.af) & 
            bulk.af <= snv.max.bulk.af & (is.na(balt.lowmq) | 
            balt.lowmq <= snv.max.bulk.alt) & (snv.allow.bulk.gt | 
            bulk.gt == "0/0")) | (muttype == "indel" & balt <= 
            indel.max.bulk.alt & bulk.dp >= indel.min.bulk.dp & 
            !is.na(bulk.af) & bulk.af <= indel.max.bulk.af & 
            (is.na(balt.lowmq) | balt.lowmq <= indel.max.bulk.alt) & 
            (indel.allow.bulk.gt | bulk.gt == "0/0"))) & talt >= 2)]
    }
    else if (mode == "legacy") {
        gatk[, `:=`(somatic.candidate, ((muttype == "snv" & balt <= 
            snv.max.bulk.alt) | (muttype == "indel" & balt <= 
            indel.max.bulk.alt)) & bulk.gt == "0/0" & talt >= 2)]
    }
    else stop("unrecognized mode")
    gatk
}


# Again, modifying 'gatk' by reference.
# Returns list of resampling auxiliary data with one entry for SNVs and one for indels
# and modifies `gatk` by reference.
gatk.resample.phased.sites = function (gatk, M = 20, seed = 0) 
{

	# Modification: Instead of properly phased snps, hsnps are chosen based on their af and depth using binomial distribution

	# Remove homozygous
	tokeep1 = which(!(gatk$balt == 0 | gatk$bulk.gt == '1/1' | gatk$bulk.gt == '0/0'))
	#gatk_filt <- gatk[!(balt == 0 | bulk.gt == '1/1' | bulk.gt == '0/0')]

	require(parallel)
	tokeep1_perc <- quantile(tokeep1, probs = seq(0, 1, by = 0.025))
	indexes <- as.numeric( sapply(tokeep1_perc, function(p) which.min(abs(tokeep1 - p))) )
	indexes = tokeep1[indexes]
	
	resL = mclapply(tokeep1, mc.cores = M, function(x){
			if(x %in% indexes){print(paste0('Choosing likely hSNPs without phasing: ', round(100 * (x / max(tokeep1))), '%'))}
			res = binom.confint(gatk[x,bulk.dp] * 0.5, gatk[x,bulk.dp], 0.99, method = 'bayes')
			interval_check = gatk[x,bulk.af] > res$lower & gatk[x,bulk.af] < res$upper
			nan_check = is.nan(res$lower) | is.nan(res$upper)
			!(nan_check) & interval_check
			})
	tokeep2 = tokeep1[which(unlist(resL))]

	# Subset the data frame to include only high-quality SNPs
	gatk[, hetSNP := .I %in% tokeep2]

	ret <- list()
	for (mt in c("snv", "indel")) {
		aux.data <- resample.germline(sites = gatk[somatic.candidate == T & muttype == mt],
						hsnps = gatk[hetSNP == T & muttype == mt], M = M, seed = seed)
		
		#is_present = paste0(gatk$chr, gatk$pos) %in% paste0(hqSNP$chr, hqSNP$pos)
		gatk[hetSNP == T & muttype == mt, `:=`(resampled.training.site, aux.data$selection$keep)]
		ret[[mt]] <- aux.data
	}
	gatk[is.na(resampled.training.site), `:=`(resampled.training.site, FALSE)]
	gatk[, hetSNP := NULL]
	return(ret)
}


# Write, bgzip and tabix index the integrated table
# 'out.tab' is TEMPORARY. The final file is out.tab.gz.
write.integrated.table <- function(inttab, out.tab,
    out.tab.gz=paste0('out.tab', '.gz'), overwrite=FALSE)
{
    if (file.exists(out.tab) & !overwrite)
        stop(paste('output file', out.tab, 'already exists, please delete it first'))

    if (file.exists(out.tab.gz) & !overwrite)
        stop(paste('output file', out.tab.gz, 'already exists, please delete it first'))

    colnames(inttab)[1] <- paste0('#', colnames(inttab)[1]) # hack to comment out the header
    data.table::fwrite(inttab, file=out.tab, sep='\t', na='NA', quote=FALSE)
    Rsamtools::bgzip(out.tab, out.tab.gz)
    # if we get here, overwrite is either TRUE or the file didn't exist,
    # so there's no danger of deleting a file without warning.
    #unlink(out.tab)
    Rsamtools::indexTabix(file=out.tab.gz, format='vcf', comment='#')
}
