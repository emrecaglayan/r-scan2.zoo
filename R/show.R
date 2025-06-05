setMethod("show", "SCAN2", function(object) {
    cat("#", is(object)[[1]], "\n")
    if (!is.null(object@region)) {
        cat("#   Region:")
        if (length(object@region) > 1) {
            cat("\n")
            print(object@region)
        } else {
            cat('',length(object@region),'intervals\n')
        }
    }
    cat("#   Genome:", object@genome.string, "\n")
    cat("#   Single cell ID:", object@single.cell, "\n")
    cat("#   Bulk ID:", object@bulk, "\n")

    show.gatk(object)
    show.abmodel(object)
    show.mut.models(object)
    show.cigar.data(object)
    show.static.filters(object)
    show.fdr.prior(object)
    show.depth.profile(object)
    show.call.mutations(object)
    show.mutburden(object)
    show.mutsig.rescue(object)
    show.spatial.sensitivity(object)
})

setGeneric("show.gatk", function(object) standardGeneric("show.gatk"))
setMethod("show.gatk", "SCAN2", function(object) {
    cat("#   GATK:")
    if (is.null(object@gatk)) {
        cat(" (no data)\n")
    } else {
        if (is.compressed(object)) {
            cat(' (compressed)\n')
        } else {
            cat('', nrow(object@gatk), "raw sites\n")
        }
    }
})


setGeneric("show.abmodel.training.sites", function(object) standardGeneric("show.abmodel.training.sites"))
setMethod("show.abmodel.training.sites", "SCAN2", function(object) {
    cat("#   AB model training hSNPs:")
    if (is.compressed(object)) {
        cat(' (table compressed)\n')
    } else {
        if (!('training.site' %in% colnames(object@gatk))) {
            cat(" (no data)\n")
        } else if (nrow(object@gatk[training.site==TRUE]) == 0) {
            cat(" (training data removed, likely minimized SCAN2 object)\n")
        } else {
            # germline indels are not used for AB model training
            #per.hap <- object@gatk[training.site==TRUE & muttype == 'snv', .N, by=phased.gt]
            #tdata <- object@gatk[training.site==TRUE & muttype=='snv']
            #cat('', nrow(tdata),
            #    sprintf("phasing: %s=%d, %s=%d\n",
            #        per.hap$phased.gt[1], per.hap$N[1],
            #        per.hap$phased.gt[2], per.hap$N[2]))
            neighbor.approx <- approx.abmodel.covariance(object, bin.breaks=10^(0:5))
            cors <- round(neighbor.approx$observed.cor, 3)
                cat('#       OBSERVED VAF correlation between neighboring hSNPs:\n')
            cat('#           <10 bp', cors[1], '<100 bp', cors[2],
                '<1000 bp', cors[3], '<10 kbp', cors[4], '<100 kbp', cors[5], '\n')
            if ('resampled.training.site' %in% colnames(object@gatk)) {
                cat('#        ', nrow(object@gatk[resampled.training.site == TRUE & muttype == 'snv']),
                    'resampled hSNPs\n')
                cat('#        ', nrow(object@gatk[resampled.training.site == TRUE & muttype == 'indel']),
                    'resampled hIndels\n')
            }
        }
    }
})
    

setGeneric("show.abmodel.params", function(object) standardGeneric("show.abmodel.params"))
setMethod("show.abmodel.params", "SCAN2", function(object) {
    cat("#   AB model parameters:")
    if (is.null(object@ab.fits)) {
        cat(" (no data)\n")
    } else {
        cat(sprintf('\n#       average (over chromosomes): a=%0.3f, b=%0.3f, c=%0.3f, d=%0.3f\n',
            mean(object@ab.fits$a),
            mean(object@ab.fits$b),
            mean(object@ab.fits$c),
            mean(object@ab.fits$d)))
    }
})


setGeneric("show.abmodel.ab.distn", function(object) standardGeneric("show.abmodel.ab.distn"))
setMethod("show.abmodel.ab.distn", "SCAN2", function(object) {
    cat("#   Allele balance:")
    if (is.compressed(object)) {
        cat(' (table compressed)\n')
    } else { 
        if (is.null(object@ab.estimates)) {
            cat(" (not computed)\n")
        } else {
            s <- summary(object@gatk$gp.sd)
            cat('\n#       mean (0 is neutral):',
                round(mean(object@gatk$gp.mu), 3), '\n')
            cat('#       uncertainty (Q25, median, Q75):',
                round(s['1st Qu.'], 3),
                round(s['Median'], 3),
                round(s['3rd Qu.'], 3), '\n')
            if ('training.site' %in% colnames(object@gatk) & nrow(object@gatk[training.site==TRUE]) > 0) {
                xs <- round(object@gatk[training.site==TRUE & muttype == 'snv',
                    .(mean=mean(gp.mu), cor=cor(af, ab, use='complete.obs'))],3)
                cat('#       mean at training hSNPs:', xs$mean, '\n')
            }
        }
    }
})


setGeneric("show.abmodel", function(object) standardGeneric("show.abmodel"))
setMethod("show.abmodel", "SCAN2", function(object) {
    show.abmodel.training.sites(object)
    show.abmodel.params(object)
    show.abmodel.ab.distn(object)
})

    
setGeneric("show.mut.models", function(object) standardGeneric("show.mut.models"))
setMethod("show.mut.models", "SCAN2", function(object) {
    cat("#   Mutation models:")
    if (is.null(object@mut.models)) {
        cat(" (not computed)\n")
    } else {
        cat(" computed\n")
    }
})


setGeneric("show.cigar.data", function(object) standardGeneric("show.cigar.data"))
setMethod("show.cigar.data", "SCAN2", function(object) {
    cat("#   CIGAR data:")
    if (is.null(object@cigar.data)) {
        cat(" (no data)\n")
    } else {
        cat(' single cell:', object@cigar.data$sc.sites, "sites, ")
        cat('bulk:', object@cigar.data$bulk.sites, "sites\n")
    }
})


setGeneric("show.static.filters", function(object) standardGeneric("show.static.filters"))
setMethod("show.static.filters", "SCAN2", function(object) {
    cat("#   Static filters: ")
    if (is.compressed(object)) {
        cat(' (table compressed)\n')
    } else {
        if (!('static.filter' %in% colnames(object@gatk))) {
            cat("(not applied)\n")
        } else {
            if ('mode' %in% names(object@static.filter.params))  # supporting older versions of SCAN2 objects
                cat(paste0('mode=', object@static.filter.params$mode))
            cat('\n')
            na.or.val <- function(x, val=0) ifelse(is.na(x), val, x)
            for (mt in c('snv', 'indel')) {
                tb <- table(object@gatk[muttype == mt, static.filter], useNA='always')
                cat(sprintf('#       %6s: %8d retained %8d removed %8d NA\n',
                    mt, na.or.val(tb['TRUE']), na.or.val(tb['FALSE']), na.or.val(tb['NA'])))
            }
        }
    }
})


setGeneric("show.fdr.prior", function(object) standardGeneric("show.fdr.prior"))
setMethod("show.fdr.prior", "SCAN2", function(object) {
    cat("#   FDR prior data: ")
    if (is.null(object@fdr.prior.data)) {
        cat("(not computed)\n")
    } else {
        if ('mode' %in% names(object@fdr.prior.data))  # supporting older versions of SCAN2 objects
            cat(paste0('mode=', object@fdr.prior.data$mode))
        cat('\n')
        #for (mt in names(object@fdr.prior.data)) {
        for (mt in c('snv', 'indel')) {
            cat(sprintf("#       %6s: %8d candidates %8d germline hets %8d max burden\n",
                mt, object@fdr.prior.data[[mt]]$candidates.used,
                object@fdr.prior.data[[mt]]$hsnps.used,
                object@fdr.prior.data[[mt]]$burden[2]))
        }
    }
})


setGeneric("show.depth.profile", function(object) standardGeneric("show.depth.profile"))
setMethod("show.depth.profile", "SCAN2", function(object) {
    cat("#   Depth profile: ")
    if (is.null(object@depth.profile)) {
        cat("(not added)\n")
    } else {
        mc=mean.coverage(object)
        cat(sprintf('\n#   mean depth: single cell: %0.2fx,   bulk: %0.2fx (max DP set to %d)\n',
            mc['single.cell'], mc['bulk'], object@depth.profile$clamp.dp))
        for (mt in c('snv', 'indel')) {
            sfp <- object@static.filter.params[[mt]]
            dptab <- object@depth.profile$dptab
            cat(sprintf("#       %6s: genome <:   min. bulk DP %0.1f%%,   min. sc DP %0.1f%%,   either %0.1f%%\n",
                mt,
                100*sum(dptab[, 1:sfp$min.bulk.dp])/sum(dptab),
                100*sum(dptab[1:sfp$min.sc.dp,])/sum(dptab),
                100*(1 - (sum(dptab[(sfp$min.sc.dp+1):nrow(dptab),(sfp$min.bulk.dp+1):ncol(dptab)])/sum(dptab)))))
        }
    }
})


setGeneric("show.call.mutations", function(object) standardGeneric("show.call.mutations"))
setMethod("show.call.mutations", "SCAN2", function(object) {
    cat("#   Somatic mutation calls: ")
    if (is.null(object@call.mutations)) {
        cat("(not called)\n")
    } else {
        cat(sprintf("target.fdr=%0.3f\n",
            object@call.mutations$target.fdr))
        for (mt in c('snv', 'indel')) {
            cat(sprintf("#       %6s: %8d called %8d resampled training calls\n",
                mt,
                as.integer(object@call.mutations[[paste0(mt, '.pass')]]),
                as.integer(object@call.mutations[[paste0(mt, '.resampled.training.pass')]])))
        }
        if (object@call.mutations$suppress.all.indels) {
            cat(sprintf("#       ALL indel calls have been suppressed (insufficient number of single cells in cross-sample filter\n"))
        } else if (object@call.mutations$suppress.shared.indels) {
            cat(sprintf("#       indel calls shared between cells have been suppressed (insufficient number of unique individuals in cross-sample filter\n"))
        }
    }
})


setGeneric("show.mutburden", function(object) standardGeneric("show.mutburden"))
setMethod("show.mutburden", "SCAN2", function(object) {
    cat("#   Somatic mutation burden: ")
    if (is.null(object@mutburden)) {
        cat("(not computed)\n")
    } else {
        cat('\n')
        for (mt in c('snv', 'indel')) {
            mb <- object@mutburden[[mt]][2,]  # row 2 is middle 50%
            cat(sprintf("#       %6s: %6d somatic,   %0.1f%% sens,   %0.3f callable Gbp,   %0.1f muts/haploid Gbp,   %0.1f muts per genome%s%s\n",
                mt, mb$ncalls, 100*mb$callable.sens, mb$callable.bp/1e9, mb$rate.per.gb, mutburden(object, muttype=mt),
                ifelse(any(mb$unsupported.filters), ' (INVALID: static.filter.params)', ''),
                ifelse(mt=='indel' & (object@call.mutations$suppress.all.indels | object@call.mutations$suppress.shared.indels), ' (INVALID: cross-sample panel insufficient)', '')
                ))
        }
    }
})


setGeneric("show.mutsig.rescue", function(object) standardGeneric("show.mutsig.rescue"))
setMethod("show.mutsig.rescue", "SCAN2", function(object) {
    cat("#   Mutation rescue by signature: ")
    if (is.null(object@mutsig.rescue)) {
        cat("(not computed)\n")
    } else {
        # XXX: assumes SNV and indel use the same FDR.  currently correct, may break later
        cat(sprintf('rescue.target.fdr=%0.3f\n',
            object@mutsig.rescue[['snv']]$rescue.target.fdr))
        for (mt in c('snv', 'indel')) {
            msr <- object@mutsig.rescue[[mt]]
            cat(sprintf("#       %6s: %6d/%d candidates rescued,   %0.1f%% rel. error,   sig. weights:  %0.3f true,   %0.3f artifact\n",
                mt, nrow(object@gatk[muttype == mt & rescue]), msr$nmuts,
                100*msr$relative.error, msr$weight.true, msr$weight.artifact))
        }
    }
})



setGeneric("show.spatial.sensitivity", function(object) standardGeneric("show.spatial.sensitivity"))
setMethod("show.spatial.sensitivity", "SCAN2", function(object) {
    cat("#   Spatial sensitivity: ")
    ss <- object@spatial.sensitivity
    if (is.null(ss)) {
        cat("(not computed)\n")
    } else {
        cat(sprintf("applies only to VAF-based calling (NOT mutsig rescue!)\n#       tiles of width %d, neighborhood size for training hSNPs %d (%d tiles)\n",
            ss$tilewidth, ss$tilewidth * (1 + 2*ss$neighborhood.tiles),
            ss$neighborhood.tiles))
        ss <- object@spatial.sensitivity$somatic.sensitivity
        cat(sprintf("#                  either allele              maj allele               min allele\n"))
        cat(sprintf("#              germ.   predicted (SD)   germ.   predicted (SD)   germ.   predicted (SD)\n"))
        for (mt in c('snv', 'indel')) {
            cat(sprintf("#       %6s: %2.1f%%   %2.1f%% +/- %2.1f   %2.1f%%   %2.1f%% +/- %2.1f   %2.1f%%   %2.1f%% +/- %2.1f\n",
                mt,
                100*ss[[mt]]['germline.training', 'both'],
                100*ss[[mt]]['predicted.somatic.mean', 'both'],
                100*ss[[mt]]['predicted.somatic.stdev', 'both'],
                100*ss[[mt]]['germline.training', 'maj'],
                100*ss[[mt]]['predicted.somatic.mean', 'maj'],
                100*ss[[mt]]['predicted.somatic.stdev', 'maj'],
                100*ss[[mt]]['germline.training', 'min'],
                100*ss[[mt]]['predicted.somatic.mean', 'min'],
                100*ss[[mt]]['predicted.somatic.stdev', 'min']))
        }

        # Don't show these yet. They're experimental and will likely only confuse users.
        #b <- object@spatial.sensitivity$burden
        #cat(sprintf("#      EXPERIMENTAL! equal-allele burden genome-wide estimates:\n"))
        #for (mt in c('snv', 'indel')) {
            #cat(sprintf("#       %6s: %0.1f (95%% CI: %0.1f-%0.1f)\n",
                #mt, b[[mt]]$median, b[[mt]]$ci95[1], b[[mt]]$ci95[2]))
        #}
    }
})
