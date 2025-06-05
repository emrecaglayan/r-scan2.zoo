# given two data frames of somatic and germline locations, annotate
# the somatic data frame with the position of the nearest germline entry.
find.nearest.germline <- function (som, germ, chrs) {
    som$nearest.het <- NA
    for (chr in chrs) {
        gpos <- germ$pos[germ$chr == chr]
        spos <- som$pos[som$chr == chr]
        gidx <- findInterval(spos, gpos)
        gidx[gidx==0] <- 1  # somatic is to the left of lowest germline, so
                            # the lowest germline is the nearest
        nearest.idx <- ifelse(abs(gpos[gidx] - spos) <= abs(gpos[gidx + 
            1] - spos), gidx, gidx + 1)
        som$nearest.het[som$chr == chr] <- gpos[nearest.idx]
    }
    som
}

# d is a vector of distances to the nearest hSNP
get.distance.distn <- function(d, min=1, max=5) {
    h <- hist(d[d >= min & d <= max], breaks=50, plot=FALSE)
    h$density <- h$density / sum(h$density)
    h
}


convert_to_na <- function(arg) {
  if (arg == "NA") {
    return(NA)
  } else {
    return(arg)
  }
}


retry_genome_string_to_tiling <- function(genome, tilewidth = 1e+07, group = "auto", max_attempts = 5) {
  attempts <- 0
  success <- FALSE
  result <- NULL
  
  while (!success && attempts < max_attempts) {
    attempts <- attempts + 1
    tryCatch({
      result <- genome.string.to.tiling(genome, tilewidth = tilewidth, group = group)
      success <- TRUE
    }, error = function(e) {
      cat("Attempt", attempts, "failed: ", e$message, "\n")
    })
  }
  
  if (!success) {
    stop("Failed to run genome.string.to.tiling after ", max_attempts, " attempts")
  }
  
  return(result)
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

get.3mer <- function(df, chr, pos, refnt, altnt,
    genome=BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5,
    verbose=FALSE)
{
    if (!missing(df)) {
        chr <- df$chr
        pos <- df$pos
        refnt <- df$refnt
        altnt <- df$altnt
    } else {
        if (missing(chr) | missing(pos) | missing(refnt) | missing(altnt))
            stop('must supply either "df" or all of chr, pos, refnt, alt')
    }
    
    muttype <- muttype.map[paste0(refnt, '>', altnt)]

    tmp <- Biostrings::getSeq(genome, GRanges(seqnames=chr,
                              ranges=IRanges(start=pos-1, end=pos+1)))

    ctx <- as.character(tmp)
    ctx.rc <- as.character(Biostrings::reverseComplement(tmp))
    type.and.ctx <- paste0(ifelse(refnt == 'C' | refnt == 'T', ctx, ctx.rc), ':', muttype)
    # sbs96 converts the string form into a factor to enforce ordering
    sbs96(type.and.ctx)
}


# Extremely simple heuristic for improving short-range phasing
# Ties caused by the abs() > abs() statement produce a bias in
# adjusted AFs.
# 10kb is probably much larger than the correlation signal
# caries, but it doesn't actually matter. 
#
# IMPORTANT: adjust.phase does NOT perform accurate long-range
# phasing. By setting the phase to the most agreeable VAF of
# the nearest hSNP, the phase for small ranges of VAFs near
# 0.5 will likely be incorrect; however, these small phase
# swaps are unlikely to affect the mutation models.
# The intention is to resolve regions with high imbalance and
# to thereby prevent AB "ping-ponging" between ABs near 0 and 1.
#
# N.B. with enough single cells, a method that harmonizes VAFs
# for multiple cells simultaneously would perform accurate long
# range phasing.
#
# 'pht' is a data.table that is modified by reference.
adjust.phase <- function(pht, dist.cutoff=1e4, quiet=FALSE) {
    if (!quiet) cat("WARNING: adjust.phase is EXPERIMENTAL!\n")

    pos <- pht$pos
    dp <- pht$phased.hap1 + pht$phased.hap2
    # copy will be overwritten as the sapply below runs
    af <- pht$phased.hap1 / dp
    cat(sprintf("adjust.phase: neighboring hSNP correlation before adjustment: %0.3f\n",
        cor(af[-1], af[-nrow(pht)], use='complete.obs')))

    # this can't be by independently evaluating each site. after evaluation of one
    # site, the next site's results may change.
    swapped <- rep(F, length(af))
    sapply(2:length(af), function(i) {
        swap.i <- pos[i] - pos[i-1] < dist.cutoff &
                !is.na(af[i-1]) & !is.na(af[i]) &
                abs(af[i-1]-af[i]) > abs(af[i-1]-(1-af[i]))
        swapped[i] <<- swap.i
        if (swap.i)
            af[i] <<- 1-af[i]
    })

    pht[, c('phased.hap1', 'phased.hap2', 'phase.adjusted') :=
        list(ifelse(swapped, phased.hap2, phased.hap1),
             ifelse(swapped, phased.hap1, phased.hap2),
             swapped)]

    if (!quiet) {
        cat(sprintf("adjust.phase: neighboring hSNP correlation after adjustment: %0.3f\n",
            cor(af[-1], af[-nrow(pht)], use='complete.obs')))
    }
}
