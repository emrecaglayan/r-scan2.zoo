# use sigprofilermatrixgenerator to classify mutations
# note that sigprofilermatrixgenerator always generates the most granular classification
# (e.g., SBS6144 for SNVs and ID415 for indels)
# mostly useful for indel classification, but also good for doing TSB on SNVs
classify.muts <- function(df, genome.string, spectype='SNV',
    sample.name='dummy', save.plot=F, auto.delete=T, verbose=FALSE)
{
    if (nrow(df) == 0)
        return(df)

    recognized.spectypes <- c('SNV', 'ID')
    if (!(spectype %in% recognized.spectypes))
        stop(sprintf("unrecognized spectype '%s', currently only supporting %s",
            spectype, paste('"', recognized.spectypes, '"', collapse=', ')))

    require(SigProfilerMatrixGeneratorR)
    spmgd <- tempfile()  # For multithreaded workflows, it is CRITICAL that library(future)
                         # supply different random seeds to child processes.
    if (file.exists(spmgd))
        stop(paste('temporary directory', spmgd, 'already exists'))

    dir.create(spmgd, recursive=TRUE)

    # Write out the VCF
    out.file <- paste0(spmgd, '/', sample.name, '.vcf')
    f <- file(out.file, "w")
    vcf.header <- c("##fileformat=VCFv4.0", "##source=SCAN2", 
"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">", 
        paste(c("#CHROM", "POS", "ID", "REF", "ALT", 
            "QUAL", "FILTER", "INFO", "FORMAT", sample.name), collapse = "\t"))
    writeLines(vcf.header, con = f)

    # SigProfilerMatrixGenerator doesn't classify duplicated mutations
    # from the same sample, it throws errors instead. It also will not
    # detect duplicates if they are not adjacent in the file.  If the
    # duplicate is not detected in this way, it causes the bug where
    # the final newdf dataframe does not match the input df.
    # To circumvent all these headaches: just remove duplicates up front.
    s <- df
    mutid <- paste(s$chr, s$pos, s$refnt, s$altnt)
    dupmut <- duplicated(mutid)
    if (verbose)
        cat("Removing", sum(dupmut), "/", nrow(s), "duplicated mutations before annotating\n")
    s <- s[!dupmut,]

    # will write, eg., position=7000000 as 7e6, which will
    # confuse sigprofilermatrixgenerator
    old.opt <- options('scipen')$scipen
    options(scipen=10000)
    writeLines(paste(s$chr, s$pos, '.', s$refnt, s$altnt, 
        ".", "PASS", ".", "GT", "0/1", sep = "\t"), con = f)
    close(f)
    options(scipen=old.opt)

    if (verbose) {
        mat <- SigProfilerMatrixGeneratorR::SigProfilerMatrixGeneratorR(sample.name, genome.string, spmgd, seqInfo=TRUE, plot=save.plot)
    } else {
        # Prevent sigprofilermatrixgenerator's output from being printed
        # XXX: should probably do some error handling here
        reticulate::py_capture_output(
            mat <- SigProfilerMatrixGeneratorR::SigProfilerMatrixGeneratorR(sample.name, genome.string, spmgd, seqInfo=TRUE, plot=save.plot))
    }

    # Read in the types
    # SigProfilerMatrixGenerator output files have chromosomes named 1..22 without
    # the 'chr' prefix, even if the prefix was present in the input vcf.
    chrs.to.read <- sub('chr', '', s$chr[!duplicated(s$chr)])
    annot.files <- paste0(spmgd, '/output/vcf_files/', spectype, '/', chrs.to.read, "_seqinfo.txt")
    if (spectype == 'ID') {
        colclasses <- c(V2='character', V5='character', V6='character')
    } else if (spectype == 'SNV') {
        colclasses <- c(V2='character')
    }

    annots <- do.call(rbind, lapply(annot.files, function(f) {
        tryCatch(x <- read.table(f, header=F, stringsAsFactors=FALSE,
                colClasses=colclasses),
            error=function(e) NULL)
    }))

    # SigProfilerMatrixGenerator always strips leading 'chr' prefixes. Add
    # them back if the input had them. Column 2 is chromosome (1 is sample name).
    if (substr(s$chr[1], 1, 3) == 'chr')
        annots[[2]] <- paste0('chr', annots[[2]])

    if (spectype == 'ID') {
        colnames(annots) <- c('sample', 'chr', 'pos', 'iclass', 'refnt', 'altnt', 'unknown')
        newdf <- plyr::join(df, annots[2:6], by = colnames(annots)[-c(1, 4, 7)])
    } else if (spectype == 'SNV') {
        colnames(annots) <- c('sample', 'chr', 'pos', 'iclass', 'unknown')
        newdf <- plyr::join(df, annots[2:4], by = colnames(annots)[-c(1, 4, 5)])
    }

    if (save.plot) {
        plotfiles <- list.files(paste0(spmgd, '/output/plots/'), full.names=T)
        file.copy(plotfiles, '.')
    }

    if (!all(df$chr == newdf$chr))
        stop('df and newdf do not perfectly correspond: df$chr != newdf$chr')
    if (!all(df$pos == newdf$pos))
        stop('df and newdf do not perfectly correspond: df$pos != newdf$pos')

    if (auto.delete)
        unlink(spmgd, recursive=TRUE)
    df$muttype <- newdf$iclass
    df
}

genome.to.spmgr.format = c("GRCh37", "GRCh38", "mm10", "rheMac10", "panTro6", "calJac4")
names(genome.to.spmgr.format) = c('hs37d', 'hg38', 'mm10', 'rheMac10', 'panTro6', 'calJac4')


# new version: just returns the vector of indel classes
# FORCES ID83 FOR NOW
# N.B.: don't provide a default genome.string; it's dangerous.
classify.indels <- function(df, genome.string, sample.name='dummy', save.plot=F, auto.delete=T, verbose=FALSE) {
    # SigProfilerMatrixGenerator returns ID415 by default, which is ID83 plus one of
    # 5 transcription strand states: B, N, Q, T, U. The format of the string is, e.g.,
    #    U:1:Del:T:1
    # Removing the first two characters "U:" leaves an ID83 type.
    id415 <- classify.muts(df=df, genome.string=genome.to.spmgr.format[genome.string],
        spectype='ID', sample.name=sample.name, save.plot=save.plot,
        auto.delete=auto.delete, verbose=verbose)$muttype
    # id83() converts strings into a factor
    id83(substr(id415, 3, nchar(id415[1])))
}
