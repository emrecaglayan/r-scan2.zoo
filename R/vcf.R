# ALL DEAD CODE --- keeping around for possible later use to generate VCF files, which users tend to want.


# Create VCF output for a single sample called by SCAN-SNV
# output.fmt may contain a single '%s' which will be replaced by
# the sample name.
scansnv.to.vcf <- function(ss.dir, output.fmt, type='somatic', muttype='snv', overwrite=FALSE) {
    if (!(type %in% c('somatic', 'mosaic', 'hsnp_spikein')))
        stop(sprintf("type must be either somatic, mosaic or hsnp_spikein, not '%s'", type))

    if (!(muttype %in% c('snv', 'indel')))
        stop(spritnf("muttype must be either 'snv' or 'indel', not %s", muttype))

    ss.config <- file.path(ss.dir, "scan.yaml")
    if (!file.exists(ss.config))
        stop(sprintf("expected SCAN-SNV config file does not exist: %s\n",
            ss.config))

    yaml <- yaml::read_yaml(ss.config)
    sc.samples <- names(yaml$sc_bams)
    for (s in sc.samples) {
        path.fmt <- "%s_genotypes.rda"
        if (muttype == 'indel' & type == 'somatic')
            path.fmt <- "%s_genotypes.pon_filter.rda"
        f <- file.path(ss.dir, muttype, s, sprintf(path.fmt, type))
        print(f)
        load(f)
        # we assume the loaded variable is called 'somatic' below
        # but the mosaic results are called 'mosaic'.
        # just stick it in a variable named somatic anyway. it doesn't matter.
        somatic <- get(ifelse(type == 'hsnp_spikein', 'spikeins', type))

        out.file <- sprintf(output.fmt, s)
        scansnv.df.to.vcf(df=somatic, out.file=out.file, yaml=yaml,
            sample.name=s, overwrite=overwrite)
    }
}



# Write out a results data frame to out.file
scansnv.df.to.vcf <- function(df, out.file, ss.config, yaml, sample.name,
    overwrite=FALSE, chrs=c(1:22, 'X', 'Y')) {
    if (!missing(yaml) & !missing(ss.config))
        stop('only one of "ss.config" or "yaml" can be specified')

    if (!missing(ss.config))
        yaml <- yaml::read_yaml(ss.config)

    if (!missing(yaml) | !missing(ss.config)) {
        ref.genome <- yaml$ref
        if (missing(chrs))
            chrs <- yaml$chrs
        # read the FASTA index for the reference
        fai <- read.table(paste0(ref.genome, '.fai'), sep='\t',
            stringsAsFactors=F)
    }

    if (!overwrite & file.exists(out.file))
        stop(sprintf("output file %s already exists, please delete it first",
            out.file))
    f <- file(out.file, 'w')
    cat(sprintf("writing to %s..\n", out.file))

    vcf.header <- c(
        '##fileformat=VCFv4.0',
        '##source=SCAN2.ZOO',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    if (!missing(ss.config) | !missing(yaml))
        vcf.header <- c(vcf.header, 
            sprintf('##reference=%s', yaml$ref),
            sprintf('##contig=<ID=%s,length=%d>', fai[,1], fai[,2]))
    vcf.header <- c(vcf.header,
            paste(c("#CHROM", "POS", 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
                'INFO', 'FORMAT', sample.name), collapse='\t'))

    writeLines(vcf.header, con=f)
    s <- df[!is.na(df$pass) & !is.na(df$chr) & df$pass,]
    s <- do.call(rbind, lapply(chrs, function(chr) {
        ss <- s[s$chr==chr,]
        ss[order(ss$pos),]
    }))
    if (nrow(s) > 0) {
        writeLines(paste(s$chr, s$pos, s$dbsnp, s$refnt, s$altnt,
            '.', 'PASS', '.', 'GT', '0/1', sep='\t'),
            con=f)
    }
    close(f)
}

