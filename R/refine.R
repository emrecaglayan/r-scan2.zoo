
region_overlap = function(query, subject){

	# Keep only chr, start, end
	query_simp = query[, c('chr', 'pos')]
	query_simp$end = query_simp$pos
	colnames(query_simp) = c('chr', 'start', 'end')
	homop = homop[, 1:3]
	colnames(homop) = c('chr', 'start', 'end')

	# Turn into granges
	query_simp_gr = makeGRangesFromDataFrame(query_simp)
	homop_gr = makeGRangesFromDataFrame(homop)

	# Overlap and return the overlaps in query
	hits = findOverlaps(query_simp_gr, homop_gr)
	return(query[queryHits(hits), ])
}


get_nxp_artifacts = function(passes, sc, meta, panel, stringency = 0.5, cnt_alt = 3, bsgenome, gb){

	# Load artifacts
	passes$id = paste0(passes$chr, '_', passes$pos, '_', passes$refnt, '_', passes$altnt)
	panel_pos = panel[panel$pos %in% passes$pos,]

	panel_pos_passesL = list()
	for(i in 1:nrow(panel_pos)){
		#panel_pos$id = paste0(panel_pos$chr, '_', panel_pos$pos, '_', panel_pos$refnt, '_', panel_pos$altnt)
		panel_pos_i = panel_pos[i,]	
		ids = paste0(panel_pos_i[,'chr'], '_', panel_pos_i[,'pos'], '_', panel_pos_i[,'refnt'], '_', unlist(strsplit(panel_pos_i[,altnt], ',')))
		
		# Extract artifacts from the panel if an ALT allele is also found in passes
		if(sum(ids %in% passes$id) >= 1){
			panel_pos_passesL[[i]] = panel_pos_i
		} else{
			panel_pos_passesL[[i]] = NULL
		}
	}

	panel_pos = as.data.frame(do.call(rbind, panel_pos_passesL))
	panel_pos_num = panel_pos[,7:(ncol(panel)-1)]

	
	meta = read.csv(meta, header = T)
	
	if(sc %in% meta$sample){

		# Get donor name
		donor = meta[meta$sample == sc, 'donor']

		# Exclude samples from the same donor
		donorS = colnames(panel_pos)[grepl(donor, colnames(panel_pos))]
		panel_pos_num_others = panel_pos_num[, !(colnames(panel_pos_num) %in% donorS)]
	} else{
		panel_pos_num_others = panel_pos_num
	}

	# Get minimum number of samples that should have the alt allele if it will be filtered
	sampN = round(ncol(panel_pos_num_others)*stringency)

	# Find number of samples (sc or bulk) with more than X ALT allele at the mutation site
	tmp = apply(panel_pos_num_others, 1, function(x){y = sum(x >= cnt_alt)})

	# Find the sites that meet the above criteria for at least X percent of the samples, depending on stringency
	panel_pos$id0 = paste0(panel_pos$chr, '_', panel_pos$pos)
	passes$id0 = paste0(passes$chr, '_', passes$pos)
	artifacts = passes[passes$id0 %in% panel_pos[which(tmp > sampN), 'id0'], ]
	artifacts
}


get_snv_mpileup_pass = function(samp_passes_snv1, outdir, samp_bulk_bam, tmpdir, tmppref = c('bulkbam', 'snv')){

	# Get sites that have realigned
	sites_filt = paste0(samp_passes_snv1$chr, ':', samp_passes_snv1$pos, '-', samp_passes_snv1$pos)

	sites_gr <- GRanges(
	  seqnames = sub(":.*", "", sites_filt),
	  ranges = IRanges(
	    start = as.integer(sub(".*:(\\d+)-\\d+", "\\1", sites_filt)),
	    end = as.integer(sub(".*-([0-9]+)$", "\\1", sites_filt))
	  )
	)

	param <- Rsamtools::ScanBamParam(which = sites_gr, what = "pos")

	# Read alignments overlapping your positions
	outbam = paste0(outdir, '/Regenotype_', tmppref[1], '_', tmppref[2], '_', samp, '.bam')
	reads <- Rsamtools::scanBam(outbam, param = param)

	# Check which positions had reads
	has_reads <- sapply(reads, function(x) length(x$pos) > 0)

	# Filter to only sites with reads
	covered_sites <- sites_filt[has_reads]
	covered_sites = gsub(':|-', '_', covered_sites)

	# Separate data frame with realigned sites
	all_sites = paste0(samp_passes_snv1$chr, '_', samp_passes_snv1$pos, '_', samp_passes_snv1$pos)
	samp_passes_snv1_realigned = samp_passes_snv1[all_sites %in% covered_sites,]

	## MPILEUP - Regenotype new passes, remove if present in bulk ##
	all_mp = remove_bulk_mpileup(passes_file = samp_passes_snv1,
				bulk_bam_file = samp_bulk_bam,
				outfile = tempfile(tmpdir = tmpdir))

	samp_passes_snv1_realigned2 = samp_passes_snv1_realigned[samp_passes_snv1_realigned$id %in% all_mp$id,]
	ralign_mp = remove_bulk_mpileup(passes_file = samp_passes_snv1_realigned2,
				bulk_bam_file = outbam,
				outfile = tempfile(tmpdir = tmpdir))

	norealign_pass = all_mp[!(all_mp$id %in% samp_passes_snv1_realigned2$id), id]
	realign_pass = ralign_mp[,id]
	bulk_pass = c(norealign_pass, realign_pass)
	return(bulk_pass)
}

regenotype_gatk = function(intervals_file, samp, outdir, fastafile, bamfile, pref){

	outvcf = paste0(outdir, '/Regenotype_', pref[1], '_', pref[2], '_', samp, '.vcf')
	outbam = paste0(outdir, '/Regenotype_', pref[1], '_', pref[2], '_', samp, '.bam')

	if('indel' %in% pref){

		gatk_args <- c(
		  "HaplotypeCaller",
		  "-R", fastafile,
		  "-I", bamfile,
		  "-O", outvcf,
		  "--minimum-mapping-quality", "0",
		  "--max-reads-per-alignment-start", "0",
		  "--dont-use-soft-clipped-bases", "true",
		  "--min-base-quality-score", "0",
		  "--emit-ref-confidence", "GVCF",
		  "--standard-min-confidence-threshold-for-calling", "0",
		  "--pcr-indel-model", "NONE",
		  "--min-pruning", "0",
		  "--interval-padding", "500",
		  "--min-dangling-branch-length", "0",
		  "--allow-non-unique-kmers-in-ref", "true",
		  "--bam-output", outbam,
		  "-L", intervals_file
		)

		# Run the command
		system2("gatk", args = gatk_args)
		indel_regen = read.vcfR(outvcf)
		indel_regen_fix = indel_regen@fix

		# Check if alt allele is present
		alts = strsplit(indel_regen_fix[, 'ALT'], ',')

		regen_ids = lapply(1:length(alts), function(k){paste(indel_regen_fix[k, 'CHROM'],
				indel_regen_fix[k, 'POS'], indel_regen_fix[k, 'REF'], alts[[k]], sep = '_')}) %>% unlist
	} else if('snv' %in% pref){

		gatk_args <- c(
		  "HaplotypeCaller",
		  "-R", fastafile,
		  "-I", bamfile,
		  "-O", outvcf,
		  "--minimum-mapping-quality", "0",
		  "--max-reads-per-alignment-start", "0",
		  "--dont-use-soft-clipped-bases", "true",
		  "--min-base-quality-score", "0",
		  "--emit-ref-confidence", "GVCF",
		  "--standard-min-confidence-threshold-for-calling", "0",
		  "--min-pruning", "0",
		  "--interval-padding", "100",
		  "--min-dangling-branch-length", "0",
		  "--allow-non-unique-kmers-in-ref", "true",
		  "--bam-output", outbam,
		  "-L", intervals_file
		)

		# Run the command
		system2("gatk", args = gatk_args)
		indel_regen = read.vcfR(outvcf)
		indel_regen_fix = indel_regen@fix

		# Check if alt allele is present
		alts = strsplit(indel_regen_fix[, 'ALT'], ',')

		regen_ids = lapply(1:length(alts), function(k){paste(indel_regen_fix[k, 'CHROM'],
				indel_regen_fix[k, 'POS'], indel_regen_fix[k, 'REF'], alts[[k]], sep = '_')}) %>% unlist
	}
	return(regen_ids)
}

decide_cigar_cutoff = function(burdens, sample_name, slope_threshold = -0.3, loess_span = 0.75){

	# This function only accepts a data frame of burdens with varying cigar cutoffs
	# This function is only used when there is a big drop in burden estimate between
	# the most lenient and most stringent cigar cutoffs.
	# The goal of the function is to find the sweet spot so that drop in the mutation burden
	# with increasing cigar stringency is accomplished without losing unnecessary amount of
	# sensitivity.

	cutoffs = burdens$cigar_cutoff
	mutation_burden = rev(burdens$burden)
	mutation_burden_sc = scale(mutation_burden)
	library(scales)
	# Smooth the curve
	loess_fit <- loess(mutation_burden_sc ~ cutoffs, span = loess_span)
	smoothed <- predict(loess_fit)

	loess_fit_plt <- loess(mutation_burden ~ cutoffs, span = loess_span)
	smoothed_plt <- predict(loess_fit_plt)

	# Calculate slope (first derivative)
	slopes <- diff(smoothed) / diff(cutoffs)

	# Define a slope threshold (flattening = near-zero slope)
	flattening_indices <- which(slopes > slope_threshold)

	# Find the first point where flattening starts
	if (length(flattening_indices) > 0) {
	  elbow_index <- flattening_indices[1]
	} else {
	  elbow_index <- which.min(slopes)  # fallback: steepest point
	}
	elbow_cutoff <- cutoffs[elbow_index]

	# Plot
	pdf(paste0('CIGAR_BY_BURDEN_', sample_name, '.pdf'))
	plot(cutoffs, mutation_burden, type = "b", pch = 19, col = "gray",
	     main = "Last Point Before Flattening", xlab = "Cutoff", ylab = "Mutation Burden Z-Score")
	lines(cutoffs, smoothed_plt, col = "blue", lwd = 2)
	points(cutoffs[elbow_index], smoothed_plt[elbow_index], col = "red", pch = 19, cex = 1.5)
	text(cutoffs[elbow_index], smoothed_plt[elbow_index],
	     labels = paste0("Cutoff: ", round(elbow_cutoff, 3)), pos = 3, col = "red")
	dev.off()

	cat("Chosen cigar cutoff:", elbow_cutoff, "\n")
	return(elbow_cutoff)
}

compute_mutburden_new = function(object, gbp.per.genome = get.gbp.by.genome(object), 
    quiet = FALSE, q1 = 0.3, q2 = 0.8) 
{
    check.slots(object, c("call.mutations", "depth.profile"))
    muttypes <- c("snv", "indel")
    object@mutburden <- setNames(lapply(muttypes, function(mt) {

        pre.geno.burden <- object@fdr.prior.data[[mt]]$burden[2]
        sfp <- object@static.filter.params[[mt]]
        g <- object@gatk[resampled.training.site == TRUE & muttype == 
            mt]
        dptab <- object@depth.profile$dptab
        dptab <- dptab[1:min(max(g$dp) + 1, nrow(dptab)), ]
        if (nrow(g) < 100) {
            warning(paste("only", nrow(g), "resampled germline", 
                mt, "sites were detected; aborting genome-wide extrapolation. Typical whole-genome experiments include ~10-100,000 germline sites"))
            ret <- data.frame(ncalls = NA, callable.sens = NA, 
                callable.bp = NA)[c(1, 1, 1), ]
        }
        else {
            s <- object@gatk[pass == TRUE & muttype == mt]
            q = 4
            qstouse <- c(1, 2, 4)
            qbreaks <- quantile(g$dp, prob = c(0, q1, (q1+q2)/2, q2, 1.0))
            if (length(unique(qbreaks)) != q + 1) {
                ncalls <- rep(NA, 3)
                callable.sens <- rep(NA, 3)
                rowqs <- rep(NA, 3)
                callable.bp <- rep(NA, 3)
                warning(paste("could not derive unique breakpoints for quartiles of sequencing depth at germline hSNPs.  this usually indicates that sequencing depth is heavily skewed toward low depths (typically DP=0)\ngot qbreaks = ", 
                  deparse(qbreaks)))
            }
            else {
                s$dpq <- cut(s$dp, qbreaks, include.lowest = T, 
                  labels = F)
                s$dpq[s$dpq == 3] <- 2
                g$dpq <- cut(g$dp, qbreaks, include.lowest = T, 
                  labels = F)
                g$dpq[g$dpq == 3] <- 2
                rowqs <- cut(0:(nrow(dptab) - 1), qbreaks, include.lowest = T, 
                  labels = F)
                rowqs[rowqs == 3] <- 2
                s <- s[dpq %in% qstouse]
                g <- g[dpq %in% qstouse]
                ncalls <- sapply(qstouse, function(q) sum(s[dpq == 
                  q]$pass, na.rm = TRUE))
                callable.bp <- sapply(split(dptab[, -(1:sfp$min.bulk.dp)], 
                  rowqs), sum)
                callable.sens <- sapply(qstouse, function(q) mean(g[bulk.dp >= 
                  sfp$min.bulk.dp & dpq == q & resampled.training.site == 
                  TRUE]$training.pass, na.rm = TRUE))
            }
            ret <- data.frame(ncalls = ncalls, callable.sens = callable.sens, 
                callable.bp = callable.bp)
        }
        ret$callable.burden <- ret$ncalls/ret$callable.sens
        ret$rate.per.gb <- ret$callable.burden/ret$callable.bp * 
            1e+09/2
        ret$burden <- ret$rate.per.gb * gbp.per.genome
        ret$somatic.sens <- ret$ncalls/ret$burden
        ret$pre.genotyping.burden <- pre.geno.burden
        ret$unsupported.filters <- sfp$max.bulk.alt > 0 | sfp$max.bulk.af > 0
        if (any(ret$unsupported.filters)) {
            warning("mutation burdens must be extrapolated without clonal mutations (i.e., max.bulk.alt=0 and max.bulk.af=0)! burdens will be estimated, but they are invalid")
        }
        ret
    }), muttypes)
    object
}


recompute_burden = function(obj, calls_pass_snv, calls_pass_indel, cigar_change = F, lysis_change = F, lysis_thresh = 0.01, cigar_cutoff_fc = 1, cov_change = F, cov_thresh = 25, quantile1 = 0.25, quantile2 = 0.75, mtypes = c('snv', 'indel')){

	# Extract gatk
	raw_gatk = obj@gatk

	# Give unique id to sites
	raw_gatk[, id := paste0(chr, '_', pos, '_', refnt, '_', altnt)]
	calls_pass_snv$id = paste0(calls_pass_snv$chr, '_', calls_pass_snv$pos, '_', calls_pass_snv$refnt, '_', calls_pass_snv$altnt)
	calls_pass_indel$id = paste0(calls_pass_indel$chr, '_', calls_pass_indel$pos, '_', calls_pass_indel$refnt, '_', calls_pass_indel$altnt)

	# Assign new mutations
	raw_gatk$pass = F
	raw_gatk[(id %in% calls_pass_snv$id) | (id %in% calls_pass_indel$id), pass := T]

	if(cigar_change == T){

	for(mtype in mtypes){

		if(mtype == 'snv'){
			# Extract score cutoffs
			id_score_q = obj@excess.cigar.scores$snv$id.score.q
			hs_score_q = obj@excess.cigar.scores$snv$hs.score.q

			# Assign more strict cigar passes
			raw_gatk = raw_gatk[, pass_strictCigar := muttype == 'snv' &
					id.score > id_score_q*cigar_cutoff_fc &
					hs.score > hs_score_q*cigar_cutoff_fc, ]
			
			# Remove somatic mutation calls and training.passes if they fail the new cigar filter
			raw_gatk[muttype == 'snv' & pass == T & pass_strictCigar == F, pass := F]
			raw_gatk[muttype == 'snv' & training.pass == T & pass_strictCigar == F, training.pass := F]

		} else if(mtype == 'indel'){
			# Extract score cutoffs
			id_score_q = obj@excess.cigar.scores$indel$id.score.q
			hs_score_q = obj@excess.cigar.scores$indel$hs.score.q

			# Assign more strict cigar passes
			raw_gatk = raw_gatk[, pass_strictCigar := muttype == 'indel' &
					id.score > id_score_q*cigar_cutoff_fc &
					hs.score > hs_score_q*cigar_cutoff_fc, ]
			
			# Remove somatic mutation calls and training.passes if they fail the new cigar filter
			raw_gatk[muttype == 'indel' & pass == T & pass_strictCigar == F, pass := F]
			raw_gatk[muttype == 'indel' & training.pass == T & pass_strictCigar == F, training.pass := F]
		}
	}
	}

	if(lysis_change == T){
		raw_gatk[pass == T & lysis.pv >= lysis_thresh, pass := F]
		raw_gatk[training.pass == T & lysis.pv >= lysis_thresh, training.pass := F]
	}

	if(cov_change == T){
		raw_gatk[pass == T & dp <= cov_thresh, pass := F]
		raw_gatk[training.pass == T & dp <= cov_thresh, training.pass := F]
	}

	obj@gatk = raw_gatk
	old_mb = obj@mutburden

	obj = compute_mutburden_new(object = obj,
				gbp.per.genome = scan2.zoo::get.gbp.by.genome(obj),
				quiet = F, q1 = quantile1, q2 = quantile2)

	return(obj)
}


find_clustered_snvs = function(obj = res, extsize = 5000, inpass){

	artpass = obj@gatk[obj@gatk$static.filter & obj@gatk$lysis.fdr < 0.9,]
	indL = list()
	for(k in 1:nrow(inpass)){

		snv = inpass[k, ]
		st = snv$pos - extsize
		end = snv$pos + extsize
		snvs = artpass[chr == snv$chr & pos > st & pos < end,]
		if(nrow(snvs) > nrow(snv)){indL[[k]] = k}
	}
	return(inpass[unlist(indL),])
}

remove_bulk_mpileup = function(passes_file, bulk_bam_file, outfile, altnt = 'altnt', refnt = 'refnt', pos = 'pos', chr = 'chr'){

	tmpL = lapply(1:nrow(passes_file), function(x){

		if(passes_file[x, 'muttype'] == 'snv'){

		# Extract the snv and its basecalls
		site = passes_file[x,]
		alt = passes_file[x, altnt]
		ref_basecall = extract_basecall_and_read(bam = bulk_bam_file, outfile = outfile,
								chr = as.character(site[,chr]),
								pos = as.character(site[,pos]))
		}
		if(x%%50==0){print(x)}

		ref_basecall$V3 = gsub("(.)\\1+", "\\1", ref_basecall$V3)
		sum(ref_basecall$V3 %in% alt) / nrow(ref_basecall)
	})

	# Fail if it returns null
	tmpL = lapply(tmpL, function(x){ifelse(length(x)==0, 0, x)})
	passes_file_bulkFilt = passes_file[!(unlist(tmpL) > 0),]
	return(passes_file_bulkFilt)
}

safe_read_table <- function(file, ...) {
  tryCatch(
    {
      read.table(file, ...)
    },
    error = function(e) {
      message("Error reading file: ", file, "\nReturning NULL.")
      return(NULL)
    }
  )
}


plot_mutprof_from_vcf = function(mutdf, genome_build, ref_genome,
				prefix = '', hg = 7, wd1 = 14, wd2 = 10,
				add_chr = F, ymax1 = 0.2, ymax2 = 1, plot = T,
				oldcols = c('chr', 'pos', 'refnt', 'altnt')){

	setnames(mutdf, old = oldcols, new = c("chr","start","REF","ALT"), skip_absent=TRUE)

	if(add_chr == T){mutdf$chr = paste0('chr', mutdf$chr)}

	mutdf$end = mutdf$start
	mutgr = makeGRangesFromDataFrame(mutdf, keep.extra.columns=T)
	GenomeInfoDb::genome(mutgr) = genome_build
	
	mut_mat <- mut_matrix(vcf_list = mutgr, ref_genome = ref_genome)

	if(plot == T){

		pdf(paste(prefix, 'MutationSignatures_TripleNucleotide.pdf', sep = '_'), width = wd1, height = hg)
		print( plot_96_profile2(mut_mat, plot_size_scale = 2, condensed = T, ymax = ymax1) )
		dev.off()
	}

	type_occurrences <- mut_type_occurrences(mutgr, ref_genome)
	
	# Combine two specific groups (e.g., "G1" and "G2")
	df = as.data.frame(t(type_occurrences))
	df$ctgroup = rownames(df)
	df <- df %>%
	  mutate(group_combined = ifelse(ctgroup %in% c("C>T at CpG", "C>T other"), "C>T", ctgroup))
	df = df[df$ctgroup != 'C>T',]
	df$relcount = df$My_sample / sum(df$My_sample)

	# Plot
	cols = c('skyblue', 'black', 'pink', 'red', 'gray', 'darkolivegreen3', 'mistyrose')

	if(plot == T){

		pdf(paste(prefix, 'MutationSignatures_SingleNucleotide.pdf', sep = '_'), width = wd2, height = hg)
		print(	

		ggplot(df, aes(x = group_combined, y = relcount, fill = ctgroup)) +
		geom_bar(stat = "identity", position = "stack") +
		theme_classic() +
		xlab('') + ylab('Relative Contribution') +
		scale_fill_manual(values = cols) +
		theme(axis.text.x = element_text(size=20),
			axis.text.y = element_text(size=20),
			axis.title = element_text(size=20)) +
		ylim(0, ymax2) +
		rotate_x_text(90)

		)
		dev.off()
	}

	return(mut_mat)
}

extract_basecall_and_read = function(bam, outfile, chr, pos){

	# Define your parameters
	bam_file <- bam
	chromosome <- chr
	position <- pos

	python_script <- "/home/ch253794/SCRIPTS/Extract_readname_and_basecall2.py"  # Use the full path to your Python script

	# Construct the command to run the Python script
	command <- sprintf("python %s %s %s %s %s", python_script, bam_file, outfile, chromosome, position)

	# Execute the command
	system(command)
	
	# Read and return as data frame
	outfile_in = safe_read_table(outfile, colClasses = "character")
	return(outfile_in)
	
}

# Change behavior of plot_96_profile in mutationalpatterns
plot_96_profile2 = function(mut_matrix, colors = NA, ymax = 0.2, condensed = FALSE, plot_size_scale = 1) 
{
    freq <- full_context <- substitution <- context <- NULL
    if (is.na(colors)) {
        colors <- get("COLORS6", envir = asNamespace("MutationalPatterns"))
    }
    if (length(colors) != 6) {
        stop("Provide colors vector with length 6", call. = FALSE)
    }
    norm_mut_matrix <- apply(mut_matrix, 2, function(x) x/sum(x))
    tb <- norm_mut_matrix %>% as.data.frame() %>% tibble::rownames_to_column("full_context") %>% 
        dplyr::mutate(substitution = stringr::str_replace(full_context, 
            "\\w\\[(.*)\\]\\w", "\\1"), context = stringr::str_replace(full_context, 
            "\\[.*\\]", "\\.")) %>% dplyr::select(-full_context) %>% 
        tidyr::pivot_longer(c(-substitution, -context), names_to = "sample", 
            values_to = "freq") %>% dplyr::mutate(sample = factor(sample, 
        levels = unique(sample)))
    if (condensed == TRUE) {
        width <- 1
        spacing <- 0
    }
    else {
        width <- 0.6
        spacing <- 0.5
    }
    plot <- ggplot(data = tb, aes(x = context, y = freq, fill = substitution, 
        width = width)) + geom_bar(stat = "identity", colour = "black", 
        size = 0.2) + scale_fill_manual(values = colors) + facet_grid(sample ~ 
        substitution) + ylab("Relative contribution") + coord_cartesian(ylim = c(0, 
        ymax)) + scale_y_continuous(breaks = seq(0, ymax, 0.1)) + 
        guides(fill = "none") + theme_bw() + theme(axis.title.y = element_text(size = 12*plot_size_scale, 
        vjust = 1), axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 12*plot_size_scale), 
        axis.text.x = element_text(size = 5*plot_size_scale, angle = 90, vjust = 0.5), 
        strip.text.x = element_text(size = 9*plot_size_scale), strip.text.y = element_text(size = 9*plot_size_scale), 
        panel.grid.major.x = element_blank(), panel.spacing.x = unit(spacing, 
            "lines"))
    return(plot)
}

