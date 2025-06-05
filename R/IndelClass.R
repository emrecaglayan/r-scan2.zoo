
getSequenceContext <- function(genome, position, chr, offsetL= 10, offsetR=50){
	  sequence <- Biostrings::getSeq(genome, 
				GRanges(seqnames=chr,ranges=IRanges(start=position-offsetL, end=position+offsetR)))
	  context_string <- as.character(sequence)
	  sequenceContext <- list(sequence=sequence,context_string=context_string)
}




attribute_sequence_contex_indel <- function(genome, in_CHROM, in_POS, 
				in_REF,in_ALT, in_verbose = FALSE, 
				in_offsetL=50,in_offsetR=50){

	diff <- abs(nchar(in_REF)-nchar(in_ALT))
	out_dat <- list()
	out_dat$CHROM <- in_CHROM
	out_dat$POS <- in_POS
	out_dat$REF <- in_REF
	out_dat$ALT <- in_ALT

	ref_length <- nchar(in_REF)
	alt_length <- nchar(in_ALT)
	if(ref_length > alt_length){
		out_dat$Type <- "Del"
	}else{
		out_dat$Type <- "Ins"
	}
	
	context_list <- getSequenceContext(genome = genome,
	                                     position = in_POS,
	                                     chr = in_CHROM,
	                                     offsetL= in_offsetL,
	                                     offsetR= in_offsetR)
	out_dat$SequenceContext <- context_list$context_string
	out_dat$Differnce <- diff
	out_dat$Change <- paste0(in_REF,">",in_ALT)
	return(out_dat)
}

rgestStartSubstr<-function(word1, word2){ 
    word1vec<-unlist(strsplit(word1, "", fixed=TRUE))
    word2vec<-unlist(strsplit(word2, "", fixed=TRUE))
    indexes<-intersect(1:nchar(word1), 1:nchar(word2))
    bools<-word1vec[indexes]==word2vec[indexes]
    if(bools[1]==FALSE){
        ""
    }else{
        lastChar<-match(1,c(0,diff(cumsum(!bools))))-1
        if(is.na(lastChar)){
            lastChar<-indexes[length(indexes)]
        }
        substr(word1, 1,lastChar)
    }
}



attribution_of_indels <-function(genome, in_CHROM, in_POS,in_REF,in_ALT) { 


	ref_length <- nchar(in_REF) 
	alt_length <- nchar(in_ALT)
	max_length <- max(ref_length,alt_length)

	in_offsetL = max(6*max_length,50)
	in_offsetR = max(6*max_length,50)

	in_dat_return <- attribute_sequence_contex_indel(genome=genome, 
					in_CHROM=in_CHROM, 
					in_POS=in_POS,
					in_REF=in_REF,
					in_ALT=in_ALT,
					in_offsetL=in_offsetL,
					in_offsetR=in_offsetR)
	sequence_string <- DNAString(in_dat_return$SequenceContext)
	motif_end_index <- ref_length + in_offsetL
	rest_start_index <- motif_end_index + 1 

	if(in_dat_return$Type == "Del"){
		motif <- substring(in_REF, ref_length-in_dat_return$Differnce+1)
	}
	if(in_dat_return$Type == "Ins"){
		motif <- substring(in_ALT, alt_length-in_dat_return$Differnce+1)
	}

	motif_length <- nchar(motif)
	rest_string_R <-subseq(sequence_string, start=rest_start_index)
	rest_string_L <-subseq(sequence_string, start=1, end=in_offsetL+1)
	match <- matchPattern(motif, rest_string_R)
	
	###classify deletion
	if(motif_length > 5){
		in_dat_return$key1=5
	}else{
		in_dat_return$key1=motif_length
	}


	####define key3
	if(motif_length ==1){
		if(as.character(motif)== "T" | as.character(motif) == "A"){
			in_dat_return$key3='T'
		}else if(as.character(motif)== "G" | as.character(motif) == "C"){
			in_dat_return$key3='C'
		}else{ 
			in_dat_return$key3='NA' }
		#####define key4 for repeats
		i=0
		rep_size = c(1,2,3,4,5,6)
		for ( irep in rep_size)	{
			if(irep %in%  start(match)){
				i = i + 1
			}else{  break  }
		}	
		repeat_size = min(i,5)
		in_dat_return$key4=repeat_size
	}else{
		###repeats
		i=0
		rep_size = c(1,2,3,4,5,6)
		for(irep in rep_size){
			match_site = 1 + i * motif_length
			if(match_site %in% start(match)){
				i = i + 1
			}else{  break  }
		}
		repeat_size = min(i,5)

		####also check for microhomology if it is a insertion
		if(in_dat_return$Type=='Del'){
			###microhomology
			homologR <- as.character(rest_string_R)
			homologL <- as.character(Biostrings::reverse(rest_string_L))
			motif_rev <- Biostrings::reverse(motif)
			match_R = rgestStartSubstr(homologR, motif)
			match_L = rgestStartSubstr(homologL, motif_rev)
			homolog_size = min(5,max(nchar(match_R),nchar(match_L)))
			
			if(repeat_size >= 1){
				in_dat_return$key3='R'
				in_dat_return$key4=repeat_size
			} else if(homolog_size > 0){
				in_dat_return$key3='M'
				in_dat_return$key4=homolog_size
			} else {
				in_dat_return$key3='R'
				in_dat_return$key4=repeat_size
			}
		}else{
			in_dat_return$key3='R'
			in_dat_return$key4=repeat_size
		}
	}	
	in_dat_return$Indelsig = paste(in_dat_return$key1,in_dat_return$Type, 
					in_dat_return$key3, in_dat_return$key4, sep = ":")
	in_dat_return$Indelsig
}

