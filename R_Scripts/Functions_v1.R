# Pepsum Class
setClass("pepsum", representation(
  data = "ANY", data_fp = "ANY", data_fplm = "ANY", data_fplmc = "ANY", 
  data_uniqpeps = "ANY", data_uniqsites = "ANY", totalheat = "matrix", 
  filetype = "character", scanvec = "ANY", modposition.protein = "ANY", 
  chargevec = "ANY", pepvec = "ANY", protvec = "ANY", filenames="ANY", 
  filepath="ANY", modposition.peptide="ANY", modindex="list", 
  modsummary="list"))

# Clean protein accession IDs
simplifyID <- function(idlist){
    proteins <- substr(idlist, start = regexpr("\\|", idlist) + 1,
                       stop = nchar(idlist))
    proteins <- substr(proteins, start = 1, 
                       stop = regexpr("\\|", proteins) - 1)
    return(proteins)
}

# Clean peptides with modification info
simplifyPep <- function(peplist){
    peptides <- as.character(list())	
    peptides <- substr(peplist, start = 3, stop = nchar(peplist) - 2)
    peptides <- gsub("(\\[|\\[-)([0-9]+)(.)([0-9]+)(]+)", 
                     replacement="", peptides)
    peptides <- gsub("^n", replacement = "", peptides)
    return(peptides)
}

# Import MS file as a dataframe
read.ms <- function(input, type){
  object <- new("pepsum")
  object@filetype <- type
  if(type == "SpectraST"){	
    object@data <- read.delim(input, header = T, stringsAsFactors = FALSE)
  }
  if(type == "PeptideProphet"){	
    object@data <- read.delim(input, header = T, stringsAsFactors = FALSE)
    scanlist <- strsplit(as.character(object@data$spectrum), split = ".",
                         fixed = T)
    scanlen <- length(scanlist)
    scanvec <- as.numeric(unlist(scanlist)[seq(2, length(scanlist) * 4, 
                                               by = 4)])
    chargevec <- as.numeric(unlist(scanlist)[seq(4, length(scanlist) * 4, 
                                                 by = 4)])
    peptides <- simplifyPep(as.character(object@data$peptide))
    proteins <- simplifyID(object@data$protein)
    object@data$pepvec <- peptides
    object@data$protvec <- proteins
    object@data$scanvec <- scanvec
    object@data$chargevec <- chargevec
  }
  if(type == "PSP"){
    object@data <- read.delim(input, header=T, skip = 3)
    object@data <- object@data[which(object@data[, "ORGANISM"] == "human"), ]
    positionvec <- gsub('[A-Z,a-z,-]', '', object@data[, "MOD_RSD"])
    protacc <- unique(object@data[, "ACC_ID"])
    for(x in protacc){
	    position <- positionvec[which(object@data[, "ACC_ID"] == x)]
      object@modsummary[[x]] <- as.numeric(position)
    }
  }
  
  if(type == "Nature2"){
	  object@data <- read.delim(input, header = T)
		study.lines <- object@data[object@data[, "Count....Studies.SUMO.2.Modified."] >= 1, ]
		protacc <- unique(as.factor(as.character(study.lines[, "Protein"])))
		for(x in protacc){
		  linematches <- object@data[which(object@data[, "Protein"] == x & object@data[,"Count....Studies.SUMO.2.Modified."] >= 1),]
		  object@modsummary[[x]] <- linematches$"Positions.within.proteins"
		}
	}

  if(type == "Nature2016"){
	  object@data <- read.delim(input, header = T)
	  protacc <- levels(object@data[, "Uniprot"])
	  for(x in protacc){
		  object@modsummary[[x]] <- object@data[which(object@data[, "Uniprot"] == x), "Position"]
	  }
  }
  
  if(type == "Uniprot-SUMO"){
	  object@data <- read.delim(input, header = T, stringsAsFactors = FALSE)
	  protacc <- object@data[, 1]
	  for(x in 1:length(protacc)){
		  crosslinks <- strsplit(object@data$"Cross.link"[x], ";")
		  crosslinkmatch <- crosslinks[[1]][grep("SUMO", crosslinks[[1]])]
		  object@modsummary[[protacc[x]]] <- as.numeric(substr(start = regexpr("\\ [[:digit:]]", crosslinkmatch) + 1, 
		                                                       stop = regexpr("[[:digit:]]\\ ", crosslinkmatch), crosslinkmatch))
		}
  }
  
  if(type == "Uniprot-Ub"){
	  object@data <- read.delim(input, header = T, stringsAsFactors = FALSE)
	  protacc <- object@data[, 1]
	  for(x in 1:length(protacc)){
		  crosslinks <- strsplit(object@data$"Cross.link"[x], ";")
		  crosslinkmatch <- crosslinks[[1]][grep("ubiquitin", crosslinks[[1]])]
		  object@modsummary[[protacc[x]]] <- as.numeric(substr(start = regexpr("\\ [[:digit:]]", crosslinkmatch) + 1, 
		                                                       stop = regexpr("[[:digit:]]\\ ", crosslinkmatch), crosslinkmatch))
	  }
  }

  object@filetype <- type
  object@filenames <- basename(input)
  return(object)
}

# Filter dataframe by a minscore probability cutoff
filterProbability <- function(object, minscore = 0.9){
	keepthese <- c()
	prob_scores <- as.character(object@data$"probability")
	object@data_fp <- object@data[which(object@data$probability > 0.9), ]
	return(object)
}

# Filter peptides for diglys and apply localization score cutoff
localizeMods <- function(object, minscore = 0.75, 
  fasta = "C:/MS-GF+/database/UniProt.human.20141017.RNFISnr.150contams.fasta"){
  require(seqinr)
  require(Biostrings)
  keepthese <- c()
  modsite.peptide <- list()
  modsite.protein <- list()
  flank_pep <- list()
  ptm_lines <- object@data_fp[as.character(object@data_fp$ptm_peptide) != 
	                              "[unavailable]", ]
  ptm_lines <- ptm_lines[grep("242|250", ptm_lines$peptide), ]
  lines <- 1
  for(i in 1:length(ptm_lines$peptide)){
    tempscore <- as.numeric(unlist(regmatches(ptm_lines$ptm_peptide[i],
      gregexpr("[[:digit:]]+\\.*[[:digit:]]*", ptm_lines$ptm_peptide[i]))))
    if(length(which(tempscore >= minscore) > 0) > 0){
      keepthese <- c(keepthese,i)
      n <- which(tempscore >= minscore)
      modsite.peptide[[lines]] <- unlist(gregexpr("[[:digit:]]+\\.*[[:digit:]]*",
				ptm_lines$ptm_peptide[i]))[n] - ((n - 1) * 7 + 2)
      lines <- lines + 1
    }
  }
  object@data_fplm <- ptm_lines[keepthese, ]
  fastaobj <- read.fasta(fasta, seqtype = "AA", as.string = TRUE)
  fasta_names <- simplifyID(names(fastaobj))
  allproteins <- object@data_fplm$protvec
  uniqproteins <- as.character(unique(allproteins))
  pepvec <- as.character(object@data_fplm$pepvec)
  pepvec.l <- length(pepvec)
	for(i in 1:pepvec.l){
    currentprot <- fastaobj[[which(fasta_names == allproteins[i])]][1]
    currentpep <- pepvec[i]
		matched <- matchPattern(currentpep, currentprot)
    seq <- paste(paste(rep("-", times = 5), collapse = ""), currentprot,
			paste(rep("-",times=5),collapse=""),sep="")
		start <- matched@ranges@start
		stop <- matched@ranges@start + matched@ranges@width
		flank_pep[[i]] <- paste(substr(seq, start - 5 + 5, start - 1 + 5), ".", 
		                        currentpep, ".", substr(seq, stop + 5, stop + 9), 
		                        sep = "")
		if(length(matched@ranges) == 0){
			modsite.protein[[i]] <- "null"
		}else{
      modsite.protein[[i]] <- matched@ranges[[1]][1] + 
			  (unlist(modsite.peptide[i]) - 1)
		}
	}
  object@data_fplm$modsite <- as.character(modsite.protein)
  object@data_fplm$flank_pep <- flank_pep
  object@data_fplmc <- object@data_fplm
  object@data_fplmc$modsite <- sapply(object@data_fplmc$modsite, as.character)
  for(n in 1:length(object@data_fplmc$modsite)){
    if(grepl("c\\(", object@data_fplmc$modsite[n])){
      n1 <- gsub("c\\((\\d*)\\, (\\d*)\\)","\\1", 
                 toString(object@data_fplmc$modsite[n]))
			n2 <- gsub("c\\((\\d*)\\, (\\d*)\\)","\\2", 
			           toString(object@data_fplmc$modsite[n]))
			rowholder <- object@data_fplmc[n, ]
			rowholder$modsite <- as.character(n2)
			object@data_fplmc$modsite[n] <- as.character(n1)
			object@data_fplmc <- rbind(object@data_fplmc, rowholder)
		}
	}
	object@data_fplmc$modsite <- sapply(object@data_fplmc$modsite, as.double)
	object@data_fplmc$modsite <- sapply(object@data_fplmc$modsite, as.factor)
	return(object)
}

# Summarize unique modification positions in modsummary table
summarizeProtPositions <- function(object){

  # Get naked accession ID
  allproteins <- object@data_fplmc$protvec
  uniqueproteins <- unique(allproteins)
  
  # Iterates through unique accessions
    # List lines matching each accession
    # Relist positions for single accession as a time (unique)
    # Relist row numbers matching positions for a single accession
    # Give full list of mods matching rows in prot.lines (nonunique)
    # Iterate though unique mods for a single accession
      # Assign an index to each line corresponging to a position match 
  uniq.mods.by.prot <- list()
  lines.position <- list()
  position.index.list <- list()
  index <- 1
  for(i in 1:length(uniqueproteins)){
    uniq.mods.by.prot[[uniqueproteins[i]]] <- unique(unlist(
      object@data_fplmc$modsite[which(allproteins == uniqueproteins[i])]))
    uniq.positions <- uniq.mods.by.prot[[uniqueproteins[i]]]
    uniq.positions.l <- length(uniq.positions)
    prot.lines <- which(allproteins == uniqueproteins[i])
    all.positions <- object@data_fplmc$modsite[prot.lines]
    for(x in uniq.positions){
      lines.position <- which(all.positions == x)
      position.index.list[prot.lines[lines.position]] <- index
      index <- index + 1
    }
  }
  
  # Adds prot.pos.list data to object1@modsummary
  # Adds position.index.list to dataframe
  object@modsummary <- uniq.mods.by.prot
  object@data_fplmc$modindex <- position.index.list
  return(object)
}

# Determine unique peptides and summarize in single lines
uniqPeps <- function(object){
	object@data_uniqpeps <- object@data_fplmc
	duplicates <- which(duplicated(object@data_uniqpeps$modsite) & 
	                    duplicated(object@data_uniqpeps$peptide) & 
	                    duplicated(object@data_uniqpeps$protein) == TRUE)
	for(n in duplicates){
		parentrow <- which(
		  object@data_uniqpeps$modsite == object@data_uniqpeps[n,]$modsite & 
		  object@data_uniqpeps$peptide == object@data_uniqpeps[n,]$peptide &
		  object@data_uniqpeps$protein == object@data_uniqpeps[n,]$protein)[1]	
		object@data_uniqpeps[parentrow,]$probability <- 
		  paste(object@data_uniqpeps[parentrow,]$probability, ",", 
		        object@data_uniqpeps[n,]$probability, sep="")
		object@data_uniqpeps[parentrow,]$spectrum <- 
		  paste(object@data_uniqpeps[parentrow,]$spectrum, ",", 
		        object@data_uniqpeps[n,]$spectrum, sep="")
		object@data_uniqpeps[parentrow,]$spectrum <- 
		  paste(object@data_uniqpeps[parentrow,]$start_scan, ",", 
		        object@data_uniqpeps[n,]$start_scan, sep="")
		}
	object@data_uniqpeps <- object@data_uniqpeps[-duplicates, ]
	return(object)
	}

# Determine unique sites and use index values to summarize peptide information
uniqSites <- function(object){
  # List all indices (ordered by row), then unique indices
  # Initialize Lists
	all.indices <- unlist(object@data_fplmc$modindex)
	uniq.indices <- unique(all.indices)
	numuniq <- length(uniq.indices)
	unique.position.singleline <- rep(0,times=numuniq)
	unique.position.singleweight <- rep(0,times=numuniq)
	stats.scans <- list()
	stats.spectrum <- list()
	stats.ptm <- list()
	stats.xpress <- list()
	stats.heavy <- list()
	stats.light <- list()

	stats.xpress.average<- list()
	stats.xpress.stdev <- list()
	stats.heavy.average <- list()
	stats.heavy.stdev <- list()
	stats.light.average <- list()
	stats.light.stdev <- list()
	
	# Iterate unique indices
	  # Uses first pep match as reference
	  # Combines all pep matches into a single list for the reference
	  # Calculates statistics for combined lists
	i <- 1
	pep.list <- list()
	for(x in uniq.indices){
		pep.list[i] <- list(which(all.indices==x))
		unique.position.singleline[i]<-pep.list[[i]][1]
		if(length(pep.list[[i]]) < 50){									
			stats.scans[[i]] <- object@data_fplmc$"start_scan"[pep.list[[i]]]
			stats.spectrum[[i]] <- object@data_fplmc$"spectrum"[pep.list[[i]]]
			stats.xpress[[i]] <- object@data_fplmc$"xpress"[pep.list[[i]]]
			stats.heavy[[i]] <- object@data_fplmc$"heavy_area"[pep.list[[i]]]
			stats.light[[i]] <- object@data_fplmc$"light_area"[pep.list[[i]]]
		}else if(length(pep.list[[i]])>50){
			stats.scans[[i]] <- "Overflow"
			stats.spectrum[[i]] <- "Overflow"
			stats.xpress[[i]] <- "Overflow"
			stats.heavy[[i]] <- "Overflow"
			stats.light[[i]] <- "Overflow"
		}
		stats.xpress.average[i] <- mean(na.omit(unlist(
			object@data_fplmc$"xpress"[pep.list[[i]]])))
		stats.xpress.stdev[i] <- sd(na.omit(unlist(
			object@data_fplmc$"xpress"[pep.list[[i]]])))
		stats.heavy.average[i] <- mean(na.omit(unlist(
			object@data_fplmc$"heavy_area"[pep.list[[i]]])))
		stats.heavy.stdev[i] <- sd(na.omit(unlist(
			object@data_fplmc$"heavy_area"[pep.list[[i]]])))
		stats.light.average[i] <- mean(na.omit(unlist(
			object@data_fplmc$"light_area"[pep.list[[i]]])))
		stats.light.stdev[i] <- sd(na.omit(unlist(
			object@data_fplmc$"light_area"[pep.list[[i]]])))
		i <- i + 1
	}
	
	# Adds columns to dataframe
	# Trims dataframe to length(unique sites)
	count.peplist <- lapply(pep.list,length)
	object@data_uniqsites<-object@data_fplmc[unique.position.singleline, ]
	object@data_uniqsites$"NumScans"<- as.character(count.peplist)
	object@data_uniqsites$"All_Spectra" <- as.character(stats.spectrum)
	object@data_uniqsites$"All_Scans" <- as.character(stats.scans)

## NOT WORKING
	object@data_uniqsites$"All_Heavy" <- as.character(stats.heavy)
	object@data_uniqsites$"Ave_Heavy" <- as.character(stats.heavy.average)
	object@data_uniqsites$"SD_Heavy" <- as.character(stats.heavy.stdev)
	object@data_uniqsites$"All_Light" <- as.character(stats.light)
	object@data_uniqsites$"Ave_Light" <- as.character(stats.light.average)
	object@data_uniqsites$"SD_Light" <- as.character(stats.light.stdev)
	object@data_uniqsites$"All_Xpress" <- as.character(stats.xpress)
	object@data_uniqsites$"Ave_Xpress" <- as.character(stats.xpress.average)
	object@data_uniqsites$"SD_Xpress" <- as.character(stats.xpress.stdev)
## NOT WORKING
	return(object)
}

# Export Tables
exportTables=function(object, input){
	require(xlsx)
	require(tools)
	wb <- createWorkbook()
	allpep <- createSheet(wb, sheetName = "AllPeptides")
	alldiglypep <- createSheet(wb, sheetName = "AllDiglyPeptides")
	uniqdiglypep <- createSheet(wb, sheetName = "UniqueDiglyPeptides")
	uniqsites <- createSheet(wb, sheetName = "UniqueSites")
	addDataFrame(object@data_fp, allpep)
	addDataFrame(object@data_fplm, alldiglypep)
	addDataFrame(object@data_uniqpeps, uniqdiglypep)
	addDataFrame(object@data_uniqsites, uniqsites)
	saveWorkbook(wb, paste(file_path_sans_ext(object@filenames),".xlsx",sep=""))
}

# Get Significant
getSignificant <- function(object, name = "filename.filtered", 
                           writetsv = T, stdev = 0.9404181){
  # AllIDs are alldata, targetIDs are summarized diglys
  allIDs <- object@data_fp
	targetIDs <- object@data_fplmc
  unique.indexes <- unique(unlist(object@data_fplmc$modindex))
	unique.indexes.len <- length(unique.indexes)
	#unique.indexes
	weighted.ratios <- list()
	linesum = 0
	for(j in 1:unique.indexes.len){
		templines <- which(position.index.list == unique.indexes[[j]])
		targetIDs[templines, ]
		if(unique.indexes[[j]] >= 1){
			### divide by temparea
			temparea = sum(targetIDs[templines, "light_area"]) + 
			  sum(targetIDs[templines, "heavy_area"])
			templines.l <- length(templines)
			tempsum <- 0
			linesum = linesum + templines.l
			weights <- rep(0, times = templines.l)
			for(i in 1:templines.l){
				weights[i] <- (targetIDs[templines[i], "light_area"] + 
				                 targetIDs[templines[i], "heavy_area"]) / temparea
				tempsum = tempsum + targetIDs[templines[i], "xpress"] * weights[i]
			}
			for(i in 1:templines.l){
				weighted.ratios[[templines[i]]]<-tempsum
			}
		}
	}
	targetIDs$weighted.ratios <- weighted.ratios
  
  # Performs log2 of xpress scores
  # Binds log data to dataframe
	logRatiosAll <- log(allIDs$xpress, base=2)
	allIDs <- cbind(allIDs, logRatiosAll)
	
	# Stats for set of all log2 scores
	meanAll <- mean(logRatiosAll)
	medianAll <- median(logRatiosAll)
	logRatiosAllnorm <- logRatiosAll - medianAll
	norm.median.all <- median(logRatiosAllnorm)

	# istribution of weighted averages for unique sites
	numuniq <- length(unique(unlist(targetIDs$modindex))) - 1
	unique.position.singleline <- rep(0, times = numuniq)
	unique.position.singleweight <- rep(0, times = numuniq)

	#! Weighted.Ratios needs to be previously assigned
	# Only a single line is taken; weighted ratio assigned to every line
	# Takes log2 of weighted ratio
	for(i in 1:numuniq){
		unique.position.singleline[i] <- which(targetIDs$modindex == i)[1]
		unique.position.singleweight[i] <- 
		  log(targetIDs[which(targetIDs$modindex == i)[1], 
		                     "weighted.ratios"],base=2)
	}
	# median(unique.position.singleweight)
	# range(unique.position.singleweight)
	# Median raw log2 subtracted from weighted ratio to give normalization
	norm.uniq.pos.singleweight <- unique.position.singleweight - medianAll

	# Calculate SD upper and lower limits
	if(length(stdev) == 1){
		upperlim <- norm.median.all + 2 * stdev
		lowerlim <- norm.median.all - 2 * stdev
		onesdup <- norm.median.all + 1 * stdev
		onesddown <- norm.median.all - 1 * stdev
		}
	if(length(stdev) == 0){
		upperlim <- norm.median.all + 2 * sd(logRatiosAllnorm)
		lowerlim <- norm.median.all - 2 * sd(logRatiosAllnorm)
		onesdup <- norm.median.all + 1 * sd(logRatiosAllnorm)
		onesddown <- norm.median.all - 1 * sd(logRatiosAllnorm)
		}
	names(norm.uniq.pos.singleweight) <- unique.position.singleline
	sort.norm.uniq.pos.singleweight <- sort(norm.uniq.pos.singleweight * -1)

	# Make new plot
	plot.new()
	if(writetsv){
		postscript(paste(name, ".ps", sep=""))
	}
	plot(x = 1:numuniq, y = sort.norm.uniq.pos.singleweight)
	abline(h = upperlim, col = "red")
	abline(h = lowerlim, col = "red")
	if(writetsv){
		dev.off()
	}
	print("one standard deviation is")
	print(onesdup)
	abline(h = onesdup, col = "red")
	abline(h = onesddown, col = "red")
	#### now flag those IDs that are outside 2 std. dev.
	idlen <- length(targetIDs@data$xpress)

	# Calculate Sigma information
	oversigma <- rep(0, times = idlen)
	count <- 0
	#### loop through and test which are outside 1 std. dev
	for(i in 1:numuniq){
		if(norm.uniq.pos.singleweight[i] >= onesdup | 
		   norm.uniq.pos.singleweight[i] <= onesddown){
			count <- count + 1
			oversigma[unique.position.singleline[i]] <- 1
		}
	}
	over2sigma <- rep(0, times = idlen)
	count <- 0
	#### loop through and test which are outside 2 std. dev
	for(i in 1:numuniq){
		if(norm.uniq.pos.singleweight[i] >= upperlim | 
		   norm.uniq.pos.singleweight[i] <= lowerlim){
			count = count + 1
			oversigma[unique.position.singleline[i]] <- 2
		}
	}
	
	# for those missing their weighted average value, put in xpress value
	for(i in 1:length(targetIDs[, "weighted.ratios"])){
		if(targetIDs[i, "weighted.ratios"] == 0){
			targetIDs[i, "weighted.ratios"] <- targetIDs[i, "xpress"]
		}
	}

	# Adds "log.ratios" and "oversigma" columns to dataframe
	targetIDs <- cbind(targetIDs, 
	                        log.ratios = log(targetIDs[, "weighted.ratios"],
	                                         base = 2) - medianAll)
	targetIDs <- cbind(targetIDs, oversigma)

	### part to output only unique significant changers
	uniquelines <- targetIDs[unique.position.singleline, ]
	changelines <- uniquelines[uniquelines[,"oversigma"] >= 1, ]
	nochangelines <- uniquelines[uniquelines[,"oversigma"] == 0, ]

	### part to output only unique non changers

	### now have new object with binary whether outside 2*sigma
	print(i)
	print("unique sites")
	print(numuniq)
	print("over 2 stdev")
	print(length(oversigma[oversigma == 2]))
	print("over 1 stdev")
	print(length(oversigma[oversigma == 1]))
	if(writetsv){
		write.table(file = paste(name, ".all.tsv",sep = ""), 
		            targetIDs, quote = F, sep = "\t", col.names = T, row.names = F)
		write.table(file = paste(name, ".changing.tsv",sep = ""), 
		            changelines, quote = F, sep = "\t", col.names=T, row.names = F)
		write.table(file = paste(name, ".nochange.tsv",sep = ""), 
		            nochangelines, quote = F,sep = "\t", col.names=T, row.names = F)
		cat("FilterFromAll->GetSignificant", 
		    file = paste(name, ".stats.tsv", sep=""), sep = "\n", append = FALSE)
		cat(paste("Unique= ", numuniq, sep=""), 
		    file = paste(name, ".stats.tsv", sep = ""), sep = "\n", append = TRUE)
		cat(paste("Over 2 stdev= ", length(oversigma[oversigma == 2]), sep=""), 
		    file = paste(name, ".stats.tsv", sep = ""), sep = "\n", append = TRUE)
		cat(paste("Over 1 stdev= ", length(oversigma[oversigma == 1]), sep=""), 
		    file = paste(name, ".stats.tsv", sep = ""), sep = "\n", append = TRUE)
	}
	return(list(allIDs, targetIDs, nochangelines, changelines))
}

match.col.iso=function(dataobject,matchobject){	
	count=0
	object.col<-rep("",times=length(dataobject[,"protein.position"]))
	dataprot<-gsub("^[^|]*[|]([^|]*)[|].*$", "\\1", dataobject[,"protein"])

	for(x in dataprot){
		pos.in.match<-which(x==names(matchobject))					
		pos.in.dataprot<-which(x==dataprot)
		sites.dataprot<-dataobject[pos.in.dataprot,"protein.position"]
		if(length(pos.in.match)>0){

			if(length(na.omit(pos.in.dataprot[match(matchobject[[x]],sites.dataprot)]))>0){

				### this give the positions that should be + in the 
				object.col[na.omit(pos.in.dataprot[match(matchobject[[x]],sites.dataprot)])]<-"+"
				}
			count=count+1
			}
		}
	return(object.col)
}
				
compare.tables=function(object1,object2,object3,natureobject,nature2016,ub,sumo,methyl,acetyl,
	uniprot,name="test.tsv",agnostic=FALSE){
	object1mods<-object1@data_uniquesites
	object2mods<-object2@modsummary
	object3mods<-object3@modsummary
	naturemods<-natureobject@modsummary
	nature<-nature2016@modsummary
	sumomods<-sumo@modsummary
	ubmods<-ub@modsummary
	methylmods<-methyl@modsummary
	acetylmods<-acetyl@modsummary
	uniprotmods<-uniprot@modsummary
	
	overlapiso<-list()
	nooverlapiso<-list()
	overlapagno<-list()
 	nooverlapagno<-list()

	nature.overlap<-list()
	sumo.overlap<-list()
	ub.overlap<-list()
	uniprot.overlap<-list()
	
	nature.col<-match.col.iso(object1mods,naturemods)
	nature2016.col<-match.col.iso(object1mods,nature2016)
	sumo.col<-match.col.iso(object1mods,sumomods)
	ub.col<-match.col.iso(object1mods,ubmods)
	methyl.col<-match.col.iso(object1mods,methylmods)
	acetyl.col<-match.col.iso(object1mods,acetylmods)
	uniprot.col<-match.col.iso(object1mods,uniprotmods)
	object2.col<-match.col.iso(object1mods,object2mods)
	object3.col<-match.col.iso(object1mods,object3mods)
	
	output.table<-cbind(object1@data_uniquesites, Object2 = object2.col, Object2=object3.col)
	output.table<-cbind(output.table, Nature.sumo = nature.col, Nature.2016 = nature2016.col)
	output.table<-cbind(output.table, PSP.sumo = sumo.col, PSP.ub = ub.col)
	output.table<-cbind(output.table, PSP.acetyl = acetyl.col, PSP.methyl = methyl.col)
	output.table<-cbind(output.table, UNIPROT.sumo = uniprot.col)
	write.table(output.table,file=name,quote=F,sep="\t",col.names=T,row.names=F)
	return(output.table)
}

compare.objects=function(object1,object2,object3,name="test.tsv"){

	datamods<-object1@data
	object2mods<-object2@modsummary
	object2.col<-match.col.iso(datamods,object2mods)
	object3mods<-object3@modsummary
	object3.col<-match.col.iso(datamods,object3mods)
	output.table<-cbind(datamods,Object2=object2.col,Object3=object3.col)
	write.table(output.table,file=name,quote=F,sep="\t",col.names=T,row.names=F)
	return(output.table)
}

compare.datasets=function(data1,data2){
	count=0
	total1=0
	total2=0
	dataset1=data1@modsummary
	dataset2=data2@modsummary
	for(i in names(dataset1)){
		if(length(which(i==names(dataset2)))>0){
			print(paste("Match for: ",i,sep=""))
			prot=i
			for(j in dataset1[[prot]]){
				print(j)
				if(j %in% dataset2[[prot]]){
					count=count+1
					}
				}
			}
		}
	for(i in names(dataset1)){
		for(j in dataset1[[i]]){
			total1=total1+1
			}
		}
	for(i in names(dataset2)){
		for(j in dataset2[[i]]){
			total2=total2+1
			}
		}
	print(paste("Total in Dataset1: ",total1,sep=""))
	print(paste("Total in Dataset2: ",total2,sep=""))
	print(paste("Count match: ",count,sep=""))

	return(count)
}

createFolders=function(subDir,rawdata){
	dir.create(file.path(tppmainDir,subDir))
	dir.create(file.path(tppmainDir,subDir,"1_Spectra"))
	dir.create(file.path(tppmainDir,subDir,"2_MS-GF+"))
	dir.create(file.path(tppmainDir,subDir,"3_TPP"))
	dir.create(file.path(tppmainDir,subDir,"4_Stats"))
	file.create(file.path(tppmainDir,subDir,paste(rawdata,"_log.txt",sep="")))
	cat(paste("CreateFolders executed for ",subDir,sep=""),file=file.path(tppmainDir,subDir,paste(rawdata,"_log.txt",sep="")),sep="\n",append=TRUE)
	}

msConvert=function(subDir,rawdata){
	cat(paste("MSConvert executed for ",rawdata,sep=""),file=file.path(tppmainDir,subDir,paste(rawdata,"_log.txt",sep="")),sep="\n",append=TRUE)
	cmd=paste('C:/Inetpub/tpp-bin/msconvert \"',file.path(tppmainDir,subDir,"1_Spectra",paste(rawdata,".raw",sep="")),'\" -v --mzXML -o \"',file.path(tppmainDir,subDir,"1_Spectra"),'\" >>',file.path(tppmainDir,subDir,paste(rawdata,"_log.txt",sep=""))," 2>&1",sep="") 
		shell(paste('C:/Inetpub/tpp-bin/msconvert \"',file.path(tppmainDir,tppsubDir,paste(rawdataname,".raw",sep="")),'\" -v --mzXML -o \"',file.path(tppmainDir,tppsubDir),'\" >>',file.path(tppmainDir,tppsubDir,paste(rawdataname,"_log.txt",sep=""))," 2>&1",sep=""), flag = "/c", intern = TRUE, wait = TRUE, translate = FALSE, mustWork = FALSE)


	finish<-shell(cmd, flag = "/c", intern = TRUE, wait = TRUE, translate = FALSE, mustWork = FALSE)
	finished<-finish[length(finish)]
}

# Enzyme: 0=unspecific 1=Trypsin 9=No cleavage 10=WaLP
# Instrument: 0=LowRes 1=HighRes 2=TOF 3=QExactive
# Fragmentation: 0=As Written 1=CID 2=ETD 3=HCD 4=Merge spectra from same precursor
msGF=function(subDir,rawdata,enz=10,inst=1,frag=0,outputname,mods){
	cat(paste("MS-GF+ executed for ",rawdata,sep=""),file=file.path(tppmainDir,subDir,paste(rawdata,"_log.txt",sep="")),sep="\n",append=TRUE)
	setwd(file.path(tppmainDir,subDir,rawdata))
	cmd=paste('C: && cd \"C:/MS-GF+/\" && \"',file.path('C:/Program Files/Java/jre1.8.0_25/bin/java.exe'),'\" -Xmx12000M -d64 -jar MSGFplus.jar -thread 8 -s \"',file.path(tppmainDir,subDir,"1_Spectra",paste(rawdata,".mzXML",sep="")),'\" -o \"',file.path(tppmainDir,subDir,"2_MS-GF+",paste(rawdata,".mzid",sep="")),'\" -e ',enz,' -ti 0,2 -tda 1 -inst ',inst,' -ntt 0 -m ',frag,' -t 10ppm -d \"database/UniProt.human.20141017.RNFISnr.150contams.fasta\" -mod diGLY.txt >>',file.path(tppmainDir,subDir,paste(rawdata,"_log.txt",sep=""))," 2>&1",sep="")
	shell(cmd, flag = "/c", intern = FALSE, wait = TRUE, translate = FALSE, mustWork = FALSE)
	if(silac = "true"){
	  cmd=paste('C: && cd \"C:/MS-GF+/\" && \"',file.path('C:/Program Files/Java/jre1.8.0_25/bin/java.exe'),'\" -Xmx12000M -d64 -jar MSGFplus.jar -thread 8 -s \"',file.path(tppmainDir,subDir,"1_Spectra",paste(rawdata,".mzXML",sep="")),'\" -o \"',file.path(tppmainDir,subDir,"2_MS-GF+",paste(rawdata,"_heavy.mzid",sep="")),'\" -e ',enz,' -ti 0,2 -tda 1 -inst ',inst,' -ntt 0 -m ',frag,' -t 10ppm -d \"database/UniProt.human.20141017.RNFISnr.150contams.fasta\" -mod diGLYh.txt >>',file.path(tppmainDir,subDir,paste(rawdata,"_log.txt",sep=""))," 2>&1",sep="")
	  shell(cmd, flag = "/c", intern = FALSE, wait = TRUE, translate = FALSE, mustWork = FALSE)
	}
}

mzidFix=function(subDir,rawdata){
	cat(paste("mzidFix executed for ",rawdata,sep=""),file=file.path(tppmainDir,subDir,paste(rawdata,"_log.txt",sep="")),sep="\n",append=TRUE)
	setwd(file.path(tppmainDir,subDir,"2_MS-GF+")) 	# set location of your mzid file
	doc <- xmlParse(paste(rawdata,".mzid",sep=""))				# parse and find peptide elements
	r = xmlRoot(doc)						#gives content of root

	# mzid has peptide IDs in the third column,so this value will always be 3
	tags <- xmlElementsByTagName(r[[3]], "Peptide") #get the peptide tags and store in tags object.

	# function to fix all tags, run all lines in order
	manipulate <- function(tag) {
    		## get 'Modification' node set
    		dMod <- tag["Modification"]
    		## get 'location' numbers
    		loc <- sapply(dMod, xmlGetAttr, "location")
    		## get the sum of 'monoisotopicMassDelta' 
    		lapply(unique(loc), function(i) {
        		if(length(d <- dMod[loc == i]) > 1) {
            		nm <- "monoisotopicMassDelta"
            		s <- sapply(d, xmlGetAttr, nm)
            		xmlAttrs(d[[1]])[nm] <- sum(as.numeric(s))
        		}
    		})
    		## remove duplicated location nodes
    		removeNodes(dMod[duplicated(loc)])
    		## return the adjusted tag
    		tag
		}

	# this command takes about 20 minutes, 
	# R appears to freeze but it does not actually, 
	# just let it work until you can interact with the GUI again
	lapply(tags, manipulate)   				# applies manipulate to the XML structure to fix tags

	saveXML(r,file.path(tppmainDir,subDir,"2_MS-GF+",paste(rawdata,"_fix.mzid",sep="")), indent=TRUE)	# saves new XML file
	heavyparse <- readLines(file.path(tppmainDir,subDir,"2_MS-GF+",paste(rawdata,"_fix.mzid",sep="")))	
	sink(file.path(tppmainDir,subDir,"2_MS-GF+",paste(rawdata,"_fix.mzid",sep="")))
	cat('<?xml version="1.0" encoding="windows-1252"?>\n')
	sink()
	write(heavyparse,file.path(tppmainDir,subDir,"2_MS-GF+",paste(rawdata,"_fix.mzid",sep="")),append=TRUE)
	}

parseMzid=function(subDir,rawdata){
	cat(paste("parseMzid executed for ",rawdata,sep=""),file=file.path(tppmainDir,subDir,paste(rawdata,"_log.txt",sep="")),sep="\n",append=TRUE)
	mzidparse <- readLines(file.path(tppmainDir,subDir,"2_MS-GF+",paste(rawdata,".mzid",sep="")))
	mzidparse <- gsub('<Enzyme missedCleavages="1000" semiSpecific="true" id="WaLP">','<Enzyme missedCleavages="1000" semiSpecific="true" id="Tryp">',mzidparse)
	mzidparse <- gsub('<userParam name="wt-aLP"/>','<cvParam accession="MS:1001251" cvRef="PSI-MS" name="Trypsin"/>',mzidparse)
	mzidparse <- gsub('post="\\?"','post="-"',mzidparse)
	writeLines(mzidparse,file.path(tppmainDir,subDir,"2_MS-GF+",paste(rawdata,"_parse.mzid",sep="")))
	}

parseMzidMalp=function(subDir,rawdata){
	cat(paste("parseMzid executed for ",rawdata,sep=""),file=file.path(tppmainDir,subDir,paste(rawdata,"_log.txt",sep="")),sep="\n",append=TRUE)
	mzidparse <- readLines(file.path(tppmainDir,subDir,"2_MS-GF+",paste(rawdata,".mzid",sep="")))
	mzidparse <- gsub('<Enzyme missedCleavages="1000" semiSpecific="true" id="MaLP">','<Enzyme missedCleavages="1000" semiSpecific="true" id="Tryp">',mzidparse)
	mzidparse <- gsub('<userParam name="M190A-aLP"/>','<cvParam accession="MS:1001251" cvRef="PSI-MS" name="Trypsin"/>',mzidparse)
	writeLines(mzidparse,file.path(tppmainDir,subDir,"2_MS-GF+",paste(rawdata,"_parse.mzid",sep="")))
	}

idConvert=function(subDir,rawdata){
	cat(paste("idConvert executed for ",rawdata,sep=""),file=file.path(tppmainDir,subDir,paste(rawdata,"_log.txt",sep="")),sep="\n",append=TRUE)
	setwd(file.path(tppmainDir,subDir,"3_TPP"))
	cmd=""
	cmd=paste('C:/idconvert/idconvert.exe \"',file.path(tppmainDir,subDir,"2_MS-GF+",paste(rawdata,"_parse.mzid",sep="")),'\" --pepXML -v -e .pep.xml -o \"',file.path(tppmainDir,subDir,"3_TPP"),'\" >>',file.path(tppmainDir,subDir,paste(rawdata,"_log.txt",sep=""))," 2>&1",sep="")
	shell(cmd, flag = "/c", intern = FALSE, wait = TRUE, translate = FALSE, mustWork = FALSE)
	}

idTrypConvert=function(subDir,rawdata){
	cat(paste("idConvert executed for ",rawdata,sep=""),file=file.path(tppmainDir,subDir,paste(rawdata,"_log.txt",sep="")),sep="\n",append=TRUE)
	setwd(file.path(tppmainDir,subDir,"3_TPP"))
	cmd=""
	cmd=paste('C:/idconvert/idconvert.exe \"',file.path(tppmainDir,subDir,"2_MS-GF+",paste(rawdata,".mzid",sep="")),'\" --pepXML -v -e .pep.xml -o \"',file.path(tppmainDir,subDir,"3_TPP"),'\" >>',file.path(tppmainDir,subDir,paste(rawdata,"_log.txt",sep=""))," 2>&1",sep="")
	shell(cmd, flag = "/c", intern = FALSE, wait = TRUE, translate = FALSE, mustWork = FALSE)
	}

pepProphet=function(subDir,rawdata){	
	cat(paste("pepProphet executed for ",rawdata,sep=""),file=file.path(tppmainDir,subDir,paste(rawdata,"_log.txt",sep="")),sep="\n",append=TRUE)
	cmd=""
	cmd=paste('cd ',file.path(tppmainDir,subDir,"3_TPP"),' && C:/Inetpub/tpp-bin/xinteract -N',paste(rawdata,"_parse_pp.pep.xml",sep=""),' -p0.05 -l7 -c2.5 -PPM -eN -PPM -OANMFPp -dXXX_ -I1 ',paste(rawdata,"_parse.pep.xml",sep=""),' >>',file.path(tppmainDir,subDir,paste(rawdata,"_log.txt",sep=""))," 2>&1",sep="") 
	shell(cmd, flag = "/c", intern = FALSE, wait = TRUE, translate = FALSE, mustWork = FALSE)
	}

pepTrypProphet=function(subDir,rawdata){	
	cat(paste("pepProphet executed for ",rawdata,sep=""),file=file.path(tppmainDir,subDir,paste(rawdata,"_log.txt",sep="")),sep="\n",append=TRUE)
	cmd=""
	cmd=paste('cd ',file.path(tppmainDir,subDir,"3_TPP"),' && C:/Inetpub/tpp-bin/xinteract -N',paste(rawdata,"_pp.pep.xml",sep=""),' -p0.05 -l7 -c2.5 -PPM -eT -PPM -OANMFPp -dXXX_ -I1 ',paste(rawdata,".pep.xml",sep=""),' >>',file.path(tppmainDir,subDir,paste(rawdata,"_log.txt",sep=""))," 2>&1",sep="") 
	shell(cmd, flag = "/c", intern = FALSE, wait = TRUE, translate = FALSE, mustWork = FALSE)
	}

parsePepxml=function(subDir,rawdata){
	cat(paste("parsePepxml executed for ",rawdata,sep=""),file=file.path(tppmainDir,subDir,paste(rawdata,"_log.txt",sep="")),sep="\n",append=TRUE)
	mzidpepparse <- readLines(file.path(tppmainDir,subDir,"3_TPP",paste(rawdata,"_parse_pp.pep.xml",sep="")))
	mzidpepparse <- gsub('<modification_info mod_nterm_mass="-16.0187240729" modified_peptide="Q','<modification_info modified_peptide="Q[111]',mzidpepparse)
	writeLines(mzidpepparse,file.path(tppmainDir,subDir,"3_TPP",paste(rawdata,"_parse_pp_parse.pep.xml",sep="")))
	}

parseTrypPepxml=function(subDir,rawdata){
	cat(paste("parsePepxml executed for ",rawdata,sep=""),file=file.path(tppmainDir,subDir,paste(rawdata,"_log.txt",sep="")),sep="\n",append=TRUE)
	mzidpepparse <- readLines(file.path(tppmainDir,subDir,"3_TPP",paste(rawdata,"_pp.pep.xml",sep="")))
	mzidpepparse <- gsub('<modification_info mod_nterm_mass="-16.0187240729" modified_peptide="Q','<modification_info modified_peptide="Q[111]',mzidpepparse)
	writeLines(mzidpepparse,file.path(tppmainDir,subDir,"3_TPP",paste(rawdata,"_pp_parse.pep.xml",sep="")))
	}

ptmProphet=function(subDir,rawdata,outputname){
	cat(paste("ptmProphet executed for ",rawdata,sep=""),file=file.path(tppmainDir,subDir,paste(rawdata,"_log.txt",sep="")),sep="\n",append=TRUE)
	setwd(file.path(tppmainDir,subDir,"3_TPP"))
	cmd=""
	cmd=paste('C:/Inetpub/tpp-bin/PTMProphetParserMod K,114.0429,Q,-17.02655,AS,42.0106,M,15.9949 MZTOL=0.1 NOUPDATE ',file.path(tppmainDir,subDir,"3_TPP",
	paste(rawdata,"_parse_pp_parse.pep.xml",sep="")),' -p0.05 -l7 -c2.5 -PPM -eN -PPM -OANMFPp -dXXX_ -I1 ',paste(outputname,"_parse_pp_parse_ptm.pep.xml",sep=""),' >>',file.path(tppmainDir,subDir,paste(rawdata,"_log.txt",sep=""))," 2>&1",sep="")
	shell(cmd, flag = "/c", intern = FALSE, wait = TRUE, translate = FALSE, mustWork = FALSE)
	}

ptmTrypProphet=function(subDir,rawdata,outputname){
	cat(paste("ptmProphet executed for ",rawdata,sep=""),file=file.path(tppmainDir,subDir,paste(rawdata,"_log.txt",sep="")),sep="\n",append=TRUE)
	setwd(file.path(tppmainDir,subDir,"3_TPP"))
	cmd=""
	cmd=paste('C:/Inetpub/tpp-bin/PTMProphetParserMod K,114.0429,Q,-17.02655,AS,42.0106,M,15.9949 MZTOL=0.1 NOUPDATE ',file.path(tppmainDir,subDir,"3_TPP",
	paste(rawdata,"_pp_parse.pep.xml",sep="")),' -p0.05 -l7 -c2.5 -PPM -eT -PPM -OANMFPp -dXXX_ -I1 ',paste(outputname,"_pp_parse_ptm.pep.xml",sep=""),' >>',file.path(tppmainDir,subDir,paste(rawdata,"_log.txt",sep=""))," 2>&1",sep="")
	shell(cmd, flag = "/c", intern = FALSE, wait = TRUE, translate = FALSE, mustWork = FALSE)
	}

combineAnalyses=function(subDir,files,ext=""){
	setwd(file.path(tppmainDir,subDir,"3_TPP"))
	cmdraw=""
	for(item in files){
		rawdata=file_path_sans_ext(item)
		cmdraw=paste(cmdraw,file.path(tppmainDir,subDir,"3_TPP",paste(rawdata,ext,"_parse_pp.pep.xml",sep=""))," ",sep="")
		}
	ablabel=""
	if(ext=="_light"){
		ablabel="A_"
	}else if(ext=="_heavy_fix"){
		ablabel="B_"
	}
	cmd=paste("cd ",file.path(tppmainDir,subDir,"3_TPP")," && c:/Inetpub/tpp-bin/InterProphetParser NONSS ",cmdraw,ablabel,"combined",ext,"_parse_pp.pep.xml", 
	" && c:/Inetpub/tpp-bin/ProteinProphet ",file.path(tppmainDir,subDir,"3_TPP",paste(ablabel,"combined",ext,"_parse_pp.pep.xml",sep=""))," ",file.path(tppmainDir,subDir,"3_TPP",paste(ablabel,"combined",ext,"_parse_pp.prot.xml",sep=""))," IPROPHET", sep="")

	shell(cmd, flag = "/c", intern = FALSE, wait = TRUE, translate = FALSE, mustWork = FALSE)
	}

combineHL=function(subDir,rawdata){
	cat(paste("combineHL executed for ",rawdata,sep=""),file=file.path(tppmainDir,subDir,paste(rawdata,"_log.txt",sep="")),sep="\n",append=TRUE)
	cmd=""
	cmd=paste('cd ',file.path(tppmainDir,subDir,"3_TPP"),' && c:/Inetpub/tpp-bin/InterProphetParser NONSS ',file.path(tppmainDir,subDir,"3_TPP",paste("A_",rawdata,"_light_parse_pp_parse_ptm.pep.xml",sep="")),' ',file.path(tppmainDir,subDir,"3_TPP",paste("B_",rawdata,"_heavy_fix_parse_pp_parse_ptm.pep.xml",sep="")),' ',rawdata,"_lightheavy_parse_pp_parse_ptm.pep.xml && ",
	"c:/Inetpub/tpp-bin/ProteinProphet ",file.path(tppmainDir,subDir,"3_TPP",paste(rawdata,"_lightheavy_parse_pp_parse_ptm.pep.xml",sep=""))," ",file.path(tppmainDir,subDir,"3_TPP",paste(rawdata,"_lightheavy_parse_pp_parse_ptm.prot.xml",sep=""))," IPROPHET >>",file.path(tppmainDir,subDir,paste(rawdata,"_log.txt",sep=""))," 2>&1",sep="")
	shell(cmd, flag = "/c", intern = FALSE, wait = TRUE, translate = FALSE, mustWork = FALSE)
	}

quantifyHL=function(subDir,rawdata){
	cat(paste("quantifyHL executed for ",rawdata,sep=""),file=file.path(tppmainDir,subDir,paste(rawdata,"_log.txt",sep="")),sep="\n",append=TRUE)
	cmd=""
	cmd=paste('cd ',file.path(tppmainDir,subDir,"3_TPP"),' && c:/Inetpub/tpp-bin/XPressPeptideParserMod ',file.path(tppmainDir,subDir,"3_TPP",paste(rawdata,"_lightheavy_parse_pp_parse_ptm.pep.xml",sep="")),' -m0.05 -nK,8.014199 -c5 -p1 >>',file.path(tppmainDir,subDir,file=paste(rawdata,"_log.txt",sep=""))," 2>&1",sep="") 
	shell(cmd, flag = "/c", intern = FALSE, wait = TRUE, translate = FALSE, mustWork = FALSE)
	}

pvalue=function(subDir,rawdata){
	cat(paste("pvalue executed for ",rawdata,sep=""),file=file.path(tppmainDir,subDir,paste(rawdata,"_log.txt",sep="")),sep="\n",append=TRUE)
	ptmfile=list.files(file.path(tppmainDir,subDir,"3_TPP"),pattern=paste("^",rawdata,".*_ptm.pep.xml$",sep=""))
	file.remove(file.path(tppmainDir,subDir,"3_TPP",list.files(file.path(tppmainDir,subDir,"3_TPP"),pattern=paste("^",rawdataname,"\\_.*\\_ptm\\.pep\\.xml\\.index$",sep=""))))
	file.remove(file.path(tppmainDir,subDir,"3_TPP",list.files(file.path(tppmainDir,subDir,"3_TPP"),pattern=paste("^",rawdataname,"\\_.*\\_ptm\\.pep\\-MODELS\\.html$",sep=""))))
	if(length(ptmfile)>1){
		print("ERROR! TOO MANY PTM FILES")
	}
	doc <- xmlParse(file.path(tppmainDir,subDir,"3_TPP",ptmfile[1]))				# parse and find peptide elements
	r = xmlRoot(doc)						#gives content of root
	tags <- xmlElementsByTagName(r[[2]][[1]][[2]], "error_point") #get the peptide tags and store in tags object.
	error=sapply(tags[22],xmlGetAttr,"error")
	prob=sapply(tags[22],xmlGetAttr,"min_prob")
	evalue=as.numeric(error[[1]])
	pvalue=as.numeric(prob[[1]])
	if(evalue==0.01){
		print(pvalue)
		cmd=paste("PepXMLViewer.cgi -B exportSpreadsheet 1 -F mPprobability ",pvalue," -C Pprobability,Mptm_peptide,Gspectrum,Gstart_scan,Iiprobability,Gions2,Gpeptide,Gprotein,Gcalc_neutral_pep_mass,Qxpress,Gassumed_charge,Gprotein_descr,GpI,Gppm,Gnum_tol_term,Qlight_area,Qheavy_area -I C:/Inetpub/wwwroot/ISB/data/",subDir,"/3_TPP/",ptmfile[1],sep="")
		shell(cmd, flag = "/c", intern = FALSE, wait = TRUE, translate = FALSE, mustWork = FALSE)
		file.rename(from=file.path(tppmainDir,subDir,"3_TPP",paste(file_path_sans_ext(ptmfile[1]),".xls",sep="")),to=file.path(tppmainDir,subDir,"4_Stats",paste(file_path_sans_ext(ptmfile[1]),".xls",sep="")))
		setwd(file.path(tppmainDir,subDir,"4_Stats"))
		xlsfile<-list.files(path=file.path(tppmainDir,subDir,"4_Stats"),pattern="\\ptm.pep.xls")
		if(length(xlsfile)>1){
			print("ERROR! TOO MANY XLS FILES")
		}
	}else{
		stop("Error")
	}
}

#GET RID OF PARSE RENAMING. Just add to a log. Then: Raw->Mzid->Pep.XML->PP.Pep.XML->PTM.Pep.XML->Xls->xlsx(filter)
tppFlow <- function(msgf='pre',setname="default",silac="false",enzyme="walp",combsplit="false"){
  enz <- vector(mode="list",length=3)
  names(enz) <- c('walp','malp','trypsin')
  enz['walp'] <- 10; enz['malp'] <- 11; enz['typsin'] <- 1;
  files=list.files(file.path(tppmainDir,"input"))
  dir.create(file.path(tppmainDir,setname))
  for(file in files){
    createFolders #same for everything
    if(msgf=='pre'){
      msconvert(setname,file) #same for everything
      msGF(setname,file,enz[enzyme],silac)
    }
    if(silac == "true"){
      mzidFix
    }
    if(enzyme == 'walp' | enzyme=='malp'){
      parseMzid(enzyme, silac)
    }
    idconvert(setname, file, silac)
    pepProphet(setname, file, silac)
  }
  if(combsplit == 'true'){
    combineAnalyses(setname, file, silac)
    file='combined'
    parsePepxml(setname, file, silac)
    ptmProphet(setname, file, silac)
    if(silac == "true"){
        combineHL(setname, file, silac)
        quantifyHL(setname, file, silac)
    }
    pvalue(setname, file, silac)
  }else{
    for(file in files){
      parsePepxml(setname, file, silac)
      ptmProphet(setname, file, silac)
      if(silac == "true"){
          combineHL(setname, file, silac)
          quantifyHL(setname, file, silac)
      }
      pvalue(setname, file, silac)
    }
  }
}

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
input <- paste("C:/Inetpub/wwwroot/ISB/data/2017-04-04/032717_25660_HG_L/4_Stats/",
            "032717_25660_HG_L_parse_pp_parse_ptm.pep.xls",sep="")
type="PeptideProphet"
obj <- read.ms(input, type)
obj@data_fp
length(obj@data$probability)
obj <- filterProbability(obj, minscore = 0.9)
length(obj@data_fp$probability)
obj <- localizeMods(obj, minscore = 0.9)
length(obj@data_fplm$probability)
obj <- summarizeProtPositions(obj)
length(obj@modsummary)
length(obj@data_fplmc$probability)
obj <- uniqPeps(obj)
length(obj@data_uniqpeps$probability)
obj <- uniqSites(obj)
obj@data_uniqsites$probability
exportTables(obj, input)
getwd()
#system.time(#functionToOptimize)
format(Sys.time(), "%a %b %d %H:%M:%S %Y")
