bigcancersimulator <- function(num_generations, generationmetastasis, num_metastasis, purity) {
  par(new)
  library(stringr)
  library(dplyr)
  library(RColorBrewer)
  library(abind)
  library(scales)
  
  ChrLengths <- c(248956422, 242193529, 198295559, 190214555, 181538259,
                  170805979, 159345973, 145138636, 138394717, 133797422,
                  135086662, 133275309, 114364328, 107043718, 101991189,
                  90338345, 83257441, 80373285, 58617616, 64444167, 
                  46709983, 50818468)
  generateCNA <- function(){
    maternal <- rpois(1, 3)
    paternal <- rpois(1, 3)
    if (maternal==paternal) {
      if (paternal ==1){
        maternal <- 0
      }
      else if (paternal == 0){
        maternal <- 1
      }
    }
    total <- maternal + paternal
    copynum <- c(maternal, paternal, total)
    return(copynum)
  }
  arrayunfolder <- function(array, z, needle) {
    index <- NA
    for (i in 1:nrow(array[,,z])) {
      if(needle %in% array[i, 1, z]:array[i, 2, z]) {
        index <- i
        break
      }
    }
    return(index)
  }
  metastasis <- function(bigtable, copynuminfo, num_generations, cnaprob, lineage) {
    ChrLengths <- c(248956422, 242193529, 198295559, 190214555, 181538259,
                    170805979, 159345973, 145138636, 138394717, 133797422,
                    135086662, 133275309, 114364328, 107043718, 101991189,
                    90338345, 83257441, 80373285, 58617616, 64444167, 
                    46709983, 50818468)
    
    generateCNA <- function(){
      maternal <- rpois(1, 3)
      paternal <- rpois(1, 3)
      if (maternal==paternal) {
        if (paternal ==1){
          maternal <- 0
        }
        else if (paternal == 0){
          maternal <- 1
        }
      }
      total <- maternal + paternal
      copynum <- c(maternal, paternal, total)
      return(copynum)
    }
    
    arrayunfolder <- function(array, z, needle) {
      index <- NA
      for (i in 1:nrow(array[,,z])) {
        if(needle %in% array[i, 1, z]:array[i, 2, z]) {
          index <- i
          break
        }
      }
      return(index)
    }
    
    #print("METASTASIS OCCURS")
    xs <- c()
    ys <- c()
    bigtable[,8] <- 1
    outtable <- matrix(0, nrow = (nrow(bigtable)+num_generations), ncol = 10)
    outtable[1:nrow(bigtable),] <- bigtable
    CNAprob <- cnaprob
    chrcopy <- c("mat", "pat")
    aberratedchrs <- copynuminfo[[1]]
    abberatedloci <- copynuminfo[[2]]
    abberatedstate <- copynuminfo[[3]]
    cluster <- 2
    tau <- runif(1)
    newclusterprob <- 1
    clusterdescriptions <- c(1,2)
    taus <- c(1, tau)
    metalineage <- c(1, 2)
    for (division in nrow(bigtable):(nrow(bigtable)+num_generations)) {
      #Part 1 - Mutate a Locus
      matorpat <- sample(1:2, 1)
      chromosome <- sample(1:22, 1)
      xs <- append(xs, (chromosome*2 - (abs(matorpat-2))))
      locus <- sample(1:ChrLengths[chromosome], 1)
      ys <- append(ys, locus/1000)
      #points(c(chromosome*2 - (abs(matorpat-2))), c(locus/1000), pch = 4, col = "red")
      ID <- ID <- paste0(c("METASTASIS", chrcopy[matorpat],"_Chr", as.character(chromosome), 
                           "_", as.character(locus), ":"), collapse="")
      outtable[division, 1:4] <- c(ID, chromosome, locus, matorpat)
      outtable[division, 7] <- 1
      if(chromosome %in% aberratedchrs){
        if(locus %in% abberatedloci[[chromosome]]) {
          row <- arrayunfolder(abberatedstate, chromosome, locus)
          outtable[division, 5:6] <- abberatedstate[row,3:4,chromosome]
        }
      }
      else {
        outtable[division, 5:6] <- c(1,1)
      }
      outtable[division,8] <- tau
      outtable[division,9] <- cluster
      outtable[division,10] <- paste0(c(LETTERS[lineage], metalineage), collapse = "")
      
      #Part 2 - If a CNA randomly occurs this division, execute it
      if (runif(1) > CNAprob) { #Copy number aberration occurs
        chromosome <- sample(1:22,1)
        abberlength <- round(rnorm(1, mean = 250000, sd = 25000))
        end <- sample(1:(ChrLengths[chromosome]), 1)
        start <- abs(end - abberlength)
        loci <- start:end
        CNA <- generateCNA()
        if (chromosome %in% outtable[,2]) {
          rows <- which(outtable[,2] %in% chromosome) 
          mutatedloci <- as.numeric(outtable[rows,3])
          if (any(mutatedloci %in% loci)) {
            #print(paste0(c("Copy Number Aberration has Occurred on Chromosome ", 
                           #chromosome, " from position ", start, " to ", 
                           #end, " and at least one SNV has been affected"), collapse = ""))
            rowsforediting <- rows[mutatedloci %in% loci]
            for (row in rowsforediting) {
              outtable[row, 5:6] <- as.numeric(outtable[row, 5:6])*CNA[1:2]
              matpat <- as.numeric(outtable[row,4])
              outtable[row, 7] <- as.numeric(outtable[row,7])*CNA[matpat]
            }
          }
        }
        #Add copynum loci to array keeping track
        abberatedchrs <- append(aberratedchrs, chromosome)
        abberatedloci[[chromosome]] <- append(abberatedloci[[chromosome]], loci)
        newrow <- array(0, dim=c(1,4,22))
        abberatedstate <- abind(abberatedstate, newrow, along=1)
        rowindex <- nrow(abberatedstate[,,chromosome])
        abberatedstate[rowindex,,chromosome] <- c(start, end, CNA[1], CNA[2])
      }
      CNAprob <- CNAprob*0.99
      if(runif(1) > newclusterprob) {
        if (cluster == 1) {
          tau <- sample(1:(10*tau), 1)/10
          cluster <- cluster + 1
          newclusterprob <- 1
          clusterdescriptions <- append(clusterdescriptions, 2)
          taus <- append(taus, tau)
          metalineage <- append(metalineage, cluster)
        }
        else {
          #decision is the new cluster a sister or a daughter?
          if (runif(1) < 0.5) { #new cluster is a sister
            currentclusterdesc <- clusterdescriptions[length(clusterdescriptions)] -1
            lastlevel <- which(clusterdescriptions %in% currentclusterdesc)
            parent <- lastlevel[length(lastlevel)]
            formertau <- taus[parent]
            #print(formertau)
            useduptau <- sum(taus[which(clusterdescriptions == clusterdescriptions[length(clusterdescriptions)])])
            #print(useduptau)
            tauchoice <- formertau-useduptau
            if(tauchoice < 0.1) {
              currentclusterdesc <- clusterdescriptions[length(clusterdescriptions)]
              formertau <- taus[match(currentclusterdesc, clusterdescriptions)]
              tau <- sample(1:(10*formertau), 1)/10
              cluster <- cluster + 1
              newclusterprob <- 1
              clusterdescriptions <- append(clusterdescriptions, (clusterdescriptions[length(clusterdescriptions)]+1))
              taus <- append(taus, tau)
              metalineage <- append(metalineage, cluster)
            }
            else {
            tau <- sample(1:(10*tauchoice), 1)/10
            cluster <- cluster + 1
            newclusterprob <- 1
            clusterdescriptions <- append(clusterdescriptions, clusterdescriptions[length(clusterdescriptions)])
            taus <- append(taus, tau)
            onestepback <- length(metalineage) - 1
            metalineage <- append(metalineage[1:onestepback], cluster)
            }
          }
          else { #newclusterisadaughter
            currentclusterdesc <- clusterdescriptions[length(clusterdescriptions)]
            formertau <- taus[match(currentclusterdesc, clusterdescriptions)]
            tau <- sample(1:(10*formertau), 1)/10
            cluster <- cluster + 1
            newclusterprob <- 1
            clusterdescriptions <- append(clusterdescriptions, (clusterdescriptions[length(clusterdescriptions)]+1))
            taus <- append(taus, tau)
            metalineage <- append(metalineage, cluster)
          }
        }
        #print(clusterdescriptions)
      }
      newclusterprob <- newclusterprob*0.99999
      incProgress(1/(num_generations+num_metastasis))
    }
    outtable <- na.omit(outtable)
    #print("METASTASIS SIMULATION COMPLETE")
    return(list(outtable, clusterdescriptions, taus, xs, ys))
  }
  
  #par(mar=c(0,0,0,0), pty = "m")
  #plot(c(1, 1), c(0, ChrLengths[1]/1000), type = "l", bty = "n", 
    #   lwd = 3, col = "Navy", xlim = c(0, 45), yaxt = "n", xaxt = "n", ylab = "",
    #   xlab = "")
  #count <- 1.01
  #cols <- c("pink", "navy")
  #for(x in 2:44) {
  #  index <- round(count)
  #  lines(c(x, x), c(0, ChrLengths[index]/1000), lwd = 3, col = cols[(x%%2)+1])
  #  count <- count + 0.5
  #}
  xs <- c()
  ys <- c()
  outtable <- matrix(NA, nrow = num_generations, ncol = 10)
  CNAprob <- 0.90
  mutationprob <- 0.001
  chrcopy <- c("mat", "pat")
  aberratedchrs <- c()
  abberatedloci <- vector("list", 22)
  abberatedstate <- array(0, dim = c(0,4,22))
  cluster <- 1
  tau <- 1
  newclusterprob <- 1
  clusterdescriptions <- c(1)
  taus <- c(1)
  lineage <- c(1)
  for (division in 1:num_generations) {
    #Does Metastatis happen this division?
    if (division == generationmetastasis) {
      branchlineage <- paste0(LETTERS[lineage], collapse = "")
      inputtable <- data.frame(outtable) %>%
        filter(X9 %in% lineage)
      metabranchpoint <- cluster
      inputtable <- as.matrix(inputtable)
      meta <- metastasis(inputtable, list(aberratedchrs, abberatedloci, abberatedstate),num_metastasis, CNAprob, lineage)
      xs <- c(xs, meta[[4]]) 
      ys <- c(ys, meta[[5]])
      newclusterprob <- 0
      }
    #Part 1 - Mutate a Locus
    matorpat <- sample(1:2, 1)
    chromosome <- sample(1:22, 1)
    locus <- sample(1:ChrLengths[chromosome], 1)
    xs <- append(xs, (chromosome*2 - (abs(matorpat-2))))
    ys <- append(ys, (locus/1000))
    #points(c(chromosome*2 - (abs(matorpat-2))), c(locus/1000), pch = 4, col = "red")
    ID <- ID <- paste0(c(chrcopy[matorpat],"_Chr", as.character(chromosome), "_",
                         as.character(locus), ":"), collapse="")
    outtable[division, 1:4] <- c(ID, chromosome, locus, matorpat)
    outtable[division, 7] <- 1
    if(chromosome %in% aberratedchrs){
      if(locus %in% abberatedloci[[chromosome]]) {
        row <- arrayunfolder(abberatedstate, chromosome, locus)
        outtable[division, 5:6] <- abberatedstate[row,3:4,chromosome]
      }
    }
    else {
      outtable[division, 5:6] <- c(1,1)
    }
    outtable[division,8] <- tau
    outtable[division,9] <- cluster
    outtable[division,10] <- paste0(LETTERS[lineage], collapse = "")
    
    #Part 2 - If a CNA randomly occurs this division, execute it
    if (runif(1) > CNAprob) { #Copy number aberration occurs
      chromosome <- sample(1:22,1)
      abberlength <- round(rnorm(1, mean = 250000, sd = 25000))
      end <- sample(1:(ChrLengths[chromosome]), 1)
      start <- abs(end - abberlength)
      loci <- start:end
      CNA <- generateCNA()
      if (chromosome %in% outtable[,2]) {
        rows <- which(outtable[,2] %in% chromosome) 
        mutatedloci <- as.numeric(outtable[rows,3])
        if (any(mutatedloci %in% loci)) {
          #print(paste0(c("Copy Number Aberration has Occurred on Chromosome ", 
                         #chromosome, " from position ", start, " to ", 
                         #end, " and at least one SNV has been affected"), collapse = ""))
          rowsforediting <- rows[mutatedloci %in% loci]
          for (row in rowsforediting) {
            outtable[row, 5:6] <- as.numeric(outtable[row, 5:6])*CNA[1:2]
            matpat <- as.numeric(outtable[row,4])
            outtable[row, 7] <- as.numeric(outtable[row,7])*CNA[matpat]
          }
        }
      }
      #Add copynum loci to array keeping track
      abberatedchrs <- append(aberratedchrs, chromosome)
      abberatedloci[[chromosome]] <- append(abberatedloci[[chromosome]], loci)
      newrow <- array(0, dim=c(1,4,22))
      abberatedstate <- abind(abberatedstate, newrow, along=1)
      rowindex <- nrow(abberatedstate[,,chromosome])
      abberatedstate[rowindex,,chromosome] <- c(start, end, CNA[1], CNA[2])
    }
    CNAprob <- CNAprob*0.99
    if(runif(1) > newclusterprob) {
      if (cluster == 1) {
        tau <- sample(1:(10*tau), 1)/10
        cluster <- cluster + 1
        newclusterprob <- 1
        clusterdescriptions <- append(clusterdescriptions, 2)
        taus <- append(taus, tau)
        lineage <- append(lineage, cluster)
      }
      else {
        #decision is the new cluster a sister or a daughter?
        if (runif(1) < 0.5) { #new cluster is a sister
          currentclusterdesc <- clusterdescriptions[length(clusterdescriptions)] -1
          lastlevel <- which(clusterdescriptions %in% currentclusterdesc)
          parent <- lastlevel[length(lastlevel)]
          formertau <- taus[parent]
          useduptau <- sum(taus[which(clusterdescriptions == clusterdescriptions[length(clusterdescriptions)])])
          tauchoice <- formertau-useduptau
          if(tauchoice < 0.1) { #not enough tau left, new cluster must be daughter
            currentclusterdesc <- clusterdescriptions[length(clusterdescriptions)]
            formertau <- taus[match(currentclusterdesc, clusterdescriptions)]
            tau <- sample(1:(10*tau), 1)/10
            cluster <- cluster + 1
            newclusterprob <- 1
            clusterdescriptions <- append(clusterdescriptions, (clusterdescriptions[length(clusterdescriptions)]+1))
            taus <- append(taus, tau)
            lineage <- append(lineage, cluster)
          }
          else {
          tau <- sample(1:(10*tauchoice), 1)/10
          cluster <- cluster + 1
          newclusterprob <- 1
          clusterdescriptions <- append(clusterdescriptions, clusterdescriptions[length(clusterdescriptions)])
          taus <- append(taus, tau)
          onestepback <- length(lineage) - 1
          lineage <- append(lineage[1:onestepback], cluster)
          }
        }
        else { #newclusterisadaughter
          currentclusterdesc <- clusterdescriptions[length(clusterdescriptions)]
          formertau <- taus[match(currentclusterdesc, clusterdescriptions)]
          tau <- sample(1:(10*tau), 1)/10
          cluster <- cluster + 1
          newclusterprob <- 1
          clusterdescriptions <- append(clusterdescriptions, (clusterdescriptions[length(clusterdescriptions)]+1))
          taus <- append(taus, tau)
          lineage <- append(lineage, cluster)
        }
      }
      #print(clusterdescriptions)
    }
    newclusterprob <- newclusterprob*0.99999
    incProgress(1/(num_generations+num_metastasis))
  }
  #chromplot <- recordPlot()
  treeplotter <- function(clusterdescriptions, taus) {
    colors <- heat.colors(15)
    par(mar = c(0,0,0,0))
    plot(0.45, (max(clusterdescriptions)+1), xlim = c(0,1), 
         ylim = c(0,(max(clusterdescriptions)+2)), pch = 19, cex = 7, col = "grey")
    prev <- c(0.45, (max(clusterdescriptions)+1))
    count <- 1
    for (i in unique(clusterdescriptions)) {
      y <- max(clusterdescriptions) - i + 1
      numatthislevel <- sum(clusterdescriptions==i)
      space <- (0.9/i) / (numatthislevel +1)
      for(j in 1:numatthislevel) {
        x <- j*space
        lines(c(prev[1], x), c(prev[2], y))
        tau <- 11 - (taus[count]*10)
        points(x, y, col = colors[tau], pch = 19, cex = 7)
        count <- count + 1
      }
      prev <- c(space, y)
      
    }
  }
  addcountsandvaf <- function(table, purity) {
    seqerror <- 0.001
    totalcountsfinder <- function(row, table, purity) {
      totalreads <- abs(round((purity*(sum(as.numeric(table[row, 5:6]))) + (1-purity)*2)*50))
      actualtotalreads <- round(rnorm(1, mean = totalreads, sd = 25))
    }
    variantcountsfinder <- function(row, table, purity, totals) {
      if (table[row, (4+as.numeric(table[row,4]))] == as.numeric(table[row, 7])) {
        n <- as.numeric(table[row, 7])
      }
      else {
        n <- 2
      }
      variantprobnumerator <- ((2*(1-purity)*seqerror) + (purity*(1-purity)*seqerror*n)+(purity*as.numeric(table[row,8])*as.numeric(table[row,7])))
      variantprobdenominator <- ((2*(1-purity)) + (purity*(1-purity)*n)+(purity*(sum(as.numeric(table[row, 5:6])))))
      variantprob <- variantprobnumerator/variantprobdenominator
      vafcount <- rbinom(1, totals[[row]], variantprob)
    }
    totals <- lapply(1:nrow(table), totalcountsfinder, table, purity)
    variants <- lapply(1:nrow(table), variantcountsfinder, table, purity, totals)
    vaf <- unlist(variants)/unlist(totals)
    newtable <- data.frame(table, total = unlist(totals), varcounts = unlist(variants), vaf = vaf)
  }
  treeplotterwithmeta <- function(clusterdescriptions, taus, metabranchpoint,metaobject, clusterlabels, metaclusterlabels){
    par(mar = c(0,0,0,0), pty = "m")
    metaclusterdescriptions <- metaobject[[2]]
    metataus <- metaobject[[3]]
    metatable <- metaobject[[1]]
    graphheight <- max(max(clusterdescriptions), max(metaclusterdescriptions)) 
    colors <- heat.colors(15)
    par(mar = c(0,0,0,0))
    plot(0.5, (graphheight+1), xlim = c(0,1), bty = "n", yaxt = "n", xaxt = "n",
         ylim = c(0,(graphheight+2)), pch = 19, cex = 7, col = "grey")
    prev <- c(0.5, (graphheight+1))
    count <- 1
    primaryxs <- c()
    primaryys <- c()
    metaxs <- c()
    metays <- c()
    for (i in unique(clusterdescriptions)) {
      y <- max(clusterdescriptions) - i + 1
      numatthislevel <- sum(clusterdescriptions==i)
      space <- (0.5) / (numatthislevel +1)
      if (count == 1) {
        space = 0.5
      }
      for(j in 1:numatthislevel) {
        x <- j*space
        lines(c(prev[1], x), c(prev[2], y))
        tau <- 11 - (taus[count]*10)
        points(x, y, col = colors[tau], pch = 19, cex = 10)
        text(x,y, clusterlabels[count])
        if (count == metabranchpoint) {
          linkcoords <- c(x,y)
        }
        primaryxs[count] <- x
        primaryys[count] <- y
        count <- count + 1
      }
      prev <- c(numatthislevel*space, y)
    }
    prev <- linkcoords
    count <- 1
    for (i in unique(metaclusterdescriptions)) {
      y <- graphheight - i + 1
      numatthislevel <- sum(metaclusterdescriptions==i)
      space <- (0.5/i) / (numatthislevel +1)
      for(j in 1:numatthislevel) {
        x <- 1 -  (j*space)
        if (count == 1) {
          lines(c(prev[1], x), c(prev[2], y), lty=2)
        }
        else {
          lines(c(prev[1], x), c(prev[2], y))
        }
        tau <- 11 - (metataus[count]*10)
        points(x, y, col = colors[tau], pch = 19, cex = 10)
        text(x,y, metaclusterlabels[count])
        metaxs[count] <- x
        metays[count] <- y
        count <- count + 1
      }
      prev <- c((1- space), y)
    }
    for (i in 1:length(metaxs)){
      tau <- 11 - (metataus[i]*10)
      points(metaxs[i], metays[i], col = colors[tau], pch = 19, cex = 10)
      text(metaxs[i], metays[i], metaclusterlabels[i])
    }
    for (i in 1:length(primaryxs)){
      tau <- 11 - (taus[i]*10)
      points(primaryxs[i], primaryys[i], col = colors[tau], pch = 19, cex = 10)
      text(primaryxs[i], primaryys[i], clusterlabels[i])
    }
    lines(c(linkcoords[1], metaxs[1]), c(linkcoords[2], metays[1]), lty = 2)
  }
  renameclusters <- function(clusterlist, branchlineage, primaryexpansions) {
    returnvector <- c()
    for (i in clusterlist) {
      numbers <- unlist(str_split(i, branchlineage))
      decomposed <- suppressWarnings(na.omit(as.numeric(unlist(str_split(numbers, "")))))
      decomposed <- decomposed + (primaryexpansions-1)
      decomposed <- na.omit(decomposed[2:length(decomposed)])
      label <- paste0(c(branchlineage, LETTERS[decomposed]), collapse = "")
      returnvector <- append(returnvector, label)
    }
    return(returnvector)
  }
  ccfplotter <- function(primary, meta) {
    primarylite <- primary %>%
      transmute(clusterlabel=X10, PrimaryCCF = X8) %>%
      distinct()
    for (i in 1:nrow(primarylite)) {
      expansions <- unlist(str_split(primarylite[i,1], ""))
      primarylite[i,1] <- expansions[length(expansions)]
    }
    metaexpansions <- (unique(unlist(str_split(meta$X10, ""))))
    primaryexpansions <- primarylite[,1]
    allexpansions <- length(unique(c(primaryexpansions, metaexpansions)))
    outtable <- matrix(0, nrow = allexpansions, ncol = 3)
    outtable[,1] <- LETTERS[1:allexpansions]
    outtable[1:nrow(primarylite),2] <- primarylite[,2]
    
    metalite <- meta %>%
      transmute(clusterlabel=X10, PrimaryCCF = X8) %>%
      distinct()
    seeds <- unlist(str_split(metalite$clusterlabel[1], ""))
    for (i in seeds) {
      row <- match(i, outtable[,1])
      outtable[row, 3] <- 1
    }
    for (i in 1:nrow(metalite)) {
      expansions <- unlist(str_split(metalite[i,1], ""))
      expansion <- expansions[length(expansions)]
      row <- match(expansion, outtable[,1])
      outtable[row,3] <- metalite[i, 2]
    }
    par(mar =c(4,4,2,2), pty = "s")
    plot(outtable[,2], outtable[,3], col = alpha(brewer.pal(nrow(outtable), "Set3"), 0.7),
         pch = 19, cex = 7, xlab = "Primary CCF", ylab = "Metastasis CCF",
         xlim = c(0,1), ylim = c(0,1))
    xcords <- as.numeric(outtable[,2])
    ycords <- as.numeric(outtable[,3])
    labs <- as.character(outtable[,1])
    text(xcords, ycords, labs)
  }
  primvsmeta_vaf <- function(primary, meta) {
    primarylite <- primary %>%
      transmute(ID = X1, clusterlabel=X9, VAF = vaf)
    metalite <- meta %>%
      transmute(ID = X1, clusterlabel=X9, VAF = vaf)
    mutations <- unique(c(primarylite$ID, metalite$ID))
    totalmutations <- length(mutations)
    mainclusters <- length(unique(primarylite$clusterlabel))
    meta <- FALSE
    count <- 1
    metalite[grepl("METASTASIS", metalite$ID), 2] <- as.numeric(metalite[grepl("METASTASIS", metalite$ID), 2]) + mainclusters -1
    outtable <- matrix(0, nrow = totalmutations, ncol = 5)
    outtable[,1] <- mutations
    for (i in 1:nrow(primarylite)) {
      ID <- primarylite[i, 1]
      row <- match(ID, mutations)
      outtable[row, 2] <- primarylite[i, 2]
      outtable[row, 3] <- primarylite[i, 3]
    }
    for (i in 1:nrow(metalite)) {
      ID <- metalite[i, 1]
      row <- match(ID, mutations)
      outtable[row, 2] <- metalite[i, 2]
      outtable[row, 4] <- metalite[i, 3]
    }
    totalclusters <- length(unique(outtable[,2]))
    par(pty="s")
    colors1 <- brewer.pal(totalclusters, "Set3")
    for (i in 1:nrow(outtable)){
      cluster <- as.numeric(outtable[i,2])
      outtable[i,5] <- colors1[cluster]
    }
    plot(outtable[,3], outtable[,4], col = outtable[,5], 
         ylab = "VAF in Metastasis Sample", xlab = "VAF in Primary Sample",
         ylim = c(0,1), xlim = c(0,1))
  }
  
  outtable <- addcountsandvaf(outtable, purity)
  meta[[1]] <- addcountsandvaf(meta[[1]], purity)
  clusterlabels <- unique(outtable$X10)
  primaryexpansions <- length(unique(unlist(str_split(clusterlabels, ""))))
  metatable <- meta[[1]]
  metaclusterlabels <- unique(renameclusters(metatable$X10, branchlineage, primaryexpansions))
  metatable$X10 <- renameclusters(metatable$X10, branchlineage, primaryexpansions)
  meta[[1]] <- metatable
  #treeplotterwithmeta(clusterdescriptions, taus, metabranchpoint, meta, clusterlabels, metaclusterlabels)
  treeinfo <- list(clusterdescriptions, taus, metabranchpoint, meta, clusterlabels, metaclusterlabels)
  #treeplot <- recordPlot()
  #ccfplotter(outtable, metatable)
  #ccfplot <- recordPlot()
  #primvsmeta_vaf(outtable, metatable)
  #vafplot <- recordPlot()
  return(list(outtable, clusterlabels, meta, treeinfo, xs, ys))
}




evofreezesimulator <- function(num_generations, generationmetastasis, num_metastasis, purity) {
  par(new)
  library(stringr)
  library(dplyr)
  library(RColorBrewer)
  library(abind)
  library(scales)
  ChrLengths <- c(248956422, 242193529, 198295559, 190214555, 181538259,
                  170805979, 159345973, 145138636, 138394717, 133797422,
                  135086662, 133275309, 114364328, 107043718, 101991189,
                  90338345, 83257441, 80373285, 58617616, 64444167, 
                  46709983, 50818468)
  generateCNA <- function(){
    maternal <- rpois(1, 3)
    paternal <- rpois(1, 3)
    if (maternal==paternal) {
      if (paternal ==1){
        maternal <- 0
      }
      else if (paternal == 0){
        maternal <- 1
      }
    }
    total <- maternal + paternal
    copynum <- c(maternal, paternal, total)
    return(copynum)
  }
  arrayunfolder <- function(array, z, needle) {
    index <- NA
    for (i in 1:nrow(array[,,z])) {
      if(needle %in% array[i, 1, z]:array[i, 2, z]) {
        index <- i
        break
      }
    }
    return(index)
  }
  metastasis <- function(bigtable, copynuminfo, num_generations, cnaprob, lineage) {
    ChrLengths <- c(248956422, 242193529, 198295559, 190214555, 181538259,
                    170805979, 159345973, 145138636, 138394717, 133797422,
                    135086662, 133275309, 114364328, 107043718, 101991189,
                    90338345, 83257441, 80373285, 58617616, 64444167, 
                    46709983, 50818468)
    
    generateCNA <- function(){
      maternal <- rpois(1, 3)
      paternal <- rpois(1, 3)
      if (maternal==paternal) {
        if (paternal ==1){
          maternal <- 0
        }
        else if (paternal == 0){
          maternal <- 1
        }
      }
      total <- maternal + paternal
      copynum <- c(maternal, paternal, total)
      return(copynum)
    }
    
    arrayunfolder <- function(array, z, needle) {
      index <- NA
      for (i in 1:nrow(array[,,z])) {
        if(needle %in% array[i, 1, z]:array[i, 2, z]) {
          index <- i
          break
        }
      }
      return(index)
    }
    
    #print("METASTASIS OCCURS")
    bigtable[,8] <- 1
    outtable <- matrix(0, nrow = (nrow(bigtable)+num_generations), ncol = 10)
    xs <- c()
    ys <- c()
    CNAprob <- cnaprob
    chrcopy <- c("mat", "pat")
    aberratedchrs <- copynuminfo[[1]]
    abberatedloci <- copynuminfo[[2]]
    abberatedstate <- copynuminfo[[3]]
    cluster <- 2
    tau <- runif(1)
    newclusterprob <- 1
    clusterdescriptions <- c(1,2)
    taus <- c(1, tau)
    metalineage <- c(1, 2)
    for (division in nrow(bigtable):(nrow(bigtable)+num_generations)) {
      #Part 1 - Mutate a Locus
      matorpat <- sample(1:2, 1)
      chromosome <- sample(1:22, 1)
      locus <- sample(1:ChrLengths[chromosome], 1)
      xs <- append(xs, (chromosome*2 - (abs(matorpat-2))))
      ys <- append(ys, (locus/1000))
      #points(c(chromosome*2 - (abs(matorpat-2))), c(locus/1000), pch = 4, col = "red")
      ID <- ID <- paste0(c("METASTASIS", chrcopy[matorpat],"_Chr", as.character(chromosome), 
                           "_", as.character(locus), ":"), collapse="")
      outtable[division, 1:4] <- c(ID, chromosome, locus, matorpat)
      outtable[division, 7] <- 1
      if(chromosome %in% aberratedchrs){
        if(locus %in% abberatedloci[[chromosome]]) {
          row <- arrayunfolder(abberatedstate, chromosome, locus)
          outtable[division, 5:6] <- abberatedstate[row,3:4,chromosome]
        }
      }
      else {
        outtable[division, 5:6] <- c(1,1)
      }
      outtable[division,8] <- tau
      outtable[division,9] <- cluster
      outtable[division,10] <- paste0(c(LETTERS[lineage], metalineage), collapse = "")
      
      #Part 2 - If a CNA randomly occurs this division, execute it
      if (runif(1) > CNAprob) { #Copy number aberration occurs
        chromosome <- sample(1:22,1)
        abberlength <- round(rnorm(1, mean = 250000, sd = 25000))
        end <- sample(1:(ChrLengths[chromosome]), 1)
        start <- abs(end - abberlength)
        loci <- start:end
        CNA <- generateCNA()
        if (chromosome %in% outtable[,2]) {
          rows <- which(outtable[,2] %in% chromosome) 
          mutatedloci <- as.numeric(outtable[rows,3])
          if (any(mutatedloci %in% loci)) {
            #print(paste0(c("Copy Number Aberration has Occurred on Chromosome ", 
                         #  chromosome, " from position ", start, " to ", 
                         #  end, " and at least one SNV has been affected"), collapse = ""))
            rowsforediting <- rows[mutatedloci %in% loci]
            for (row in rowsforediting) {
              outtable[row, 5:6] <- as.numeric(outtable[row, 5:6])*CNA[1:2]
              matpat <- as.numeric(outtable[row,4])
              outtable[row, 7] <- as.numeric(outtable[row,7])*CNA[matpat]
            }
          }
        }
        #Add copynum loci to array keeping track
        abberatedchrs <- append(aberratedchrs, chromosome)
        abberatedloci[[chromosome]] <- append(abberatedloci[[chromosome]], loci)
        newrow <- array(0, dim=c(1,4,22))
        abberatedstate <- abind(abberatedstate, newrow, along=1)
        rowindex <- nrow(abberatedstate[,,chromosome])
        abberatedstate[rowindex,,chromosome] <- c(start, end, CNA[1], CNA[2])
      }
      CNAprob <- CNAprob*0.99
      if(runif(1) > newclusterprob) {
        if (cluster == 1) {
          tau <- sample(1:(10*tau), 1)/10
          cluster <- cluster + 1
          newclusterprob <- 1
          clusterdescriptions <- append(clusterdescriptions, 2)
          taus <- append(taus, tau)
          metalineage <- append(metalineage, cluster)
        }
        else {
          #decision is the new cluster a sister or a daughter?
          if (runif(1) < 0.5) { #new cluster is a sister
            currentclusterdesc <- clusterdescriptions[length(clusterdescriptions)] -1
            lastlevel <- which(clusterdescriptions %in% currentclusterdesc)
            parent <- lastlevel[length(lastlevel)]
            formertau <- taus[parent]
            #print(formertau)
            useduptau <- sum(taus[which(clusterdescriptions == clusterdescriptions[length(clusterdescriptions)])])
            #print(useduptau)
            tauchoice <- formertau-useduptau
            #print(tauchoice)
            tau <- sample(1:(10*tauchoice), 1)/10
            cluster <- cluster + 1
            newclusterprob <- 1
            clusterdescriptions <- append(clusterdescriptions, clusterdescriptions[length(clusterdescriptions)])
            taus <- append(taus, tau)
            onestepback <- length(metalineage) - 1
            metalineage <- append(metalineage[1:onestepback], cluster)
          }
          else { #newclusterisadaughter
            currentclusterdesc <- clusterdescriptions[length(clusterdescriptions)]
            formertau <- taus[match(currentclusterdesc, clusterdescriptions)]
            tau <- sample(1:(10*formertau), 1)/10
            cluster <- cluster + 1
            newclusterprob <- 1
            clusterdescriptions <- append(clusterdescriptions, (clusterdescriptions[length(clusterdescriptions)]+1))
            taus <- append(taus, tau)
            metalineage <- append(metalineage, cluster)
          }
        }
        #print(clusterdescriptions)
      }
      newclusterprob <- newclusterprob*0.999999
    }
    outtable[1:nrow(bigtable),] <- bigtable
    outtable <- na.omit(outtable)
    #print("METASTASIS SIMULATION COMPLETE")
    return(list(outtable, clusterdescriptions, taus, xs, ys))
  }
  
  #par(mar=c(0,0,0,0), pty = "m")
  #plot(c(1, 1), c(0, ChrLengths[1]/1000), type = "l", bty = "n", 
  #     lwd = 3, col = "Navy", xlim = c(0, 45), yaxt = "n", xaxt = "n", ylab = "",
  #     xlab = "")
  #count <- 1.01
  #cols <- c("pink", "navy")
  #for(x in 2:44) {
  #  index <- round(count)
  #  lines(c(x, x), c(0, ChrLengths[index]/1000), lwd = 3, col = cols[(x%%2)+1])
  #  count <- count + 0.5
  #}
  xs <- c()
  ys <- c()
  outtable <- matrix(NA, nrow = num_generations, ncol = 10)
  CNAprob <- 0.40
  mutationprob <- 0.001
  chrcopy <- c("mat", "pat")
  aberratedchrs <- c()
  abberatedloci <- vector("list", 22)
  abberatedstate <- array(0, dim = c(0,4,22))
  cluster <- 1
  tau <- 1
  newclusterprob <- 1
  clusterdescriptions <- c(1)
  taus <- c(1)
  lineage <- c(1)
  for (division in 1:num_generations) {
    #Does Metastatis happen this division?
    if (division == generationmetastasis) {
      branchlineage <- paste0(LETTERS[lineage], collapse = "")
      inputtable <- data.frame(outtable) %>%
        filter(X9 %in% lineage)
      metabranchpoint <- cluster
      inputtable <- as.matrix(inputtable)
      meta <- metastasis(inputtable, list(aberratedchrs, abberatedloci, abberatedstate),num_metastasis, CNAprob, lineage)
      xs <- c(xs, meta[[4]])
      ys <- c(ys, meta[[5]])
      newclusterprob <- 0
      }
    #Part 1 - Mutate a Locus
    matorpat <- sample(1:2, 1)
    chromosome <- sample(1:22, 1)
    locus <- sample(1:ChrLengths[chromosome], 1)
    xs <- append(xs, (chromosome*2 - (abs(matorpat-2))))
    ys <- append(ys, (locus/1000))
    #points(c(chromosome*2 - (abs(matorpat-2))), c(locus/1000), pch = 4, col = "red")
    ID <- ID <- paste0(c(chrcopy[matorpat],"_Chr", as.character(chromosome), "_",
                         as.character(locus), ":"), collapse="")
    outtable[division, 1:4] <- c(ID, chromosome, locus, matorpat)
    outtable[division, 7] <- 1
    if(chromosome %in% aberratedchrs){
      if(locus %in% abberatedloci[[chromosome]]) {
        row <- arrayunfolder(abberatedstate, chromosome, locus)
        outtable[division, 5:6] <- abberatedstate[row,3:4,chromosome]
      }
    }
    else {
      outtable[division, 5:6] <- c(1,1)
    }
    outtable[division,8] <- tau
    outtable[division,9] <- cluster
    outtable[division,10] <- paste0(LETTERS[lineage], collapse = "")
    
    #Part 2 - If a CNA randomly occurs this division, execute it
    if (runif(1) > CNAprob) { #Copy number aberration occurs
      chromosome <- sample(1:22,1)
      abberlength <- round(rnorm(1, mean = 250000, sd = 25000))
      end <- sample(1:(ChrLengths[chromosome]), 1)
      start <- abs(end - abberlength)
      loci <- start:end
      CNA <- generateCNA()
      if (chromosome %in% outtable[,2]) {
        rows <- which(outtable[,2] %in% chromosome) 
        mutatedloci <- as.numeric(outtable[rows,3])
        if (any(mutatedloci %in% loci)) {
          #print(paste0(c("Copy Number Aberration has Occurred on Chromosome ", 
                     #    chromosome, " from position ", start, " to ", 
                     #    end, " and at least one SNV has been affected"), collapse = ""))
          rowsforediting <- rows[mutatedloci %in% loci]
          for (row in rowsforediting) {
            outtable[row, 5:6] <- as.numeric(outtable[row, 5:6])*CNA[1:2]
            matpat <- as.numeric(outtable[row,4])
            outtable[row, 7] <- as.numeric(outtable[row,7])*CNA[matpat]
          }
        }
      }
      #Add copynum loci to array keeping track
      abberatedchrs <- append(aberratedchrs, chromosome)
      abberatedloci[[chromosome]] <- append(abberatedloci[[chromosome]], loci)
      newrow <- array(0, dim=c(1,4,22))
      abberatedstate <- abind(abberatedstate, newrow, along=1)
      rowindex <- nrow(abberatedstate[,,chromosome])
      abberatedstate[rowindex,,chromosome] <- c(start, end, CNA[1], CNA[2])
    }
    CNAprob <- CNAprob*0.99
    if(runif(1) > newclusterprob) {
      if (cluster == 1) {
        tau <- sample(1:(10*tau), 1)/10
        cluster <- cluster + 1
        newclusterprob <- 1
        clusterdescriptions <- append(clusterdescriptions, 2)
        taus <- append(taus, tau)
        lineage <- append(lineage, cluster)
      }
      else {
        #decision is the new cluster a sister or a daughter?
        if (runif(1) < 0.5) { #new cluster is a sister
          currentclusterdesc <- clusterdescriptions[length(clusterdescriptions)] -1
          lastlevel <- which(clusterdescriptions %in% currentclusterdesc)
          parent <- lastlevel[length(lastlevel)]
          formertau <- taus[parent]
          useduptau <- sum(taus[which(clusterdescriptions == clusterdescriptions[length(clusterdescriptions)])])
          tauchoice <- formertau-useduptau
          tau <- sample(1:(10*tauchoice), 1)/10
          cluster <- cluster + 1
          newclusterprob <- 1
          clusterdescriptions <- append(clusterdescriptions, clusterdescriptions[length(clusterdescriptions)])
          taus <- append(taus, tau)
          onestepback <- length(lineage) - 1
          lineage <- append(lineage[1:onestepback], cluster)
        }
        else { #newclusterisadaughter
          currentclusterdesc <- clusterdescriptions[length(clusterdescriptions)]
          formertau <- taus[match(currentclusterdesc, clusterdescriptions)]
          tau <- sample(1:(10*tau), 1)/10
          cluster <- cluster + 1
          newclusterprob <- 1
          clusterdescriptions <- append(clusterdescriptions, (clusterdescriptions[length(clusterdescriptions)]+1))
          taus <- append(taus, tau)
          lineage <- append(lineage, cluster)
        }
      }
      #print(clusterdescriptions)
    }
    newclusterprob <- newclusterprob*0.99999
    incProgress(1/num_generations)
  }
  treeplotter <- function(clusterdescriptions, taus) {
    colors <- heat.colors(15)
    par(mar = c(0,0,0,0))
    plot(0.45, (max(clusterdescriptions)+1), xlim = c(0,1), 
         ylim = c(0,(max(clusterdescriptions)+2)), pch = 19, cex = 7, col = "grey")
    prev <- c(0.45, (max(clusterdescriptions)+1))
    count <- 1
    for (i in unique(clusterdescriptions)) {
      y <- max(clusterdescriptions) - i + 1
      numatthislevel <- sum(clusterdescriptions==i)
      space <- (0.9/i) / (numatthislevel +1)
      for(j in 1:numatthislevel) {
        x <- j*space
        lines(c(prev[1], x), c(prev[2], y))
        tau <- 11 - (taus[count]*10)
        points(x, y, col = colors[tau], pch = 19, cex = 7)
        count <- count + 1
      }
      prev <- c(space, y)
      
    }
  }
  addcountsandvaf <- function(table, purity) {
    seqerror <- 0.001
    totalcountsfinder <- function(row, table, purity) {
      totalreads <- abs(round((purity*(sum(as.numeric(table[row, 5:6]))) + (1-purity)*2)*500))
      actualtotalreads <- round(rnorm(1, mean = totalreads, sd = 25))
    }
    variantcountsfinder <- function(row, table, purity, totals) {
      if (table[row, (4+as.numeric(table[row,4]))] == as.numeric(table[row, 7])) {
        n <- as.numeric(table[row, 7])
      }
      else {
        n <- 2
      }
      variantprobnumerator <- ((2*(1-purity)*seqerror) + (purity*(1-purity)*seqerror*n)+(purity*as.numeric(table[row,8])*as.numeric(table[row,7])))
      variantprobdenominator <- ((2*(1-purity)) + (purity*(1-purity)*n)+(purity*(sum(as.numeric(table[row, 5:6])))))
      variantprob <- variantprobnumerator/variantprobdenominator
      vafcount <- rbinom(1, totals[[row]], variantprob)
    }
    totals <- lapply(1:nrow(table), totalcountsfinder, table, purity)
    variants <- lapply(1:nrow(table), variantcountsfinder, table, purity, totals)
    vaf <- unlist(variants)/unlist(totals)
    newtable <- data.frame(table, total = unlist(totals), varcounts = unlist(variants), vaf = vaf)
  }
  treeplotterwithmeta <- function(clusterdescriptions, taus, metabranchpoint,metaobject, clusterlabels, metaclusterlabels){
    par(mar = c(0,0,0,0), pty = "m")
    metaclusterdescriptions <- metaobject[[2]]
    metataus <- metaobject[[3]]
    metatable <- metaobject[[1]]
    graphheight <- max(max(clusterdescriptions), max(metaclusterdescriptions)) 
    colors <- heat.colors(15)
    par(mar = c(0,0,0,0))
    plot(0.5, (graphheight+1), xlim = c(0,1), 
         ylim = c(0,(graphheight+2)), pch = 19, cex = 7, col = "grey")
    prev <- c(0.5, (graphheight+1))
    count <- 1
    primaryxs <- c()
    primaryys <- c()
    metaxs <- c()
    metays <- c()
    for (i in unique(clusterdescriptions)) {
      y <- max(clusterdescriptions) - i + 1
      numatthislevel <- sum(clusterdescriptions==i)
      space <- (0.5) / (numatthislevel +1)
      if (count == 1) {
        space = 0.5
      }
      for(j in 1:numatthislevel) {
        x <- j*space
        lines(c(prev[1], x), c(prev[2], y))
        tau <- 11 - (taus[count]*10)
        points(x, y, col = colors[tau], pch = 19, cex = 10)
        text(x,y, clusterlabels[count])
        if (count == metabranchpoint) {
          linkcoords <- c(x,y)
        }
        primaryxs[count] <- x
        primaryys[count] <- y
        count <- count + 1
      }
      prev <- c(numatthislevel*space, y)
    }
    prev <- linkcoords
    count <- 1
    for (i in unique(metaclusterdescriptions)) {
      y <- graphheight - i + 1
      numatthislevel <- sum(metaclusterdescriptions==i)
      space <- (0.5/i) / (numatthislevel +1)
      for(j in 1:numatthislevel) {
        x <- 1 -  (j*space)
        if (count == 1) {
          lines(c(prev[1], x), c(prev[2], y), lty=2)
        }
        else {
          lines(c(prev[1], x), c(prev[2], y))
        }
        tau <- 11 - (metataus[count]*10)
        points(x, y, col = colors[tau], pch = 19, cex = 10)
        text(x,y, metaclusterlabels[count])
        metaxs[count] <- x
        metays[count] <- y
        count <- count + 1
      }
      prev <- c((1- space), y)
    }
    for (i in 1:length(metaxs)){
      tau <- 11 - (metataus[i]*10)
      points(metaxs[i], metays[i], col = colors[tau], pch = 19, cex = 10)
      text(metaxs[i], metays[i], metaclusterlabels[i])
    }
    for (i in 1:length(primaryxs)){
      tau <- 11 - (taus[i]*10)
      points(primaryxs[i], primaryys[i], col = colors[tau], pch = 19, cex = 10)
      text(primaryxs[i], primaryys[i], clusterlabels[i])
    }
    lines(c(linkcoords[1], metaxs[1]), c(linkcoords[2], metays[1]), lty = 2)
  }
  
  renameclusters <- function(clusterlist, branchlineage, primaryexpansions) {
    returnvector <- c()
    for (i in clusterlist) {
      numbers <- unlist(str_split(i, branchlineage))
      decomposed <- na.omit(as.numeric(unlist(str_split(numbers, ""))))
      decomposed <- decomposed + (primaryexpansions-1)
      decomposed <- na.omit(decomposed[2:length(decomposed)])
      label <- paste0(c(branchlineage, LETTERS[decomposed]), collapse = "")
      returnvector <- append(returnvector, label)
    }
    return(returnvector)
  }
  ccfplotter <- function(primary, meta) {
    primarylite <- primary %>%
      transmute(clusterlabel=X10, PrimaryCCF = X8) %>%
      distinct()
    for (i in 1:nrow(primarylite)) {
      expansions <- unlist(str_split(primarylite[i,1], ""))
      primarylite[i,1] <- expansions[length(expansions)]
    }
    metaexpansions <- (unique(unlist(str_split(meta$X10, ""))))
    primaryexpansions <- primarylite[,1]
    allexpansions <- length(unique(c(primaryexpansions, metaexpansions)))
    outtable <- matrix(0, nrow = allexpansions, ncol = 3)
    outtable[,1] <- LETTERS[1:allexpansions]
    outtable[1:nrow(primarylite),2] <- primarylite[,2]
    
    metalite <- meta %>%
      transmute(clusterlabel=X10, PrimaryCCF = X8) %>%
      distinct()
    seeds <- unlist(str_split(metalite$clusterlabel[1], ""))
    for (i in seeds) {
      row <- match(i, outtable[,1])
      outtable[row, 3] <- 1
    }
    for (i in 1:nrow(metalite)) {
      expansions <- unlist(str_split(metalite[i,1], ""))
      expansion <- expansions[length(expansions)]
      row <- match(expansion, outtable[,1])
      outtable[row,3] <- metalite[i, 2]
    }
    par(mar =c(4,4,2,2), pty = "s")
    plot(outtable[,2], outtable[,3], col = alpha(brewer.pal(nrow(outtable), "Set3"), 0.7),
         pch = 19, cex = 7, xlab = "Primary CCF", ylab = "Metastasis CCF",
         xlim = c(0,1), ylim = c(0,1))
    xcords <- as.numeric(outtable[,2])
    ycords <- as.numeric(outtable[,3])
    labs <- as.character(outtable[,1])
    text(xcords, ycords, labs)
  }
  primvsmeta_vaf <- function(primary, meta) {
    primarylite <- primary %>%
      transmute(ID = X1, clusterlabel=X9, VAF = vaf)
    metalite <- meta %>%
      transmute(ID = X1, clusterlabel=X9, VAF = vaf)
    mutations <- unique(c(primarylite$ID, metalite$ID))
    totalmutations <- length(mutations)
    mainclusters <- length(unique(primarylite$clusterlabel))
    meta <- FALSE
    count <- 1
    metalite[grepl("METASTASIS", metalite$ID), 2] <- as.numeric(metalite[grepl("METASTASIS", metalite$ID), 2]) + mainclusters -1
    outtable <- matrix(0, nrow = totalmutations, ncol = 5)
    outtable[,1] <- mutations
    for (i in 1:nrow(primarylite)) {
      ID <- primarylite[i, 1]
      row <- match(ID, mutations)
      outtable[row, 2] <- primarylite[i, 2]
      outtable[row, 3] <- primarylite[i, 3]
    }
    for (i in 1:nrow(metalite)) {
      ID <- metalite[i, 1]
      row <- match(ID, mutations)
      outtable[row, 2] <- metalite[i, 2]
      outtable[row, 4] <- metalite[i, 3]
    }
    totalclusters <- length(unique(outtable[,2]))
    par(pty="s")
    colors1 <- brewer.pal(totalclusters, "Set3")
    for (i in 1:nrow(outtable)){
      cluster <- as.numeric(outtable[i,2])
      outtable[i,5] <- colors1[cluster]
    }
    plot(outtable[,3], outtable[,4], col = outtable[,5], 
         ylab = "VAF in Metastasis Sample", xlab = "VAF in Primary Sample",
         ylim = c(0,1), xlim = c(0,1))
  }
  
  chromplot <- recordPlot()
  outtable <- addcountsandvaf(outtable, purity)
  meta[[1]] <- addcountsandvaf(meta[[1]], purity)
  clusterlabels <- unique(outtable$X10)
  primaryexpansions <- length(unique(unlist(str_split(clusterlabels, ""))))
  metatable <- meta[[1]]
  metaclusterlabels <- unique(renameclusters(metatable$X10, branchlineage, primaryexpansions))
  metatable$X10 <- renameclusters(metatable$X10, branchlineage, primaryexpansions)
  meta[[1]] <- metatable
  treeplotterwithmeta(clusterdescriptions, taus, metabranchpoint, meta, clusterlabels, metaclusterlabels)
  treeinfo <- list(clusterdescriptions, taus, metabranchpoint, meta, clusterlabels, metaclusterlabels)
  #ccfplotter(outtable, metatable)
  #ccfplot <- recordPlot()
  #primvsmeta_vaf(outtable, metatable)
  #vafplot <- recordPlot()
  return(list(outtable, clusterlabels, meta, treeinfo, xs, ys))
}

primvsmeta_vaf <- function(robject) {
  primary <- robject[[1]]
  meta <- robject[[3]]
  meta <- meta[[1]]
  primarylite <- primary %>%
    transmute(ID = X1, clusterlabel=X9, VAF = vaf)
  metalite <- meta %>%
    transmute(ID = X1, clusterlabel=X9, VAF = vaf)
  mutations <- unique(c(primarylite$ID, metalite$ID))
  totalmutations <- length(mutations)
  mainclusters <- length(unique(primarylite$clusterlabel))
  meta <- FALSE
  count <- 1
  metalite[grepl("METASTASIS", metalite$ID), 2] <- as.numeric(metalite[grepl("METASTASIS", metalite$ID), 2]) + mainclusters -1
  outtable <- matrix(0, nrow = totalmutations, ncol = 5)
  outtable[,1] <- mutations
  for (i in 1:nrow(primarylite)) {
    ID <- primarylite[i, 1]
    row <- match(ID, mutations)
    outtable[row, 2] <- primarylite[i, 2]
    outtable[row, 3] <- primarylite[i, 3]
  }
  for (i in 1:nrow(metalite)) {
    ID <- metalite[i, 1]
    row <- match(ID, mutations)
    outtable[row, 2] <- metalite[i, 2]
    outtable[row, 4] <- metalite[i, 3]
  }
  totalclusters <- length(unique(outtable[,2]))
  par(pty="s")
  colors1 <- brewer.pal(totalclusters, "Set3")
  for (i in 1:nrow(outtable)){
    cluster <- as.numeric(outtable[i,2])
    outtable[i,5] <- colors1[cluster]
  }
  plot(outtable[,3], outtable[,4], col = outtable[,5], 
       ylab = "VAF in Metastasis Sample", xlab = "VAF in Primary Sample",
       ylim = c(0,1), xlim = c(0,1))
}


ccfplotter <- function(robject) {
  primary <- robject[[1]]
  meta <- robject[[3]]
  meta <- meta[[1]]
  primarylite <- primary %>%
    transmute(clusterlabel=X10, PrimaryCCF = X8) %>%
    distinct()
  for (i in 1:nrow(primarylite)) {
    expansions <- unlist(str_split(primarylite[i,1], ""))
    primarylite[i,1] <- expansions[length(expansions)]
  }
  metaexpansions <- (unique(unlist(str_split(meta$X10, ""))))
  primaryexpansions <- primarylite[,1]
  allexpansions <- length(unique(c(primaryexpansions, metaexpansions)))
  outtable <- matrix(0, nrow = allexpansions, ncol = 3)
  outtable[,1] <- LETTERS[1:allexpansions]
  outtable[1:nrow(primarylite),2] <- primarylite[,2]
  
  metalite <- meta %>%
    transmute(clusterlabel=X10, PrimaryCCF = X8) %>%
    distinct()
  seeds <- unlist(str_split(metalite$clusterlabel[1], ""))
  for (i in seeds) {
    row <- match(i, outtable[,1])
    outtable[row, 3] <- 1
  }
  for (i in 1:nrow(metalite)) {
    expansions <- unlist(str_split(metalite[i,1], ""))
    expansion <- expansions[length(expansions)]
    row <- match(expansion, outtable[,1])
    outtable[row,3] <- metalite[i, 2]
  }
  par(mar =c(4,4,2,2), pty = "s")
  plot(outtable[,2], outtable[,3], col = alpha(brewer.pal(nrow(outtable), "Set3"), 0.7),
       pch = 19, cex = 7, xlab = "Primary CCF", ylab = "Metastasis CCF",
       xlim = c(0,1), ylim = c(0,1))
  xcords <- as.numeric(outtable[,2])
  ycords <- as.numeric(outtable[,3])
  labs <- as.character(outtable[,1])
  text(xcords, ycords, labs)
}



treeplotterwithmeta <- function(clusterdescriptions, taus, metabranchpoint,metaobject, clusterlabels, metaclusterlabels){
  par(mar = c(0,0,0,0), pty = "m")
  metaclusterdescriptions <- metaobject[[2]]
  metataus <- metaobject[[3]]
  metatable <- metaobject[[1]]
  graphheight <- max(max(clusterdescriptions), max(metaclusterdescriptions)) 
  colors <- heat.colors(15)
  par(mar = c(0,0,0,0))
  plot(0.5, (graphheight+1), xlim = c(0,1), bty = "n", yaxt = "n", xaxt = "n",
       ylim = c(0,(graphheight+2)), pch = 19, cex = 7, col = "grey")
  prev <- c(0.5, (graphheight+1))
  count <- 1
  primaryxs <- c()
  primaryys <- c()
  metaxs <- c()
  metays <- c()
  for (i in unique(clusterdescriptions)) {
    y <- max(clusterdescriptions) - i + 1
    numatthislevel <- sum(clusterdescriptions==i)
    space <- (0.5) / (numatthislevel +1)
    if (count == 1) {
      space = 0.5
    }
    for(j in 1:numatthislevel) {
      x <- j*space
      lines(c(prev[1], x), c(prev[2], y))
      tau <- 11 - (taus[count]*10)
      points(x, y, col = colors[tau], pch = 19, cex = 10)
      text(x,y, clusterlabels[count])
      if (count == metabranchpoint) {
        linkcoords <- c(x,y)
      }
      primaryxs[count] <- x
      primaryys[count] <- y
      count <- count + 1
    }
    prev <- c(numatthislevel*space, y)
  }
  prev <- linkcoords
  count <- 1
  for (i in unique(metaclusterdescriptions)) {
    y <- graphheight - i + 1
    numatthislevel <- sum(metaclusterdescriptions==i)
    space <- (0.5/i) / (numatthislevel +1)
    for(j in 1:numatthislevel) {
      x <- 1 -  (j*space)
      if (count == 1) {
        lines(c(prev[1], x), c(prev[2], y), lty=2)
      }
      else {
        lines(c(prev[1], x), c(prev[2], y))
      }
      tau <- 11 - (metataus[count]*10)
      points(x, y, col = colors[tau], pch = 19, cex = 10)
      text(x,y, metaclusterlabels[count])
      metaxs[count] <- x
      metays[count] <- y
      count <- count + 1
    }
    prev <- c((1- space), y)
  }
  for (i in 1:length(metaxs)){
    tau <- 11 - (metataus[i]*10)
    points(metaxs[i], metays[i], col = colors[tau], pch = 19, cex = 10)
    text(metaxs[i], metays[i], metaclusterlabels[i])
  }
  for (i in 1:length(primaryxs)){
    tau <- 11 - (taus[i]*10)
    points(primaryxs[i], primaryys[i], col = colors[tau], pch = 19, cex = 10)
    text(primaryxs[i], primaryys[i], clusterlabels[i])
  }
  lines(c(linkcoords[1], metaxs[1]), c(linkcoords[2], metays[1]), lty = 2)
}




karyotypeplotter <- function(xs, ys) {
  ChrLengths <- c(248956422, 242193529, 198295559, 190214555, 181538259,
                  170805979, 159345973, 145138636, 138394717, 133797422,
                  135086662, 133275309, 114364328, 107043718, 101991189,
                  90338345, 83257441, 80373285, 58617616, 64444167, 
                  46709983, 50818468)
  par(mar=c(0,0,0,0), pty = "m")
  plot(c(1, 1), c(0, ChrLengths[1]/1000), type = "l", bty = "n", 
     lwd = 3, col = "deepskyblue", xlim = c(-1, 45), yaxt = "n", xaxt = "n", ylab = "",
     xlab = "")
  count <- 1.01
  cols <- c("pink", "deepskyblue")
  for(x in 2:44) {
    index <- round(count)
    lines(c(x, x), c(0, ChrLengths[index]/1000), lwd = 3, col = cols[(x%%2)+1])
    count <- count + 0.5
  }
  for (i in 1:length(xs)){
    points(c(xs[i]), c(ys[i]), pch = 4, col = "red")
  }
}





