#' Constructs the cleavage windows images
#'
#' This is an example for how the windows aka cleavage pictures might be created.
#' feel free to use your own pictures. However if using the CNN model we designed it is recommended
#' To adopt this script to your data as other pictures might not be recognized by the model.
#' @param gffTrans - a dataframe of your imported gff with V10 = transcript ID see  for example extendGffTrans in makeWindowPictureFunc for example
#' @param aliFilesPattern1 - regular expression defining the pattern of your first degradome library
#' @param aliFilesPattern2 - regular expression defining the pattern of your second degradome library
#' @param dirO - output path
#' @param cleavageData - must be a dataframe containing at least columns genesT (the target genes) and posT (target position in the transcript)
#' @param edgesExtend1 - how many positions from the cleavage site that are included in the pictures. Default is c(1,2),
#' defining 1 pos upstream and 21 positions downstream. As the degradome reads are 20 nt each this defines
#' 1 position upstream and downstream of the degradome reads. Only one, to try and exclude as much background
#' as possible but still catch the caracteristic "cleavage tower".
#' @param ylim1 - defines the minimum plotted height of the y-axis. In practise this means a lower "cleavage tower"
#' will not be recognized by the CNN. This to exclude potential false positives caused by noice
#' @paramonlyOneTest - sometimes the script is missing one cleavage picture for unknown reasons. Hence,
#' not reporting that all files are created but still not creating any more. If this is the case
#' one can test for this by running with "onlyOneTest = T", default is F.
#' @param addToDir - if one wants to add the pictures to the directory. Default is T. If F the user also need to
#' remove the has in the function script, this as a precaution not to remove any pictures by accident
#' @param aliFilesPath -path to the bamfiles  for the user defined degradome. The easiest is to use transcriptome
#' generated bam files. If you do so the path must end with "bamTranscriptome/". However, this option cannot check
#' for cleavages outside the transcripts (which is totally  fine if you have a well annotated genome).
#' If you have trascripts and a genome but no gff file you can create a gff with the function
#' makeGFF (reqires spaln installed on your system).
#' Also included is the function readGFF which is a basic function to import your GFF.
#' IMPORTANT!!! is that your gffTrans has a 10th column with your the transcript ID. Often this must be
#' extracted from the 9th column (V9). However, most GFF files I've been exposed to are different.
#' In the attached function extendGffTrans you can see how this was achieved from the Jupe dataset (Jupe et al. 2013)
#' I was using.
#' @keywords cleavageWindows
#' @export
#' @examples
#' cleavageWindows(dirO = paste0("pathOut/"),
#' cleavageData = cleavageDataDataset,
#' aliFilesPath = "path/bamTranscriptome/",
#' aliFilesPattern1 = "pattern1.sorted.bam$",
#' aliFilesPattern2 = "pattern2.sorted.bam$",
#' ylim1 = 5,
#' addToDir = T,
#' onlyOneTest = F,
#' edgesExtend1 = c(1,21),
#' gffTrans = gffTrans
#' )


cleavageWindows <- function(dirO = paste0("pathOut/"),
                            cleavageData = cleavageDataDataset,
                            aliFilesPath = "path/bamTranscriptome/",
                            aliFilesPattern1 = "pattern1.sorted.bam$",
                            aliFilesPattern2 = "pattern2.sorted.bam$",
                            ylim1 = 5,
                            addToDir = T,
                            onlyOneTest = F,
                            edgesExtend1 = c(1,21),
                            gffTrans = gffTrans,
                            savePics = T,
                            jpegWidHei = c(480,480),
                            qual = 75,
                            pz = 12
){
  transAndPosToBam <- function(
    PotInPotvsNoInf4 = data.frame(genesT="PGSC0003DMT400019853", posT=(974)),
    edgesExtend1 = 100,
    aliFilesPattern1 = "B_GF(.*?).sorted.bam$",
    aliFilesPath = "bamTranscript/",
    aliFilesPatternCtrl1 = "",
    individual = T,
    ylim1 = 5,
    gffTrans
  ){
    readDepthBam <- function(posT,
                             aliFilesPath = "/home/kristianh/p/Degradome/analysis/bamPot/",
                             aliFilesPattern = "B_GF(.*?).sorted.bam$",
                             edgesExtend,
                             chr1,
                             aliFilesPatternCtrl = "",
                             printPage = "",
                             individual = F,
                             ylim1 = 10,
                             strand3){




      library(Rsamtools)
      #par(mfrow=c(2,2), mar=c(0,0,0,0))

      aliFilesPath2 <- list.files(path = aliFilesPath,pattern = aliFilesPattern,full.names = T)
      printP("Deg:",aliFilesPath2)
      bamFile <- BamFile(file = aliFilesPath2[1])
      seqinfoBam <- seqinfo(bamFile)
      bamChr <- grep(pattern = chr1, x = names(seqinfoBam))
      seqMax <- seqinfoBam[names(seqinfoBam)[bamChr],]
      seqMax2 <- seqlengths(seqMax)

      posT <- facToNum(posT)
      chr1 <- as.character(names(seqMax2))

      printP("Transcript",chr1)
      if(length(edgesExtend) == 1){
        edgesExtend <- c(edgesExtend,edgesExtend)
      }

      if(strand3 == "+"){
        intStart <- posT - edgesExtend[1]
        intEnd <- posT + edgesExtend[2]
      }else if (strand3 == "-"){
        intStart <- posT - edgesExtend[2]
        intEnd <- posT + edgesExtend[1]
      }


      if(intStart < 1){
        intStart <- 1
      }
      if(intEnd > seqMax2){
        intEnd <- seqMax2
      }

      mygene <- GRanges(chr1, ranges=IRanges(intStart, intEnd))
      z1 <- GRanges(chr1,IRanges(intStart,intEnd),strand="+" )
      z2 <- GRanges(chr1,IRanges(intStart,intEnd),strand="-" )

      plotBamGene <- function(aliFile,
                              mygene,
                              linePlot = F,
                              edgesExtend = 100,
                              posT,
                              ylim1 = ylim1,
                              strand3){

        mygene.reads <- readGAlignments(file=aliFile,
                                        param=ScanBamParam(which=mygene,
                                                           what=c("seq","mapq","flag","qual","isize")
                                        )
        )

        mygene.readsP <- mygene.reads[which(strand(mygene.reads) == "+")]
        mygene.readsM <- mygene.reads[which(strand(mygene.reads) == "-")]

        covP <- coverage(mygene.readsP)
        numP <- as.numeric(covP[[chr1]][ranges(z1)])

        covM <- coverage(mygene.readsM)
        numM <- as.numeric(covM[[chr1]][ranges(z2)])

        lim1 <- max(numM,numP)
        if(lim1 < ylim1){
          lim1 <- ylim1
        }
        if(strand3 == "+"){
          xlim1 <- posT - edgesExtend[1]
          xlim2 <- posT + edgesExtend[2]
          posTend <- posT + (edgesExtend[2]-edgesExtend[1])
        }else if (strand3 == "-"){
          xlim1 <- posT - edgesExtend[2]
          xlim2 <- posT + edgesExtend[1]
          posTend <- posT - (edgesExtend[2]-edgesExtend[1])
        }

        if(linePlot == T){
          plot(numP, type="l", col="blue", lwd=2,ylim=c(-lim1,lim1))
          lines(-numM, col="red", lwd=2)
          abline(v = length(numM)/2)
        }else{
          basenam <- basename(aliFile)

          pos <- data.frame(Chr=rep(chr1, length(numP)), Strand=rep("+", length(numP)), Position= intStart:intEnd, Coverage= numP,
                            bam = rep(basenam, length(numP)))
          neg <- data.frame(Chr=rep(chr1, length(numM)), Strand=rep("-", length(numM)), Position= intStart:intEnd, Coverage= -numM,
                            bam = rep(basenam, length(numP)))

          covdf <- rbind(neg,pos)
          covdf$Strand <- factor(covdf$Strand,levels = c("+", "-"))
          p <- ggplot(covdf, aes(Position, Coverage, fill=Strand)) +
            geom_bar(stat="identity", position="identity")+ # heading facet_wrap(~bam)+
            geom_vline(xintercept = posT, color ="black", size = .3)+
            geom_vline(xintercept = posTend, color ="black", size = .3)+
            geom_hline(yintercept = 0, color ="black", size = .3) +
            ylim(-lim1, lim1) + xlim(xlim1-1, xlim2+1)  +
            theme_classic() +theme(legend.position = "none")
          return(p)
        }



      }

      if(individual == T){
        p <- lapply(aliFilesPath2, function(x){
          plotBamGene(aliFile = x,
                      mygene = mygene,
                      edgesExtend = edgesExtend,
                      posT = posT,
                      ylim1 = ylim1,
                      strand3 = strand3)
        }
        )
        if(aliFilesPatternCtrl != ""){
          aliFilesPathC <- list.files(path = aliFilesPath,pattern = aliFilesPatternCtrl,full.names = T)
          p1 <- lapply(aliFilesPathC, function(x){
            plotBamGene(aliFile = x,mygene, edgesExtend = edgesExtend, posT = posT, ylim1 = ylim1)
          }
          )
          p <-c(p,p1)
        }
        print(p)
      }else{

        p <- lapply(aliFilesPath2, function(x){
          plotBamGene(aliFile = x, edgesExtend = edgesExtend, posT = posT, ylim1 = ylim1)
        }
        )
        if(aliFilesPatternCtrl != ""){
          print("Ctrl found")
          aliFilesPathC <- list.files(path = aliFilesPath,pattern = aliFilesPatternCtrl,full.names = T)
          p1 <- lapply(aliFilesPathC, function(x){
            plotBamGene(aliFile = x,mygene,
                        edgesExtend = edgesExtend,
                        posT = posT,
                        ylim1 = ylim1,
                        strand3 = strand3)
          }
          )
          p2 <- c(p,p1)
          library(gridExtra)
          lenP <- length(p2)

          lay1 <- rbind(c(1,1),
                        cbind(2:((lenP/2)+1),
                              (lenP/2+1): lenP+1))

          mytheme <- gridExtra::ttheme_default(
            core = list(fg_params=list(cex = 0.2)),
            colhead = list(fg_params=list(cex = 0.2)),
            rowhead = list(fg_params=list(cex = 0.2)))

          ncolP <- ncol(printPage)

          gridTab <- rbind(names(printPage[,1:ceiling(ncolP/3)]),as.character(printPage[,1:ceiling(ncolP/3)]),
                           names(printPage[,(ceiling(ncolP/3)+1):(ceiling(ncolP/3*2))]), as.character(printPage[,(ceiling(ncolP/3)+1):(ceiling(ncolP/3*2))]),
                           names(printPage[,(ceiling(ncolP/3*2)+1):ncolP]),as.character(printPage[,(ceiling(ncolP/3*2)+1):ncolP]))

          gridTab2 <-list(tableGrob(gridTab,theme = mytheme))
          print(do.call(grid.arrange, c(c(gridTab2,p2),
                                        layout_matrix=list(lay1)
          )))


        }else{
          library(gridExtra)
          print("No Ctrl found")
          print(do.call(grid.arrange, c(p, ncol=2, layout_matrix=list(cbind(c(1,3), c(2,4))))))
        }
      }


    }
    print("tra chrom pos path")
    if(length(PotInPotvsNoInf4$genesT) == 1){
      ChromPosT <- rbind(as.character(PotInPotvsNoInf4$posT), as.character(PotInPotvsNoInf4$genesT))
      strand3 <- "+"

    }else{
      ChromPosT <- sapply(1:length(PotInPotvsNoInf4$genesT), function(x){
        print(x)
        if(x > length(PotInPotvsNoInf4$genesT)){
          stop()
        }
      }
      )
    }
    if(any(sapply(ChromPosT,is.null))){
      null <-  which(sapply(ChromPosT,is.null))
      printP(null," is NULL")
      return()
    }
    ChromPosT <- t(unique(t(ChromPosT)))
    if(ncol(ChromPosT) == 1){
      ChromPosT <- t(ChromPosT)
    }
    if(is.list(ChromPosT)){
      ChromPosT <- do.call(rbind,ChromPosT)
      chromPos1 <- paste0(ChromPosT[,2])
    }else{
      chromPos1 <- paste0(ChromPosT[,2])
    }

    strand3 <- "+"

    if(individual != T){
      p <- sapply(1:nrow(ChromPosT), function(x){
        printP(x, "plot, cleavage pos: ", ChromPosT[x,1])
        readDepthBam(posT = ChromPosT[1,x],
                     aliFilesPath = aliFilesPath,
                     aliFilesPattern = aliFilesPattern1,
                     aliFilesPatternCtrl = aliFilesPatternCtrl1,
                     chr1 = chromPos1[x],
                     edgesExtend = edgesExtend1,
                     printPage = PotInPotvsNoInf4[x,],
                     individual = individual,
                     ylim1 = ylim1,
                     strand3 = strand3
        )
      }
      )
    }else if(individual == T){
      p <- sapply(1:nrow(ChromPosT), function(x){
        printP(x, "plot, cleavage pos: ", ChromPosT[x,1])
        readDepthBam(posT = ChromPosT[1,x],
                     aliFilesPath = aliFilesPath,
                     aliFilesPattern = aliFilesPattern1,
                     aliFilesPatternCtrl = "",
                     chr1 = chromPos1[x],
                     edgesExtend = edgesExtend1,
                     printPage = PotInPotvsNoInf4[x,],
                     individual = T,
                     ylim1 = ylim1,
                     strand3 = strand3
        )
      }
      )
    }

  }
  if(!is.null(cleavageData$NFA)){
    cleavageData = cleavageData[order(cleavageData$NFA, decreasing = T),]
  }
  run1 <- unique(cleavageData[,c("genesT","posT")])
  if(addToDir == F){
    #Careful!!, if you want to delete old pictures in your output folder,
    #please remove the hash of the unlink command and run the script again
    stop("remove hash under to unlink")
    #unlink(dirO, recursive = TRUE)
  }
  dir.create(path = dirO, recursive = T)
  existingFiles <- list.files(dirO)
  if(length(existingFiles) != 0){

    substrToChrFromEnd <- function(text, pattern1,nrPatFromEnd = 0){
      #substr from the end of the string to the character "pattern1" assigned
      #nrPatFromEnd defines the Nth occurrence of pattern1 from the end

      staText <- gregexpr(text = text,pattern = pattern1, fixed = T)
      nrPatLast <- sapply(staText, length)
      staText <-sapply(seq_along(staText), function(x){
        sapply(staText[x], "[[", (nrPatLast[x]-nrPatFromEnd))
      }
      )
      return(substr(text,1,staText- 1))
    }

    existingFiles2 <- substrToChrFromEnd(text = existingFiles,pattern1 = "_", nrPatFromEnd = 0)
    run2 <- facToCharDf(run1)
    run3 <- paste(run2[,1], run2[,2],sep = "_")

    if(onlyOneTest == T){
      onlyOne <- which(table(existingFiles2) == 1)
      onlyOne2 <- existingFiles2[onlyOne]
      run1 <- run1[which(!is.na(match(run3,onlyOne2))),]
    }else if(length(which(is.na(match(run3,existingFiles2)))) == 0){
      stop("no files to add")
    }else{
      run1 <- run1[which(is.na(match(run3,existingFiles2))),]
    }
  }
  for(x in(1:nrow(run1))){
    printP(x, " of ", nrow(run1))
    if(savePics == T){
      jpeg(file = paste0(dirO,"/",run1$genesT[x],"_",run1$posT[x],"_R1.jpg"),
           width = jpegWidHei[1], height = jpegWidHei[2], quality = qual, pointsize = pz)
    }

    transAndPosToBam(PotInPotvsNoInf4 = run1[x,],
                     edgesExtend1 = edgesExtend1,
                     aliFilesPattern1 = aliFilesPattern1,
                     aliFilesPath = aliFilesPath,
                     aliFilesPatternCtrl1 = "",
                     individual = T,
                     ylim1 = ylim1,
                     gffTrans = gffTrans
    )
    if(savePics == T){
      dev.off()
      jpeg(file = paste0(dirO,"/",run1$genesT[x],"_",run1$posT[x],"_R2.jpg"),
           width = jpegWidHei[1], height = jpegWidHei[2], quality = qual, pointsize = pz)
    }
    transAndPosToBam(PotInPotvsNoInf4 = run1[x,],
                     edgesExtend1 = edgesExtend1,
                     aliFilesPattern1 = aliFilesPattern2,
                     aliFilesPath = aliFilesPath,
                     aliFilesPatternCtrl1 = "",
                     individual = T,
                     ylim1 = ylim1,
                     gffTrans = gffTrans
    )
    if(savePics == T){
      dev.off()
    }
  }
}

#' Constructs a gff object
#'
#'Reqires spaln installed on your system, a genome and transcripts
#' @param genome Path to genome
#' @param transcript Path to transcripts
#' @param outDir Path to dir where the gff will be written
#' @keywords makeGFF
#' @export
#' @examples
#' makeGFF(genome,transcript, outDir)
makeGFF <- function(genome,transcript, outDir){
  system(paste0("spaln ", genome, " ", transcript, " -o S ", outDir))
}

#' Read a gff object int R
#'
#' @param path Path to the gff file
#' @keywords readGFF
#' @export
#' @examples
#' readGFF(path)
readGFF <- function(path){
  miRNAgff <- read.delim(path, header=F, comment.char="#")
  return(miRNAgff)
}

#' An example how to adopt the gff file to be used with cleavageWindows
#' @param gffTrans Path to the gff file
#' @keywords extendGffTrans
#' @export
#' @examples
#' extendGffTrans(gffTrans)

extendGffTrans <- function(gffTrans = gffTrans){
  findNthChar <- function(pattern,string,n){
    stringBack <- unlist(gregexpr(pattern, string))[n]
  }
  gffTrans2 <- substr(gffTrans,9,10)
  j1 <- sapply(gffTrans, pattern = "_",findNthChar, 1)
  j2 <- sapply(gffTrans, pattern = "_",findNthChar, 2)
  j3 <- sapply(gffTrans, pattern = "_",findNthChar, 3)
  j4 <- sapply(gffTrans, pattern = "_",findNthChar, 4)
  j2[which(is.na(j2))] <- nchar(gffTrans)[is.na(j2)]
  gffTrans3 <- as.data.frame(cbind(V1 = gffTrans2,
                                   V2 = "ownDesign",
                                   V3 = "mRNA",
                                   V4 = substr(gffTrans,j2+1,j3-1),
                                   V5 = substr(gffTrans,j3+1,j4-1),
                                   V6 = ".",
                                   V7 = ".",
                                   V8 = ".",
                                   V9 = as.character(gffTrans),
                                   V10 = substr(gffTrans,j1+1,j2-1)
  )
  )
  colnames(gffTrans3)[1] <- "chrom"
  colnames(gffTrans3)[4] <- "pos"

  colnames(gffTrans3)[4] <- "posStart"
  colnames(gffTrans3)[5] <- "pos"
  return(gffTrans3)
}

