# smartPARE

This is an R package for sRNA cleavage confirmation on degradome data, that can be applied on sRNA cleavage predictions.
# I used PAREsnip2 which is to date the most efficient miRNA cleavage analysis tool according to their paper 
# (Thody et al. 2018). This extension module was developed also to identify other sRNA cleavages because 
# our data was rich on noise, possibly caused by degraded mRNA.

# The analysis is separated in 3 parts. "Preparation of cleavage windows", "Cleavage picture training" and
# Cleavage confirmations. 

# If you decide to use your own pictures, the first step "Preparation of cleavage pictures"
# might not be necessary (but hopefully helpful).
# If you go with our CNN model, the second step "Cleavage window training" can be skipped.

# We hope to help you to make analysis of your degradome more efficient. Please site our article if you are using the package.

#Set up smartPare:
#install_github('smartPare','github_username')
#library(smartPare)

# If you have any questions or encounter any code related problems, don't hesitate to ask or inform.



# Dependencies - The following packages must be installed. Please see our publication to see what 
# versions I used, if struggling. 
library("reticulate")
reticulate::use_condaenv("hdf5=1.10.5", required = TRUE) # This might not be necessary for you but it helped me. 
library("EBImage")
library("fftwtools")
library("keras")
library("tensorflow")
library("kerasR")
library("mcparallelDo")
library("cowplot")
library("reshape2")
library("ggplot2")
library("rBayesianOptimization")
library("zoo")


######### Preparation of cleavage windows
# This is just an example for how the windows aka cleavage pictures might be created.
# feel free to use your own pictures. However if using the CNN model we designed it is recommended
# To adopt this script to your data as other pictures might not be recognized by the model. 

###
# cleavageData - must be a dataframe containing at least columns genesT (the target genes) and posT (target position in the transcript)

# edgesExtend1 - how many positions from the cleavage site that are included in the pictures. Default is c(1,2),
# defining 1 pos upstream and 21 positions downstream. As the degradome reads are 20 nt each this defines 
# 1 position upstream and downstream of the degradome reads. Only one, to try and exclude as much background 
# as possible but still catch the caracteristic "cleavage tower".

# ylim1 - defines the minimum plotted height of the y-axis. In practise this means a lower "cleavage tower"
# will not be recognized by the CNN. This to exclude potential false positives caused by noice

# onlyOneTest - sometimes the script is missing one cleavage picture for unknown reasons. Hence, 
# not reporting that all files are created but still not creating any more. If this is the case 
# one can test for this by running with "onlyOneTest = T", default is F.

# addToDir - if one wants to add the pictures to the directory. Default is T. If F the user also need to 
# remove the has in the function script, this as a precaution not to remove any pictures by accident  

# aliFilesPath -path to the bamfiles  for the user defined degradome. The easiest is to use transcriptome   
# generated bam files. If you do so the path must end with "bamTranscriptome/". However, this option cannot check  
# for cleavages outside the transcripts (which is totally  fine if you have a well annotated genome).
# If you have trascripts and a genome but no gff file you can create a gff with the function 
# makeGFF (reqires spaln installed on your system).
# Also in makeWindowPictureFunc you can find the function readGFF which is a basic function to import your GFF. 
# IMPORTANT!!! is that your gffTrans has a 10th column with your the transcript ID. Often this must be 
# extracted from the 9th column (V9). However, most GFF files I've been exposed to are different. 
# In the attached function extendGffTrans you can see how this was achieved from the Jupe dataset (Jupe et al. 2013) 
# I was using. 

# gffTrans - a dataframe of your imported gff with V10 = transcript ID see  for example extendGffTrans 
# in makeWindowPictureFunc for example

# aliFilesPattern1 - regular expression defining the pattern of your first degradome library

# aliFilesPattern2 - regular expression defining the pattern of your second degradome library

# dirO - output path 


cleavageWindows(dirO = paste0("pathOut/"), 
                 cleavageData = cleavageDataDataset,
                 aliFilesPath = "path/bamTranscriptome/",
                 aliFilesPattern1 = "pattern1.sorted.bam$",
                 aliFilesPattern2 = "pattern2.sorted.bam$",
                 ylim1 = 5,
                 edgesExtend1 = c(1,21),
                 gffTrans = gffTrans
)






##### Cleavage picture training - If you deside to create your own model, otherwise skip this part

#1. Manually put pictures of true cleavages in subdirs to 
#homePath/train/goodUp (true cleavages on 5' strand)
#homePath/train/goodDown (true cleavages on 3' strand)
#homePath/train/bad (false cleavages)
#make sure to get as much variation in teh pictures as possible

#2 
# Create the training dataset  
homePath1 = "keras/3categories/1_21_5_2/"
kerasCreateDataset_2d(homePath = homePath1 ,pixels = 28)

#3
# Tune the cyclical learning rate
tuneCLR(batch_size2 = 64,
        epochs_find_LR = 20,
        lr_max = 0.1, 
        optimizer2 = optimizer_sgd(lr=lr_max, decay=0), #optimizer_rmsprop(lr=lr_max, decay=0),
        validation_split2 = 0.2,
        rollmeanSplit = 3
)

#4 
#double check the assignment of the Learning_rates
#If you need to change the Learning_rate_l or Learning_rate_h manually
#The algorithm is a bit shaky for non-smooth curves
#Learning_rate_l should be at minimum of the curve and Learning_rate_h at max
#Learning_rate_l = 5e-04 #1e-02
#Learning_rate_h = 0.03
#Learning_rate_h = 1*10^-3

rm1 <- 20
plot(rollmean(accDforig$lr, rm1),
     rollmean(accDforig$acc, rm1),
     log="x", type="l", pch=16, cex=0.3,
     xlab="learning rate", ylab="accuracy: rollmean(100)")
abline(v=Learning_rate_l, col="blue")
abline(v=Learning_rate_h, col="red")

#5
#Define the bounds of the variables for the Bayesian optimization 
search_bound <- list(unitPower2 = c(0,4), 
                     epochs2 = c(100, 300),
                     batch_size2 = c(32,128) ,
                     dropout2 = c(0, 0.3),
                     validation_split2 = c(0.1,0.4),
                     NOfilters2 = c(1,4),
                     NO_pooling2 = c(1,2)
)
#6
#Define initiation grid, ie. the start values of the variables 
search_grid <- data.frame(unitPower2 = c(1,2), 
                          epochs2 = c(5,5),
                          batch_size2 = c(32,64),
                          dropout2 = c(0.05,0.1),
                          validation_split2 = c(0.1,0.2),
                          NOfilters2 = c(1,2),
                          NO_pooling2 = c(0,1)
)

#7
#Initiate the model counting, it is run in the global environment
#and run the Bayesian optimization
#Makes sure to set the other standards of the input variables of your 
#choice in the function defitition of runCLR 
#for instance the pathOut - dir of your output
count1 <-1
bayes_ucb <-
  BayesianOptimization(FUN = runCLR, 
                       bounds = search_bound, 
                       init_grid_dt = search_grid, 
                       init_points = 0,
                       n_iter = 100,
                       acq =  "ucb" #"ei" "ucb"
  )

#8
#Check which model you prefer
order(bayes_ucb$Pred, decreasing = T)
which(bayes_ucb$Pred == max(bayes_ucb$Pred))
1/max(bayes_ucb$Pred)





######### Cleavage confirmations

#1
#Load your preferred model
#The model file names are written so that they start with the number of the model,
#then accuracy (0-1), then loss etc

#a if you designed your own model
model <- load_model_hdf5("keras/3categories/bayesmodels/modelNumberAndContinousName.h5")
model %>% summary()
#b if you are using our model
model <- load_model_hdf5("CNNmodel.h5")
model %>% summary()

#2
#Define your directories you want to examine
extDirs <- unique(dirname(list.files(homePath1,rec=T)))
#exclude the training dirs
extDirs[-which(startsWith(prefix = "train",extDirs))]
extDirs <- extDirs[-which(startsWith(prefix = "train",extDirs))]
rootExt <- paste0(homePath1,extDirs)

#3
#Examine the cleavages of your directories
examineCleavages(examinePath = rootExt, model = model,pixels = pixs)

#4 Constructs a list of the true cleavages that can be
# used to filer away false cleavages from the original dataset 
###important that all final dirs end with _mpd ex athA_mpd or gpotInPinf_mpd
kerasListTrue(pathToTrue = homePath1)
