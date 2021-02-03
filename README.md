# smartPARE

#This is an R package for sRNA cleavage confirmation of degradome data, that can be applied on sRNA cleavage prediction data.  
#I used PAREsnip2 to predict the sRNA cleavages, which is to date the most efficient miRNA cleavage analysis tool according to their paper  
#(Thody et al. 2018). smartPARE was then run to evaluate the PAREsnip2 predicted cleavages. Our data was rich on noise, possibly caused by degraded mRNA, causing a lot of false positives in the PAREsnip2 generated data when including all types of sRNA.  smartPARE helped in identifying the false positives in the PAREsnip2 data. 

#The analysis is separated into 3 parts. "Preparation of cleavage windows", "Cleavage window training" and Cleavage confirmations.  

#If you decide to use your own images, the first step "Preparation of cleavage windows"  
#might not be necessary (but hopefully helpful).  
#If you go with our CNN model, the second step "Cleavage window training" can be skipped.  

#We hope to help you to make analysis of your degradome more efficient. Please site our article if you are using the package (still unpublished).  

#Dependencies - The following packages are required: 
#reticulate, EBImage, fftwtools, keras, tensorflow, kerasR, mcparallelDo, cowplot, reshape2, ggplot2, rBayesianOptimization, zoo, GenomicRanges, GenomicAlignments, Rsamtools,  gridExtra, data.table, reshape2, generics, IRanges, BiocGenerics, magrittr 

#Installing keras can be a challange depending on your system. My system was complaining about missing the hdf5=1.10.5 version why the following steps worked for me but I recommend  first trying a standard installation of keras.   
library(reticulate)  
reticulate::conda_install(c("hdf5=1.10.5"), pip = TRUE)  
reticulate::use_condaenv("hdf5=1.10.5", required = TRUE)  
install.packages("tensorflow")  
install.packages("keras")  
install.packages("kerasR")  
library("keras")  
install_keras()  

#Set up smartPare:  
#devtools::install_github('kristianHoden/smartPARE')  
#library(smartPare)  

#If you have any questions or encounter any code related problems, don't hesitate to ask or inform.  

# Preparation of cleavage windows
#This is just an example for how the windows aka cleavage pictures might be created.  
#feel free to use your own pictures. However if using the CNN model we designed it is recommended  
#To adopt this script to your data as other pictures might not be recognized by the model.  

cleavageWindows(dirO = paste0("pathOut/"),  
                cleavageData = cleavageDataDataset,  
                aliFilesPath = "path/bamTranscriptome/",  
                aliFilesPattern1 = "pattern1.sorted.bam$",  
                aliFilesPattern2 = "pattern2.sorted.bam$",  
                ylim1 = 5,  
                edgesExtend1 = c(1,21),  
                gffTrans = gffTrans,
                savePics = T,
                jpegWidHei = c(480,480),
                qual = 75,
                pz = 12
)  

# Cleavage window training - If you deside to create your own model, otherwise skip this part  

#1.  
#Manually put pictures of true cleavages in subdirs to  
#homePath/train/goodUp (true cleavages on 5' strand)  
#homePath/train/goodDown (true cleavages on 3' strand)  
#homePath/train/bad (false cleavages)  
#make sure to get as much variation in the pictures as possible  

#2  
#Create the training dataset  
homePath1 = "example/"  
kerasCreateDataset_2d(homePath = homePath1 ,pixels = 28)  

#3  
#Tune the cyclical learning rate  
tuneCLR(batch_size2 = 64,  
        epochs_find_LR = 20,  
        lr_max = 0.1,  
        optimizer2 = keras::optimizer_sgd(lr=lr_max, decay=0), #optimizer_rmsprop(lr=lr_max, decay=0),  
        validation_split2 = 0.2,  
        rollmeanSplit = 3  
)  

#4  
#Double check the assignment of the Learning_rates  
#If you need to change the Learning_rate_l or Learning_rate_h manually  
#The algorithm is a bit shaky for non-smooth curves  
#Learning_rate_l should be at minimum of the curve and Learning_rate_h at max  
#Learning_rate_l = 5e-04 #1e-02  
#Learning_rate_h = 1*10^-3  

rm1 <- 20  
plot(zoo::rollmean(accDforig$lr, rm1),  
     zoo::rollmean(accDforig$acc, rm1),  
     log="x", type="l", pch=16, cex=0.3,  
     xlab="learning rate", ylab="accuracy: rollmean(100)")  
abline(v=Learning_rate_l, col="blue")  
abline(v=Learning_rate_h, col="red")  

#5  
#Define the bounds of the variables for the Bayesian optimization   
search_bound <- list(unitPower2 = c(0,4),  
                     epochs2 = c(100, 300),  
                     batch_size2 = c(32,128),  
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
  rBayesianOptimization::BayesianOptimization(FUN = runCLR,   
                                              bounds = search_bound,   
                                              init_grid_dt = search_grid,   
                                              init_points = 0,  
                                              n_iter = 3,  
                                              acq =  "ucb" #"ei" "ucb"  
  )  

#8  
#Check which model you prefer  
order(bayes_ucb$Pred, decreasing = T)  
which(bayes_ucb$Pred == max(bayes_ucb$Pred))  
1/max(bayes_ucb$Pred)  

# Cleavage confirmations  

#1  
#Load your preferred model  
#The model file names are written so that they start with the number of the model,  
#then accuracy (0-1), then loss etc  

#a if you designed your own model (your models are saved in your homeDir/bayesmodels
#the performance of the model can be interpreted from the model name should bayes_ucb$Pred
#have been overwritten.  
model <- keras::load_model_hdf5("example/bayesmodels/modelNumberAndContinousName.h5")#  
model %>% summary()  

#b if you are using our model  
model <- keras::load_model_hdf5("data/model/CNNmodel.h5")  
model %>% summary()  

#2
#Define the directories you want to examine  
extDirs <- unique(dirname(list.files(homePath1,rec=T)))  
#exclude the training dirs  
extDirs[-which(startsWith(prefix = "train",extDirs))]  
extDirs <- extDirs[-which(startsWith(prefix = "train",extDirs))]  
rootExt <- paste0(homePath1,extDirs)  

#3
#Examine the cleavages of your directories  
examineCleavages(examinePath = rootExt, model = model,pixels = pixs)  

#4 
#Constructs a list of the true cleavages that can be  
#used to filer away false cleavages from the original dataset    
###important that all final dirs end with _mpd ex athA_mpd or exampleDir_mpd  
kerasListTrue(pathToTrue = homePath1)  
