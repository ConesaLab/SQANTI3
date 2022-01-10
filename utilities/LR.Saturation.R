### Rarefraction functions Ana


# For calculating saturation of Pacbio LR data we use to functions
# LR.rarefraction takes a eSet object where several samples area possible and performs the simulation of counts
# for different sequencing depths. The parameters are:
  # input --> eSet. Make sure you have factors and possibly biotype and category defined too
  # ticks --> how many values of sequencing depth  (SD) will be sampled. If NULL this is ten, you can provide specific values/
  # rep --> how many times to repleat the sample at each SD. The times for computation is proportional
  # to the number of samples, ticks and reps. If there are many samples you can use less reps. Better not to use less ticks.
  # samples. You indicate the samples you want by number. The funcion ALWAYS makes a simulation for a global sample.
  # this is a sample that sums the gene counts across samples.

# plot.rarefraction takes the result of LR.rarefraction and displays the saturation graphs. Parameters are
  # rarefract --> result of LR.rarefraction
  # k = 1 --> filter the number of counts, i.e. if k = 3, the saturation shows the number of transcripts that have at least 3 counts and each SD.
  # sample --> which sample to plot. If you put multiple samples, and have break.by != "none", then specify must be only one type of transcript, i.e. "FSM" or "coding"
  # deph.increase --> expansion of SD to project the saturation curve. If = 2 means that the curve goes up to 2x the real SD.
  # break.by --> if you want results by transcript qualifier. Possible values "none", "category", "biotype"
  # specify = "all transcripts"  --> Will give you which transcripts types for break.by. For example if break.by = "category", then specify could be c("FSM", "ISM"). See restrictions when samples > 1
  # what --> type of graph to plot. Could be either saturation ("sat") or transcript increase/ million additional reads ("diff"). Or both.

# you additionally need the package NOIseq (not really sure about this) and the newdata function for "NOISEq"



##### Functions needed
#######################


library("NOIseq")
source("./LR.rarefraction.R")
source ("./plot.rarefaction.R")
source("./NewReadData.R")

##########################################
##########  Mouse data 2016 ##########
##########################################
# Data
##############

# counts
mycounts <- read.delim("~/Florida/Work/ProjectsUF/ProjectIdeas/SaturatioPlot/mycounts.txt", as.is = TRUE, sep = " ")

# biotype
mybiotype = sample(c(rep("coding" , round(nrow(mycounts)*0.8)), 
                     rep("non-coding", round(nrow(mycounts)*0.2))),
                   nrow(mycounts))
mybiotype = data.frame(biotype = mybiotype) ; rownames(mybiotype)  = rownames(mycounts)
head(mybiotype); length(mybiotype[,1]) == nrow(mycounts)

## category
mycategory = sample(c(rep("FSM", round(nrow(mycounts)*0.6)+1), 
                      rep("ISM", round(nrow(mycounts)*0.1)),
                      rep("NIC", round(nrow(mycounts)*0.25)),
                      rep("NNC", round(nrow(mycounts)*0.05))),
                    nrow(mycounts))
mycategory = data.frame(category= mycategory) ; rownames(mycategory)  = rownames(mycounts)
head(mycategory); length(mycategory[,1]) == nrow(mycounts)

# factor
myfactor <- read.delim("myfactor.txt", as.is = TRUE, sep = " ")

# create data
library(NOISeq)
source ("NewReadDataFunctionForNOISeq.R")
mydata = readData(data = mycounts, factors = myfactor, biotype = mybiotype, category = mycategory)


# Analysis
###########
rarefract <- LR.rarefraction (mydata, samples = c(1:15), rep = 3)
layout(matrix(c(1:4),2,2))

plot.rarefaction(rarefract, sample = 1, k = 1, depth.increase = 2)
plot.rarefaction(rarefract, sample = 16, k = 1, depth.increase = 2)
plot.rarefaction(rarefract, sample = c(1:16), k = 1, depth.increase = 2)
plot.rarefaction(rarefract, sample = c(1:16), k = 5, depth.increase = 2)

layout(matrix(c(1:4),2,2))
plot.rarefaction(rarefract, sample = c(16), k = 1, depth.increase = 2, break.by ="none", specify = "FSM")
plot.rarefaction(rarefract, sample = c(1,16), k = 1, depth.increase = 2, break.by ="category", specify = "ISM")
plot.rarefaction(rarefract, sample = c(1,16), k = 1, depth.increase = 2, break.by ="category", specify = "NIC")
plot.rarefaction(rarefract, sample = c(1,16), k = 1, depth.increase = 2, break.by ="category", specify = "NNC")







##########################################
##########  Alzheimer data 2016 ##########
##########################################
# Data
######
Alzheimer2016 <- read.delim("/Users/Dr.Conesa/Florida/Work/ProjectsUF/Pending/R21_SQANTI-QC/RAnalyses/IsoSeq_Alzheimer_2016edition_polished.confident.unimapped.abundance.txt",
                            as.is = T, skip  = 14) ; head(Alzheimer2016)

Alzheimer2016.counts <- cbind (Alzheimer2016[,2], Alzheimer2016[,2]+Alzheimer2016[,3])
rownames(Alzheimer2016.counts) <- Alzheimer2016 [,1]
colnames(Alzheimer2016.counts) <- c("Alzheimer2016_FL", "Alzheimer2016_FL+nFL")
head(Alzheimer2016.counts)
myfactors <- data.frame(sample = c("Alz", "Alz"))
rownames(myfactors) = colnames(Alzheimer2016.counts)
mydata.Alz = readData(data = Alzheimer2016.counts, factors = myfactors)

# Analysis
# Analysis
##########
rarefract.Alz <- LR.rarefraction (mydata.Alz , samples = c(1:2))
layout(matrix(c(1:4),2,2))
plot.rarefaction(rarefract.Alz, sample = 1, k = 1, depth.increase = 2)
plot.rarefaction(rarefract.Alz, sample = 1, k = 10, depth.increase = 2)
plot.rarefaction(rarefract.Alz, sample = c(1,2), k = 10, depth.increase = 2)
plot.rarefaction(rarefract.Alz, sample = 2, k = 5, depth.increase = 2)




##########################################
##########  Alzheimer data 2019 ##########
##########################################
# Data
######
library("NOISeq")
source("./LR.rarefaction.R")
source ("./plot.rarefaction.R")
source("./NewReadData.R")

Alzheimer2019 <- read.delim("/Users/Dr.Conesa/Florida/Work/ProjectsUF/Pending/R21_SQANTI-QC/RAnalyses/out.abundance_Alzheimer.txt",
                            as.is = T, skip  = 7) ; head(Alzheimer2019)


#### el *.counts es una matriz de una unica columna que es la de FL (sumando las FL de distintas muestras si las hubiera)
#### el myfactors puede ser coding/non-coding
#### el category puede ser las strcutural categories

Alzheimer2019.counts <- cbind (Alzheimer2019[,2])
rownames(Alzheimer2019.counts) <- Alzheimer2019[,1]
colnames(Alzheimer2019.counts) <- "Alzheimer2019_FL"
head(Alzheimer2019.counts)
myfactors <- data.frame(sample = c("Alz"))
rownames(myfactors) = colnames(Alzheimer2019.counts)
#mybiotype=as.matrix(data.class$coding)
#rownames(mybiotype)=data.class$isoform
#mycategory=as.matrix(data.class$structural_category)
#rownames(mycategory)=data.class$isoform
#mydata.Alz19 = readData(data = Alzheimer2019.counts, factors = myfactors, biotype = mybiotype, category=mycategory)

# Analysis
##########
rarefract.Alz19 <- LR.rarefraction(mydata.Alz19 , samples = 1)
layout(matrix(c(1:4),2,2))
plot.rarefraction(rarefract.Alz19, sample = 1, k = 1, depth.increase = 2)
plot.rarefraction(rarefract.Alz19, sample = 1, k = 2, depth.increase = 2)
plot.rarefraction(rarefract.Alz19, sample = 1, k = 3, depth.increase = 2)
plot.rarefraction(rarefract.Alz19, sample = 1, k = 5, depth.increase = 2)
Alheimer_2016 <- read.delim()
