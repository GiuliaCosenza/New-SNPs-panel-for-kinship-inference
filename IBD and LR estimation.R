```R
############################################################################################
###----------Pedigree creation from genotype data with pedbuildr------------###
############################################################################################

################################################################
##Before starting, make sure you have two csv files ready:    #
##(1)A file containing genetic information of individuals,     #
##   following the following schema (header= TRUE)             #   
##                                                             #
##   id fid mid sex  <rs1> <rs2> <rs3> <rs4> <rs5>             #
#     1   0   0   1   1/2   1/2   1/1   1/2   1/2              #
#     2   0   0   1   1/1   1/2   1/2   2/2   1/2              #
#     3   0   0   1   1/1   1/2   1/2   2/2   1/1              #
##                                                             #
##(2)A file containing information (rs and allele frequency)   #
##   related to the markers,                                   #
##   following the following schema (header= FALSE)            # 
##                                                             #
##    rs1  NA                                                  #
##     1   0.5                                                 #
##     2   0.5                                                 #
##    rs2  NA                                                  #
##     1   0.75      Place all diallelic markers first         #
##     2   0.25     followed by triallelic ones                #
##    rs3  NA                                                  #
##     1   0.33     etc...                                     #
##     2   0.33                                                #
##     3   0.33                                                #
##                                                             #
##                                                             #
################################################################

#------------------------------------------------#
#|        Load essential packages               |#
#------------------------------------------------#
library('pedtools')
library("pedbuildr")
library('usethis')
library("devtools")
library('forrel')


#------------------------------------------------#
#|  Build the ped file with genotype info       |#
#|  of the individuals in the cohort             |#
#------------------------------------------------#

##Load the file containing genetic information of 18000 individuals (dataframe 18000*10381)
file<-read.csv("genotypes.csv", sep = ',', header = FALSE)
file<-file[,-c(2,4,6,8,10,11)]

##Build a new dataframe from the first one with correctly renamed columns 
##and eliminating parent information (mid and fid always =0)
ped<-as.data.frame(file[,1:2])
ped<-cbind(ped, rep(0,length(ped[,1])))
ped<-cbind(ped, rep(0,length(ped[,1])))
ped<-cbind(ped, file[,5])
names(ped)<-c("famid", "id", "fid", "mid", "sex")

genotypes<-file[,6:length(file[1,])]

vect<-c(1:length(genotypes))*2
vect<-vect[1:(length(vect)/2)]
for ( i in vect){
  ped<-cbind(ped, as.data.frame(paste(genotypes[,i], genotypes[,i-1], sep = '/')))
}
names<-read.csv("Header_familias.csv", sep = ';')
names(ped)[6:length(ped[1,])]<-names(names)

#--------Environment cleaning-------------#
rm(list = ls()[!(ls() %in% "ped")])

#--------------------------------------------------------------#
#|        Remove markers in linkage Disequilibrium             |#
#--------------------------------------------------------------#
removed<-read.csv('rem.csv', header = TRUE, sep = ',')
ped<-ped[,-c(removed[,2])]

#------------------------------------------------#
#|  Build object with marker info               |#
#|  of the individuals in the cohort            |#
#------------------------------------------------#

##ATTENTION: do not run the commands automatically, adapt them to your files##

#Split the dataframe into two, one with biallelic loci and one with triallelic ones
markers<-read.csv("merkers_no_Linkage_Disequilibrium.csv", header = FALSE, sep = ';')
markers<-markers[,1:2]
marker2<-markers[1:14472,]
marker3<-markers[14473:14576,]

loc=as.list(NULL)

#Create a list with all the loci where each element of the list represents a locus
#Each element of the list is characterized by 3 variables:
#$name (the identifier of the marker), $alleles (the name of the allelic alternatives), $afreq (the frequency of each allele)
for (i in 0:4823){
  a=(3*i)+1
  b=a+1
  c=b+1
  namel=marker2[a,1]
  allele=c(2,1)
  afq=c(marker2[b,2],marker2[c,2])
  loc[[i+1]]= list(name = namel , alleles = allele , afreq = afq)
  
}

for (i in 0:25){
  a=(4*i)+1
  b=a+1
  c=b+1
  d=c+1
  namel=marker3[a,1]
  allele=c(3,1,2)
  afq=c(marker3[b,2],marker3[c,2], marker3[d,2])
  loc[[4823+i+1]]= list(name = namel , alleles = allele , afreq = afq)
  
}


#--------Environment cleaning-------------#
rm(list = ls()[!(ls() %in% c("ped", "loc"))])

#############################################################

  #------------------------------------------------#  
  #|  Build a list of 1000 elements               |#
  #|  each representing a family.                 |#
  #|  The initial df is then fragmented into      |#
  #|  a single list with multiple elements        |#
  #------------------------------------------------#

# The pedlist object is a list of ped elements, each representing a family of 18 individuals.  
# The data.frames containing the genotype information of the individuals composing each family must 
# be necessarily converted into ped elements so that they can be the argument of the ibdEstimate function 
# of the next step.                                                                                            

##See Ped.R file for explanation of as.ped command  

fam<-data.frame()
pedlist<-list()

for (j in 1:max(ped$famid)){   
  for (i in 1:length(ped[,1])){
    if (ped[i,1] %in% j){
      fam<-rbind(fam, ped[i,])
    }
  }
  names(fam)<-names(ped)
  pedlist[[j]]<-as.ped(fam, locusAttributes = loc)
  fam<-data.frame()
}
  

#--------Environment cleaning-------------#
rm(list = ls()[!(ls() %in% "pedlist")])

#---------------------------------------------------------------------------------------------------------------#
#|                                    IBD and LR estimation and triangle plot                                  |#  
#|                                                                                                             |#
#| Considering all possible pairs of individuals within a family composed of 18 individuals                    |#
#| a loop is built to calculate the values k0, k1, k2 of the relationship between each pair                    |#
#| and the associated LR value                                                                                   |#
#---------------------------------------------------------------------------------------------------------------#
checkpairwise_result <- data.frame()
for (i in 1:length(pedlist)){
  checkpairwise_result <- rbind(checkpairwise_result, checkPairwise(pedlist[[i]]))
}

#---------------------------------------------------------------------------------------------------------------#
#|                                                         RMSE                                                |#  
#|  The Root Mean Square Error (RMSE) is a measure of the discrepancy between estimated and observed values.   |#
#|  In this context, it is used to associate the values of k0, k1, and k2 obtained from checkpairwise          |#
#|  with the corresponding relationship label.                                                                 |#
#---------------------------------------------------------------------------------------------------------------#
kinshipIBD <- data.frame("kinship" = c("UN", "PC", "FS", "HS,ZI,GP", "CO,PP,GG", "1C1R,PPP","C2"), 
                         "k0" = c(1,0,0.25,0.5,0.75,0.875,0.9375), 
                         "k1" = c(0,1,0.5,0.5,0.25,0.125,0.0625), 
                         "k2" = c(0,0,0.25,0,0,0,0))

for (i in 1:nrow(checkpairwise_result)) {
  RMSE = 10
  kin = ""
  for (j in 1:nrow(kinshipIBD)) {
    if ((((checkpairwise_result[i,4]-kinshipIBD[j,2])^2 + (checkpairwise_result[i,5]-kinshipIBD[j,3])^2 + (checkpairwise_result[i,6]-kinshipIBD[j,4])^2)/3)^0.5 < RMSE)
    {
      RMSE = (((checkpairwise_result[i,4]-kinshipIBD[j,2])^2 + (checkpairwise_result[i,5]-kinshipIBD[j,3])^2 + (checkpairwise_result[i,6]-kinshipIBD[j,4])^2)/3)^0.5
      kin = kinshipIBD[j,1]
    }
  }
  checkpairwise_result$kinship[i] <- kin
}

#--------Environment cleaning-------------#
rm(list = ls()[!(ls() %in% c("pedlist", "kinshipIBD", "checkpairwise_result"))])
```
