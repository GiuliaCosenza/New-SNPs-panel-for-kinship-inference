############################################################################################
###----------creazione di pedigree a partire da dati genotipici con pedbuildr------------###
############################################################################################


################################################################
##prima di iniziare assicurarsi di avere pronti due file csv:  #
##(1)Un file contenete le info genetiche degli individui,      #
##   secondo il seguente schema (header= TRUE)                 #   
##                                                             #
##   id fid mid sex  <rs1> <rs2> <rs3> <rs4> <rs5>             #
#     1   0   0   1   1/2   1/2   1/1   1/2   1/2              #
#     2   0   0   1   1/1   1/2   1/2   2/2   1/2              #
#     3   0   0   1   1/1   1/2   1/2   2/2   1/1              #
##                                                             #
##(2)Un file con le info (rs e frequenza allelica)             #
##   relative ai marcatori,                                    #
##   secondo il seguente schema (header= FALSE)                # 
##                                                             #
##    rs1  NA                                                  #
##     1   0.5                                                 #
##     2   0.5                                                 #
##    rs2  NA                                                  #
##     1   0.75      Mettere prima tutti                       #
##     2   0.25     i marcatori diallelici                     #
##    rs3  NA       seguiti dai triallelici                    #
##     1   0.33     ecc...                                     #
##     2   0.33                                                #
##     3   0.33                                                #
##                                                             #
##                                                             #
################################################################


#------------------------------------------------#
#|    caricamento dei pacchetti essenziali      |#
#------------------------------------------------#

#install.packages("pedbuildr")
#install.packages("devtools")
#install.packages("shiny")
#install.packages("stringi")
#install.packages("Rtools")
#devtools::install_github("magnusdv/pedbuildr")
#devtools::install_github("magnusdv/forrel")
library('pedtools')
library("pedbuildr")
library('usethis')
library("devtools")
library('forrel')


#------------------------------------------------#
#|  costruzione del ped file con le info        |#
#|  genotipiche degli individui della coorte    |#
#------------------------------------------------#


##carico il file contenente le informazioni genotipiche dei 18000 individui (dataframe 18000*10381)
file<-read.csv("myped-out-from-ped2-10000repl.csv", sep = ',', header = FALSE)
file<-file[,-c(2,4,6,8,10,11)]

##costruisco un nuovo dataframe a partire dal primo con le colonne rinominate correttamente 
##ed eliminando le informazioni relative ai genitori (mid e fid sempre =0)
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

#--------Pulizia ambiente-------------#
rm(list = ls()[!(ls() %in% "ped")])

#--------------------------------------------------------------#
#|        Rimuovo i marcatori in linkage Disequilibrium       |#
#--------------------------------------------------------------#
removed<-read.csv('rem.csv', header = TRUE, sep = ',')
ped<-ped[,-c(removed[,2])]

#------------------------------------------------#
#|  costruzione dell'oggetto con le info        |#
#|  sui marcatori degli individui della coorte  |#
#------------------------------------------------#

##ATTENZIONE: non lanciare i comandi automaticamente, adattali ai tuoi file##

#divido il df in due, uno con i loci bialleleici e uno con quelli triallelici
markers<-read.csv("inputfamilias_DNA_db_num_no_Linkage_Disequilibrium.csv", header = FALSE, sep = ';')
markers<-markers[,1:2]
marker2<-markers[1:14472,]
marker3<-markers[14473:14576,]

loc=as.list(NULL)

#creo la lista con tutti i loci per cui ogni elemento della lista rappresenta un locus
#ogni elemento della lista è caratterizzato da 3 variabili:
#$name (il codice identificativo del marcatore), $alleles (il nome delle alternative alleliche), $afreq (la frequenza di ciascun allele)
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


#--------Pulizia ambiente-------------#
rm(list = ls()[!(ls() %in% c("ped", "loc"))])

#############################################################

  #------------------------------------------------#  
  #|  costruzione di una lista di 1000 elementi   |#
  #|  ciascuno dei quali rappresenta una famiglia.|#
  #|  Il df iniziale viene quindi frammentato in  |#
  #|  una singola lista a più elementi            |#
  #------------------------------------------------#

# L'oggetto pedlist è una lista di elementi ped che rappresentano ciascuno una famiglia di 18 individui.  
# I data.frame contenenti le informazioni genotipiche degli individui che compongono ciascuna famiglia devono 
# essere obbligatoriamente convertiti in elementi ped affinchè pssano essere l'argomento della funzione ibdEstimate 
# dello step successivo.                                                                                            

##Vedi file Ped.R per la spiegazione del comando as.ped  

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
  

#--------Pulizia ambiente-------------#
rm(list = ls()[!(ls() %in% "pedlist")])

#---------------------------------------------------------------------------------------------------------------#
#|                                    IBD and LR estimation and triangle plot                                  |#  
#|                                                                                                             |#
#| Considerando tutte le possibili coppie di individui all'interno di una famiglia composta da 18 individui    |#
#| viene costruito un ciclo che calcoli i valori k0,k1,k2 della relazione tra ciascuna coppia                  |#
#| e il valore di LR associato                                                                                 |#
#---------------------------------------------------------------------------------------------------------------#
checkpairwise_result <- data.frame()
for (i in 1:length(pedlist)){
  checkpairwise_result <- rbind(checkpairwise_result, checkPairwise(pedlist[[i]]))
}

#---------------------------------------------------------------------------------------------------------------#
#|                                                         RMSE                                                |#  
#|                                                                                                             |#
#|                                                                                                             |#
#---------------------------------------------------------------------------------------------------------------#
kinshipIBD <- data.frame("kinship" = c("UN", "PC", "FS", "HS,ZI,GP", "CO,PP,GG", "1C1R,PPP","C2"), "k0" = c(1,0,0.25,0.5,0.75,0.875,0.9375), "k1" = c(0,1,0.5,0.5,0.25,0.125,0.0625), "k2" = c(0,0,0.25,0,0,0,0))

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

#--------Pulizia ambiente-------------#
rm(list = ls()[!(ls() %in% c("pedlist", "kinshipIBD", "checkpairwise_result"))])
