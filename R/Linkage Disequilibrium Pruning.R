#######################################################
## Selecting independent variants by evaluating      ##
##             LINKAGE DISEQUILIBRIUM                ##
#######################################################

#--------------------------------------------------------#
#|                   Loading essential                  |#
#|                        packages                      |#
#--------------------------------------------------------#
library(LDlinkR)


#----------------------------------------------#
#|   Load the file containing the variants    |#
#|   to be analyzed and potentially filtered  |#
#----------------------------------------------#

markers <- read.csv("markers_order_chr.csv", header = FALSE, sep = ';')

#------------------------------------------------#
#|   Group the variants by chromosome           |#
#------------------------------------------------#

ordered_markers <- list()

for (i in 1:22) {
  ordered_markers[[i]] <- data.frame()
  for (k in 1:length(markers[,1])) {
    if (markers[k,1] %in% i) {
      ordered_markers[[i]] <- rbind(ordered_markers[[i]], markers[k,])
    }
  }
}

#----------------------------------------------------------------------#
#|                   Calculate Linkage Disequilibrium                  |#
#----------------------------------------------------------------------#

SNP_LD <- data.frame()
for (q in 1:22) {
  SNP_LD <- rbind(SNP_LD, SNPclip(ordered_markers[[q]][,2], 
                                  pop = "ALL",
                                  r2_threshold = "0.1", 
                                  maf_threshold = "0.01", 
                                  token = "token", 
                                  file = FALSE,
                                  genome_build = "grch37"
  ))
}


# Filter out removed variants
removed <- data.frame()
for (p in 1:length(SNP_LD[,1])) {
  if (!(SNP_LD[p,4] %in% 'Variant kept.')) {
    removed <- rbind(removed, SNP_LD[p,])
  }
}
