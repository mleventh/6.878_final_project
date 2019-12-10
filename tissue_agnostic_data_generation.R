# Tissue agnostic significant differences between mutants and wild-types dependencies
# @achill: cellular dependencies scores ('Achilles_gene_effect.csv')
# @metadata: CCLE Cell line metadata ('sample_info.csv')
# @mut: CCLE mutational data ('CCLE_mutations.csv')

library(data.table)

arguments <- commandArgs(TRUE)
achill<-fread(arguments[1],sep=',',header=T,stringsAsFactors=F)
metadata<-fread(arguments[2],sep=',',header=T,stringsAsFactors=F)
mut<-fread(arguments[3],sep=',',header=T,stringsAsFactors=F)

# Pre-processing the data
colnames(achill)[1]<-"DepMap_ID"
metadat.sub<-metadata[,c("DepMap_ID","lineage","lineage_subtype","disease_subtype",
                         "sex","age")]

achill.meta<-merge(achill,metadat.sub,by='DepMap_ID')
colnames(achill.meta)<-gsub(" \\([0-9]+\\)","",colnames(achill.meta))
mut.meta<-merge(mut,metadat.sub,by="DepMap_ID")

# Initialization of the training and test data
train.delta <- data.frame(matrix(ncol=2, nrow = 1), stringsAsFactors = F)
colnames(train.delta) <- c("Gene", "Delta")

test.delta <- data.frame(matrix(ncol=2, nrow = 1), stringsAsFactors = F)
colnames(test.delta) <- c("Gene", "Delta")

i=1
j=1

# Genearte the training and test dataset 
for(gene in unique(mut.meta$Hugo_Symbol)){
  
  col.idx<-which(colnames(achill.meta)==gene)
  
  mutants <- achill.meta[achill.meta$DepMap_ID %in% mut.meta$DepMap_ID[mut.meta$Hugo_Symbol %in% gene],
                         col.idx,with=FALSE]
  
  wt <- achill.meta[achill.meta$DepMap_ID %in% mut.meta$DepMap_ID[mut.meta$Hugo_Symbol %in% gene==FALSE],
                    col.idx,with=FALSE]
  
  # ONLY genes with at least 5 mutant and wt dependencies are added to the training dataset, otherwise test dataset
  # Calculate the difference between the mutant and wt dependency scores
  if(nrow(mutants) >= 5 & nrow(wt) >= 5){
    print(paste(gene, ": train", sep = ""))
    
    if (i == 1) {
      train.delta[i,] <- c(gene,median(mutants[[gene]])-median(wt[[gene]]))
      i <- i + 1
    } else{
      train.delta <- rbind(train.delta, c(gene, median(mutants[[gene]])-median(wt[[gene]])))
    }
  } else{
    if(nrow(mutants) != 0 & nrow(wt) != 0){
      if (j == 1) {
        print(paste(gene, ": test, ", sep = ""))
        
        test.delta[j,] <- c(gene,median(mutants[[gene]])-median(wt[[gene]]))
        j <- j + 1
      } else{
        test.delta <- rbind(test.delta, c(gene, median(mutants[[gene]])-median(wt[[gene]])))
      }
    }
  }
}

# Remove genes with 'NA' deltas
train.delta.nona <- train.delta[is.na(train.delta$Delta)==FALSE,]
test.delta.nona <- test.delta[is.na(test.delta$Delta)==FALSE,]

# Export the training and test datasets to a file
write.table(train.delta.nona,"train_genes_deltas.txt", sep='\t', quote=F, row.names=F)
write.table(test.delta.nona,"test_genes_deltas.txt", sep='\t', quote=F, row.names=F)

train.delta.nona$Delta <- as.numeric(train.delta.nona$Delta)

# Compute FDR-corrected normal p-value of change in CERES scores (training dataset)
train.delta.nona[,"normp"] <- pnorm(train.delta.nona$Delta,mean(train.delta.nona$Delta),sd=sqrt(var(train.delta.nona$Delta)))
train.delta.nona[,"q"] <- p.adjust(train.delta.nona$normp,'fdr')
delta.sig <- train.delta.nona[train.delta.nona$q < 0.1,]

# Export the training data FDR-corrected normal p-value of change in CERES scores to a file
write.table(delta.sig,"train_genes_delta_fdr.txt", sep='\t', quote=F, row.names=F)




