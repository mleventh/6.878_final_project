library(data.table)

achill<-fread('Achilles_gene_effect.csv',sep=',',header=T,stringsAsFactors=F)
metadata<-fread('sample_info.csv',sep=',',header=T,stringsAsFactors=F)
mut<-fread('CCLE_mutations.csv',sep=',',header=T,stringsAsFactors=F)

colnames(achill)[1]<-"DepMap_ID"

metadat.sub<-metadata[,c("DepMap_ID","lineage","lineage_subtype","disease_subtype",
                         "sex","age")]

achill.meta<-merge(achill,metadat.sub,by='DepMap_ID')

colnames(achill.meta)<-gsub(" \\([0-9]+\\)","",colnames(achill.meta))

mut.meta<-merge(mut,metadat.sub,by="DepMap_ID")
mut.meta.sub<-mut.meta[mut.meta$lineage%in%achill.meta$lineage&
                         mut.meta$Variant_annotation%in%c("damaging","other non-conserving"),]

i<-0
j<-0

#train<-data.table()
#test<-data.table()

for(tissue in unique(mut.meta.sub$lineage)){
  for(gene in unique(mut.meta.sub$Hugo_Symbol)){
    
    col.idx<-which(colnames(achill.meta)==gene)
    
    mutants<-achill.meta[achill.meta$lineage==tissue&
            achill.meta$DepMap_ID%in%mut.meta.sub$DepMap_ID[
              mut.meta.sub$Hugo_Symbol%in%gene&mut.meta.sub$lineage==tissue]
            ,col.idx,with=FALSE]
    wt<-achill.meta[achill.meta$lineage==tissue&
                      achill.meta$DepMap_ID%in%mut.meta.sub$DepMap_ID[
                        mut.meta.sub$Hugo_Symbol%in%gene==FALSE&mut.meta.sub$lineage==tissue],
                    col.idx,with=FALSE]
    
    if(nrow(mutants)>=5&nrow(wt)>=5){
      p.value<-t.test(mutants[,gene,with=FALSE],wt[,gene,with=FALSE],alternative="less")$p.value
      #q<-p.adjust(p.value,n=nrow(mutants)+nrow(wt),method='fdr')
      print(c(tissue,gene,p.value))
      
      if(i==0){
        train<-data.frame(tissue=tissue,Hugo_Symbol=gene,p=p.value,stringsAsFactors=F)
        #print(train)
        i=i+1
      }else{
        #print(data.frame(tissue=tissue,Hugo_Symbol=gene,p=p.value,stringsAsFactors=F))
        train<-rbind(train,data.frame(tissue=tissue,Hugo_Symbol=gene,p=p.value,stringsAsFactors=F))
      }
    }else{
      if(j==0){
        test<-data.frame(tissue=tissue,Hugo_Symbol=gene,stringsAsFactors=F)
        j=j+1
      }else{
        test<-rbind(test,data.frame(tissue=tissue,Hugo_Symbol=gene,stringsAsFactors=F))
      }
    }
  }
}

write.table(train,"train_genes.txt",sep='\t',quote=F,row.names=F)
write.table(test,"test_genes.txt",sep='\t',quote=F,row.names=F)

train.sig<-train[train$p<0.01,]
train$q<-p.adjust(train$p,method='fdr')
train.sig.q<-train[train$q<0.1,]

mut.meta.sub$gene_tissue_anno<-paste(mut.meta.sub$Hugo_Symbol,mut.meta.sub$lineage,sep=";")
train$gene_tissue_anno<-paste(train$Hugo_Symbol,train$tissue,sep=";")

mut.meta.train<-mut.meta.sub[mut.meta.sub$gene_tissue_anno%in%train$gene_tissue_anno,]
mut.meta.train$sig_dependency<-0
mut.meta.train$sig_dependency[
  mut.meta.train$gene_tissue_anno%in%train$gene_tissue_anno[train$q<0.1]]<-1


#add deltas as a feature
delta<-list()
i=1
for(tissue in unique(mut.meta.train$lineage)){
  for(gene in unique(mut.meta.train$Hugo_Symbol)){
    
    col.idx<-which(colnames(achill.meta)==gene)
    
    mutants<-achill.meta[achill.meta$lineage==tissue&
                           achill.meta$DepMap_ID%in%mut.meta.train$DepMap_ID[
                             mut.meta.train$Hugo_Symbol%in%gene&mut.meta.train$lineage==tissue]
                         ,col.idx,with=FALSE]
    wt<-achill.meta[achill.meta$lineage==tissue&
                      achill.meta$DepMap_ID%in%mut.meta.train$DepMap_ID[
                        mut.meta.train$Hugo_Symbol%in%gene==FALSE&mut.meta.train$lineage==tissue],
                    col.idx,with=FALSE]
    
    delta[[i]]<-data.frame(tissue=tissue,gene=gene,diff=
      median(mutants[[gene]])-median(wt[[gene]]),stringsAsFactors=F)
    i=i+1
  }
}

delta.df<-data.table(matrix(unlist(delta), nrow=length(delta), byrow=T))
colnames(delta.df)<-c("tissue","Hugo_Symbol","delta")
delta.df$gene_tissue_anno<-paste(delta.df$Hugo_Symbol,delta.df$tissue,sep=";")

delta.df.nona<-delta.df[is.na(delta.df$delta)==FALSE,]
delta.df.nona$delta<-as.numeric(delta.df.nona$delta)
delta.df.nona$normp<-pnorm(delta.df.nona$delta,
                      mean(delta.df.nona$delta),sd=sqrt(var(delta.df.nona$delta)))

delta.df.nona$q<-p.adjust(delta.df.nona$normp,'fdr')
delta.sig<-delta.df.nona[delta.df.nona$q<0.1,]

write.table(delta.df,"train_delta.txt",sep='\t',quote=F,row.names=F)

### annotate training data with features that can be learned by a model, converting ###
### categorical data to numerical data ###
mut.meta.train$sig_dep<-0
mut.meta.train$sig_dep[mut.meta.train$gene_tissue_anno%in%delta.sig$gene_tissue_anno]<-1

train.formatted<-mut.meta.train

train.formatted$isDeleterious[train.formatted$isDeleterious==TRUE]<-1
train.formatted$isDeleterious[train.formatted$isDeleterious==FALSE]<-0

train.formatted$isTCGAhotspot[train.formatted$isTCGAhotspot==TRUE]<-1
train.formatted$isTCGAhotspot[train.formatted$isTCGAhotspot==FALSE]<-0

train.formatted$isCOSMIChotspot[train.formatted$isCOSMIChotspot==TRUE]<-1
train.formatted$isCOSMIChotspot[train.formatted$isCOSMIChotspot==FALSE]<-0

#annotate protein changes with blosum scores
blosum<-fread("blosum62.txt",stringsAsFactors=F)

rownames(blosum)<-blosum$V1
train.formatted$blosum_score<-0

for(i in 1:nrow(train.formatted)){
  train.formatted$blosum_score[i]<-blosum[
    rownames(blosum)%in%substring(train.formatted$Protein_Change[i],3,3),
    colnames(blosum)%in%substring(train.formatted$Protein_Change[i],
                                  nchar(train.formatted$Protein_Change[i]),
                                  nchar(train.formatted$Protein_Change[i])),with=FALSE]
}

#annotate with expression data
express<-fread("CCLE_expression.csv",stringsAsFactors=F)
colnames(express)[1]<-"DepMap_ID"
colnames(express)<-gsub(" \\([0-9]+\\)","",colnames(express))
train.formatted$expression<-0

for(i in 1:nrow(train.formatted)){
  train.formatted$expression[i]<-tryCatch(express[express$DepMap_ID%in%train.formatted$DepMap_ID[i],
                                         train.formatted$Hugo_Symbol[i],with=FALSE],
                                         error=function(e){0})
}

write.table(train.formatted.df,"training_partially_formatted_11_17.txt",sep='\t',quote=F,row.names=F)
train.formatted.df$damaging<-0
train.formatted.df$non_damaging<-0
train.formatted.df$male<-0
train.formatted.df$female<-0

train.formatted.df$damaging[train.formatted$Variant_annotation=="damaging"]<-1
train.formatted.df$non_damaging[train.formatted$Variant_annotation=="other non-conserving"]<-1

train.formatted.df$male[train.formatted.df$sex=="Male"]<-1
train.formatted.df$female[train.formatted$sex=="Female"]<-1

train.formatted.df[,"n_damage":=sum(damaging),by="gene_tissue_anno"]
train.formatted.df[,"n_non_damage":=sum(non_damaging),by="gene_tissue_anno"]
train.formatted.df[,"mean_exp":=mean(expression,na.rm=TRUE),by="gene_tissue_anno"]
train.formatted.df[,"n_male":=sum(male),by="gene_tissue_anno"]
train.formatted.df[,"n_female":=sum(female),by="gene_tissue_anno"]
train.formatted.df[,"n_mutations":=.N,by="gene_tissue_anno"]
train.formatted.df[,"n_COSMIC":=mean(isCOSMIChotspot),by="gene_tissue_anno"]
train.formatted.df[,"n_TCGA":=mean(isTCGAhotspot),by="gene_tissue_anno"]
train.formatted.df[,"n_deleterious":=mean(isDeleterious),by="gene_tissue_anno"]

train.formatted.df$damaging_ratio<-train.formatted.df$n_damage/(
  train.formatted.df$n_non_damage+train.formatted.df$n_damage)

train.formatted.df$non_damaging_ratio<-train.formatted.df$n_non_damage/(
  train.formatted.df$n_non_damage+train.formatted.df$n_damage)

train.formatted.df$n_male<-train.formatted.df$n_male/train.formatted.df$n_mutations
train.formatted.df$n_female<-train.formatted.df$n_female/train.formatted.df$n_mutations

train.noDup<-train.formatted.df[!duplicated(train.formatted.df$gene_tissue_anno),]

write.table(train.noDup,"training_formatted_11_18.txt",sep='\t',quote=F,row.names=F)

ccle.meta<-metadata[,c("DepMap_ID","CCLE Name")]
rr<-fread("CCLE_RRBS_tss_CpG_clusters_20181022.txt",stringsAsFactors=F)

train.ccle<-merge(train.formatted.df,ccle.meta,by="DepMap_ID")

train.cols<-unique(train.ccle$`CCLE Name`)[unique(train.ccle$`CCLE Name`)%in%colnames(rr)]
rr.train<-rr[,c("cluster_id","avg_coverage",train.cols),with=FALSE]

rr.train$cluster_id<-gsub("_*","",rr.train$cluster_id)
rr.train.concat<-rr.train[,lapply(.SD, sum, na.rm=TRUE),by="cluster_id"]

train.ccle$cpg_dens<-0

for(i in 1:nrow(train.ccle)){
  if(train.ccle$`CCLE Name`[i]%in%colnames(rr.train.concat)){
    train.ccle$ccle_dens[i]<-rr.train.concat[
      rr.train.concat$cluster_id==train.ccle$Hugo_Symbol[i],
      train.ccle$`CCLE Name`[i],with=FALSE]
  }
}
