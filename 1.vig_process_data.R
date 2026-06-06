
## This script reads in vignette datasets for the mediation project
## limit to only g level and 2 outcomes

#clear environment
rm(list = ls())

#set working directory
# dir_wd<-
# setwd(dir_wd)

#load libraries
library("stringr")
library("dplyr")
library("ape")
library("readr")

###### functions #######
keep_species <-function(dat_taxa){
  temp = dat_taxa[grepl("s__",rownames(dat_taxa)) & (!grepl("t__",rownames(dat_taxa))) & (!grepl("unclassified$",rownames(dat_taxa))),]
  rownames(temp) = gsub(".*g__.*\\|","",rownames(temp))
  return(temp)
}

is_common <- function(otu_df,cutoff=0.1){
  num_row = dim(otu_df)[1]
  num_col = dim(otu_df)[2]
  otu_barcode = otu_df > 0.0001
  common_index = apply(otu_barcode,1,mean) > cutoff
  return(common_index)
}

top_abun <- function(df,num=10){
  # select top 25 abundunt species
  index = sort(apply(df,1,mean),decreasing=TRUE,index.return=TRUE)$ix[1:num]
  return(1:dim(df)[1] %in% index)
}

top_abun_median <- function(df,num=25){
  # select top 25 abundunt species
  index = sort(apply(df,1,median),decreasing=TRUE,index.return=TRUE)$ix[1:num]
  return(1:dim(df)[1] %in% index)
}


keep_pathways <-function(dat_path){
  temp = dat_path[!grepl("\\|",rownames(dat_path)),]
  # rownames(temp) = gsub(":.*","",rownames(temp))
  return(temp)
}

select_pathways <-function(dat_path,dat_taxa){
  taxa.names = rownames(dat_taxa)
  pathway.names = rownames(dat_path)
  index = rep(FALSE,length(pathway.names))
  for(taxa in taxa.names){
    index = index | grepl(taxa,pathway.names)
  }
  pathway.names.select = pathway.names[index]
  pathway.names.select = unique(gsub("\\|.*","",pathway.names.select))
  dat_path = dat_path[pathway.names.select,]
  # rownames(dat_path) = gsub(":.*","",rownames(dat_path))
  return(dat_path)
}

median_matrix <-function(dat_path,dat_taxa){
  taxa.names = rownames(dat_taxa)
  pathway.names = rownames(dat_path)
  #print(dim(dat_path))
  #print( sum(apply(dat_path>0,1,sum) <= 0 ) )
  
  index = rep(FALSE,length(pathway.names))
  for(taxa in taxa.names){
    index = index | grepl(taxa,pathway.names)
  }
  pathway.names = pathway.names[index]
  dat_path = dat_path[pathway.names,]
  #print(dim(dat_path))
  #print( sum(apply(dat_path>0,1,sum) <= 0 ) )
  
  v_median_all = apply(dat_path,1,median)
  pathway.names.id = gsub("\\|.*","",rownames(dat_path))
  pathway.names.id.unique = unique(pathway.names.id)
  
  out = matrix(0,nrow=length(taxa.names),ncol=length(pathway.names.id.unique))
  
  rownames(out) = taxa.names
  colnames(out) = pathway.names.id.unique
  
  for(taxa in taxa.names){
    for(pathway in pathway.names.id.unique){
      m = v_median_all[grepl(taxa,pathway.names) & pathway == pathway.names.id]
      if(length(m) ==1){
        out[taxa,pathway] = m
      }
      if(length(m) >1){
        print(length(m))
      }
    }
  }
  #out = out[,apply(out,2,median)>0]
  out = out[,names(sort(apply(out,2,sum),decreasing = TRUE)[1:50])]
  return(out)
}

same_median_matrix <- function(dat_median,dat_path){
  taxa.names = rownames(dat_median)
  pathway.names.id.unique = colnames(dat_median)
  
  out = matrix(0,nrow=length(taxa.names),ncol=length(pathway.names.id.unique))
  
  rownames(out) = taxa.names
  colnames(out) = pathway.names.id.unique
  
  v_median_all = apply(dat_path,1,median)
  
  pathway.names = rownames(dat_path)
  pathway.names.id = gsub("\\|.*","",rownames(dat_path))
  
  for(taxa in taxa.names){
    for(pathway in pathway.names.id.unique){
      m = v_median_all[grepl(taxa,pathway.names) & pathway == pathway.names.id]
      if(length(m) ==1){
        out[taxa,pathway] = m
      }
      if(length(m) >1){
        print(length(m))
      }
    }
  }
  return(out)
}

#Abundance filtering - filter out bugs (columns) that are not >0.1% in 5% of samples
abund_filter<-function(bugs_count){
  bugs<-bugs_count/apply(bugs_count,1,sum)
  col_keep_sum<-rep(NA,ncol(bugs))
  col_keep<-rep(TRUE,ncol(bugs))
  for(i in 1:ncol(bugs)){
    bugsum<-sum(bugs[,i]>=(0.1/100),na.rm=TRUE)
    #must be in at least 5% of samples
    if(bugsum<floor(nrow(bugs)*0.05)){
      col_keep[i]<-FALSE
    }
  }
  other_bugs<-!col_keep
  other<-rep(0,nrow(bugs))
  
  for(i in 1:ncol(bugs_count)){
    if(other_bugs[i]){
      other<-other+bugs_count[,i]
    }
  }
  bugs_count<-bugs_count[,col_keep]
  bugs_count<-cbind(bugs_count,other)
  
  
  return(bugs_count)
}

featurewise_pseudocount<-function(bugs){
  #by feature pseudocount
  bugs[,which(apply(bugs,2,sum)==0)]<-1
  for(i in 1:ncol(bugs)){
    col_i<-bugs[,i]
    min_val<-min(col_i[col_i>0])
    zero_vals<-which(col_i==0)
    if(length(zero_vals)>0){
      bugs[zero_vals,i]<-min_val/2
    }
  }
  return(bugs)
}

####################

#MLVS data
#Not normalized abundance data
bugs<-as.data.frame(t(read.csv(paste(dir_wd,"MLVS/bugs_dna_929_unFilt.tsv",sep=""),
                               header=TRUE,sep="\t",check.names=FALSE,row.names=1)))

df<-as.data.frame(read.csv(paste(dir_wd,"MLVS/meta_biomarker_468.csv",sep=""),
                           header=TRUE,sep=",",check.names=FALSE,row.names=1))
df_1<-df[which(df$cvisit==1),]
df_2<-df[which(df$cvisit==2),]
rm(df)
df<-rbind(cbind(df_1,"link"=rep("05",nrow(df_1))),
          cbind(df_1,"link"=rep("06",nrow(df_1))),
          cbind(df_2,"link"=rep("07",nrow(df_2))),
          cbind(df_2,"link"=rep("08",nrow(df_2))))
df$sample<-paste0(df$SubjectID,"_SF",df$link)
row.names(df)<-df$sample
rm(df_1,df_2)
df<-df[df$cvisit==2,] #subset to second time point


df<-df[intersect(df$sample,row.names(bugs)),]
bugs<-bugs[intersect(df$sample,row.names(bugs)),]
df<-df[,-which(colnames(df)=="link")]
df<-df[str_sort(row.names(df),decreasing=FALSE,numeric=TRUE),]
bugs<-bugs[str_sort(row.names(bugs),decreasing=FALSE,numeric=TRUE),]

#if subjects have multiple samples -> average them
sub<-df[,2,drop=F]
double<-unique(data.frame("SubjectID"=sub[sub$SubjectID %in% sub[duplicated(sub$SubjectID),"SubjectID"],]))
single<-data.frame("SubjectID"=sub$SubjectID[which(!(sub$SubjectID %in% double$SubjectID))])
sub$SampleID<-row.names(sub)
double<-merge(x=double,y=sub,by="SubjectID",all.x=TRUE)
single<-merge(x=single,y=sub,by="SubjectID",all.x=TRUE)
rm(sub)

bugs_single<-bugs[which(row.names(bugs) %in% single$SampleID),]
bugs_double<-bugs[which(row.names(bugs) %in% double$SampleID),]

bugs_merged<-data.frame(matrix(0,nrow=nrow(bugs_double)/2,ncol=ncol(bugs_double)))
colnames(bugs_merged)<-colnames(bugs_double)

row.names(double)<-double$SampleID
double<-double[str_sort(as.character(double$SampleID),decreasing=FALSE,numeric=TRUE),]
bugs_double<-bugs_double[str_sort(row.names(bugs_double),decreasing=FALSE,numeric=TRUE),]
row.names(bugs_merged)<-unique(double$SubjectID)

for(i in 1:nrow(bugs_merged)){
  bugs_merged[i,]<-apply(bugs_double[(i*2-1):(i*2),],2,mean)
}

row.names(bugs_single)<-sapply(strsplit(row.names(bugs_single),"_"), `[`, 1)
bugs<-rbind(bugs_single,bugs_merged)
rm(bugs_double,bugs_merged,bugs_single,double,single)
bugs<-bugs[str_sort(row.names(bugs),decreasing=FALSE,numeric=TRUE),]
# total of 256 subjects

df<-df[,-which(colnames(df)=="sample")]
df<-unique(df)
row.names(df)<-df$SubjectID
df<-df[str_sort(row.names(df),decreasing=FALSE,numeric=TRUE),]

#Make new metadata variables
df$treatment<-ifelse(df$emed122ch>median(df$emed122ch),1,0)
df$mscore_norm<-(df$mscore-mean(df$mscore))/sd(df$mscore)
meta.vars<-c("calor122cn", "totMETs_paq", "probio_2mo_qu", "stool_type", 
             "age_fecal", "pril12", "metfo12", "ant_12mo_qu")
meta.vars.names<-c('Total energy intake','Physical activity',
                   'Probiotics use','Bristol stool scale','Age',
                   'Proton pump inhibitors','Metformin','Antibiotics') 

# check confounding using residualized outcome
fit <- lm(mscore_norm ~ calor122cn + totMETs_paq + probio_2mo_qu + stool_type + age_fecal + pril12 + metfo12 + ant_12mo_qu,data=df) # fit the model
df$mscore_norm_predicted <- predict(fit)   # Save the predicted values
df$mscore_norm_residuals <- residuals(fit) # Save the residual values
save(df,file=paste0(dir_wd,"MLVS/metadata.RData"))

sum_phage<-bugs[,grep("k__Viruses$",colnames(bugs),perl=TRUE)]
bugs<-bugs[,-grep("k__Viruses",colnames(bugs))]
bugs<-bugs[,grep("g__",colnames(bugs))]
bugs<-bugs[,-grep("t__",colnames(bugs))]
col_s<-grep("s__",colnames(bugs))
col_g<-setdiff(seq(1,ncol(bugs)),col_s)
col_f<-col_g[grep("unclassified",colnames(bugs)[col_g])]

names_col_s<-sapply(strsplit(colnames(bugs)[col_s],"s__"), `[`, 2)
names_col_f<-sapply(strsplit(colnames(bugs)[col_f],"g__"), `[`, 2)
names_long<-c(colnames(bugs)[col_s],colnames(bugs)[col_f])
bugs<-bugs[,c(col_s,col_f)]
colnames(bugs)<-c(names_col_s,names_col_f)
# sort((apply(bugs,1,sum)+sum_phage)) #adds up to 100 except for rounding error
bugs<-bugs[str_sort(row.names(bugs),decreasing=FALSE,numeric=TRUE),]

names_long<-data.frame("long"=names_long,stringsAsFactors = FALSE)
names_long$kingdom<-sapply(strsplit(names_long$long,"|",fixed=TRUE), `[`, 1)
names_long$phylum<-sapply(strsplit(names_long$long,"|",fixed=TRUE), `[`, 2)
names_long$class<-sapply(strsplit(names_long$long,"|",fixed=TRUE), `[`, 3)
names_long$order<-sapply(strsplit(names_long$long,"|",fixed=TRUE), `[`, 4)
names_long$family<-sapply(strsplit(names_long$long,"|",fixed=TRUE), `[`, 5)
names_long$genus<-sapply(strsplit(names_long$long,"|",fixed=TRUE), `[`, 6)
names_long$species<-sapply(strsplit(names_long$long,"|",fixed=TRUE), `[`, 7)
for(i in 2:8){
  names_long[,i]<-sapply(strsplit(names_long[,i],"__",fixed=TRUE), `[`, 2)
}
g_rows<-which(is.na(names_long$species))
names_long$species[g_rows]<-names_long$genus[g_rows]

bugs_t<-data.frame(t(bugs))
bugs_t<-merge(x=names_long,y=bugs_t,by.x="species",by.y="row.names",all.x=TRUE, all.y=TRUE)
sample_cols<-seq(grep("genus",colnames(bugs_t))+1,ncol(bugs_t),1)
bugs_p<-aggregate(x=bugs_t[,sample_cols],by=list(grp=factor(bugs_t$phylum)),FUN=sum)
bugs_g<-aggregate(x=bugs_t[,sample_cols],by=list(grp=factor(bugs_t$genus)),FUN=sum)

row.names(bugs_p)<-bugs_p$grp
bugs_p<-bugs_p[,-1]
row.names(bugs_g)<-bugs_g$grp
bugs_g<-bugs_g[,-1]
bugs_p<-data.frame(t(bugs_p))
bugs_g<-data.frame(t(bugs_g))
rm(bugs_t)

# all species level
# all normalized by sample
# no categorized into other --> bugs_all
# with PC --> bugs_s_pc
# no PC --> bugs_s
bugs_all<-bugs # species level
bugs_all<-bugs_all/(apply(bugs_all,1,sum)) # normalize by sample
bugs<-abund_filter(bugs_count=bugs)
bugs_pc<-featurewise_pseudocount(bugs)
bugs<-bugs/apply(bugs,1,sum)
bugs_pc<-bugs_pc/apply(bugs_pc,1,sum)
bugs_s<-bugs
bugs_s_pc<-bugs_pc
save(bugs_s,file=paste0(dir_wd,"MLVS/bugs_s.RData"))
save(bugs_s_pc,file=paste0(dir_wd,"MLVS/bugs_s_pc.RData"))
save(bugs_all,file=paste0(dir_wd,"MLVS/bugs_all.RData"))

# genus level
bugs_g<-abund_filter(bugs_count=bugs_g)
bugs_g_pc<-featurewise_pseudocount(bugs_g)
bugs_g<-bugs_g/apply(bugs_g,1,sum)
bugs_g_pc<-bugs_g_pc/apply(bugs_g_pc,1,sum)
save(bugs_g,file=paste0(dir_wd,"MLVS/bugs_g.RData"))
save(bugs_g_pc,file=paste0(dir_wd,"MLVS/bugs_g_pc.RData"))

# phylum level
bugs_p<-abund_filter(bugs_count=bugs_p)
bugs_p_pc<-featurewise_pseudocount(bugs_p)
bugs_p<-bugs_p/apply(bugs_p,1,sum)
bugs_p_pc<-bugs_p_pc/apply(bugs_p_pc,1,sum)
save(bugs_p,file=paste0(dir_wd,"MLVS/bugs_p.RData"))
save(bugs_p_pc,file=paste0(dir_wd,"MLVS/bugs_p_pc.RData"))

names_long[nrow(names_long)+1,]<-rep("other",ncol(names_long))
save(names_long,file=paste0(dir_wd,"MLVS/names.RData"))
####################
####################

####################
####################
rm(list = setdiff(ls(),c("dir_wd","abund_filter","featurewise_pseudocount")))
#blueberry data

#Bring in the metadata
meta.dat = data.frame(read.csv(paste0(dir_wd,"blueberry/metadata_emma.csv"), header = T, row.names = 1))
meta.dat$sample<-row.names(meta.dat)

#exclude samples with low read count and no taxonomic information
#These are the individuals who did not stick to the study design/had poor DNA extraction
meta.dat = subset(meta.dat, subject != "s52" & subject != "s70" & subject != "s81" & group != "")

#Make metadata easier to work with: The lazy way
meta.dat$timepoint = as.character(meta.dat$timepoint)
meta.dat$metabolicsyndromecriteria = as.character(meta.dat$metabolicsyndromecriteria)
meta.dat$statin = as.character(meta.dat$statin)
meta.dat$group = as.character(meta.dat$group)
meta.dat$subject = as.character(meta.dat$subject)
meta.dat$bp_med = as.character(meta.dat$bp_med)
# str(meta.dat)

#Read in the taxonomy data
bugs = data.frame(t(read.csv(paste0(dir_wd,"blueberry/metaphlan2_taxonomic_profiles_emma.csv"),
                             header=TRUE,sep=",",check.names=FALSE,row.names=1)))

sum_phage<-bugs[,grep("k__Viruses$",colnames(bugs),perl=TRUE)]
bugs<-bugs[,-grep("k__Viruses",colnames(bugs))]
bugs<-bugs[,grep("g__",colnames(bugs))]
bugs<-bugs[,-grep("t__",colnames(bugs))]
col_s<-grep("s__",colnames(bugs))
col_g<-setdiff(seq(1,ncol(bugs)),col_s)
col_f<-col_g[grep("unclassified",colnames(bugs)[col_g])]

names_col_s<-sapply(strsplit(colnames(bugs)[col_s],"s__"), `[`, 2)
names_col_f<-sapply(strsplit(colnames(bugs)[col_f],"g__"), `[`, 2)
names_long<-c(colnames(bugs)[col_s],colnames(bugs)[col_f])

bugs<-bugs[,c(col_s,col_f)]
colnames(bugs)<-c(names_col_s,names_col_f)
# sort((apply(bugs,1,sum)+sum_phage)) #adds up to 100 except for rounding error
bugs<-bugs[intersect(meta.dat$sample,row.names(bugs)),]
bugs<-bugs[str_sort(row.names(bugs),decreasing=FALSE,numeric=TRUE),]

names_long<-data.frame("long"=names_long,stringsAsFactors = FALSE)
names_long$kingdom<-sapply(strsplit(names_long$long,".",fixed=TRUE), `[`, 1)
names_long$phylum<-sapply(strsplit(names_long$long,".",fixed=TRUE), `[`, 2)
names_long$class<-sapply(strsplit(names_long$long,".",fixed=TRUE), `[`, 3)
names_long$order<-sapply(strsplit(names_long$long,".",fixed=TRUE), `[`, 4)
names_long$family<-sapply(strsplit(names_long$long,".",fixed=TRUE), `[`, 5)
names_long$genus<-sapply(strsplit(names_long$long,".",fixed=TRUE), `[`, 6)
names_long$species<-sapply(strsplit(names_long$long,".",fixed=TRUE), `[`, 7)
for(i in 2:8){
  names_long[,i]<-sapply(strsplit(names_long[,i],"__",fixed=TRUE), `[`, 2)
}
g_rows<-which(is.na(names_long$species))
names_long$species[g_rows]<-names_long$genus[g_rows]

bugs_t<-data.frame(t(bugs))
bugs_t<-merge(x=names_long,y=bugs_t,by.x="species",by.y="row.names",all.x=TRUE, all.y=TRUE)
sample_cols<-seq(grep("genus",colnames(bugs_t))+1,ncol(bugs_t),1)
bugs_p<-aggregate(x=bugs_t[,sample_cols],by=list(grp=factor(bugs_t$phylum)),FUN=sum)
bugs_g<-aggregate(x=bugs_t[,sample_cols],by=list(grp=factor(bugs_t$genus)),FUN=sum)

# genus & phylum level
row.names(bugs_p)<-bugs_p$grp
bugs_p<-bugs_p[,-1]
row.names(bugs_g)<-bugs_g$grp
bugs_g<-bugs_g[,-1]
bugs_p<-data.frame(t(bugs_p))
bugs_g<-data.frame(t(bugs_g))
rm(bugs_t)

meta.dat<-meta.dat[str_sort(row.names(meta.dat),decreasing=FALSE,numeric=TRUE),]
meta_dat_6M = subset(meta.dat, timepoint == "6")
meta_dat_0M = subset(meta.dat, timepoint == "0")
identical(meta_dat_0M$subject, meta_dat_6M$subject) # TRUE
any(is.na(meta_dat_0M)) # TRUE
any(is.na(meta_dat_6M)) # TRUE

###################################
#This are the major variables that Aedin et al. have indicated are important. FMD is the one they piroirtize the most. 
#We have also found interesting results with the 2 or more category. (Subject is listed as a responder in two or more the "important" clinical endpoints)
meta_compare = data.frame(meta_dat_0M$subject, meta_dat_0M$group)
names(meta_compare) = c("subject", "group")
meta_compare$FMD_change = meta_dat_6M$FMD - meta_dat_0M$FMD
meta_compare$HDLC_change = meta_dat_6M$HDLC - meta_dat_0M$HDLC
meta_compare$HDLP_change = meta_dat_6M$HDLP - meta_dat_0M$HDLP
meta_compare$cGMP_change = meta_dat_6M$cGMPpmolmL - meta_dat_0M$cGMPpmolmL
meta_compare$AIx_change = meta_dat_6M$AIx75bpm - meta_dat_0M$AIx75bpm
meta_compare$HOMA_change = meta_dat_6M$HOMA_IR_-meta_dat_0M$HOMA_IR_
meta_compare$PWV_change = meta_dat_6M$PWV_-meta_dat_0M$PWV_
meta_compare$SBP_change = meta_dat_6M$SBP_-meta_dat_0M$SBP_
meta_compare$DBP_change = meta_dat_6M$DBP_-meta_dat_0M$DBP_
meta_compare$insulin_change = meta_dat_6M$insulin_-meta_dat_0M$insulin_
meta_compare$glucose_change = meta_dat_6M$glucose-meta_dat_0M$glucose
meta_compare$HB1Ac_change = meta_dat_6M$HB1Ac_-meta_dat_0M$HB1Ac_
meta_compare$totalNO_change = meta_dat_6M$totalNO-meta_dat_0M$totalNO
meta_compare$totalCholesterol_change = meta_dat_6M$totalCholesterol-meta_dat_0M$totalCholesterol
meta_compare$Trigl_change = meta_dat_6M$Trigl_-meta_dat_0M$Trigl_
meta_compare$framingham_CV_score_change = meta_dat_6M$framingham_CV_score-meta_dat_0M$framingham_CV_score
meta_compare$framingham_CV_score=meta_dat_6M$framingham_CV_score

meta_compare$sample<-meta_dat_6M$sample

meta_compare$cGMP_change_norm<-(meta_compare$cGMP_change-mean(meta_compare$cGMP_change))/sd(meta_compare$cGMP_change)

meta_compare$subject = as.character(meta_compare$subject)
meta_dat_6M <- meta_dat_6M[meta_dat_6M$subject %in% meta_compare$subject,]
identical(meta_dat_6M$subject, meta_compare$subject) # TRUE

grep("total_urine_derived_metabolites",colnames(meta_dat_6M))
grep("_24h",colnames(meta_dat_6M))

grep("total_serum_derived_metabolites",colnames(meta_dat_6M))
grep("_serum",colnames(meta_dat_6M))

df_urine<-meta_dat_6M[,c(grep("total_urine_derived_metabolites",colnames(meta_dat_6M)),
                         grep("_24h",colnames(meta_dat_6M)))]
df_serum<-meta_dat_6M[,c(grep("total_serum_derived_metabolites",colnames(meta_dat_6M)),
                         grep("_serum",colnames(meta_dat_6M)))]
which(apply(df_urine,2,sum)==0)
which(apply(df_serum,2,sum)==0)

df_urine<-df_urine[-which(apply(is.na(df_urine),1,sum)==ncol(df_urine)),]
df_serum<-df_serum[-which(apply(is.na(df_serum),1,sum)==ncol(df_serum)),]

df_urine<-df_urine[,-which(apply(is.na(df_urine),2,sum)>0)]
df_serum<-df_serum[,-which(apply(is.na(df_serum),2,sum)>0)]

pcs<-prcomp(df_urine,center=TRUE,scale=TRUE)
load_urine<-data.frame(pcs$x[,c(1,2,3)],stringsAsFactors = FALSE)
colnames(load_urine)<-paste0(colnames(load_urine),"_urine")
load_urine$sample<-row.names(load_urine)

pcs<-prcomp(df_serum,center=TRUE,scale=TRUE)
load_serum<-data.frame(pcs$x[,c(1,2,3)],stringsAsFactors = FALSE)
colnames(load_serum)<-paste0(colnames(load_serum),"_serum")
load_serum$sample<-row.names(load_serum)

meta_compare<-merge(x=meta_compare,y=load_urine,by.x="sample",all.x=TRUE)
meta_compare<-merge(x=meta_compare,y=load_serum,by.x="sample",all.x=TRUE)
row.names(meta_compare)<-meta_compare$sample
meta_compare<-meta_compare[str_sort(row.names(meta_compare),decreasing=FALSE,numeric=TRUE),]

meta_compare$treatment_1<-ifelse(meta_compare$group=="Y",1,0)
meta_compare$treatment_2<-ifelse(meta_compare$group=="Y"|meta_compare$group=="X",1,0)

save(meta_dat_6M,file=paste0(dir_wd,"blueberry/meta_dat_6M.RData"))
save(meta_compare,file=paste0(dir_wd,"blueberry/meta_compare.RData"))
df<-meta_compare
save(df,file=paste0(dir_wd,"blueberry/metadata.RData"))

###################################
# species level
# only at 6 months saved
# subjects that have NA values have been removed
bugs_all<-bugs
bugs_all<-bugs_all/(apply(bugs_all,1,sum))
bugs<-abund_filter(bugs_count=bugs)
bugs_pc<-featurewise_pseudocount(bugs)
bugs<-bugs/apply(bugs,1,sum)
bugs_pc<-bugs_pc/apply(bugs_pc,1,sum)
bugs_s<-bugs[meta_dat_6M$sample,]
bugs_s_pc<-bugs_pc[meta_dat_6M$sample,]

save(bugs_s,file=paste0(dir_wd,"blueberry/bugs_s.RData"))
save(bugs_s_pc,file=paste0(dir_wd,"blueberry/bugs_s_pc.RData"))
save(bugs_all,file=paste0(dir_wd,"blueberry/bugs_all.RData"))

# genus level
bugs_g<-abund_filter(bugs_count=bugs_g)
bugs_g_pc<-featurewise_pseudocount(bugs_g)
bugs_g<-bugs_g/apply(bugs_g,1,sum)
bugs_g_pc<-bugs_g_pc/apply(bugs_g_pc,1,sum)
bugs_g<-bugs_g[meta_dat_6M$sample,]
bugs_g_pc<-bugs_g_pc[meta_dat_6M$sample,]

save(bugs_g,file=paste0(dir_wd,"blueberry/bugs_g.RData"))
save(bugs_g_pc,file=paste0(dir_wd,"blueberry/bugs_g_pc.RData"))

# phylum level
bugs_p<-abund_filter(bugs_count=bugs_p)
bugs_p_pc<-featurewise_pseudocount(bugs_p)
bugs_p<-bugs_p/apply(bugs_p,1,sum)
bugs_p_pc<-bugs_p_pc/apply(bugs_p_pc,1,sum)
bugs_p<-bugs_p[meta_dat_6M$sample,]
bugs_p_pc<-bugs_p_pc[meta_dat_6M$sample,]

save(bugs_p,file=paste0(dir_wd,"blueberry/bugs_p.RData"))
save(bugs_p_pc,file=paste0(dir_wd,"blueberry/bugs_p_pc.RData"))

names_long[nrow(names_long)+1,]<-rep("other",ncol(names_long))
save(names_long,file=paste0(dir_wd,"blueberry/names.RData"))


##########
## parameter to be used for method running 
outcomes_mlvs<-c("mscore_norm","mscore_norm_residuals")
outcomes_blueberry <- c("cGMP_change_norm")

#no SparseMCMM & MedTest for now
# these two methods have separate scripts
df_vig<-tidyr::crossing(method=c("ccmm","hdma","hima_gaussian","hima_compositional","naive","pcr"), #no SparseMCMM for now
                 microbiome_feature_level=c("g"),
                 dataset=c("MLVS","blueberry"),
                 outcome=c(outcomes_mlvs,outcomes_blueberry))
df_vig<-df_vig[-which(df_vig$dataset=="MLVS" & df_vig$outcome %in% outcomes_blueberry),]
df_vig<-df_vig[-which(df_vig$dataset=="blueberry" & df_vig$outcome %in% outcomes_mlvs),]
df_vig$index<-seq(1,nrow(df_vig),1)
df_vig$exposure_var<-NA
df_vig$exposure_var[which(df_vig$dataset=="MLVS")]<-"treatment"
df_vig$exposure_var[which(df_vig$dataset=="blueberry")]<-"treatment_1"
df_vig$exposure_var[which(df_vig$outcome=="cGMP_change"|df_vig$outcome=="cGMP_change_norm")]<-"treatment_2"
save(df_vig, file = paste0(dir_wd, "results/df_vig.RData"))

## only assess treatment 2 for blueberry RCT
## which is dichotomized treatment (as long as have blueberry --> treatment2=1)


