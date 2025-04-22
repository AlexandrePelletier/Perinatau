
#proteo transcripto integr
out<-'outputs/TauHFDRegulome'
dir.create(out)
source('../utils/r_utils.R') #https://github.com/AlexandrePelletier/alexandre-utils/blob/main/r_utils.R
source('../utils/visualisation.R')#https://github.com/AlexandrePelletier/alexandre-utils/blob/main/visualisation.R

#install.packages('BiocManager')
#BiocManager::install("GENIE3")
#install.packages('R.utils')
#BiocManager::install("RcisTarget")

library(GENIE3)
library(RcisTarget)
library(fgsea)
library(edgeR)
library(DESeq2)
# Sys.setenv(LIBARROW_MINIMAL = "false") 
# Sys.setenv(ARROW_WITH_ZSTD = "ON")
# install.packages('arrow')
# Trscriptomic change asso to a specific TF? 

#load ranking
motifRankings<-importRankings('ref-data/cistarget/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather',
                              indexCol = 'motifs')
motifRankings@rankings
#load motif annotations
data("motifAnnotations_mgi_v9")
trans<-fread('ref-data/Normal/raw_expr/raw_rna_counts_male_4M.csv')[,.(gene_id,gene_name)]
tfmotifs<-merge(motifAnnotations_mgi_v9[,gene_name:=TF],
                             trans,by='gene_name')
tfs<-unique(tfmotifs$TF)
length(tfs)#1523

#O) #data quality check####
# get the full matrix of transcripto /proteo

#create metadata
#from samples names
samples_proteo<-setdiff(union(colnames(fread('ref-data/Tau/proteo_f_imput.csv')),
                      colnames(fread('ref-data/Tau/proteo_m_imput.csv'))),c('gene_id',paste0('V',1:100)))
samples_proteo<-data.table(prot_id=samples_proteo,
                           mouse=str_extract(samples_proteo,'[0-9]+$'))
samples_proteo[,sample_id:=prot_id][,
                                  sex:=str_extract(sample_id,'^F|M')][,
                                                                      genotype:=str_extract(sample_id,'T|W')][,
                                                                                                              diet:=str_extract(sample_id,'H|C')][,
                                                                                                                                                  age:=str_extract(sample_id,'7|4')]
fwrite(samples_proteo,'ref-data/Tau/metadata_proteo.csv')
mtd_proteo<-fread('ref-data/Tau/metadata_proteo.csv')

samples_rna<-colnames(fread('ref-data/transfer_8413177_files_667f6b10/perinatau_salmon.merged.gene_counts.tsv'))[-c(1,2)]
                             
samples_rna
samples_rna<-data.table(rna_id=samples_rna,
                        mouse=str_extract(samples_rna,'[0-9]+'))


samples_comm<-merge(samples_proteo,samples_rna)
samples_comm#41

#which samples not common?
setdiff(samples_rna$mouse,samples_proteo$mouse)
# [1] "241" "250" "269" "273" "278" "281" "286" "288" "290" "354" "356" "363" "364" "365" "368" "373" "387" "390" "391" "810"
# [21] "82"  "841" "845" "846" "850" "870" "871" "877" "882" "884" "887" "890" "891" "895"

setdiff(samples_proteo$mouse,samples_rna$mouse)
#[1] "809"  "84"   "385"  "282"  "294"  "849"  "3365"
fwrite(samples_comm,'ref-data/Tau/metadata_common.csv')
mtd_comm<-fread('ref-data/Tau/metadata_common.csv')

mtd_rna<-data.table(rna_id=samples_rna$rna_id,
                     sex=str_extract(samples_rna$rna_id,'(F[0-9]|M[0-9])$')|>str_remove('[0-9]'),
                     genotype=str_extract(samples_rna$rna_id,'T|W'),
                     diet=str_extract(samples_rna$rna_id,'H|C'),
                     age=str_extract(samples_rna$rna_id,'7$|4$')
)
fwrite(mtd_rna,'ref-data/Tau/metadata_rna.csv')
mtd_rna<-fread('ref-data/Tau/metadata_rna.csv')

#rna####
rna_files<-list.files('ref-data/Tau',pattern ='RNAseq.*\\.csv',full.names = T)

# genes<-intersect(fread(rna_files[1])[order(gene_id)]$gene_id,
#                  fread(rna_files[2])[order(gene_id)]$gene_id)


# rna_mat<-Reduce(function(x,y)merge(x,y,by=c('gene_id')),
#                 lapply(rna_files, function(f)fread(f)))[,.SD,.SDcols=c('gene_id',samples_comm)][genes,on='gene_id']

rna_mat<-fread('ref-data/transfer_8413177_files_667f6b10/perinatau_salmon.merged.gene_counts.tsv')[,.SD,.SDcols=c('gene_id',
                                                                                                                  'gene_name',
                                                                                                                  mtd_comm$rna_id)]
#rm duplicates
# rna_mat[duplicated(gene_id)]#113
# 
# rna_matf<-unique(rna_mat[,avg.expr:=rowMeans(.SD),.SDcols=is.numeric][order(gene_id,-avg.expr)],by='gene_id')[,-'avg.expr']
# 
# fwrite(rna_matf,'ref-data/Tau/raw_rna_counts_prot_match_clean.csv.gz')
# 

## Exclude Gm and Rik and apply 'size factors all'
annotGm=fread('ref-data/transfer_8413177_files_667f6b10/annotGm.csv')
annotRik=fread('ref-data/transfer_8413177_files_667f6b10/annotRik.csv')

finalbadGm=annotGm[which(annotGm[,2]==""),1]
finalbadRik=annotRik[which(annotRik[,2]==""),1]

rna_matf <- rna_mat[!gene_name %in% c(finalbadGm, finalbadRik)][,-'gene_name']
dim(rna_matf)
rna_matf[1:5,1:5]


#filter for expressed genes
#BiocManager::install("edgeR")
library(edgeR)
isexpr <- rowSums(cpm(data.frame(rna_matf,
                                 row.names = 'gene_id'))>1) >= 0.1 * length(mtd_comm$sample_id)
sum(isexpr) #16.5k

rna_matf <- rna_mat[isexpr,]
fwrite(rna_matf,fp(out,'raw_rna_counts_prot_match_clean_filtered.csv.gz'))
rna_matf<-fread(fp(out,'raw_rna_counts_prot_match_clean_filtered.csv.gz'))

plot(density(as.matrix(rna_matf[,-'gene_id']))) #not norm

#proteo####
prot_files<-list.files('ref-data/Tau',pattern ='proteo.*\\_imput.csv',full.names = T)
prot_mat<-Reduce(function(x,y)merge(x,y,by=c('gene_id')),
                lapply(prot_files, function(f)fread(f)))[,.SD,.SDcols=c('gene_id',
                                                                        mtd_proteo$sample_id)][str_detect(gene_id,'ENSMUSG')]


plot(density(as.matrix(prot_mat[,-'gene_id']))) #norm
fwrite(prot_mat,'ref-data/Tau/log2_prot_imput_abund.csv.gz')


proteo_mat<-fread('ref-data/Tau/log2_prot_imput_abund_rna_match_clean.csv.gz')
proteo_mat # 3289 prot
proteo_matf<-proteo_mat[,.SD,.SDcols=c('gene_id',
                          mtd_comm$sample_id)]

fwrite(prot_matf,'ref-data/Tau/log2_prot_imput_abund_rna_match_clean.csv.gz')

#filter for TFs prot
tf_mat<-proteo_matf[gene_id%in%tfmotifs$gene_id]
tf_mat#139

#merge tf abund and trscripto
tf_mat[,mol:=paste0(gene_id,'_prot')]
rna_matf[,mol:=paste0(gene_id,'_trs')]

#annot duplicates
tf_mat<-tf_mat[,avg.expr:=rowMeans(.SD),.SDcols=is.numeric][order(mol,-avg.expr)][,-'avg.expr']
tf_mat[duplicated(mol),mol:=paste0(mol,1:.N),by='mol']
#rename rna
setnames(rna_matf,
         mtd_comm$rna_id,
         mtd_comm$prot_id)
mols<-rbind(tf_mat,
            rna_matf)
dim(mols ) #43 common samples

fwrite(mols,fp(out,'tf_trscripto_matrix.csv'))
mols<-fread(fp(out,'tf_trscripto_matrix.csv'))


# Trscriptomic change asso to a specific TF?####
#bulk refined SCENIC approach: TF-gene weight, > module identif > module filtering and gene filtering based on TF motif enrichment > merging all filtered gene module by TF

# 1) TF-gene expr weight using GENIE3

#perform GENIE coregulatory network analysis
# tfs<-mols[str_detect(mol,'_prot$')]$mol

weight_mat<-GENIE3(as.matrix(data.frame(mols[,-c('gene_id')],row.names = 'mol')),
                   regulators =tf_mat$mol,
                   targets = rna_matf$mol,
                   nCores = 16)
dim(weight_mat)
weight_mat[1:10,1:10]

## Get ranking of edges
linkList <- getLinkList(weight_mat)
head(linkList)
plot(density(linkList$weight))
linkList<-data.table(linkList)

fwrite(linkList,fp(out,'tf_gene_weight.csv'))

table(linkList[weight>0.02]$regulatoryGene)

#add sens of correlation using spearman correl
linkList<-fread(fp(out,'tf_gene_weight.csv'))
mols<-fread(fp(out,'tf_trscripto_matrix.csv'))
samples<-colnames(mols)[str_detect(colnames(mols),'^F|M')]

linkList[,r:=cor(unlist(mols[mol==regulatoryGene,.SD,.SDcols=samples]),
                 unlist(mols[mol==targetGene,.SD,.SDcols=samples]),
                 method = 'pearson'),by=.(regulatoryGene,targetGene)]

fwrite(linkList,fp(out,'tf_gene_weight_and_spearman.csv.gz'))
linkList<-fread(fp(out,'tf_gene_weight_and_spearman.csv.gz'))

# cors<-sapply(unique(linkList$regulatoryGene),
#              function(tf)sapply(unique(linkList[regulatoryGene==tf]$targetGene),
#                                 function(g)cor(unlist(mols[mol==tf,.SD,.SDcols=samples]),
#                                                unlist(mols[mol==g,.SD,.SDcols=samples]),
#                                                method = 'spearman'
#                                                )
#                                 )
#              )



#2) module identif: use differnt methods
#In all the diff methods, only the links with IM > 0.001 were taken into account.
#+ module by pos or neg correlation
linkList[,reg:=ifelse(r>0,'up','down')]
plot(density(linkList$weight))
abline(v=0.02)


unique(linkList$regulatoryGene)
plot(density(linkList[regulatoryGene=='ENSMUSG00000060373_prot']$weight))
plot(density(linkList[regulatoryGene=='ENSMUSG00000026456_prot']$weight))
plot(density(linkList[regulatoryGene=='ENSMUSG00000024268_prot']$weight))
plot(density(linkList[regulatoryGene=='ENSMUSG00000032294_prot1']$weight))
plot(density(linkList[regulatoryGene=='ENSMUSG00000025142_prot']$weight))
plot(density(linkList[regulatoryGene=='Crtc1_prot']$weight))
plot(density(linkList[regulatoryGene=='Ctnnb1_prot']$weight))
plot(density(linkList[regulatoryGene=='Ddx5_prot']$weight))
abline(v=0.01)
linkListf<-linkList[weight>0.01]

#(i) setting several IM thresholds (IM > 0.02 and IM > 0.05),
#(ii) taking the 50 targets with highest IM for each TF

tfs<-unique(linkListf$regulatoryGene)
mods_groups<-c('up','down')

mods<-lapply(mods_groups,function(s){
  mod_tfs<-lapply(tfs,function(tf){
    return(list(thr0.02=linkListf[reg==s&regulatoryGene==tf][weight>0.02]$targetGene,
                thr0.05=linkListf[reg==s&regulatoryGene==tf][weight>0.05]$targetGene,
                top50=linkListf[reg==s&regulatoryGene==tf][head(order(-weight),50)]$targetGene))
    
  })
  names(mod_tfs)<-tfs
  return(mod_tfs)
})
names(mods)<-mods_groups

#(iii) keeping only the top 5, and 10  TFs for each target gene (then, split by TF).
mods2<-list(tartop5=5,
            tartop10=10)
mods2dt<-lapply(mods2,function(t){
  return(linkListf[,.SD[head(order(-weight),t)],by=.(targetGene,reg)])
})

for(s in names(mods)){
  
  for(tf in names(mods[[s]])){
    
    for(t in names(mods2dt)){
      mods[[s]][[tf]][[t]]<-mods2dt[[t]][reg==s&regulatoryGene==tf]$targetGene
      
    }
    
  }
}
unlist(lapply(mods$up,function(tfmods)lapply(tfmods,function(x)length(x))))

saveRDS(mods,fp(out,'tf_modules_before_filtering.rds'))
mods<-readRDS(fp(out,'tf_modules_before_filtering.rds'))

#3) module filtering and gene filtering based on TF motif enrichment
#rm '_trs' to have similar gene names
mods<-lapply(mods, function(lsens)lapply(lsens,function(lTF)lapply(lTF,function(lTHR)str_remove(lTHR,'_trs'))))
#rm '_prot' to have similar prot names
mods<-lapply(mods, function(lsens)setNames(lsens,str_remove(names(lsens),'_prot')))



#Run enrichment with fgsea####
#get motifs ranks
tfs<-str_remove(tfs,'_prot')
motifAnnof<-tfmotifs[gene_id%in%tfs]
motifAnnof
motifs_ranks<-melt(data.table(data.frame(motifRankings@rankings,row.names = 'motifs'),
                              keep.rownames = 'motif'),id.vars = 'motif',variable.name = 'gene_name',value.name = 'rank')

#annotate with ensembl ID of the TF
motifAnnof[,TF_ensembl:=gene_id]
motifs_ranks<-merge(motifs_ranks,motifAnnof[,-c('description','gene_id','gene_name')],by=c('motif'))

#annotate with ensembl ID for the genes target

motifs_ranks<-merge(motifs_ranks,
                    trans,by='gene_name',allow.cartesian = T)


#filter for expressed genes
rna_matf<-fread(fp(out,'raw_rna_counts_prot_match_clean_filtered.csv.gz'))

genes<-rna_matf$gene_id
length(genes) #16.5k
motifs_ranksf<-motifs_ranks[gene_id%in%genes]

motifs<-unique(motifs_ranksf$motif)
length(motifs)#104

#rm non unique motif-TF ranking
motifs_ranksff<-unique(motifs_ranksf[order(motif,TF,-directAnnotation,-inferred_Orthology,rank)],by=c('motif','TF','gene_name'))
unique(motifs_ranksff$TF)#42

#test motif enrichment for tested tfs 
res_up<-motifs_ranksff[,fgsea(pathways = mods[['up']][[TF_ensembl[1]]],
                             stats = setNames(rank,gene_id),
                             minSize = 5,scoreType = "pos"),by=.(motif,TF_ensembl)]
res_up[pval<0.01] #
res_up[,padj:=p.adjust(pval,n = length(unique(motifs_ranksff$TF))),by=.(motif,TF_ensembl)]
unique(res_up[pval<0.01]$TF) #7
unique(res_up[padj<0.1]$TF) #6

res_dn<-motifs_ranksf[,fgsea(pathways = mods[['down']][[TF_ensembl[1]]],
                             stats = setNames(rank,gene_id),
                             minSize = 5,scoreType = "pos"),by=.(motif,TF_ensembl)]
res_dn[pval<0.01] #
res_dn[,padj:=p.adjust(pval,n = length(unique(motifs_ranksff$TF))),by=.(motif,TF_ensembl)]
unique(res_dn[pval<0.01]$TF) #7
unique(res_dn[padj<0.1]$TF) #6

#merge 
res<-rbind(res_up[,sens:='upregulation'],res_dn[,sens:='downregulation'])
res<-merge(res,unique(motifAnnof[,.(TF,TF_ensembl)]))

fwrite(res,fp(out,'res_modules_tf_motif_enrichment.csv.gz'))

#(i) module filtering: 
#add NES threshold?
plot(density(res$NES))
#thr: pval<0.01
res_sig<-res[pval<0.01]

#(ii) gene filtering
#get all Leading Edges of signif modules by TFs
TFregs<-res_sig[,Reduce(union,leadingEdge),by=c('TF','sens')]
setnames(TFregs,'V1','gene_id')

TFregs<-merge(TFregs,trans)


#stats
head(TFregs[order(TF,gene_name)],100)
TFregs[,size.regulon:=.N,by=c('TF','sens')]
TFregs[,regulon:=paste(TF,str_remove(sens,'regulation'),sep='_')]
table(TFregs$regulon)
  # Bclaf1_up    Ctbp1_up Ctnnb1_down    Hdac2_up    Mecp2_up     Nono_up Prdm16_down 
  #       892        4521        1808           4        2927        1592          23 
  # Sfpq_down  Smarca4_up Trim28_down    Ube2k_up 
  #      4122          54         220          26 
fwrite(TFregs,fp(out,'TF_regulons.csv'))
linkList
linkList[,gene_id:=str_extract(regulatoryGene,'ENSMUSG[0-9]+')]
linkList=merge(linkList,unique(TFregs[,.(gene_id,gene_name)]))
setnames(linkList,'gene_id','regulatoryGeneid')
setnames(linkList,'gene_name','regulatoryGenename')

TFregs[,avg.link:=mean(linkList[targetGene%in%gene_id&regulatoryGenename==TF[1]]$weight,na.rm=T),by='regulon']
ggplot(TFregs)+geom_bar(aes(x=regulon,fill=avg.link))+theme_bw()+
  scale_x_discrete(guide = guide_axis(angle=60))+scale_y_log10()


#Are the regulons enriched in HFD vs Chow DEs genes? #### 
TFregs<-fread(fp(out,'TF_regulons.csv'))
res_all<-fread(fp(out,'res_deseq2_hfdvschow.csv.gz'))
#filter for expressed genes and tested genes
res_def<-res_all[gene%in%genes&!is.na(padj)]
# res_def[,contrast:=paste0(test,'vs',ref)]
#rm duplicate
res_def<-unique(res_def[order(pvalue)],by = c('gene','sex','genotype','age'))
unique(res_def,by='gene') #14k
res_def[,gene_stat:=sign(log2FoldChange)*-log10(pvalue),by=c('sex','genotype','age')]

res_reg<-res_def[,fgsea(pathways = split(TFregs$gene_id,TFregs$regulon),
                        stats = setNames(gene_stat,gene)),by=c('sex','genotype','age')]

table(res_reg[padj<0.01][,.(paste(genotype,sex,age))])

fwrite(res_reg,fp(out,'res_fgsea_HFDvsChow_in_regulons.csv.gz'))
res_reg<-fread(fp(out,'res_fgsea_HFDvsChow_in_regulons.csv.gz'))

ggplot(res_reg[pathway%in%pathway[padj<0.01]][genotype=='Tau'])+
  geom_col(aes(x=-log10(padj),y=pathway,fill=NES))+theme_bw()+
  scale_fill_gradient2(high = 'red',low = 'blue')+facet_wrap(sex~age)+labs(y='regulon')+
  geom_vline(xintercept = -log10(0.05),color='grey',linetype='dashed')
ggsave(fp(out,'barplot_regulons_enriched_Tau_hfdvschow.pdf'),width = 7,height = 6)

#same LE F and M ?
res_reg[padj<0.01]$pathway|>unique()

leenr<-res_reg[pathway%in%pathway[padj<0.01]&genotype=='Tau',OR3(strsplit(leadingEdge[sex=='F'&age=='M7'],
                                                                   '\\|')[[1]],
                                                          list('male4'=strsplit(leadingEdge[sex=='M'&age=='M4'],'\\|')[[1]]),
                                                          background =TFregs[regulon==pathway[1]]$gene_id),by='pathway']
leenr[padj<0.05][,-'genes.overlap']
# pathway   term term.size n.query n.overlap pct.query.overlap precision pct.term.overlap background_size
# <char> <char>     <int>   <int>     <int>             <num>     <num>            <num>           <int>
#   1:   Bclaf1_up  male4       314     486       205         0.4218107 0.4218107        0.6528662           16588
# 2:    Ctbp1_up  male4      1619    1929       746         0.3867289 0.3867289        0.4607783           16588
# 3: Ctnnb1_down  male4       720     732       371         0.5068306 0.5068306        0.5152778           16588
# 4:    Hdac2_up  male4         1       1         0         0.0000000 0.0000000        0.0000000           16588
# 5:    Mecp2_up  male4      1230    1069       504         0.4714687 0.4714687        0.4097561           16588
# 6:     Nono_up  male4       735     635       313         0.4929134 0.4929134        0.4258503           16588
# 7: Prdm16_down  male4         1      11         0         0.0000000 0.0000000        0.0000000           16588
# 8:   Sfpq_down  male4      1541    1554       686         0.4414414 0.4414414        0.4451655           16588
# 9:  Smarca4_up  male4        15      16         6         0.3750000 0.3750000        0.4000000           16588
# 10: Trim28_down  male4        76      67        20         0.2985075 0.2985075        0.2631579           16588
# pct.term.background          pval          padj fold.enrichment
# <num>         <num>         <num>           <num>
#   1:        1.892935e-02 1.267182e-250 1.267182e-250       22.283426
# 2:        9.760068e-02 7.139994e-310 7.139994e-310        3.962359
# 3:        4.340487e-02  0.000000e+00  0.000000e+00       11.676814
# 4:        6.028454e-05  1.000000e+00  1.000000e+00        0.000000
# 5:        7.414999e-02 1.842130e-311 1.842130e-311        6.358311
# 6:        4.430914e-02 2.783122e-271 2.783122e-271       11.124418
# 7:        6.028454e-05  1.000000e+00  1.000000e+00        0.000000
# 8:        9.289848e-02  0.000000e+00  0.000000e+00        4.751869
# 9:        9.042681e-04  1.379965e-15  1.379965e-15      414.700000
# 10:        4.581625e-03  5.369330e-32  5.369330e-32       65.153181
leenr[,padj:=padj+1e-310]
ggplot(leenr)+geom_col(aes(y=pathway,x=pct.query.overlap*100,fill=-log10(pval)))+
  labs(x='% of leading edges genes in 4M male also in 7M female')+theme_bw()

#comp to WT


ggplot(res_reg[pathway%in%pathway[padj<0.01]][genotype=='WT'])+
  geom_col(aes(x=-log10(padj),y=pathway,fill=NES))+theme_bw()+
  scale_fill_gradient2(high = 'red',low = 'blue')+facet_wrap(sex~age)+labs(y='regulon')+
  geom_vline(xintercept = -log10(0.05),color='grey',linetype='dashed')

#boxplot PC1 of regulons
mat<-fread(fp(out,'vst_norm_counts.csv.gz'))|>data.frame(row.names = 'gene_id')|>as.matrix()
res_pca<-rbindlist(lapply(unique(TFregs$regulon),function(reg){
  regulon<-TFregs[regulon==reg]$gene_id
  
  pca<-prcomp(t(mat[rna_matf[gene_id%in%regulon]$gene_id,]),scale.=TRUE)
  pca_dt<-merge(data.table(pca$x[,1,drop=F],keep.rownames = 'rna_id'),mtd_rna)
  return(pca_dt[,regulon:=reg])
}),fill = T)


#PCsignCorr: allign the sign of the PC coordinate with the expression matrix. ie. a sample with a high expression will have the higher PC coord
#x: pc coords with sample names
#mat= scaled normalized expression matrix used for PC. if not for correct the PC1, should be filtered to have only the feature contributing tho the PC to be corrected
PCsignCorr<-function(x,mat,pc_col='PC1',sample_col='sample_id'){
  if(colSums(mat)[which.max(x)]-colSums(mat)[which.min(x)]<0){
    message('top and bottom sample have been inverted, correcting')
    return(-x)
  }else{
    message('PC and expression matrix well aligned')
    return(x)
  }
}

matsc<-t(scale(t(mat)))
setnames(res_pca,'regulon','reg')     
res_pca[,PC1corr:=PCsignCorr(setNames(PC1,rna_id),matsc[TFregs[regulon==reg[1]]$gene_id,rna_id]),by='reg']
setnames(res_pca,'reg','regulon')    

regs=res_reg[padj<0.001][genotype=='Tau']$pathway
ggplot(res_pca[regulon%in%regs][genotype=='T'],aes(x=paste(age,sex),y=PC1corr,fill=diet))+
  geom_boxplot()+
  facet_wrap('regulon',scales = 'free_y')+theme_bw()+
  scale_x_discrete(guide = guide_axis(angle = 60))+labs(y='Regulon Eigengene')
ggsave(fp(out,'boxplot_pc1_regulons_enriched_Tau_hfdvschow.pdf'),width = 7,height = 7)


#stats hfd vs
ggplot(res_pca)+
  geom_density(aes(x=PC1,col=regulon))+
  theme_bw()+facet_wrap('regulon')
reg_changes<-res_pca[genotype=='T'][,.(pval=wilcox.test.test(PC1[diet=='H'],
                                           PC1[diet=='C'])$p.value),
                                    by=c('regulon','sex','genotype','age')]
reg_changes[pval<0.11]
res_pca[genotype=='T',data.table(summary(lm(PC1~age*diet))$coefficients,
                                 keep.rownames = 'cov'),
        by=c('regulon','sex')][`Pr(>|t|)`<0.1]

table(res_pca[genotype=='T'][age==7][sex=='M'][,.(regulon,diet)])

#Smarca4
ggplot(res_pca[regulon%in%regs][regulon=='Smarca4_up'][genotype=='T'],
       aes(x=paste(age,sex),y=PC1corr,col=diet))+
  geom_boxplot()+
  geom_point(position = position_dodge(width = 0.7))+
  theme_bw()+
  scale_x_discrete(guide = guide_axis(angle = 60))+labs(y='Regulon Eigengene')


genes_smarca4<-res_reg[padj<0.001][genotype=='Tau'][pathway=='Smarca4_up']$leadingEdge|>strsplit('\\|')|>unlist()

setnames(res_all,'gene','gene_id')
res_all<-merge(res_all,rna_matf[,.(gene_id,gene_name)])

res_all[gene_id%in%genes_smarca4][order(pvalue)][1:20]
table(res_all[padj<0.05][,.(sex,genotype,age)])
pdf(fp(out,'heatmap_HFDvsChow_logFC_Smarca4_regulon_leadingedge_only.pdf'),width = 5,height = 6)
CompDEGs(res_all[gene_id%in%genes_smarca4],group.by = c('sex','genotype','age'),
         pval_column = 'pvalue',gene_column = 'gene_name',col_range = c(-1.5,1.5))
dev.off()

pdf(fp(out,'heatmap_HFDvsChow_logFC_Smarca4_regulon.pdf'),width = 5,height = 7)
CompDEGs(res_all[gene_id%in%TFregs[regulon=='Smarca4_up']$gene_id],group.by = c('sex','genotype','age'),
         pval_column = 'pvalue',gene_column = 'gene_name',col_range = c(-2,2))
dev.off()


#see prot change
prot<-fread('ref-data/Tau/log2_prot_imput_abund.csv.gz')
prot[gene_id=='ENSMUSG00000032187']
smar=data.table(sample_id=colnames(prot)[-1],abund=2^prot[gene_id=='ENSMUSG00000032187'][,-'gene_id']|>as.vector()|>as.numeric())
smar<-merge(smar,mtd_proteo)
ggplot(smar[genotype=='T'])+geom_boxplot(aes(x=paste(age,sex),y=abund,fill=diet))+theme_bw()
plot(density(smar$abund))

#OR Smarca4
msigdb<-fread('~/adpelle1@bu.edu - Google Drive/My Drive/P.Alexandre_Pelletier/Common_resources/MSigDB/all_CPandGOs_gene_and_genesets_mouse.csv.gz')
msigdb$pathway|>unique()|>length()
res_or<-OR3(TFregs[regulon=='Smarca4_up']$gene_id,terms_list = split(msigdb$ensembl_gene,msigdb$pathway),background = genes)
res_or[,mlog10padj:=-log10(padj)]
emmaplot(res_or[padj<0.05],col.var = 'mlog10padj',min_edge = 0.5)

#pathways in all regulons
res_ora<-TFregs[,OR3(gene_id,
                     terms_list = split(msigdb$ensembl_gene,msigdb$pathway),
                     background = genes),by='regulon']

res_ora[,mlog10padj:=-log10(padj)]
res_ora[,top:=rank(pval)<=12&padj<0.25,by='regulon']
res_ora[,log2FE:=log2(fold.enrichment+1)]

pdf(fp(out,'heatmap_top12_padj0.25_pathway_enriched_regulons.pdf'),
    width = 12,height = 8)
CompPathways(res_ora[term%in%term[(top)]],
             group.by = 'regulon',
             pathw_col = 'term',
             effect_col = 'fold.enrichment',
             max_color = 4,colors = c('white','red'))
dev.off()
fwrite(res_ora,fp(out,'res_pathway_enriched_regulons.csv.gz'))
fwrite(res_ora[padj<0.25],fp(out,'res_pathway_enriched_regulons_padj0.25.csv'))

#emmaplot per regulons

for(reg in regs){
  res_oraf<-res_ora[regulon==reg][padj<0.05]|>head(80)
  if(nrow(res_oraf)>1){
    emmaplot(res_oraf,
             col.var = 'mlog10padj',
             min_edge = 0.5)
    ggsave(fp(out,ps('emmaplot_',reg,'.pdf')),width = 10,height = 7)
  }
 
}

#TF abundance boxplot
res_ora<-fread(fp(out,'res_pathway_enriched_regulons_padj0.25.csv'))
tfs<-res_ora$regulon|>unique()|>str_remove('_down|_up')
abund
res_prot<-tfs
prot_mat<-fread('ref-data/Tau/log2_prot_imput_abund.csv.gz')
trans<-fread('ref-data/Normal/raw_expr/raw_rna_counts_male_4M.csv')[,.(gene_id,gene_name)]
mtd_proteo<-fread('ref-data/Tau/metadata_proteo.csv')
table(mtd_proteo[,.(genotype,age,sex,diet)])
prot_mat<-merge(prot_mat,trans)
prot_dt<-melt(prot_mat,id.vars = c('gene_id','gene_name'),
              variable.name = 'prot_id',value.name = 'log2abund')
prot_dt<-merge(prot_dt,mtd_proteo)
ggplot(prot_dt[genotype=='T'][gene_name%in%tfs])+geom_boxplot(aes(x=paste(age,sex),
                                                y=log2abund,fill=diet))+theme_bw()+
  facet_wrap('gene_name',scales = 'free_y')
ggsave(fp(out,'boxplot_TFabundance_for_regulons_enriched_Tau_hfdvschow.pdf'),
       width = 9,height = 9)

tf_changes<-prot_dt[genotype=='T'][gene_name%in%tfs][,.(pval=wilcox.test(log2abund[diet=='H'],
                                                             log2abund[diet=='C'])$p.value),
                                    by=c('gene_name','sex','genotype','age')]
tf_changes[pval<0.11]






