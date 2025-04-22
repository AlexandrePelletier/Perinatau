out<-'outputs/figures_perinatau'
dir.create(out)
source('../utils/r_utils.R') #https://github.com/AlexandrePelletier/alexandre-utils/blob/main/r_utils.R
source('../utils/visualisation.R')#https://github.com/AlexandrePelletier/alexandre-utils/blob/main/visualisation.R

# A/ régulons et significativité selon âge et sexe
res_reg<-fread('outputs/TauHFDRegulome/res_fgsea_HFDvsChow_in_regulons.csv.gz')
res_reg[,age:=ifelse(age=='M4','4-month-old','7-month-old')]
res_reg[,sex:=ifelse(sex=='M','♂','♀')]
res_reg[,group:=paste(age,sex,sep=' - ')]
ggplot(res_reg[pathway%in%pathway[padj<0.01]][genotype=='Tau'])+
  geom_col(aes(x=-log10(padj),y=pathway,fill=NES))+theme_bw()+
  scale_fill_gradient2(high = 'red',low = 'blue')+
  facet_wrap('group')+labs(y='regulon',x='-log10(adjusted p-value)')+
  geom_vline(xintercept = -log10(0.05),color='grey',linetype='dashed')

ggsave(fp(out,'barplot_regulons_enriched_Tau_hfdvschow.pdf'),width = 7,height = 6)


# B/ enrichissement des gènes associés aux différents régulons.
# Top5 voies par régulon 

res_ora<-fread('outputs/TauHFDRegulome/res_pathway_enriched_regulons.csv.gz')
res_ora[,top:=rank(pval)<=5&padj<0.25,by='regulon']
res_ora[,pathway:=str_remove(term,'GOBP|GOMF|GOCC|')|>str_replace_all('_',' ')|>str_to_sentence()]
res_ora[,pathway:=str_replace(pathway,'(r|R)na','RNA')]

pdf(fp(out,'heatmap_top5_padj0.25_pathway_enriched_regulons.pdf'),width = 8,height = 6)
CompPathways(res_ora[regulon%in%res_reg[padj<0.01]$pathway][term%in%term[(top)]],
             group.by = 'regulon',
             pathw_col = 'pathway',
             effect_col = 'fold.enrichment',
             max_color = 4,colors = c('white','red'))
dev.off()
