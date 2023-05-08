library(plotly)
library(reshape2)
metadata = read.table('cornell-bladder-metadata-final.txt',header=TRUE,as.is=TRUE,sep='\t',fill=T)

sbs_data = read.table('sbs96.HighConf.cosmicv3.2.benign.counts.txt',header=TRUE, as.is=TRUE, sep="\t")
dbs_data = read.table('dbs78.HighConf.cosmicv3.2.benign.counts.txt',header=TRUE, as.is=TRUE, sep="\t")

postchemo = sbs_data[sbs_data$Sample %in% metadata[metadata$platinum_chemotherapy == 'Post-chemo',]$tumor,]
postchemo = merge(postchemo[,c('Sample','SBS31','SBS35')],metadata[,c('tumor','patient','months_from_chemo_start_to_sample_collection')], by.x='Sample',by.y='tumor')
postchemo$Chemotherapy = (postchemo$SBS31+postchemo$SBS35)/postchemo$months_from_chemo_start_to_sample_collection
postchemo = aggregate(postchemo, by=list(postchemo$patient), FUN=mean)
postchemo$Chemotherapy[postchemo$Chemotherapy==0] = NA

sbs_data = merge(sbs_data[,c('Sample','SBS1','SBS2','SBS13','SBS31','SBS35')],metadata[,c('tumor','patient','months_from_birth_to_sample_collection','months_from_diagnosis_to_sample_collection','months_from_chemo_start_to_sample_collection')],by.x='Sample',by.y='tumor')
sbs_data$Aging = sbs_data$SBS1/sbs_data$months_from_birth_to_sample_collection
sbs_data$APOBEC = (sbs_data$SBS2+sbs_data$SBS13)/(sbs_data$months_from_diagnosis_to_sample_collection+120)
sbs_data$Chemotherapy = (sbs_data$SBS31+sbs_data$SBS35)/sbs_data$months_from_chemo_start_to_sample_collection
sbs_data = aggregate(sbs_data, by=list(sbs_data$patient), FUN=mean,na.rm=T)

sbs_data$Aging[sbs_data$Aging==0] = NA
sbs_data$APOBEC[sbs_data$APOBEC==0] = NA

postchemo_dbs = dbs_data[dbs_data$Sample %in% metadata[metadata$platinum_chemotherapy == 'Post-chemo',]$tumor,]
postchemo_dbs = merge(postchemo_dbs[,c('Sample','DBS5')],metadata[,c('tumor','patient','months_from_chemo_start_to_sample_collection')], by.x='Sample',by.y='tumor')
postchemo_dbs$Chemotherapy = (postchemo_dbs$DBS5)/postchemo_dbs$months_from_chemo_start_to_sample_collection
postchemo_dbs = aggregate(postchemo_dbs, by=list(postchemo_dbs$patient), FUN=mean)
postchemo_dbs$Chemotherapy[postchemo_dbs$Chemotherapy==0] = NA

dbs_data = merge(dbs_data[,c('Sample','DBS5','DBS11')],metadata[,c('tumor','patient','months_from_diagnosis_to_sample_collection','months_from_chemo_start_to_sample_collection')], by.x='Sample',by.y='tumor')
dbs_data$APOBEC = dbs_data$DBS11/(dbs_data$months_from_diagnosis_to_sample_collection + 120)
dbs_data$Chemotherapy = (dbs_data$DBS5)/dbs_data$months_from_chemo_start_to_sample_collection
dbs_data = aggregate(dbs_data, by=list(dbs_data$patient), FUN=mean,na.rm=T)
dbs_data$APOBEC[dbs_data$APOBEC==0] = NA

fig <- plot_ly(type = "box")
fig <- fig %>% add_trace(y = postchemo$Chemotherapy, name='Chemotherapy (SBS31/35)',boxpoints = 'all',pointpos = 0,
                         marker = list(color = '#b2182b'),
                         line = list(color = '#b2182b'),
                         fillcolor = list(color = '#e31a1c'))
fig <- fig %>% add_trace(y = postchemo_dbs$Chemotherapy, name='Chemotherapy (DBS5)', boxpoints = 'all',pointpos = 0 ,
                         marker = list(color = '#ef8a62'),
                         line = list(color = '#ef8a62'),
                         fillcolor = list(color = '#e31a1c'))
fig <- fig %>% add_trace(y = sbs_data$APOBEC, name='APOBEC (SBS2/13)',boxpoints = 'all',pointpos = 0,
                        marker = list(color = '#2166ac'),
                        line = list(color = '#377eb8'),
                        fillcolor = list(color = '#4daf4a'))
fig <- fig %>% add_trace(y = dbs_data$APOBEC, name='APOBEC (DBS11)', boxpoints = 'all',pointpos = 0 ,
                         marker = list(color = '#67a9cf'),
                         line = list(color = '#67a9cf'),
                         fillcolor = list(color = '#4daf4a'))
fig <- fig %>% add_trace(y = sbs_data$Aging, name='Aging (SBS1)', boxpoints = 'all',pointpos = 0 ,
                         marker = list(color = '#1b7837'),
                         line = list(color = '#1b7837'),
                         fillcolor = list(color = '#4daf4a'))
fig <- layout(fig, yaxis = list(type = "log", title='Mutations per month of exposure'))

fig

wilcox.test(sbs_data$APOBEC, sbs_data$Aging)$p.value
wilcox.test(sbs_data$APOBEC, postchemo$Chemotherapy)$p.value
wilcox.test(postchemo$Chemotherapy, sbs_data$Aging)$p.value
wilcox.test(dbs_data$APOBEC, postchemo_dbs$Chemotherapy)$p.value

sbs_data = read.table('sbs96.HighConf.cosmicv3.2.benign.counts.txt',header=TRUE, as.is=TRUE, sep="\t")
dbs_data = read.table('dbs78.HighConf.cosmicv3.2.benign.counts.txt',header=TRUE, as.is=TRUE, sep="\t")
metadata = read.table('cornell-bladder-metadata-final.txt',header=TRUE,as.is=TRUE,sep='\t',fill=T)
sbs_data = merge(sbs_data[,c('Sample','SBS1','SBS2','SBS13')],metadata[,c('tumor','patient','months_from_birth_to_sample_collection','months_from_diagnosis_to_sample_collection')],by.x='Sample',by.y='tumor')

sbs_data$APOBEC_10years = (sbs_data$SBS2+sbs_data$SBS13)/(sbs_data$months_from_diagnosis_to_sample_collection+120)
sbs_data$APOBEC = (sbs_data$SBS2+sbs_data$SBS13)/(sbs_data$months_from_diagnosis_to_sample_collection)
sbs_data$APOBEC_birth = (sbs_data$SBS2+sbs_data$SBS13)/(sbs_data$months_from_birth_to_sample_collection)
sbs_data = aggregate(sbs_data, by=list(sbs_data$patient), FUN=mean,na.rm=T)
sbs_data$APOBEC_10years[sbs_data$APOBEC_10years==0] = NA
sbs_data$APOBEC[sbs_data$APOBEC==0] = NA
sbs_data$APOBEC_birth[sbs_data$APOBEC_birth==0] = NA

dbs_data = merge(dbs_data[,c('Sample','DBS11')],metadata[,c('tumor','patient','months_from_diagnosis_to_sample_collection','months_from_birth_to_sample_collection')], by.x='Sample',by.y='tumor')
dbs_data$APOBEC_10years = dbs_data$DBS11/(dbs_data$months_from_diagnosis_to_sample_collection + 120)
dbs_data$APOBEC = dbs_data$DBS11/(dbs_data$months_from_diagnosis_to_sample_collection)
dbs_data$APOBEC_birth = dbs_data$DBS11/(dbs_data$months_from_birth_to_sample_collection)
dbs_data = aggregate(dbs_data, by=list(dbs_data$patient), FUN=mean,na.rm=T)
dbs_data$APOBEC_10years[dbs_data$APOBEC_10years==0] = NA
dbs_data$APOBEC[dbs_data$APOBEC==0] = NA
dbs_data$APOBEC_birth[dbs_data$APOBEC_birth==0] = NA

fig <- plot_ly(type = "box")
fig <- fig %>% add_trace(y = postchemo$Chemotherapy, name='Chemotherapy (SBS31/35)',boxpoints = 'all',pointpos = 0,
                         marker = list(color = '#b2182b'),
                         line = list(color = '#b2182b'),
                         fillcolor = list(color = '#e31a1c'))
fig <- fig %>% add_trace(y = sbs_data$APOBEC, name='APOBEC (SBS2/13)',boxpoints = 'all',pointpos = 0,
                         marker = list(color = '#2166ac'),
                         line = list(color = '#377eb8'),
                         fillcolor = list(color = '#4daf4a'))
fig <- fig %>% add_trace(y = sbs_data$APOBEC_10years, name='APOBEC_10years (SBS2/13)',boxpoints = 'all',pointpos = 0,
                         marker = list(color = '#2166ac'),
                         line = list(color = '#377eb8'),
                         fillcolor = list(color = '#4daf4a'))
fig <- fig %>% add_trace(y = sbs_data$APOBEC_birth, name='APOBEC_birth (SBS2/13)',boxpoints = 'all',pointpos = 0,
                         marker = list(color = '#2166ac'),
                         line = list(color = '#377eb8'),
                         fillcolor = list(color = '#4daf4a'))
fig <- fig %>% add_trace(y = postchemo_dbs$Chemotherapy, name='Chemotherapy (DBS5)', boxpoints = 'all',pointpos = 0 ,
                         marker = list(color = '#ef8a62'),
                         line = list(color = '#ef8a62'),
                         fillcolor = list(color = '#e31a1c'))
fig <- fig %>% add_trace(y = dbs_data$APOBEC, name='APOBEC (DBS11)', boxpoints = 'all',pointpos = 0 ,
                         marker = list(color = '#67a9cf'),
                         line = list(color = '#67a9cf'),
                         fillcolor = list(color = '#4daf4a'))
fig <- fig %>% add_trace(y = dbs_data$APOBEC_10years, name='APOBEC_10years (DBS11)', boxpoints = 'all',pointpos = 0 ,
                         marker = list(color = '#67a9cf'),
                         line = list(color = '#67a9cf'),
                         fillcolor = list(color = '#4daf4a'))
fig <- fig %>% add_trace(y = dbs_data$APOBEC_birth, name='APOBEC_birth (DBS11)', boxpoints = 'all',pointpos = 0 ,
                         marker = list(color = '#67a9cf'),
                         line = list(color = '#67a9cf'),
                         fillcolor = list(color = '#4daf4a'))
fig <- layout(fig, yaxis = list(type = "log", title='Mutations per month of exposure'))

fig


APOBEC_sbs_data = sbs_data[,c('APOBEC_10years','APOBEC','APOBEC_birth')]
APOBEC_dbs_data = dbs_data[,c('APOBEC_10years','APOBEC','APOBEC_birth')]

wilcox.test(postchemo$Chemotherapy, APOBEC_sbs_data$APOBEC)$p.value
wilcox.test(postchemo$Chemotherapy, APOBEC_sbs_data$APOBEC_10years)$p.value
wilcox.test(postchemo$Chemotherapy, APOBEC_sbs_data$APOBEC_birth)$p.value

wilcox.test(postchemo_dbs$Chemotherapy, APOBEC_dbs_data$APOBEC)$p.value
wilcox.test(postchemo_dbs$Chemotherapy, APOBEC_dbs_data$APOBEC_10years)$p.value
wilcox.test(postchemo_dbs$Chemotherapy, APOBEC_dbs_data$APOBEC_birth)$p.value



######################################
#SBS92 by smoking status

metadata = read.table('cornell-bladder-metadata-final.txt',header=TRUE,as.is=TRUE,sep='\t',fill=T)
sbs_data = read.table('sbs96.HighConf.cosmicv3.2.benign.output',header=TRUE, as.is=TRUE, sep="\t")
dbs_data = read.table('dbs78.HighConf.cosmicv3.2.benign.counts.txt',header=TRUE, as.is=TRUE, sep="\t")

never = sbs_data[sbs_data$Sample %in% metadata[metadata$smoking_status == 'Never',]$tumor,]
current = sbs_data[sbs_data$Sample %in% metadata[metadata$smoking_status == 'Current',]$tumor,]
former = sbs_data[sbs_data$Sample %in% metadata[metadata$smoking_status == 'Former',]$tumor,]
fig <- plot_ly(type = "box")
fig <- fig %>% add_trace(y = never$SBS92, name='Never',boxpoints = 'all',pointpos = 0,
                         marker = list(color = '#b2182b'),
                         line = list(color = '#b2182b'),
                         fillcolor = list(color = '#e31a1c'))
fig <- fig %>% add_trace(y = former$SBS92, name='Former',boxpoints = 'all',pointpos = 0,
                         marker = list(color = 'orange'),
                         line = list(color = 'orange'),
                         fillcolor = list(color = 'orange'))
fig <- fig %>% add_trace(y = current$SBS92, name='Current',boxpoints = 'all',pointpos = 0,
                         marker = list(color = 'green'),
                         line = list(color = 'green'),
                         fillcolor = list(color = 'green'))
fig <- layout(fig, yaxis = list(title='SBS92 by Smoking Status'))

fig
wilcox.test(never$SBS92, former$SBS92)$p.value
wilcox.test(never$SBS92, current$SBS92)$p.value
wilcox.test(former$SBS92, current$SBS92)$p.value
