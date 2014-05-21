rm(list=ls())
##PD randomizations for each community
library(picante) #to compute pd
library(parallel)
library(raster) #to generate output rasters
#=======Functions==================

#In this function m refers to number of randomizations, n to number of species in a community, sn to the names of species in the global pool and filogenia to the phylogeny
give_complete_commununity_better= function(m=10000, n, sn,phylo) {
  #Generate community matrix
  v.names <- paste("random_community",1:m,sep="") #row names
  random_communities<-matrix(0,nrow=m,ncol=length(sn),dimnames=list(v.names,sn))
  for (k in 1:m) {
      sampling_names <- sample(sn, n)
      random_communities[k, sampling_names] <- 1
    }
  pd.result<-pd(random_communities,phylo,include.root=TRUE)
  return(pd.result)
}


#===========================================================

setwd("/home/andreap/Andrea_Filodiversidad/pd_zooregiones/zona_andes")
##Load phylogeny
filogenia<-read.nexus("consenso_andes.nex")
##Load community compostition matrix (each row represents one pixel and each column a species, the last column is the observed PD)
comunidades<-read.table("pd_aves_andes.txt",h=T,row.names=1)
species_per_pixel <- unique(apply(comunidades[,1:length(comunidades[1,]) -1],1,sum))
species_names<-names(comunidades)[1:(length(comunidades[1,])-1)]

b<-mclapply(species_per_pixel,function(x) { give_complete_commununity_better(n=x, sn= species_names,phylo=filogenia) }, mc.cores=16 )
save(b,file="aleat_aves_andes")
#Create randomization maps
  #First compare observed PD with expected PD
comunidades1<-comunidades
comunidades1$numero_sp<-apply(comunidades[,1:length(comunidades[1,]) -1],1,sum)
comunidades1$valor_p<-NA
comunidades1$valor_p_inf<-NA
comunidades1$valor_p_sup<-NA

for (i in 1:length(species_per_pixel))
{
  subset_per_count <- subset(comunidades1, comunidades1$numero_sp == species_per_pixel[i])
  match_rows<-(rownames(subset_per_count))
  valores_observados<-comunidades1[match_rows,"pd"]
  number_match<-length(valores_observados)
  valores_p_sup<-vector()
  valores_p_inf<-vector()
  for (j in 1:number_match){
    cond<-lapply(b,function(x) x[2][,1][1]==species_per_pixel[i])
    valores_p_sup[j]<-((sum(b[unlist(cond)][[1]][,1]>=valores_observados[j]))/(10000+1))*2 
    valores_p_inf[j]<-((sum(b[unlist(cond)][[1]][,1]<=valores_observados[j]))/(10000+1))*2
  }

  comunidades1[match_rows,"valor_p_inf"]<-valores_p_inf
  comunidades1[match_rows,"valor_p_sup"]<-valores_p_sup
  
}

#Chose P-value betwwwen the inferior and superior value.
a<-ifelse(comunidades1$valor_p_inf < comunidades1$valor_p_sup, comunidades1$valor_p <- comunidades1$valor_p_inf,  comunidades1$valor_p <- comunidades1$valor_p_sup)
comunidades1$valor_p<-a

#Create a column indicating whether the observed PD is Higher or Lower than expected.
comunidades1$desviacion<-NA
comunidades1$desviacion[comunidades1$valor_p==comunidades1$valor_p_sup]<-"Superior"
comunidades1$desviacion[comunidades1$valor_p==comunidades1$valor_p_inf]<-"Inferior"

#Subsetting the data frame to obtain one per category (one higher one lower)

comunidades_sup<-subset(comunidades1,desviacion=="Superior")
comunidades_inf<-subset(comunidades1,desviacion=="Inferior")

#Create rasters and plots for all significantly different than expected pixels, all signifantly higher pixels and all signifcantly lower pixels 

#1-Create an empty raster using a base to ensure correct extent, size and resolution
r<-raster("Accipiter_striatus.asc")
values(r)<-NA #eresa all values from base map
p_value_ras<-r #1
p_value_sup_ras<-r #2
p_value_inf_ras<-r #3




#2- Assign P-values to pixels
p_value_ras[as.integer(rownames(comunidades1))]<-comunidades1$valor_p#1
p_value_sup_ras[as.integer(rownames(comunidades_sup))]<-comunidades_sup$valor_p #2
p_value_inf_ras[as.integer(rownames(comunidades_inf))]<-comunidades_inf$valor_p #3

#3- Save rasters to file
writeRaster(p_value_ras,"P_value_randomization_birds_andes.asc")
writeRaster(p_value_sup_ras,"P_value_sup_randomization_birds_andes.asc")
writeRaster(p_value_inf_ras,"P_value_inf_randomization_birds_andes.asc")

#4-Plot maps in R
par(mfrow=c(1,3))
plot(p_value_ras)
plot(p_value_sup_ras)
plot(p_value_inf_ras)







