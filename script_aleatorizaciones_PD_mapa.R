rm(list=ls())
##PD randomizations for each community
library(picante) #to compute pd
library(parallel) #this code uses parallel computing
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
phylogeny<-read.nexus("consenso_andes.nex")
##Load community compostition matrix (each row represents one pixel and each column a species, the last column is the observed PD)
##The first colum must contain the pixel numbers, these will become the rownames
communities<-read.table("pd_aves_andes.txt",h=T,row.names=1)
species_per_pixel <- unique(apply(communities[,1:length(communities[1,]) -1],1,sum))
species_names<-names(communities)[1:(length(communities[1,])-1)]

b<-mclapply(species_per_pixel,function(x) { give_complete_commununity_better(n=x, sn= species_names,phylo=phylogeny) }, mc.cores=16 )
save(b,file="aleat_aves_andes")
#Create randomization maps
  #First compare observed PD with expected PD
communities1<-communities
communities1$sp_number<-apply(comunidades[,1:length(comunidades[1,]) -1],1,sum)
communities1$p_values<-NA
communities1$p_values_lower<-NA
communities1$p_values_higher<-NA

for (i in 1:length(species_per_pixel))
{
  subset_per_count <- subset(communities1, communities1$sp_number == species_per_pixel[i])
  match_rows<-(rownames(subset_per_count))
  observed_values<-communities1[match_rows,"pd"]
  number_match<-length(observed_values)
  p_values_higher<-vector()
  valores_p_inf<-vector()
  for (j in 1:number_match){
    cond<-lapply(b,function(x) x[2][,1][1]==species_per_pixel[i])
    p_values_higher[j]<-((sum(b[unlist(cond)][[1]][,1]>=observed_values[j]))/(10000+1))*2 
    p_values_lower[j]<-((sum(b[unlist(cond)][[1]][,1]<=observed_values[j]))/(10000+1))*2
  }

  communities1[match_rows,"p_values_lower"]<-p_values_lower
  communities1[match_rows,"p_values_higher"]<-p_values_higher
  
}

#Chose P-value betwwwen the inferior and superior value.
a<-ifelse(communities1$p_values_lower < communities1$p_values_higher, communities1$p_values <- communities1$p_values_lower,  communities1$p_values <- communities1$p_values_higher)
communities1$p_values<-a

#Create a column indicating whether the observed PD is Higher or Lower than expected.
communities1$desviacion<-NA
communities1$desviacion[communities1$p_values==communities1$p_values_higher]<-"Superior"
communities1$desviacion[communities1$p_values==communities1$p_values_lower]<-"Inferior"

#Subsetting the data frame to obtain one per category (one higher one lower)

comunidades_sup<-subset(communities1,desviacion=="Superior")
comunidades_inf<-subset(communities1,desviacion=="Inferior")

#Create rasters and plots for all significantly different than expected pixels, all signifantly higher pixels and all signifcantly lower pixels 

#1-Create an empty raster using a base to ensure correct extent, size and resolution
r<-raster("Accipiter_striatus.asc")
values(r)<-NA #eresa all values from base map
p_value_ras<-r #1
p_value_sup_ras<-r #2
p_value_inf_ras<-r #3




#2- Assign P-values to pixels
p_value_ras[as.integer(rownames(communities1))]<-communities1$valor_p#1
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

#5 Reclassify rasters to only show significant pixels
  ##assuming alfa of 0.05
  ##Create reclassifying matrix. Pixels significantly higher will have a value of 1 and pixels significantly lower a value of -1
  m<-c(0,0.05,1,0.05,1,NA)
  rclmat<-matrix(m,ncol=3,byrow=TRUE)
  p_value_sup<-reclassify(p_value_sup_ras, rclmat)
   m<-c(0,0.05,-1,0.05,1,NA)
  rclmat<-matrix(m,ncol=3,byrow=TRUE)
  p_value_inf<-reclassify(p_value_inf_ras, rclmat)
  
 significant_pixels<-mosaic(p_value_sup,p_value_inf, fun=mean)
writeRaster(significant_pixels,"Significant_pixels.asc")
  
  






