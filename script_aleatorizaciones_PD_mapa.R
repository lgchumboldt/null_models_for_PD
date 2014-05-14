rm(list=ls())
##Aleatorizaciones de PD en cada comunidad.
library(picante) #para calcular PD
#library(multicore)
library(parallel)
#=======Functions==================


give_complete_commununity_better= function(m=10000, n, sn,filogenia) {
  #Generar matriz de comunidades
  v.nombres <- paste("comunidad_aleatoria",1:m,sep="") #nombre de filas
  comunidades_aleatorias<-matrix(0,nrow=m,ncol=length(sn),dimnames=list(v.nombres,sn))
  for (k in 1:m) {
      sampling_names <- sample(sn, n)
      comunidades_aleatorias[k, sampling_names] <- 1
    }
  pd.result<-pd(comunidades_aleatorias,filogenia,include.root=TRUE)
  return(pd.result)
}


#===========================================================

setwd("/home/andreap/Andrea_Filodiversidad/pd_zooregiones/zona_andes")
##Cargar la filogenia
filogenia<-read.nexus("consenso_andes.nex")
##Cargar matriz con composición de comunidades
comunidades<-read.table("pd_aves_andes.txt",h=T,row.names=1)
species_per_pixel <- unique(apply(comunidades[,1:length(comunidades[1,]) -1],1,sum))
species_names<-names(comunidades)[1:(length(comunidades[1,])-1)]

b<-mclapply(species_per_pixel,function(x) { give_complete_commununity_better(n=x, sn= species_names,filogenia=filogenia) }, mc.cores=16 )
save(b,file="aleat_aves_andes")
#Hacer mapas de aleatorizaciones
  #Primero comparar
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

#Escoger entre inferior y superior cuál es el valor P.
a<-ifelse(comunidades1$valor_p_inf < comunidades1$valor_p_sup, comunidades1$valor_p <- comunidades1$valor_p_inf,  comunidades1$valor_p <- comunidades1$valor_p_sup)
comunidades1$valor_p<-a

#Crear una columna que indique si el PD es superior o inferior a lo esperado.
comunidades1$desviacion<-NA
comunidades1$desviacion[comunidades1$valor_p==comunidades1$valor_p_sup]<-"Superior"
comunidades1$desviacion[comunidades1$valor_p==comunidades1$valor_p_inf]<-"Inferior"

#Crear marcos de datos que contengan únicamente (1) Valor-P para valores superiores al esperado (2) Valor-P para valores inferiores al esperado

comunidades_sup<-subset(comunidades1,desviacion=="Superior")
comunidades_inf<-subset(comunidades1,desviacion=="Inferior")
#Crear rasters y plot para (1) Todos los pixeles significativos (2) Todos los pixeles significativamente superiores (3) Todos los pixeles significativamente inferiores

#1-Crear los rasters vacios usando una base para tener resolucion y tamaño correctos

library(raster)
r<-raster("Accipiter_striatus.asc")
values(r)<-NA #se eliminan todos los valores del mapa de base
p_value_ras<-r #1
p_value_sup_ras<-r #2
p_value_inf_ras<-r #3




#2- Asignar al raster los valores de PD que corresponden a cada pixel
p_value_ras[as.integer(rownames(comunidades1))]<-comunidades1$valor_p#1
p_value_sup_ras[as.integer(rownames(comunidades_sup))]<-comunidades_sup$valor_p #2
p_value_inf_ras[as.integer(rownames(comunidades_inf))]<-comunidades_inf$valor_p #3

#3- Guardar el raster en un archivo 
writeRaster(p_value_ras,"P_value_randomization_birds_andes.asc")
writeRaster(p_value_sup_ras,"P_value_sup_randomization_birds_andes.asc")
writeRaster(p_value_inf_ras,"P_value_inf_randomization_birds_andes.asc")

#4-Opcional ver el mapa en R
par(mfrow=c(1,3))
plot(p_value_ras)
plot(p_value_sup_ras)
plot(p_value_inf_ras)







