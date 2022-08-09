#######
#
### Setting up of a standardized database in community ecology and 
#
### how to make an unconstrained ordination based on a latent variable model 
#
### using boral package (FKC Hui, 2016; 2021)
#
### Author: Adrien Chevallier PhD
#
### Email: adchevallier@yahoo.fr
#
### Year: 2021

#___________________________________________________________________

rm(list=ls())

### packages

##basic
library(plyr)   # better to load plyr before dplyr to have access to "mapvalues" function (conflict between plyr and dplyr)
library(dplyr)
# very nice tuto of some dplyr functions: http://cmdlinetips.com/2018/01/6-most-useful-dplyr-verbs-to-manipulate-a-data-frame-in-r/
library(reshape)
#library(lubridate)
library(ggplot2)

##specific
library(boral)
library(fishMod)



path_data<- '~/Doctorado/Parallel/Archaeofish/R/Data'
setwd(path_data)

#abundance of invertebrates across archaeological eras
ref_tab <- read.csv("moluscos_taltal_v1.csv", header = T, sep = ",", quote = "",   
                    dec = ".", fill = T, stringsAsFactors = F) # NEW TAB !

#volumes excavatef for each archaeological era
volumes  <- read.csv("Malacologia_Taltal_m3_periodo_sitio.csv", header = T, sep = ",", quote = "",   
                     dec = ".", fill = T, stringsAsFactors = F) # NEW TAB !

#taxonomic and ecological information about each taxon
lookout <-read.csv("shell_lookup_biology_changolab_v1.csv", header = T, sep = ",", quote = "",   
                   dec = ".", fill = T, stringsAsFactors = F)

### merge abundance data and information about each taxon

#First, check for dissimilarities between two vectors
#setdiff(unique(ref_tab$Especie),unique(lookout$GENUS_SPECIES)) #ok
ref_tab <- merge(ref_tab,lookout,by.x="Especie", by.y="GENUS_SPECIES",all.x=T)

### Archaic times -> turn them to abbreviated numbers
ref_tab$Momento_Arcaico<-replace(ref_tab$Momento_Arcaico,ref_tab$Momento_Arcaico == "Arcaico 1",1)
ref_tab$Momento_Arcaico<-replace(ref_tab$Momento_Arcaico,ref_tab$Momento_Arcaico == "Arcaico 2",2)
ref_tab$Momento_Arcaico<-replace(ref_tab$Momento_Arcaico,ref_tab$Momento_Arcaico == "Arcaico 3",3)
ref_tab$Momento_Arcaico<-replace(ref_tab$Momento_Arcaico,ref_tab$Momento_Arcaico == "Arcaico 4",4)
ref_tab$Momento_Arcaico<-replace(ref_tab$Momento_Arcaico,ref_tab$Momento_Arcaico == "Arcaico 4-5","4-5")
ref_tab$Momento_Arcaico<-replace(ref_tab$Momento_Arcaico,ref_tab$Momento_Arcaico == "Arcaico 5",5)
ref_tab$Momento_Arcaico<-replace(ref_tab$Momento_Arcaico,ref_tab$Momento_Arcaico == "Arcaico 6",6)

###Abundance Table

#Sum MNI (abundance) by species for each Archaic moment
ecology <- aggregate(MNI ~ Especie + Momento_Arcaico, ref_tab, sum)
### long to wide format 
#Momento_Arcaico long #Especie wide #and MNI as the values
ecology <- reshape(ecology, idvar = "Momento_Arcaico", timevar = "Especie", direction = "wide")
colnames(ecology) <- gsub( "MNI.", "", as.character(colnames(ecology)))   #Erase some artefact characters in a row or a column
row.names(ecology)<-ecology$Momento_Arcaico
ecology <- ecology[, !(colnames(ecology) %in% c("Momento_Arcaico"))]
ecology[is.na(ecology)] <- 0
#sort #order
ecology <- ecology[,order(colnames(ecology),decreasing = T)]
ecology <- ecology[order(row.names(ecology)),]
#remove moment 4-5 #useless
ecology <- ecology[-5,]

### standardization of MNI (Minimal Number of Individuals) by excavated volume (m3)
volumes <- volumes[,2:7]
archaic_moment <- c("1", "2", "3", "4", "5", "6")
volume <- c(colSums(volumes, na.rm=TRUE))
volumes <- data.frame(archaic_moment, volume)
main_species <- round(ecology/volumes$volume,1)  
#remove rare taxa
totalMNI <- as.data.frame(c(colSums(main_species, na.rm=TRUE)))
totalMNI$species <- rownames(totalMNI)
hist(totalMNI$`c(colSums(main_species, na.rm = TRUE))`) # there are some very abundant taxa
commonMNI <- totalMNI[which(totalMNI$`c(colSums(main_species, na.rm = TRUE))` >= 10), ]
#tab with only rare taxa
rare_species <- main_species[, -which(colnames(main_species) %in% commonMNI$species)] #>= 10 MNI/m3
#tab with common taxa
main_species <- main_species[, colnames(main_species) %in% commonMNI$species] #>= 10 MNI/m3


### LVM modeling

system.time({ 
  main_species_tweedie <- boral(main_species, family = "tweedie", lv.control = list(num.lv = 2),
                                   row.eff = "random", save.model=TRUE)
})


## Residuals analysis

par(mfrow = c(2,2))
plot.boral(main_species_tweedie)
# residuals ok
par(mfrow = c(1,1))

## Model Plot (sites and species)
getplot <- lvsplot(main_species_tweedie, return.vals = TRUE)




#___________________________________________________________________

# © 2021 Adrien Chevallier

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
