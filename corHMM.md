# corHMM

### 1) Installation: 

#Preparing data for CorHMM
```
#install.packages("remotes")
#install.packages("readxl")
#install.packages("geiger") 
#install.packages("foreach")
#install.packages("doParallel")
#install.packages("dplyr")
#install.packages("corHMM")

library(phytools) 
library(geiger)
library(readxl)
library(foreach)
library(doParallel)
library(dplyr)
library(corHMM)
library(ape)

setwd("/Users/danielagarciacobos/Library/CloudStorage/OneDrive-AMNH/Chapter_1/ASR")
squamate.tree <- read.tree("squamates_Title_Science2024_ultrametric_constrained.tre")
snake.data_total <- read_excel("/Users/danielagarciacobos/Library/CloudStorage/OneDrive-AMNH/Chapter_1/ASR/Total_dataset_snakes_updated.xlsx")

# Cropping the data to only obtain the column "Title_phylogeny_2023" and "Ancestral_1_state". I am also cleaning the data.frame (example, eliminating any values with NA)
snake.data <- snake.data_total[ , c("Title_phylogeny_2023", "Substrate_poly", "Reproduction_poly")]
snake.data <- as.data.frame(snake.data)

#Eliminating rows with NA's
anyNA.data.frame(snake.data)
snake.data <- na.omit(snake.data)

#changing spaces for _
snake.data$Title_phylogeny_2023 <- gsub(" ", "_", snake.data$Title_phylogeny_2023)

#Changing the row names to have the species name
rownames(snake.data) <- snake.data$Title_phylogeny_2023

#changing the sign "+" for "&" which is the sign that CorHMM reads for polymorphic states
snake.data$Substrate_poly <- gsub("\\+", "&", snake.data$Substrate_poly)
snake.data$Reproduction_poly <- gsub("\\+", "&", snake.data$Reproduction_poly)

#Make sure the substrate and reproduction data is written as a factor
str(snake.data)

snake.data$Substrate_poly <- factor(snake.data$Substrate_poly, levels = c("Aq","Aq&T","Ar","Ar&T","F","F&T","M","M&Aq","T"))
snake.data$habitats <- as.integer(snake.data$Substrate_poly)
snake.data$Reproduction_poly <- factor(snake.data$Reproduction_poly, levels = c("O","V","O&V"))
snake.data$repro_mode <- as.integer(snake.data$Reproduction_poly)

which(row.names(snake.data)=="Opheodrys_vernalis")
snake.data[414, ]

#Keep the species names as rows and the habitat and reproductive information as columns
snake.data <- snake.data[ , c("Substrate_poly", "Reproduction_poly", "habitats", "repro_mode")]
anyNA.data.frame(snake.data)
sum(is.na(snake.data))
colSums(is.na(snake.data))
which(rowSums(is.na(snake.data)) > 0)
