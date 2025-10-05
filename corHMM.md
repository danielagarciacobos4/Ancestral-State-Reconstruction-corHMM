# corHMM

### 1) Preparing data for CorHMM

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

#I will eliminate the columns containing NA's to maps the states of reproduction and habitat in the tip of the tree. However, I will keep these unknown states for the CorHMM analysis 

#snake.data_poly_tree <- na.omit(snake.data_poly)

```
### 3) Tree plotting habitat and reproductive states at tips

```
###############################################################################
# PIPELINE: Match tree ↔ data, build pie-slices for polymorphic states, plot
# INPUTS assumed:
#   pruned_tree_poly      : an 'ape' phylo object already pruned to your taxa
#   corhmm_data_poly      : data.frame with at least columns Species,
#                           Substrate_poly, Reproduction_poly
# OUTPUT:
#   final_curated_habitat_repro_tree.pdf (tree with two pies per tip)
#   NOTE: Our original matrix has 2092 species in the tree, however 20 of these are NA's. We are eliminating those NA's in this tree that helps us visualize the states at the tip labels. But we will use the NAs in the models. 
###############################################################################

##Prune the tree

#Sample 500 species
set.seed(123)
sampled_snakes_habitat <- snake.data %>% sample_n(500)
sampled_species_names_habitat <- rownames(sampled_snakes_habitat)

# Prune the tree to only include sampled species
pruned_tree_habitat <- drop.tip(squamate.tree, setdiff(squamate.tree$tip.label, sampled_species_names_habitat))

# Optional: Plot to verify
plot(pruned_tree_habitat, cex = 0.5)

# Extraer los dos traits
substrate_data_habitat <- sampled_snakes_habitat$Substrate_poly
reproductive_data_habitat <- sampled_snakes_habitat$Reproduction_poly

# Asegurar que los nombres estén bien asignados
names(substrate_data_habitat) <- rownames(sampled_snakes_habitat)
names(reproductive_data_habitat) <- rownames(sampled_snakes_habitat)

# Crear un data.frame compatible con 
sampled_snakes_habitat$Species <- rownames(sampled_snakes_habitat)
rownames(sampled_snakes_habitat) <- NULL

corhmm_data_habitat <- sampled_snakes_poly[, c("Species", "Substrate_poly", "Reproduction_poly")]
corhmm_data_poly

head(corhmm_data_habitat)
str(corhmm_data_habitat)


## 1) Normalize species names in both sources (tree & table) -------------------
#    trimws() removes stray leading/trailing spaces that break matching.
tips_tree <- trimws(pruned_tree_habitat$tip.label)
sp_table  <- trimws(corhmm_data_habitat$Species)


## 2) Match species sets and diagnose mismatches -------------------------------
#    - matched: taxa present in BOTH tree and table
#    - missing_in_table: tips present in the tree but absent from the table
#    - missing_in_tree : taxa present in the table but absent from the tree

matched <- intersect(tips_tree, sp_table)
missing_in_table <- setdiff(tips_tree, sp_table)
missing_in_tree  <- setdiff(sp_table, tips_tree)

# Summary counts and a peek at the first few mismatches (if any)
message("Tips in tree: ", length(tips_tree))
message("Species in table: ", length(sp_table))
message("Matched (tree ∩ table): ", length(matched))
if (length(missing_in_table)) {
  message("First tips missing from table: ", paste(head(missing_in_table, 10), collapse=", "))
}
if (length(missing_in_tree)) {
  message("First species missing from tree: ", paste(head(missing_in_tree, 10), collapse=", "))
}

## 3) Keep only the common set to avoid zero pies / misalignment ---------------
#    Drop any tree tips that are not in the matched set.
tree_use <- drop.tip(pruned_tree_poly, setdiff(tips_tree, matched))

## 4) Subset the data table to only matched species ----------------------------
dat_use <- corhmm_data_poly[corhmm_data_poly$Species %in% matched, ]

# Quick lookups: named vectors keyed by Species
x_hab <- setNames(as.character(dat_use$Substrate_poly),    dat_use$Species)
x_rep <- setNames(as.character(dat_use$Reproduction_poly), dat_use$Species)

## 5) Rebuild matrices in the *tree tip order* ---------------------------------
#    This guarantees rows of the pie matrices align to tree_use$tip.label.
tips <- tree_use$tip.label

# Helper: convert a polymorphic state string (e.g., "M&Aq" or "V+O")
#   into an equal-weights vector over the allowed levels.
to_weights <- function(state_string, levels) {
  if (is.na(state_string) || state_string == "") return(rep(0, length(levels)))
  toks <- unlist(strsplit(state_string, "[+&]"))
  toks <- trimws(toks)
  toks <- toks[toks %in% levels]
  w <- rep(0, length(levels))
  if (length(toks)) w[match(toks, levels)] <- 1/length(toks)
  w
}

# Build the pie-weight matrices (rows = tips, cols = state levels)
hab_levels <- c("T","F","Ar","Aq","M")
rep_levels <- c("O","V")

hab_pies <- t(vapply(tips, function(tp) to_weights(x_hab[tp], hab_levels),
                     FUN.VALUE = setNames(numeric(length(hab_levels)), hab_levels)))
rep_pies <- t(vapply(tips, function(tp) to_weights(x_rep[tp], rep_levels),
                     FUN.VALUE = setNames(numeric(length(rep_levels)), rep_levels)))

## 6) Clean any stray NAs/Inf produced above -----------------------------------
hab_pies[!is.finite(hab_pies)] <- 0
rep_pies[!is.finite(rep_pies)] <- 0

## 7) Sanity checks: how many tips will actually draw a slice? ------------------
message("Tips with any HABITAT slice >0: ", sum(rowSums(hab_pies) > 0))
message("Tips with any REPRO slice  >0: ", sum(rowSums(rep_pies) > 0))

## 8) Colors aligned to column order (bullet-proof against reordering) ----------
hab_cols <- c("T"="#FF9933","F"="#9966CC","Ar"="#66CC66","Aq"="#3399FF","M"="#003366")
rep_cols <- c("O"="#FFD700","V"="#E41A1C")
hab_cols <- hab_cols[colnames(hab_pies)]
rep_cols <- rep_cols[colnames(rep_pies)]

## 9) Plot: vertical phylogram with two pies per tip ----------------------------
#    - First pie (offset ~0.02) = Habitat
#    - Second pie (offset ~0.8)  = Reproduction
#    Adjust 'height' and 'width' to your tree size; margins give space for labels.

pdf(file = "final_curated_habitat_repro_tree.pdf", height = 180, width = 15)
par(mar = c(2, 2, 2, 10))  # wider right margin for pies+labels

ape::plot.phylo(
  tree_use,
  type = "phylogram",
  direction = "rightwards",
  no.margin = TRUE,
  show.tip.label = TRUE,
  cex = 0.35,
  label.offset = 1.2
)

# IMPORTANT: Draw both sets of pies on the same plotted tree (no replot in between)
tiplabels(pie = hab_pies, piecol = unname(hab_cols), cex = 0.12, offset = 0.02)
tiplabels(pie = rep_pies, piecol = unname(rep_cols), cex = 0.12, offset = 0.8)

legend("topleft",  title = "Habitat",
       legend = names(hab_cols), pt.cex = 1.3, pch = 21, pt.bg = unname(hab_cols), bty = "n")
legend("topright", title = "Reproduction",
       legend = names(rep_cols), pt.cex = 1.3, pch = 21, pt.bg = unname(rep_cols), bty = "n")

dev.off()
```
## 3.1) Prepare corHMM-formatted data only habitat

```{r}

##Prune the tree

#Sample 500 species
set.seed(123)
sampled_snakes_habitat <- snake.data #%>% sample_n(400)
sampled_species_names_habitat <- rownames(sampled_snakes_habitat)

# Prune the tree to only include sampled species
pruned_tree_habitat <- drop.tip(squamate.tree, setdiff(squamate.tree$tip.label, sampled_species_names_habitat))

# Optional: Plot to verify
plot(pruned_tree_habitat, cex = 0.5)

# Extraer los dos traits
substrate_data_habitat <- sampled_snakes_habitat$Substrate_poly

# Asegurar que los nombres estén bien asignados
names(substrate_data_habitat) <- rownames(sampled_snakes_habitat)

# Crear un data.frame compatible con 
sampled_snakes_habitat$Species <- rownames(sampled_snakes_habitat)
rownames(sampled_snakes_habitat) <- NULL

corhmm_data_habitat <- sampled_snakes_habitat[, c("Species", "Substrate_poly")]
corhmm_data_habitat

head(corhmm_data_habitat)
str(corhmm_data_habitat)

anyNA.data.frame(corhmm_data_habitat)

cols <- c("Substrate_poly")

corhmm_data_habitat[cols] <- lapply(corhmm_data_habitat[cols], function(x) {
  x <- trimws(as.character(x))   # ensure character & strip spaces
  x[is.na(x) | x == ""] <- "?"   # turn NA or empty strings into "?"
  x
})

```
## 3.2) Running models: only habitat
```{r,}

fit_ER_1.1.h   <- corHMM(pruned_tree_habitat, corhmm_data_habitat, rate.cat = 1, model = "ER", node.states = "marginal",
                     nstarts = 48, n.cores = 48, get.tip.states = FALSE)
fit_ER_1.2.h   <- corHMM(pruned_tree_habitat, corhmm_data_habitat, rate.cat = 1, model = "ER", node.states = "marginal",
                     nstarts = 48, n.cores = 48, get.tip.states = TRUE)

fit_ER_2.1.h   <- corHMM(pruned_tree_habitat, corhmm_data_habitat, rate.cat = 2, model = "ER", node.states = "marginal", 
                     nstarts = 48, n.cores = 48, get.tip.states = TRUE)
fit_ER_2.2.h   <- corHMM(pruned_tree_habitat, corhmm_data_habitat, rate.cat = 2, model = "ER", node.states = "marginal", 
                     nstarts = 48, n.cores = 48, get.tip.states = TRUE)

fit_ER_3.1.h   <- corHMM(pruned_tree_habitat, corhmm_data_habitat, rate.cat = 3, model = "ER", node.states = "marginal", 
                     nstarts = 48, n.cores = 48, get.tip.states = TRUE)
fit_ER_3.2.h   <- corHMM(pruned_tree_habitat, corhmm_data_habitat, rate.cat = 3, model = "ER", node.states = "marginal", 
                     nstarts = 48, n.cores = 48, get.tip.states = TRUE)

fit_ER_4.1.h   <- corHMM(pruned_tree_habitat, corhmm_data_habitat, rate.cat = 4, model = "ER", node.states = "marginal", 
                     nstarts = 48, n.cores = 48, get.tip.states = TRUE)
fit_ER_4.2.h   <- corHMM(pruned_tree_habitat, corhmm_data_habitat, rate.cat = 4, model = "ER", node.states = "marginal", 
                     nstarts = 48, n.cores = 48, get.tip.states = TRUE)


fit_SYM_1.1.h   <- corHMM(pruned_tree_habitat, corhmm_data_habitat, rate.cat = 1, model = "SYM", node.states = "marginal",
                     nstarts = 48, n.cores = 48, get.tip.states = TRUE)
fit_SYM_1.2.h   <- corHMM(pruned_tree_habitat, corhmm_data_habitat, rate.cat = 1, model = "SYM", node.states = "marginal",
                     nstarts = 48, n.cores = 48, get.tip.states = TRUE)

fit_SYM_2.1.h   <- corHMM(pruned_tree_habitat, corhmm_data_habitat, rate.cat = 2, model = "SYM", node.states = "marginal",
                     nstarts = 48, n.cores = 48, get.tip.states = TRUE)
fit_SYM_2.2.h   <- corHMM(pruned_tree_habitat, corhmm_data_habitat, rate.cat = 2, model = "SYM", node.states = "marginal",
                     nstarts = 48, n.cores = 48, get.tip.states = TRUE)


fit_SYM_3.1.h   <- corHMM(pruned_tree_habitat, corhmm_data_habitat, rate.cat = 3, model = "SYM", node.states = "marginal",
                     nstarts = 48, n.cores = 48, get.tip.states = TRUE)


fit_SYM_4.1.h   <- corHMM(pruned_tree_habitat, corhmm_data_habitat, rate.cat = 4, model = "SYM", node.states = "marginal",
                     nstarts = 48, n.cores = 48, get.tip.states = TRUE)
fit_SYM_4.2.h   <- corHMM(pruned_tree_habitat, corhmm_data_habitat, rate.cat = 4, model = "SYM", node.states = "marginal",
                     nstarts = 48, n.cores = 48, get.tip.states = TRUE)


fit_ARD_1.1.h <- corHMM(pruned_tree_habitat, corhmm_data_habitat, rate.cat = 1, model = "ARD", node.states = "marginal",
                     nstarts = 5, n.cores = 5, get.tip.states = TRUE)
fit_ARD_1.2.h   <- corHMM(pruned_tree_habitat, corhmm_data_habitat, rate.cat = 1, model = "ARD", node.states = "marginal",
                     nstarts = 50, n.cores = 25, get.tip.states = TRUE)


fit_ARD_2.1.h   <- corHMM(pruned_tree_habitat, corhmm_data_habitat, rate.cat = 2, model = "ARD", node.states = "marginal",
                     nstarts = 50, n.cores = 25, get.tip.states = TRUE)
fit_ARD_2.2.h   <- corHMM(pruned_tree_habitat, corhmm_data_habitat, rate.cat = 2, model = "ARD", node.states = "marginal",
                     nstarts = 50, n.cores = 25, get.tip.states = TRUE)


fit_ARD_3.1.h   <- corHMM(pruned_tree_habitat, corhmm_data_habitat, rate.cat = 3, model = "ARD", node.states = "marginal",
                     nstarts = 50, n.cores = 25, get.tip.states = TRUE)
fit_ARD_3.2.h   <- corHMM(pruned_tree_habitat, corhmm_data_habitat, rate.cat = 3, model = "ARD", node.states = "marginal",
                     nstarts = 50, n.cores = 25, get.tip.states = TRUE)


fit_ARD_4.1.h   <- corHMM(pruned_tree_habitat, corhmm_data_habitat, rate.cat = 4, model = "ARD", node.states = "none",
                     nstarts = 5, n.cores = 5, get.tip.states = TRUE)
fit_ARD_4.2.h   <- corHMM(pruned_tree_habitat, corhmm_data_habitat, rate.cat = 4, model = "ARD", node.states = "marginal",
                     nstarts = 50, n.cores = 25, get.tip.states = TRUE)

```

## 3.3) choosing the best model habitat only

```
fit_results_habitat <- list(
 ER_1 = fit_ER_1.1.h, ER_2 = fit_ER_2.1.h, ER_3 = fit_ER_3.1.h, ER_4 = fit_ER_4.1.h,
 SYM_1 = fit_SYM_1.1.h, SYM_2 = fit_SYM_2.1.h, SYM_3 = fit_SYM_3.1.h, SYM_4 = fit_SYM_4.1.h,
 ARD_1 = fit_ARD_1.1.h, ARD_2 = fit_ARD_2.1.h, ARD_3 = fit_ARD_3.1.h, ARD_4 = fit_ARD_4.1.h)

include_models_habitat <- c(ER_1 = TRUE, ER_2 = TRUE, ER_3 = TRUE, ER_4 = TRUE, SYM_1 = TRUE, SYM_2 = TRUE, SYM_3 = TRUE, SYM_4 = TRUE, ARD_1 = TRUE, ARD_2 = TRUE, ARD_3 = TRUE, ARD_4 = TRUE)


included_fits_habitat <- fit_results_habitat#[include_models]

aic_vals_habitat <- sapply(included_fits_habitat, function(x) x$AICc)
log_liks_habitat <- sapply(included_fits_habitat, function(x) x$loglik)
min_aic_habitat <- min(aic_vals_habitat)
delta_aic_haitat <- aic_vals_habitat - min_aic_habitat
aic_weights_habitat <- exp(-0.5 * delta_aic_haitat) / sum(exp(-0.5 * delta_aic_haitat))

model_comparison_habitat <- data.frame(
  Model = names(included_fits_habitat),
  LogLik = log_liks_habitat,
  AICc = aic_vals_habitat,
  deltaAICc = delta_aic_haitat,
  AICw = round(aic_weights_habitat, 4)
)

model_comparison_habitat <- model_comparison_habitat[order(model_comparison_habitat$AICc), ]
cat("\n=== Model comparison — habitat-only ===\n")
print(model_comparison_habitat)

# Select best model and visualize
best_model_habitat_name <- model_comparison_habitat$Model[1]
best_fit_habitat <- included_fits_habitat[[best_model_habitat_name]]
print(best_fit_habitat)

```

<img width="893" height="288" alt="Screenshot 2025-10-05 at 7 11 55 PM" src="https://github.com/user-attachments/assets/3bf90ef8-d6b4-4a14-a06a-c45dcb95d48d" />



