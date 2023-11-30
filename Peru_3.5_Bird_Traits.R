library(dplyr)

load("data/BIRDS_3_COUNTS_ALL_SiteANNUALLY-X-REPLICATE-per-Species_ARRAY.Rdata")
sp_count = rownames(ssPnt)

trait_BdLf = read.csv("data/ELEData/AVONET1_BirdLife.csv")
trait_ebrd = readxl::read_xlsx(path = "data/ELEData/AVONET2_eBird.xlsx", sheet = 2)
trait_BdTr = readxl::read_xlsx(path = "data/ELEData/AVONET3_BirdTree.xlsx", sheet = 2)

names(trait_BdLf)
trait_BdLf = select(trait_BdLf, "Species1", "Family1", "Order1", "Habitat", "Habitat.Density", "Migration", "Trophic.Level", "Trophic.Niche", "Primary.Lifestyle", "Min.Latitude", "Max.Latitude", "Centroid.Latitude", "Centroid.Longitude", "Range.Size")
names(trait_BdLf)[1:3] = c("species", "family", "order")

names(trait_ebrd)
trait_ebrd = select(trait_ebrd, "Species2", "Family2", "Order2", "Habitat", "Habitat.Density", "Migration", "Trophic.Level", "Trophic.Niche", "Primary.Lifestyle")
names(trait_ebrd)[1:3] = c("species", "family", "order")
trait_ebrd[,6] = lapply(trait_ebrd[,6], FUN = as.numeric)

names(trait_BdTr)
trait_BdTr = select(trait_BdTr, "Species3", "Family3", "Order3", "Habitat", "Habitat.Density", "Migration", "Trophic.Level", "Trophic.Niche", "Primary.Lifestyle", "Min.Latitude", "Max.Latitude", "Centroid.Latitude", "Centroid.Longitude", "Range.Size")
names(trait_BdTr)[1:3] = c("species", "family", "order")
trait_BdTr[,c(5:6, 10:14)] = lapply(trait_BdTr[,c(5:6, 10:14)], FUN = as.numeric)

traits = bind_rows(trait_BdLf, trait_BdTr, trait_ebrd)
traits$species = gsub(" ", "_", traits$species)

trait = traits %>% distinct(species, .keep_all = TRUE)
sp_trait = unique(trait$species)


length(which(!sp_count %in% sp_trait))
missing = sp_count[which(!sp_count %in% sp_trait)]


#####################################
# Compare similar spellings of species

u = sort(append(missing, sp_trait))
# compare all strings against each other
d <- adist(u)
# Do not list combinations of similar words twice
d[lower.tri(d)] <- NA
# Say your threshold below which you want to consider strings similar is 4 edits:
a <- which(d > 0 & d < 2, arr.ind = TRUE)
pairs <- cbind(u[a[,1]], u[a[,2]])
pairs
sort(missing)
## Will have to curate the results yourself to avoid accidental equalization of unequal factors.
## Do this reproducably by using a named vector as a translation dictionary. For example:
dict <- c(
  # incorrect spellings          correct spellings
  # -------------------------    ----------------------------
  "Psaracolius_angustifrons"  =  "Psarocolius_angustifrons",
  "Aratinga_weddelli"         =  "Aratinga_weddellii",
  "Mealy_parrot"              =  "Amazona_farinosa",
  "Psitticara_mitratus"       =  "Psittacara_mitratus",
  "Chaetosercus_mulsant"      =  "Chaetocercus_mulsant",
  "Psitticara_wagleri"        =  "Psittacara_wagleri",
  "Mecocerculus_poecilucercus" = "Mecocerculus_poecilocercus",
  "Pteroglossus_beauharnaesii" = "Pteroglossus_beauharnaisii",
  "Pheugopedis_genibarbis"    =  "Pheugopedius_genibarbis",
  "Turdus_hauxwellii"         =  "Turdus_hauxwelli",
  "Frederickena_unduligera"   =  "Frederickena_unduliger",
  "Fornicivora_rufa"          =  "Formicivora_rufa",
  "Hylexeastes_stresemanni"   =  "Hylexetastes_stresemanni",
  "Pyriglena_leuconata"       =  "Pyriglena_leuconota",
  "Lodopleura_isabellae"      =  "Iodopleura_isabellae",
  "Cymbilaimus_sanctaemariea" =  "Cymbilaimus_sanctaemariae",
  "Myrmotherula_myotherinus"  =  "Myrmoborus_myotherinus",
  "Sporagra_magellanica"      =  "Spinus_magellanicus",
  "Tangara_chrystosis"        =  "Tangara_chrysotis",
  "Glyphorynchus_guttatus"    =  "Xiphorhynchus_guttatus",
  "Dendrocincla_longipennis"  =  "Dendrocincla_merula"
  )

# The correct levels need to be included, to
dict <- c(dict, setNames(u,u))

# Then convert factor column to character by using as.character and apply the dictionary on the original character vector df2$name_Species:
rownames(ssPnt) <- dict[rownames(ssPnt)]

missing = rownames(ssPnt)[!rownames(ssPnt) %in% traits$species]

write.csv(trait, "data/BIRDS_3.5_Bird_TRAITS.csv", row.names = F)

