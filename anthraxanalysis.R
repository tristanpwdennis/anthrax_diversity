library(tidyverse)
library(reshape2)
library(sf)
library("rnaturalearth")
library("rnaturalearthdata")
library(RColorBrewer)
library("ggrepel")
library(ggtree)
library(treeio)
library(phytools)
library(svglite)
library(cowplot)
library(ggspatial)
library(viridis)
library(rgeos)
library(data.table)


setwd('~/Projects/anthrax_diversity/anthrax_diversity/')
#read in data
metadata <- read.csv('data/metadata.csv', stringsAsFactors=FALSE, na.strings = c("", "NA"))
#read in distance matrix and convert to pairwise obcservations vs matrix
mat <- read.delim('data/MSA_Banthracis_75samples_withGapsAndNs.distmatrix.txt', check.names=FALSE)

#force lat long to numeric
metadata$long <- as.numeric(metadata$long)
metadata$lat <- as.numeric(metadata$lat)


#############################################################################################
### Turning distance matrix into within/between scale comparisons for plotting and
### analysis
#############################################################################################

#################
#create indexed observations from matrix top triangle
test <- as.matrix(mat, labels=T)
#define rownames
rownames(test) <- mat[,1]
#drop the column that became rownames
test <- test[,-1]
#generate every combination of samples (in the colnames). paired
xy <- t(combn(colnames(test), 2))
#create dataframe by retrieving every value for each pair
inddist <- data.frame(xy, dist=test[xy])
#coerce snp distances to numeric
inddist$dist <- as.numeric(inddist$dist)

#################
#rename columns and join metadata to lhs and rhs
t1 <- inddist %>% 
  rename(idx = `X1`, idy = X2, ntdiff = dist) %>% 
  left_join(., metadata, by = c("idx" = "sample_id")) %>% 
  left_join(., metadata, by = c("idy" = "sample_id"))

#################
#Create four columns - one for each scale
t2 <-
  t1 %>% 
  select(idx, idy, ntdiff, carcass.x, carcass.y, Cluster.x, Cluster.y, species.x, species.y, geog.x, geog.y, lat.x, long.x, lat.y, long.y) %>% 
  mutate(., epi = 'epi') %>%
  mutate(., geog = 'geog') %>% 
  mutate(., carcass = 'carcass') %>% 
  mutate(., species = 'species') 

#################
#Stack the four scale columns into one, duplicating the df four times
#then define the different scales of each comparison based on whether each level (carcass, epi, geog group and species) are the same or 
#different
t3 <- t2 %>% 
  pivot_longer(c(carcass, geog, epi, species), names_to = "ascale") %>% 
  mutate(withinorbetween = case_when(
    (ascale == 'carcass' & carcass.x == carcass.y) ~ 'Same Carcass',
    (ascale == 'carcass' & carcass.x != carcass.y) ~ 'Different Carcasses',
    (ascale == 'geog' & geog.x == geog.y) ~ 'Same Cluster', 
    (ascale == 'geog' & geog.x != geog.y) ~ 'Different Cluster',
    (ascale == 'species' & species.x == species.y) ~ 'Same Species', 
    (ascale == 'species' & species.x != species.y) ~ 'Different Species',
    (ascale == 'epi' & Cluster.x == Cluster.y & Cluster.x != 'Single case' & Cluster.y != 'Single Case') ~ 'Linked Cases',
    (ascale == 'epi' & Cluster.x != Cluster.y | Cluster.x == 'Single case' | Cluster.y == 'Single case') ~ 'dUnlinked Cases'))

################
#filter out same:same observations
filt_t3 <-t3 %>% filter(!(idx==idy))
#remove observations from epi clusters that are from the same carcass
filt_t4 <- filt_t3 %>% filter(!(ascale == 'epi' & carcass.x == carcass.y))


lu <- filt_t4 %>%   filter(., ascale == 'Epidemiological Cluster' & withinorbetween == "Linked Cases")

################
#generate summary stats for each group
summarystatsfiltt4 <-filt_t4 %>% 
  select(withinorbetween, ntdiff) %>% 
  group_by(withinorbetween) %>% 
  summarise(median(ntdiff), mean(ntdiff), min(ntdiff), max(ntdiff) )

################
#write to file - unfiltered df and filtered df, and summary stats
write_csv(t1, 'data/pairwise-observations-unfiltered.csv')
write_csv(filt_t4, 'data/pairwise-observations-filtered.csv')
write_csv(summarystatsfiltt4, 'data/summarystatstfromfiltered.csv')

#############################################################################################
### Prepare final df for plotting
### 
#############################################################################################

#define ascale as factor level for plotting (capitalising letters, etc)
filt_t4$ascale <- factor(filt_t4$ascale, levels = c("geog", "epi", "carcass", "species"), labels = c("Geographic Cluster",  "Epidemiological Cluster", "Carcass", "Species"))

#check withinorbetween factor elvels
levels(as.factor(filt_t4$withinorbetween))
#define palette based on order of above
colorset = c('#33a02c','#33a02c','#33a02c', '#ff7f00', '#ff7f00', '#ff7f00')

#########
#violin plots of distance within and between scales
filt_t4 %>% filter((ascale != 'Species')) %>% 
  ggplot(., aes(x=withinorbetween, y=ntdiff, fill = withinorbetween)) + 
  scale_fill_manual(values=colorset) +
  geom_violin(trim=FALSE, width = 1, show.legend = FALSE ) +
  facet_grid(. ~ ascale, scales = "free_x") +
  geom_boxplot(width=0.04, color="#252525", alpha=1.5, show.legend = FALSE ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  xlab("Difference in Scale") +
  ylab("Pairwise distance (# nucleotides)") +
  ylim(0, 80) +
  theme_minimal(base_size = 22) 

ggsave("violinfig1_cleaned_amended", plot = last_plot(), device = 'png', path = "figures/")

#########
#histogram of disdtance within carcass
filt_t4 %>% filter(ascale == 'Carcass' & withinorbetween == 'Same Carcass') %>% 
  ggplot(aes(x=ntdiff)) +
  geom_histogram(binsize =1) +
  theme_minimal()

#########
#histogram of disdtance within NCA whole
filt_t4 %>% 
  select(idx, idy, ntdiff) %>% 
  distinct() %>% 
  ggplot(aes(x=ntdiff)) +
  geom_histogram(binsize =1) +
  theme_minimal()

#violin plots of distance within and between scales FOR SPECIES ONLY
filt_t4 %>% filter((ascale == 'Species')) %>% 
  ggplot(., aes(x=withinorbetween, y=ntdiff, fill = withinorbetween)) + 
  scale_fill_manual(values=c('#33a02c','#ff7f00')) +
  geom_violin(trim=FALSE, width = 1, show.legend = FALSE ) +
  facet_grid(. ~ ascale, scales = "free_x") +
  geom_boxplot(width=0.04, color="#252525", alpha=1.5, show.legend = FALSE ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  xlab("Difference in Scale") +
  ylab("Pairwise distance (# nucleotides)") +
  ylim(0, 80) +
  theme_minimal(base_size = 22) 

#############################################################################################
###IS PAIRWISE DISTANCE A FUNCTION OF SCALE?
###i.e. do sequences tend to be more divergent the further away they are?
#############################################################################################

#get pairwise lat and long for pairs of samples
#remove imcomplete rows
lalo <- filt_t4 %>% 
  select(idx, idy, ntdiff, lat.x, long.x, lat.y, long.y) %>% 
  distinct()
c<-lalo[complete.cases(lalo), ] 

#get distances between sampling points
for(i in 1:nrow(c)) {
t<-raster::pointDistance(c(c$lat.x[i], c$long.x[i]), c(c$lat.y[i], c$long.y[i]), lonlat = T)
c[i, 8] <- t
}

#plot distance over ntdiff
ggplot(c, aes(x = `...8`, y = ntdiff)) +
  xlab("Distance in metres") +
  ylab("Distance (nucleotides)") +
  geom_point() +
  theme_minimal()

#############################################################################################
###Plot on map
#############################################################################################
#read in Tom and Grant's shapefiles for the Serengeti region
aoi <- st_read('data/V4_Serengeti_Ecosystem/v4_serengeti_ecosystem.shp', )

#get slimmed down metadata
mapmetadata <- read.csv('data/metadata-slim.csv')

#specify NCA geom coordinates
Ngorongoro = aoi[aoi$NAME == "Ngorongoro",]

#force lat long to numeric
mapmetadata$long <- as.numeric(as.character(mapmetadata$long))
mapmetadata$lat <- as.numeric(as.character(mapmetadata$lat))

#download earth data from rnaturalrearth
world <- ne_countries(scale = 'medium', returnclass = "sf", continent = 'africa')
lakes <- ne_download(scale = 'medium', type = 'lakes', category = 'physical', returnclass = "sf")

#crop the world map to contain NE Tz and SW Kenya
eastafricarop <- st_crop(world,
                         xmin = 34.5,
                         xmax = 36.5,
                         ymin = -1,
                         ymax = -4)

eastafricabbox <- st_bbox(eastafricarop)




#plot eas africa
eastafrica <- ggplot(data = eastafricarop) +
  geom_sf(fill = "#fcfcf5") +
  geom_sf(data = Ngorongoro, fill = NA, color = "red", size = 1) +
  #geom_sf(data = lakes, colour = '#a6bddb', fill = '#a6bddb')
  theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "aliceblue"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) 

eastafrica


#plot nca with points
nca <- ggplot(data = world) +
  geom_sf(fill = "#fcfcf5") +
  geom_sf(data = Ngorongoro, fill = NA) +
  geom_point(data = mapmetadata, aes(x = mapmetadata$long, y =mapmetadata$lat, color = Cluster), size = 1.5) +
  scale_color_manual(values=c("#b3b3b3", "#d7191c", "#fdae61", "#abd9e9", "#2c7bb6")) +
  geom_text_repel(data = subset(mapmetadata, Cluster == "C1" | Cluster == "C2" | Cluster == "C3" | Cluster == "C4"), aes(x = long, y = lat, label = carcass)) + 
  geom_sf(data = lakes, colour = '#a6bddb', fill = '#a6bddb') +
  coord_sf(xlim = c(34.5, 36), ylim =c(-2.4, -3.8) , expand = FALSE) +
  xlab("Longitude") + ylab("Latitude") 


nca

#draw main plot
ggplot() +
  draw_plot(nca) +
  draw_plot(eastafrica, x = 0.1, y = 0.7, width = 0.3, height = 0.3)

#plot densely annotated supp map
supp <- ggplot(data = world) +
  geom_sf(fill = "#fcfcf5") +
  geom_sf(data = Ngorongoro, fill = NA) +
  geom_point(data = mapmetadata, aes(x = mapmetadata$long, y =mapmetadata$lat, color = Cluster), size = 1.5) +
  scale_color_manual(values=c("#b3b3b3", "#d7191c", "#fdae61", "#abd9e9", "#2c7bb6")) +
  geom_text_repel(data = mapmetadata, aes(x = long, y = lat, label = carcass)) + 
  geom_sf(data = lakes, colour = '#a6bddb', fill = '#a6bddb') +
  coord_sf(xlim = c(34.5, 36), ylim =c(-2.4, -3.8) , expand = FALSE) +
  #cale_fill_viridis_c(trans = "sqrt", alpha = .4) +
  annotation_scale(location = "bl", width_hint = 0.4) +
  ggtitle("Sampling Sites", subtitle = "The Serengeti Conservation Area") +
  xlab("Longitude") + ylab("Latitude") 

#plot
supp

#############################################################################################
###Fiddle with tree
#############################################################################################
x <- read.tree('data/MSA_Banthracis_75samples_clean_shortened.afa.treefile')
metadata <- read.csv('data/metadata.csv', stringsAsFactors=FALSE, na.strings = c("", "NA"))
#midpoint root tree so is easier to look at
tree <- root(x, 'NC_007530')
#get only metadata with tips in tree
treedata <- subset(metadata, sample_id %in% tree$tip.label)



treedata

#create ggtree object and attach tree data
#plot tree with annotations
#some of the annotations mapped on funny so I ended up manually editing them in Illustrator because 
#I am too lazy to figure out what's going wrong. A false economy I will no doubt come to regret
#in the fullness of time.

treedata<- treedata %>% select(sample_id, geog, Cluster)
colnames(treedata) <- c("sample_id", "Geographical_Cluster", "Cluster")

p<-ggtree(tree)
p <- p %<+% treedata 

#main tree
 p + 
  geom_tippoint(aes(x=x+0.000015, color=Geographical_Cluster,shape = Cluster, order = Geographical_Cluster), size=4, alpha=.75) +
  geom_treescale() +
   geom_rootpoint() +
  scale_color_manual(values=c("#29ABE2", "#7AC943", "#FF931E", "#FF1D25", "#B2B2B2"), na.translate = F) +
  theme(legend.title = element_text(size = 10), 
          legend.text = element_text(size = 10),
          legend.position = c(0.75, 0.5)) +
  scale_shape_manual(values = c(15, 17, 18, 25, 16), na.translate = F) 

 
 #draw main plot
 ggplot() +
   draw_plot(nca) +
   draw_plot(eastafrica, x = 0.1, y = 0.7, width = 0.3, height = 0.3)
 
 
 
 
 
 
 
 
 
 
 
 
maintree

ggsave("maintree.svg", plot = maintree, device = "svg", path = "figures/")



#supp tree
supptree <- p + 
  geom_tippoint(aes(x=x+0.000015, color=Geographical_Cluster,shape = Cluster, order = Geographical_Cluster), size=4, alpha=.75) +
  geom_tiplab(aes(x=x+0.00003)) + 
  geom_treescale() +
  scale_color_manual(values=c("#29ABE2", "#7AC943", "#FF931E", "#FF1D25", "#B2B2B2"), na.translate = F) +
  theme(legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10),
        legend.position = c(0.75, 0.5)) +
  scale_shape_manual(values = c(15, 17, 18, 25, 16), na.translate = F) 
ggsave("supptree.svg", plot = supptree, device = "svg", path = "figures/")

supptree

hist(pairwise$ntdiff,
     main="Pairwise SNP Differences Between NCA Isolates",
     xlab="Number of nucleotide differences"
)
hist(withinsp$ntdiff,
     main="Pairwise SNP Differences Between NCA Isolates",
     xlab="Number of nucleotide differences"
)

