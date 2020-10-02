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


devtools::install_github("ropensci/rnaturalearthhires") #set working directory to source file directory

#read in data
metadata <- read.csv('data/metadata.csv', stringsAsFactors=FALSE, na.strings = c("", "NA"))
#read in distance matrix and convert to pairwise obcservations vs matrix
mat <- read.delim('data/MSA_Banthracis_75samples_clean.distmatrix.txt', check.names=FALSE)
pairwise <- setNames(melt(mat), c('rows', 'vars', 'values'))
colnames <- c('idx', 'idy', 'ntdiff')
names(pairwise) <- colnames


pairwise <- filter(pairwise, idx != 'NC_007530' & idy != 'NC_007530')

#force lat long to numeric
metadata$long <- as.numeric(metadata$long)
metadata$lat <- as.numeric(metadata$lat)

#join metadata for both samples in each comparison
f <- left_join(pairwise, metadata, by = c("idx" = "sample_id"))
t <- left_join(f, metadata, by = c("idy" = "sample_id"))
#create combines dataframe
t$id <- paste(t$idx, t$idy)

#############################################################################################
###pariwise distance over nonlinear scales
###boxplots to illustrate whether pairwise distance between isolates is very different
#depending on whether you are looking within or between a given scale (i.e. village, host species, etc)
#############################################################################################

#let's select village data and get rid of nas as metadata are incomplete
c <- t %>% select(id, ntdiff, carcass.x, carcass.y, Cluster.x, Cluster.y, geog.x, geog.y) %>% .[complete.cases(.), ]


c
#for geographic cluster
#carcass
#epi cluster
c$carcass <- paste("carcass")
c$geog <- paste("geog")
c$epi <- paste("epi")

c <- c %>% 
  pivot_longer(c(carcass, geog, epi), names_to = "ascale") %>% 
  mutate(withinorbetween = case_when(
    (ascale == 'carcass' & carcass.x == carcass.y) ~ 'Same Carcass',
    (ascale == 'carcass' & carcass.x != carcass.y) ~ 'Different Carcasses',
    (ascale == 'geog' & geog.x == geog.y) ~ 'Same Cluster', 
    (ascale == 'geog' & geog.x != geog.y) ~ 'Different Cluster',
    (ascale == 'epi' & Cluster.x == Cluster.y) ~ 'Linked Cases',
    (ascale == 'epi' & Cluster.x != Cluster.y ~'dUnlinked Cases')))




c$ascale <- factor(c$ascale, levels = c("geog", "epi", "carcass"), labels = c("Geographic Cluster",  "Epidemiological Cluster", "Carcass"))

'#ff7f00'
levels(as.factor(c$withinorbetween))

colorset = c('#33a02c','#33a02c','#33a02c', '#ff7f00', '#ff7f00', '#ff7f00')
#########
#violin plots of distance within and between scales
c %>% ggplot(., aes(x=withinorbetween, y=ntdiff, fill = withinorbetween)) + 
  scale_fill_manual(values=colorset) +
  geom_violin(trim=FALSE, width = 1, show.legend = FALSE ) +
  facet_grid(. ~ ascale, scales = "free_x") +
  geom_boxplot(width=0.04, color="#252525", alpha=1.5, show.legend = FALSE ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  xlab("Difference in Scale") +
  ylab("Pairwise distance (# nucleotides)") +
  ylim(0, 50) +
  theme_minimal(base_size = 22) 

ggsave("violinfig1_cleaned_amended", plot = last_plot(), device = 'png', path = "figures/")

#########
#histogram of disdtance within carcass
c %>% filter(ascale == 'carcass' & withinorbetween == 'Same Carcass') %>% 
  ggplot(aes(x=ntdiff)) +
  geom_density() +
  theme_minimal()

#########
#histogram of disdtance within NCA whole
c %>% 
  select(id, ntdiff) %>% 
  distinct() %>% 
  ggplot(aes(x=ntdiff)) +
  geom_density() +
  theme_minimal()




#############################################################################################
###IS PAIRWISE DISTANCE A FUNCTION OF SCALE?
###i.e. do sequences tend to be more divergent the further away they are?
#############################################################################################

#get pairwise lat and long for pairs of samples
#remove imcomplete rows
lalo <- t %>% select(id, ntdiff, lat.x, long.x, lat.y, long.y)
c<-lalo[complete.cases(lalo), ] 

#get distances between sampling points
for(i in 1:nrow(c)) {
t<-raster::pointDistance(c(c$lat.x[i], c$long.x[i]), c(c$lat.y[i], c$long.y[i]), lonlat = T)
c[i, 7] <- t
}

#plot distance over ntdiff
ggplot(c, aes(x = V7, y = ntdiff)) +
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

eastafricabbox <- 




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

