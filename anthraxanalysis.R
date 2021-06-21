##############################################################################################
#R script for anlaysis and generating figures in support of XX paper, Forde et al, 2020
#Tristan Dennis, October 2020
##############################################################################################

###############################
#import packages
pkg = c("tidyverse", "sf", "rnaturalearth", "rnaturalearthdata", "RColorBrewer", "ggtree", "treeio", "cowplot", "adegenet")
#install.packages(pkg) #install packages if you need them and load
new.packages <- pkg[!(pkg %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pkg, require, character.only = TRUE)#

#setwd
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
t2 <- t1 %>% 
  dplyr::select(idx, idy, ntdiff, carcass.x, carcass.y, Cluster.x, Cluster.y, species.x, species.y, geog.x, geog.y) %>% 
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
    #subset to 'carcass' and evaluate whether an observation is being made between samples from the same, or different, carcasses
    (ascale == 'carcass' & carcass.x == carcass.y) ~ 'Same Carcass',
    (ascale == 'carcass' & carcass.x != carcass.y) ~ 'Different Carcasses',
    #ditto for geographic group
    (ascale == 'geog' & geog.x == geog.y) ~ 'Same Group', 
    (ascale == 'geog' & geog.x != geog.y) ~ 'Different Group',
    #ditto for species
    (ascale == 'species' & species.x == species.y) ~ 'Same Species', 
    (ascale == 'species' & species.x != species.y) ~ 'Different Species',
    #ditto for epi cluster - note that for samples not in a cluster 'single case'
    #anything between a single case:single case, clsutered case:single case, etc is 'dUnlinked'
    #the d is because ggplot facets alphabetically which messes with my plotting - hacky, manual fix
    (ascale == 'epi' & Cluster.x == Cluster.y & Cluster.x != 'Single case' & Cluster.y != 'Single Case') ~ 'Same Cluster',
    (ascale == 'epi' & Cluster.x != Cluster.y | Cluster.x == 'Single case' | Cluster.y == 'Single case') ~ 'Different Cluster'))

################
#Now we need to apply some filters
#filter out same:same observations - this isn't necessary as `combn` above only outputs 
filt_t3 <-t3 # %>% filter(!(idx==idy))

t3 %>% filter(withinorbetween == 'Same Carcass') %>% ggplot(aes(x=ntdiff))+
  geom_histogram(binwidth = 1) + 
  theme_minimal()+
  labs(x='No. Nucleotide Differences', y='Count') +
  geom_vline(xintercept=2, colour = '#de2d26')

?geom_vline

s = t3 %>% filter(withinorbetween == 'Same Carcass')

#remove observations from epi clusters that are from the same carcass
filt_t4 <- filt_t3 %>% filter(!(ascale == 'epi' & carcass.x == carcass.y))
#unknown species:unknown species evaluates as spurious 'same species' - remove
filt_t4 <- filt_t4 %>% filter(!(ascale == 'species' & species.x == 'unknown' | species.y == 'unknown'))

################
#generate summary stats as three tables
#for each scale, across the nca, and distinct lineages from the same carcass

#for each 'scale' - within or between species, geog, epi, carcass
summarystatsfiltt4 <-filt_t4 %>% 
  dplyr::select(withinorbetween, ntdiff) %>% 
  group_by(withinorbetween) %>% 
  summarise(Median = median(ntdiff), Mean = mean(ntdiff), Min = min(ntdiff), Max = max(ntdiff), Lower_Quartile = quantile(ntdiff, 0.25), Upper_Quartile = quantile(ntdiff, 0.75), IQR = IQR(ntdiff))

#across the whole isolate set
ncawide <- inddist %>% 
  summarise(Median = median(dist), Mean = mean(dist), Min = min(dist), Max = max(dist), Lower_Quartile = quantile(dist, 0.25), Upper_Quartile = quantile(dist, 0.75), IQR = IQR(dist)) %>% 
  mutate(withinorbetween = paste('NCA-wide comparisons'))

#distinct lineages (non-identical isolates) from the same carcass
distinct_carc <- filt_t4 %>% 
  filter(withinorbetween == 'Same Carcass' & ntdiff > 0) %>% 
  summarise(Median = median(ntdiff), Mean = mean(ntdiff), Min = min(ntdiff), Max = max(ntdiff), Lower_Quartile = quantile(ntdiff, 0.25), Upper_Quartile = quantile(ntdiff, 0.75), IQR = IQR(ntdiff)) %>% 
  mutate(withinorbetween = paste('Distinct Isolates from Same Carcass'))

s = filt_t4 %>% 
  filter(withinorbetween == 'Same Carcass')

hist(s$ntdiff, breaks = 100)
#bind all 3 together
finalsummarystatstable <- rbind(summarystatsfiltt4, ncawide, distinct_carc)
s %>% filter(sample_id == 'LNA_D1'  )

################
#write to file - unfiltered df and filtered df, and summary stats
write_csv(t1, 'data/pairwise-observations-unfiltered.csv')
write_csv(filt_t4, 'data/pairwise-observations-filtered.csv')
write_csv(finalsummarystatstable, 'data/summarystatstfromfiltered.csv')

#############################################################################################
### Prepare final df for plotting
### 
#############################################################################################

#define ascale as factor level for plotting (capitalising letters, etc)
filt_t4$ascale <- factor(filt_t4$ascale, levels = c("geog", "epi", "carcass", "species"), labels = c("Geographic Group",  "Epidemiological Cluster", "Carcass", "Species"))

#check withinorbetween factor elvels
levels(as.factor(filt_t4$withinorbetween))

#define palette based on order of above
colorset = c('#33a02c','#33a02c','#33a02c', '#ff7f00', '#ff7f00', '#ff7f00')

#########
#violin plots of distance within and between scales
violplot <-
filt_t4 %>% filter((ascale != 'Species')) %>% 
  ggplot(., aes(x=withinorbetween, y=ntdiff, fill = withinorbetween)) + 
  scale_fill_manual(values=colorset) +
  geom_violin(trim=F, width = 1, show.legend = FALSE ) +
  facet_grid(. ~ ascale, scales = "free_x") +
  geom_boxplot(width=0.04, color="#252525", alpha=1.5, show.legend = FALSE ) +
  ylab("Pairwise distance (# nucleotides)") +
  xlab("")+
  ylim(0, 80) +
  theme_minimal(base_size = 22) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
violplot
ggsave("violinfig1.pdf", plot = violplot, device = 'pdf', path = "figures/", scale = 2)

#########
#histogram of disdtance within carcass
filt_t4
carchist <- filt_t4 %>% filter(ascale == 'Carcass' & withinorbetween == 'Same Carcass') %>% 
  ggplot(aes(x=ntdiff)) +
  geom_histogram(binwidth =1) +
  theme_minimal() +
  xlab("Nucleotide differences") +
  ylab("Count")
ggsave("histogram_carcass_ntdiff.pdf", plot = carchist, device = 'pdf', path = "figures/")
carchist
#########
#histogram of disdtance within NCA whole
ncahist <- filt_t4 %>% 
  select(idx, idy, ntdiff) %>% 
  distinct() %>% 
  ggplot(aes(x=ntdiff)) +
  geom_histogram(binwidth =1) +
  theme_minimal()+
  xlab("Nucleotide differences") +
  ylab("Count")
ggsave("histogram_ncawide_ntdiff.pdf", plot = ncahist, device = 'pdf', path = "figures/")


#########
#violin plots of distance within and between scales FOR SPECIES ONLY
specviolin <- filt_t4 %>% filter((ascale == 'Species')) %>% 
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
ggsave("species_violin.pdf", plot = specviolin, device = 'pdf', path = "figures/")
specviolin
#############################################################################################
###IS PAIRWISE DISTANCE A FUNCTION OF SCALE?
###i.e. do sequences tend to be more divergent the further away they are?
#############################################################################################
#get lat, long, idx, idt, ntdiff from uninflated indexed observations (t1)
#drop rows missing lat.long values
lalo <- t1 %>% 
  dplyr::select(idx, idy, ntdiff, lat.x, long.x, lat.y, long.y, clade.x, clade.y) %>% 
  distinct() %>% 
  drop_na()

#get distances between sampling points
for(i in 1:nrow(lalo)) {
t<-raster::pointDistance(c(lalo$lat.x[i], lalo$long.x[i]), c(lalo$lat.y[i], lalo$long.y[i]), lonlat = T)
lalo[i, 10] <- t
}

#plot distance over ntdiff
laloplot <- 
  ggplot(lalo, aes(x = V10, y = ntdiff)) +
  xlab("Distance in metres") +
  ylab("Distance (nucleotides)") +
  geom_point() +
  theme_minimal()
ggsave("allsampleslaloplot.pdf", plot = laloplot, device = 'pdf', path = "figures/")
laloplot



clade1lalo = lalo %>% filter(clade.x == 1 & clade.y == 1)
hist(clade1lalo$ntdiff)
write_csv(clade1lalo,'~/Projects/anthrax_diversity/anthrax_diversity/clade1lalo.csv')


m1 = glm(ntdiff ~ V10, family = 'gaussian', data=lalo)
summary(m1)

#plot distance over ntdiff - for clade 1 only
clade1laloplot <-
lalo %>% filter(clade.x == 1 & clade.y == 1) %>% 
  ggplot(aes(x = V10, y = ntdiff)) +
  xlab("Distance in metres") +
  ylab("Distance (nucleotides)") +
  geom_point() +
  theme_minimal()

ggsave("clade1laloplot.pdf", plot = clade1laloplot, device = 'pdf', path = "figures/")
clade1laloplot
#############################################################################################
###Make figure tree
#############################################################################################
x <- read.tree('data/MSA_Banthracis_75samples_clean_shortened.afa.treefile')
metadata <- read.csv('data/metadata.csv', stringsAsFactors=FALSE, na.strings = c("", "NA"))

#midpoint root tree so is easier to look at
tree <- root(x, 'NC_007530')
#get only metadata with tips in tree
treedata <- subset(metadata, sample_id %in% tree$tip.label)
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
maintree <- p + 
  geom_tippoint(aes(x=x+0.000015, color=Geographical_Cluster,shape = Cluster, order = Geographical_Cluster), size=4, alpha=.75) +
  geom_treescale() +
  geom_rootpoint() +
  scale_color_manual(values=c("#29ABE2", "#7AC943", "#FF931E", "#FF1D25", "#B2B2B2"), na.translate = F) +
  theme(legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10),
        legend.position = c(0.75, 0.5)) +
  scale_shape_manual(values = c(15, 17, 18, 25, 16), na.translate = F) 
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

#############
#IBD test 
#restrict to dominant clade #1 in lalo
#select ntdiff
s = lalo %>% filter(clade.x ==1) %>% filter(clade.y ==1) %>% dplyr::select(idx, idy, ntdiff) 
#select geodist
t = lalo %>% filter(clade.x ==1) %>% filter(clade.y ==1) %>% dplyr::select(idx, idy, V10) 
#fake 'bottom half'
x = s
#makie colnames for fake half the same
colnames(x) = c('idy', 'idx', 'ntdiff')
#bind into final
s = rbind(x, s)
#ditto for geodist matrix
y = t
colnames(y) = c('idy', 'idx', 'V10')
t = rbind(t, y)

#pivot into wider format to prepare for matrix coercion
s = pivot_wider(s, values_from = ntdiff, names_from = idx) 
t = pivot_wider(t, values_from = V10, names_from = idx) 

#NA to 0
s[is.na(s)] <- 0
t[is.na(t)] <- 0

#coerce to matrices
s = as.matrix(s)
t = as.matrix(t)
#make rownames
rownames(s) = s[,1]
rownames(t) = t[,1]
#remove old rowname cols
s = s[,-1]
t = t[,-1]
#mantel test
ibd = mantel.randtest(as.dist(s), as.dist(t), nrepet = 999)
#plot
plot(ibd)
#show output
ibd


