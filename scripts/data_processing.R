# Codes creating by Julie Vercelloni in February 2024
# Contact: j.vercelloni@aims.gov.au 

############################
############################ HERON ISLAND SURVEYS - DATA PROCESSING
############################

# 1. Collating table from CoralNet and ReefCloud  
# 2. Creating sub-sites of 100m area 
# 3. Extract geomorphic zones associated with images 
# 4. Save table 

rm(list=ls())

# Loading packages

source("R/packages.R")

# Read gemorphic map 
heron.morph <- st_read("data/hr_HeronManuscriptZones.shp")

############################
############################ PROCESS CORALNET 

# Read surveys table from CoralNet (2002 - 2018)
net_surveys <- read.csv("data/Surveys_2002_2018.csv")

# Load lookup table to match labels between the two programs 
lookup_label <- read.csv("data/Labelset_lookup.csv") %>%
  rename(Comm_benthic = Code_CoralNet)

# Add new labels in the table 
net_surveys_ready <- net_surveys %>%
  filter(Comm_benthic %in% lookup_label$Comm_benthic) %>% # may need to modify the lookup labelsets table to add the bleached corals (actually present in the CoralNet table)
  left_join(lookup_label) %>% # 39 class
  dplyr::select(id, year,lng, lat, DESCRIPTION, Proportion, Total, Origin)

colnames(net_surveys_ready) <- c("Name", "year", "Longitude", "Latitude", "Comm_benthic", "Cover", "Total", "Origin")

# Look at the cover changes 

heron_CoralNet <- net_surveys_ready %>%
  group_by(Comm_benthic, year) %>%
  summarize(Mean_prop = mean(Cover), 
            sd_prop = sd(Cover),
            n = n()) %>%
  mutate(SE  = sd_prop / sqrt(n)) %>%
  st_drop_geometry()

ggplot(heron_CoralNet, aes(x= year, y = Mean_prop, col = Comm_benthic)) +
  geom_point() + geom_line() + 
  geom_errorbar(aes(ymin = Mean_prop - 1.96*SE, ymax = Mean_prop + 1.96*SE), width = 0.2) + 
  theme_bw() + xlab("Year") + ylab("Cover")

############################
############################ PROCESS REEFCLOUD

# Read surveys table from ReefCloud (2020 - 2023)
xls_data <- read_excel("data/Point_Summary_Heron_RC.xls") 
# Load associated metadata
csv_data <- read.csv("data/PointsMerged.csv")

# Ensure photos are not duplicated in RC 
duplicated_image <- xls_data %>% group_by(Name) %>% tally() %>% filter(n>1) # yes 460 images are duplicated 
dat_clean <- xls_data %>% filter(!Name %in% duplicated_image$Name) 

# Merge the data based on the "Name" column
merged_data <- merge(xls_data, csv_data, by = "Name", all.x = TRUE)

# Remove rows with NA values
merged_data <- na.omit(merged_data) # check why there are NAs?? 

# Remove images with not 50 classified points 
RC_surveys_ready <- merged_data %>% filter(!total !=50) %>%
  gather(key = Comm_benthic, Proportion, 11:73) %>%
  dplyr::select(Name, year, total, Latitude, Longitude, Comm_benthic, Proportion) %>%
  mutate(Total = 50) %>%
  mutate(Origin = "ReefCloud") %>%
  filter(Comm_benthic %in% net_surveys_ready$Comm_benthic) %>% 
  mutate(Cover = as.numeric(Proportion)/100) %>%
  dplyr::select(Name, year,Longitude, Latitude, Comm_benthic, Cover, total, Origin)

colnames(RC_surveys_ready) <- c("Name", "year", "Longitude", "Latitude", "Comm_benthic", "Cover", "Total", "Origin")

# Look at the cover changes 

heron_RC <- RC_surveys_ready  %>%
  group_by(Comm_benthic, year) %>%
  summarize(Mean_prop = mean(Cover), 
            sd_prop = sd(Cover),
            n = n()) %>%
  mutate(SE  = sd_prop / sqrt(n)) %>%
  st_drop_geometry()

ggplot(heron_RC, aes(x= year, y = Mean_prop, col = Comm_benthic)) +
  geom_point() + geom_line() + 
  geom_errorbar(aes(ymin = Mean_prop - 1.96*SE, ymax = Mean_prop + 1.96*SE), width = 0.2) + 
  theme_bw() + xlab("Year") + ylab("Cover")

############################
############################ MERGING THE TWO TABLE 

all_data <- rbind(net_surveys_ready, RC_surveys_ready) # year 2019 is missing 

# remove big data tables from the environment 
rm(net_surveys, net_surveys_ready, RC_surveys_ready)

############################
############################ MAKING SUBSITES  

# beginning

# Read sites location and add it in the table 
heron.site <- st_read("data/Heron_Data_Aggregation_v2.shp")

unique_location_sf <- all_data %>%
  group_by(Name) %>%
  slice(1) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = st_crs(4326)) %>%
  st_transform(crs = st_crs(heron.site)) %>% # transform in utm coordinates 
  st_join(heron.site)

# Select needed variables 
dat_loc <- unique_location_sf %>%
  dplyr::select(Name, year, NewName) %>%
  mutate(UTMX = st_coordinates(.)[,1],
         UTMY  = st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  filter(!is.na(NewName)) %>%
  arrange(NewName, year) 

tree_loc <- dat_loc %>%
  dplyr::select(NewName, UTMX, UTMY) %>%
  split(.$NewName) 

# Clustering

clust_list <-
  tree_loc %>%
  map(.f = function(x){
    dist(x[,c("UTMX","UTMY")])
  }) %>%
  map(.f=hclust)

s_split <- unique(dat_loc$NewName)

cut.tree.list <- vector("list", length(s_split))

for (i in 1:length(s_split)){
  cut.tree.list[[i]] <-   filter(dat_loc, NewName == s_split[i] ) %>%
    mutate("100" = cutree( clust_list[[i]], h=50))%>% 
    mutate("200" = cutree( clust_list[[i]], h=100))%>% # in case of further testing 
    mutate("300" = cutree( clust_list[[i]], h=150))%>% # in case of further testing 
    gather(key = aggregation_distance, value = group,
           6:8)
}

# Keep 100m group only 
cut.tree.collapsed <- do.call(what = rbind, args = cut.tree.list)%>%
  filter(aggregation_distance==100)%>%
  mutate(splitting.var = interaction(NewName, aggregation_distance, group, drop=T)) 

## Keep groups with more than 3 images per year per group 

Rec.clean<-cut.tree.collapsed%>% 
  group_by(splitting.var,year) %>% 
  nest()%>% 
  mutate(m = map(data, function(d) nrow(d)))%>%
  filter(!m<3)%>%dplyr::select(year,splitting.var)

Rec.mod0<-inner_join(cut.tree.collapsed,Rec.clean)%>%
  arrange(aggregation_distance,group)%>%data.frame()

## At least 2 replicate through time
Rec.clean2<-Rec.mod0 %>% 
  group_by(splitting.var) %>% mutate(m2=length(unique(year)))%>% 
  filter(!m2<3)%>%dplyr::select(splitting.var)%>%
  distinct()

Rec.units<-inner_join(Rec.mod0,Rec.clean2)%>%
  arrange(aggregation_distance,group)%>%data.frame()%>%droplevels()

Rec_tally<-Rec.units%>%group_by(splitting.var,year)%>%tally
wrong_group<-names(which(table(Rec_tally$splitting.var)<=1))
# end 

############################
############################ AVERAGING BY GROUP AND COMM_BENTHIC NOW 

all_data_grouped <- all_data %>% 
  inner_join(Rec.units ,by="Name") %>%
  group_by(splitting.var, year.x, Comm_benthic) %>%
  mutate(Cover_grouped=mean(Cover))%>%
  filter(row_number()==1) %>%
  dplyr::select(NewName, year.x, Longitude, Latitude, Comm_benthic, Cover_grouped, Total, Origin)

colnames(all_data_grouped) <- c("Group", "Site", "year", "Longitude", "Latitude", "Comm_benthic", 
                                "Cover_grouped", "Total", "Origin")

############################
############################ JOIN GEOMORPHIC ZONES - need to add new variable to distinct ReefSlope North from ReefSlope South, check NAs

all_data_grouped_sf <- all_data_grouped %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = st_crs(4326))%>%
  st_transform(crs = st_crs(heron.morph)) %>%
  st_join(heron.morph)

############################
############################ PLOTS

# Plot locations over time 
ggplot() + 
geom_sf(data = heron.morph, fill="transparent") + 
geom_sf(data = all_data_grouped_sf, aes(col = as.factor(year))) + 
theme_bw()

# Plot benthic communities over time 

# At Heron scale 

heron_comm <- all_data_grouped_sf %>%
   mutate(Prop = Cover_grouped*100) %>%
   group_by(Comm_benthic, year) %>%
   summarize(Mean_prop = mean(Prop), 
   sd_prop = sd(Prop),
   n = n()) %>%
   mutate(SE  = sd_prop / sqrt(n)) %>%
   st_drop_geometry()

ggplot(heron_comm, aes(x= year, y = Mean_prop, col = Comm_benthic)) +
  geom_point() + geom_line() + 
  geom_errorbar(aes(ymin = Mean_prop - 1.96*SE, ymax = Mean_prop + 1.96*SE), width = 0.2) + 
  theme_bw() + xlab("Year") + ylab("Cover")
