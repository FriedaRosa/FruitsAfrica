library(viridis); library(scales)
## Map: fruit width (angiosperms)
tdwg_final <- read.csv(file.path(data.dir,"tdwg_all_plantsMammals_env.csv"))
table(tdwg_final$accRealm)
tdwg_final2 <- tdwg_final

tdwg_final2$accRealm[tdwg_final2$LEVEL_NAME=="Baleares" & is.na(tdwg_final2$accRealm)] <- "Asia"
tdwg_final2$accRealm[tdwg_final2$LEVEL_NAME=="Cook Is." & is.na(tdwg_final2$accRealm)] <- "Pacific"
tdwg_final2$accRealm[tdwg_final2$LEVEL_NAME=="Caroline Is." & is.na(tdwg_final2$accRealm)] <- "Pacific"
tdwg_final2$accRealm[tdwg_final2$LEVEL_NAME=="Fiji" & is.na(tdwg_final2$accRealm)] <- "Pacific"
tdwg_final2$accRealm[tdwg_final2$LEVEL_NAME=="France" & is.na(tdwg_final2$accRealm)] <- "Asia"
tdwg_final2$accRealm[tdwg_final2$LEVEL_NAME=="Italy" & is.na(tdwg_final2$accRealm)] <- "Asia"
tdwg_final2$accRealm[tdwg_final2$LEVEL_NAME=="Hawaii" & is.na(tdwg_final2$accRealm)] <- "Pacific"
tdwg_final2$accRealm[tdwg_final2$LEVEL_NAME=="Kriti" & is.na(tdwg_final2$accRealm)] <- "Asia"
tdwg_final2$accRealm[tdwg_final2$LEVEL_NAME=="Marianas" & is.na(tdwg_final2$accRealm)] <- "Pacific"
tdwg_final2$accRealm[tdwg_final2$LEVEL_NAME=="Marquesas" & is.na(tdwg_final2$accRealm)] <- "Pacific"
tdwg_final2$accRealm[tdwg_final2$LEVEL_NAME=="New Caledonia" & is.na(tdwg_final2$accRealm)] <- "Pacific"
tdwg_final2$accRealm[tdwg_final2$LEVEL_NAME=="Portugal" & is.na(tdwg_final2$accRealm)] <- "Asia"
tdwg_final2$accRealm[tdwg_final2$LEVEL_NAME=="Sicilia" & is.na(tdwg_final2$accRealm)] <- "Asia"
tdwg_final2$accRealm[tdwg_final2$LEVEL_NAME=="Sardegna" & is.na(tdwg_final2$accRealm)] <- "Asia"
tdwg_final2$accRealm[tdwg_final2$LEVEL_NAME=="Spain" & is.na(tdwg_final2$accRealm)] <- "Asia"
tdwg_final2$accRealm[tdwg_final2$LEVEL_NAME=="Samoa" & is.na(tdwg_final2$accRealm)] <- "Pacific"
tdwg_final2$accRealm[tdwg_final2$LEVEL_NAME=="Santa Cruz Is." & is.na(tdwg_final2$accRealm)] <- "Pacific"
tdwg_final2$accRealm[tdwg_final2$LEVEL_NAME=="Tokelau-Manihiki" & is.na(tdwg_final2$accRealm)] <- "Pacific"
tdwg_final2$accRealm[tdwg_final2$LEVEL_NAME=="Tonga" & is.na(tdwg_final2$accRealm)] <- "Pacific"
tdwg_final2$accRealm[tdwg_final2$LEVEL_NAME=="Tuamotu" & is.na(tdwg_final2$accRealm)] <- "Pacific"
tdwg_final2$accRealm[tdwg_final2$LEVEL_NAME=="Vanuatu" & is.na(tdwg_final2$accRealm)] <- "Pacific"
tdwg_final2$accRealm[tdwg_final2$LEVEL_NAME=="Wallis-Futuna Is." & is.na(tdwg_final2$accRealm)] <- "Pacific"

table(tdwg_final2$accRealm)

colnames(tdwg_final2)[colnames(tdwg_final2)      # Rename two variable names
                      %in% c("bio1_mean", "bio4_mean", "bio12_mean", "bio15_mean")] <- c("MAT", "TempSeas", "AnnualPrec", "PrecSeas")

names(tdwg_final2)

tdwg_final3 <- tdwg_final2 %>% dplyr::select("max95FL_palms", "max95FW_angios", "max95SL_angios",
                                             "curr_max95BodySize","presNat_max95BodySize", 
                                             
                                             "medianSL_angios", "medianFW_angios", 
                                             "curr_medianBodySize", "presNat_medianBodySize",
                                             
                                             "MAT", "TempSeas", "AnnualPrec", "PrecSeas",
                                             "CH_Mean", "CH_PercCover",
                                             
                                             "LEVEL_3_CO", "sp_richness_palms", "sp_richness_angios", "sp_richness_PresNat", "sp_richness_curr", "LEVEL_NAME", 
                                             "LAT", "LONG", "X_COORD", "Y_COORD", "ISISLAND", "accRealm", "REALM_LONG", "CONTINENT")


names(tdwg_final3)
tdwg_final3 <- tdwg_final3 %>% filter(sp_richness_palms >3)

dd <- tdwg_final3


# For scale in the map: 
q_values_FS <- c(1.339301,  2.648429,  3.500000,  5.000000,  5.981220,  7.500000,  8.965212, 10.831859, 14.506237)
a <- c(18,21)
q_values <- append(q_values_FS, a)
size_v <- seq(1,3, length.out = 8)
b <- c(3.33, 3.7)
size_v <- append(size_v, b)

# Fruit Width mapped ========
target.col <- c("LAT", "LONG", "LEVEL_3_CO", "max95FW_angios")
tdwg_maxFW <- melt(dd[,target.col], id.vars = c("LAT", "LONG", "LEVEL_3_CO"), measure.vars = "max95FW_angios", na.rm=T)
tdwg_maxFW$value <- tdwg_maxFW$value/10
max_q_mdpt <- vector()#

my_colors <- c("#440154", "#46327e", "#365c8d", "#277f8e", "#1fa187", "#4ac16d", "#a0da39", "#fde725", "#f48849", "#db5c68" )

for(i in 1:(length(q_values)-1)){ max_q_mdpt[i] <- (q_values[i] + q_values[i+1])/2 }
tdwg_maxFW$maxFS_q <- cut(tdwg_maxFW$value, q_values, include.lowest = T)
levels(tdwg_maxFW$maxFS_q) <- gsub(levels(tdwg_maxFW$maxFS_q), pattern = ",", replacement = " - ")
levels(tdwg_maxFW$maxFS_q) <- gsub(levels(tdwg_maxFW$maxFS_q), pattern = "\\(|\\[|\\]", replacement = "")

maxFSsim_p_new <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), fill = "grey80", size = 0.1, colour = "white", data = shp.df2) +
  geom_point(aes(x = LONG, y = LAT, fill = maxFS_q, size = maxFS_q), data =tdwg_maxFW, alpha = 0.9, colour = "black", pch = 21) +
  coord_fixed() +
  map_theme + theme(plot.title = element_blank()) +
  guides(size="none", fill = guide_legend(title = "Maximum\n(95th percentile)\nfruit width [cm]",  
                                          override.aes = list(size = size_v))) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_size_manual(values =size_v) +
  scale_color_manual(values = my_colors, aesthetics = c("colour", "fill"))+
  #                   breaks = levels(tdwg_maxFW$maxFS_q), labels = levels(tdwg_maxFW$maxFS_q)) +
  theme(legend.position = "right", text = element_text(size = 13), legend.title = element_text(size = 13), legend.text = element_text(size = 13))

maxFSsim_p_new
