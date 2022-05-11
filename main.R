library(data.table)
library(spatialprobit)
library(ggplot2)
#load data
if(file.exists("./sei_mar.csv")){
  sei_mar <- fread("./sei_mar.csv")
} else{
  load("dasee_database.RData")
  dasee <- data.table(dasee)
  cols_keep <- c("ANNEE", "CODGEO", "PMUN", "EMPLT", "DGF_COM", "CFT_COM", "DIT_COM", "DETTE_COM", "IND_PRIX", "DEPTOT_COM", "EVPOP", "EVEMP", "EMPHAB", "DIPLMIN", "CAPBEP", "BAC", "SUP", "ACTOCCUP1564", "CHOM1564", "POP1824", "REVUC", "denspopemp", "SUPERFICIE")
  `%!in%` <- Negate(`%in%`)

  dasee[, which(colnames(dasee) %!in% cols_keep == TRUE) := NULL]
  colSums(is.na(dasee))
  ## Taking paris department region
  sei_mar <- dasee[grep("^77", CODGEO), ]
  sei_mar <- sei_mar[ANNEE == 2014, ] ## selecting the year 2014
  dasee[, c("REVUC", "ANNEE") := NULL]
  fwrite(sei_mar, "./sei_mar.csv")
}

if(file.exists("./SHP/Seine-et-Marne.shp")){
  shp <- rgdal::readOGR("./SHP/", "Seine-et-Marne")
} else{
  shp <- rgdal::readOGR("./SHP/", "communes-20210101")
  shp <- shp[grep("^77", shp@data$insee), ]
  rgdal::writeOGR(shp, "./SHP/", "Seine-et-Marne", "ESRI Shapefile")
}
plot(shp) # plotting the shapefile
#Finding the common communes between the shp and the data
common <- intersect(sei_mar$CODGEO, shp@data$insee)
sei_mar <- sei_mar[which(CODGEO %in% common == TRUE), ]
shp <- shp[shp@data$insee %in% common, ]

## Creating the weight matrix
nb_obj <- poly2nb(shp)
# which(card(nb_obj) == 0) # to see if any polygons have 0 neighbours
listw <- nb2listw(nb_obj)
weights_mat <- splm::listw2dgCMatrix(listw = listw)
frm <- DGF_COM ~ PMUN + EMPLT + DEPTOT_COM
# Create median variables
colz <- c("DGF_COM", "CFT_COM", "DIT_COM", "DETTE_COM")
dei_mar <- sei_mar[, c("CODGEO", ..colz)]
fwrite(dei_mar, "dei_mar.csv")
sei_mar[, (colz) := lapply(.SD, function(x){(x >= median(x))*1}), .SDcols = colz]

shp@data <- merge(shp@data, sei_mar, by.x = "insee", by.y = "CODGEO", all.X = TRUE)
rgdal::writeOGR(shp, "SHP/", layer = "Seine-et-Marne-data", driver = "ESRI Shapefile")

model <- sarprobit(DETTE_COM ~ PMUN + DIPLMIN, weights_mat, sei_mar, burn.in = 40000, showProgress = TRUE)

model1 <- sarprobit(DETTE_COM ~ PMUN + DIPLMIN, weights_mat, sei_mar, burn.in = 15000, showProgress = TRUE)
model2 <- sarprobit(DETTE_COM ~ PMUN + DIPLMIN - 1, weights_mat, sei_mar, burn.in = 15000, showProgress = TRUE)
model3 <- sarprobit(DETTE_COM ~ PMUN + DIPLMIN + CFT_COM, weights_mat, sei_mar, burn.in = 15000, showProgress = TRUE)
model4 <- sarprobit(DETTE_COM ~ PMUN + DIPLMIN + CFT_COM - 1, weights_mat, sei_mar, burn.in = 15000, showProgress = TRUE)
model5 <- sarprobit(DETTE_COM ~ 1, weights_mat, sei_mar, burn.in = 15000, showProgress = TRUE)
model6 <- sarprobit(DETTE_COM ~ 1 + CFT_COM, weights_mat, sei_mar, burn.in = 15000, showProgress = TRUE)
model7 <- sarprobit(DETTE_COM ~ CFT_COM + PMUN, weights_mat, sei_mar, burn.in = 15000, showProgress = TRUE)
model8 <- sarprobit(DETTE_COM ~ CFT_COM + DIPLMIN, weights_mat, sei_mar, burn.in = 15000, showProgress = TRUE)
model9 <- sarprobit(DETTE_COM ~ DIPLMIN, weights_mat, sei_mar, burn.in = 15000, showProgress = TRUE)

saveRDS(model1, "model1.RDS")
saveRDS(model2, "model2.RDS")
saveRDS(model3, "model3.RDS")
saveRDS(model4, "model4.RDS")
saveRDS(model5, "model5.RDS")
saveRDS(model6, "model6.RDS")
saveRDS(model7, "model7.RDS")
saveRDS(model8, "model8.RDS")

#model1 <- readRDS("model1.RDS")
#model2 <- readRDS("model2.RDS")
#model3 <- readRDS("model3.RDS")
#model4 <- readRDS("model4.RDS")
#model5 <- readRDS("model5.RDS")
#model6 <- readRDS("model6.RDS")
#model7 <- readRDS("model7.RDS")
#model8 <- readRDS("model8.RDS")

plot(density(model3$bdraw[, 4]), main = "Posterior Distribution of CFT_COM")
abline(v = 0, col = "red")
plot(density(model4$bdraw[, 3]), main = "Posterior Distribution of CFT_COM")
abline(v = 0, col = "red")
plot(density(model6$bdraw[, 4]), main = "Posterior Distribution of CFT_COM")
abline(v = 0, col = "red")
plot(density(model7$bdraw[, 4]), main = "Posterior Distribution of CFT_COM")
abline(v = 0, col = "red")
plot(density(model8$bdraw[, 4]), main = "Posterior Distribution of CFT_COM")
abline(v = 0, col = "red")


ggplot(as.data.frame(model3$bdraw, col), aes(x = V4))+
  geom_density()+
  geom_vline(xintercept = 0, colour = "red")+
  xlab("")+
  ylab("Density")+
  labs(title = "Posterior Distribution of CFT_COM", caption = "Model 4")+
  theme_classic()