require('INLA')
require('raster')
require('plotly')
require('lattice')
require('grid')
require('gridExtra')
require('raster')
require('plotly')
require('dplyr')
require('tidyr')
require('tidyverse')
require('sp')
require('rgdal')
library(devtools)
library(INLAutils)
library(inlabru)
library(sf)
library(scales)
library(inlabru)
library(tidyverse)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

#data 
efb_data <- read.csv("FullVegData11_21.csv")
#remove weird NA's at end of spreadsheet. hope that hasn't affected models? 
efb_data <- drop_na(efb_data)

#model objects 
pc_fixed <- get(load("~/R/EFBDrivers2011_2021/AllTerms_PCfixed.RData")) 

#####----stacks and formulas for reference----##### 
#create beta (y) amd Bernoulli (z) response vars (from Juanmi's code)
z <- as.numeric(efb_data$hyd_bin)
y <- ifelse(efb_data$EFB_cover > 0, yes=efb_data$EFB_cover,no = NA)

y_noNA <- y[!is.na(y)]
#no zeros, but a few 1's 
#need censoring value to deal with 1 (100 % cover) and 
#small observed 0.5 % cover
#represent these as exactly 0 and 1, from inla.doc('^beta$')
#changed due to weird issue with censor value in inla() hurdle model
#instead just shave off value from 1's 
cens <- 0.001
#y_noNA[y_noNA <= cens] <- 0
y_noNA[y_noNA >= 1-cens] <- 1-cens

#apply to real y 
y[y >= 1-cens] <- 1-cens

# formula_all <- alldata ~ 0 + b0Y + b0Z +
#   f(s.index_mY,model=spde) +
#   f(s.index_mZ, copy = "s.index_mY", hyper = list(beta = list(fixed = FALSE))) +
#   f(idY, depth, hyper = FALSE) +
#   f(idZ, depth, hyper = FALSE) +
#   f(idY2, typha, hyper = FALSE) +
#   f(idZ2, typha, hyper = FALSE) +
#   f(idY3, boats, hyper = FALSE) +
#   f(idZ3, boats, hyper = FALSE) +
#   f(idY4, fetch, hyper = FALSE) +
#   f(idZ4, fetch, hyper = FALSE)

#####-------- mesh plot ------##### 
######-----mesh building----#####
GLbuffer <- shapefile("Great_Lakes_Buffer.shp")
proj4string(GLbuffer) <- "+proj=longlat +datum=WGS84 +no_defs"
GL_utm <- spTransform(GLbuffer, CRS("+proj=utm +zone=16,17 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

plot(GL_utm)
points(efb_data[,c(43,44)], col = "red", pch = 16)
#lake polys overlap all sites

#check for replicates
nrow(unique(efb_data[,c(43,44)])) 

# detect coordinates duplicated
dup <- efb_data[duplicated(efb_data[,c(43,44)]),c(43,44)]
dup
#none (as it should be!)
# full_veg[c(1068:1070),]
# # discard duplicate points, different year?
# full_veg <- full_veg[-as.numeric(row.names(dup)), ]

# distance between points
dis <- as.matrix(dist(coordinates(efb_data[, c(43,44)])))
diag(dis) <- NA
dis2 <- dis[lower.tri(dis, diag = FALSE)] #returns lower triangle of matrix distance values 
range(dis2) #very small min distance (transects ~0.1 apart)
hist(dis2)
cutoff <- 8.049944e-02

# # extract the coordinates to build the mesh
# coordsTra <- cbind(data$LON, data$LAT)
# create the study area polygon as a boundary in the mesh creator
#boundary <- inla.nonconvex.hull(points = coordsTra,convex = 0.5, concave = 0.4)

# build the mesh
mesh1 <- inla.mesh.2d(boundary=GL_utm,loc=efb_data[c(43,44)],cutoff=cutoff, max.edge=c(222000,888000),
                      crs = crs(GL_utm)) # may need to convert to m, max.edge = c(2,8)
# plot the mesh and sampled locations
plot(mesh1);points(efb_data[,c(43,44)], col = "red", pch = 16)

efb_coords <- efb_data[,c("ï..OID_", "UTM_X","UTM_Y")]
coordinates(efb_coords) <- ~UTM_X + UTM_Y
proj4string(efb_coords) <- "+proj=utm +zone=16,17 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

efb_scale <- as.data.frame(efb_coords)

library(stringr)
library(ggsn)

ggplot() +
  gg(mesh1) + 
  gg(efb_coords) + 
  labs(x = "UTM Easting (m)") + 
  labs(y = "UTM Northing (m)") + 
  north(scale = 1, location= "bottomright", symbol = 12, x.min = 1500000, x.max=1600000,
        y.min = 4300000,y.max=4700000) +
  scalebar(location = "bottomright", dist = 250000, dist_unit = "m", st.size=3, st.dist = 0.2, transform = FALSE,
           x.min = 1300000, x.max = 1650000, y.min = 4400000, y.max = 4600000)
ggsave("Figure2.jpg",width=8,height=6,units=c("in"),dpi=600)


######-----plotting posteriors----#### 
#####----posterior marginals for beta observations (shared spatial effect)----#### 
plot(inla.tmarginal(function(x) x, pc_clog$marginals.hyperpar$`Beta for s.index_mZ`), type ='l',
     main = 'Beta Obs. for Shared Spatial Random Effect')
#same thing with internal marginal hyperparams 
#plot(inla.tmarginal(function(x) x, pc_clog$internal.marginals.hyperpar$`Beta_intern for s.index_mZ`), type ='l')
#variance of beta obser. 
plot(inla.tmarginal(function(x) 1/x,
                    pc_clog$marginals.hyperpar$`precision parameter for the beta observations`), type = 'l',
     main = 'Marg. Variance for Spatial Random Effect') 


# y is cover 
plot(inla.tmarginal(function(x) x, pc_clog$marginals.fixed$b0Y), type ='l', 
     main = 'Posterior Marginal of % Cover (PC)')
# beta intercept still whacky with new data... 

# z is occurrence 
# still whacky without spatial terms (11-18 data)
# and crazy distrbution? tails beyond +/- 100 
plot(inla.tmarginal(function(x) x, pc_clog$marginals.fixed$b0Z), type ='l',
     main = 'Posterior Marginal of Occurence (PC cloglog)')
#differnt comparisons, but z looks much better with new dataset (at least with ZIB)


plot(inla.tmarginal(function(x) x, pc_clog$marginals.random$idY$index.1), type ='l',
     main = "Depth on Cover (PC Cloglog)")
plot(inla.tmarginal(function(x) x, pc_clog$marginals.random$idZ$index.1), type ='l',
     main = "Depth on Occ. (PC Cloglog)")
plot(inla.tmarginal(function(x) x, pc_clog$marginals.random$idY2$index.1), 
     type ='l', 
     main = "Typha on Cover (PC Cloglog)")
plot(inla.tmarginal(function(x) x, pc_clog$marginals.random$idZ2$index.1), 
     type ='l', 
     main = "Typha on Occ. (PC Cloglog)")
plot(inla.tmarginal(function(x) x, pc_clog$marginals.random$idY3$index.1), 
     type ='l',
     main = "Boats on Cover (PC Cloglog)")
plot(inla.tmarginal(function(x) x, pc_clog$marginals.random$idZ3$index.1), 
     type ='l', 
     main = "Boats on Occ. (PC Cloglog)")
plot(inla.tmarginal(function(x) x, pc_clog$marginals.random$idY4$index.1), 
     type ='l',
     main = "Fetch on Cover (PC Cloglog)")
plot(inla.tmarginal(function(x) x, pc_clog$marginals.random$idZ4$index.1), 
     type ='l', 
     main = "Fetch on Occ. (PC Cloglog)")



#####----spatial params & residuals----##### 
#issue in post modeling: 
#rds or rdata objects are NOT inla class, thus it is difficult to obtain spatial params
#outside of script in cluster... ASK Juanmi! 

# #to obtain range of spatial error terms across the nodes
length(EFB.hurdlemodel.inla$summary.random$s.index_mY$mean)
spatial_error <- range(EFB.hurdlemodel.inla$summary.random$s.index_mY$mean)
print("Range of Spatial Errors"); spatial_error
# #
# # # project the spatial random effect
gproj <- inla.mesh.projector(mesh1)
g.mean <- inla.mesh.project(gproj, EFB.hurdlemodel.inla$summary.random$s.index_mY$mean)
g.sd <- inla.mesh.project(gproj, EFB.hurdlemodel.inla$summary.random$s.index_mY$sd)
grid.arrange(levelplot(g.mean, scales=list(draw=F), xlab='', ylab='', main='mean',col.regions=terrain.colors(16)),
             levelplot(g.sd, scal=list(draw=F), xla='', yla='', main='sd',col.regions=terrain.colors(16)),nrow=1)
par(mfrow = c(1, 2))
#
# #not sure what's happening?
#
# # get the spatial parametres of the spatial random effect
spatial.parameters <- inla.spde2.result(inla = pc_fixed, name = "s.index_mY",
                                        spde = spde,
                                        do.transform = T)
# # nominal variance (the sill)
sill <- inla.emarginal(function(x) x, spatial.parameters$marginals.variance.nominal[[1]])
sill
# # plot posterior distribution
plot(spatial.parameters$marginals.variance.nominal[[1]],type='l',main='Sill')
#
# # range
range <- inla.emarginal(function(x) x, spatial.parameters$marginals.range.nominal[[1]])
range
# # plot posterior distribution
plot(spatial.parameters$marginals.range.nominal[[1]],type='l',main='Range')
#
# # # nugget, doesn't seem to work
nugget <- inla.emarginal(function(x) 1/x, EFB.hurdlemodel.inla$marginals.hyperpar$`precision parameter for the beta observations`)
nugget
# # # plot posterior distribution
plot(inla.tmarginal(function(x) 1/x, EFB.hurdlemodel.inla$marginals.hyperpar$`precision parameter for the beta observations`),
     type='l', main='Nugget')
#

# # plot model residuals of cover 
fitted <- pc_clog$"summary.fitted.values"[,1][1:length(efb_data$hyd_bin)][is.na(y) == F]
cov_sbst <- efb_data[which(efb_data$EFB_cover > 0),]
residuals <- (fitted - cov_sbst$EFB_cover) 

plot(fitted,residuals,main=('Residuals vs Fitted (Including spatial term)')); abline(h=0,lty='dashed',col='red')
#not so much randomly scatted? suggests that we may need quads in this model... 

#need to figure out how to calculate r2 for hurdle model
observed <- cov_sbst$EFB_cover
meanOb <- mean(cov_sbst$EFB_cover)
numerator <- sum((fitted - meanOb)^2)
denominator <- sum((observed - meanOb)^2)
r2 <- (numerator) / (denominator)
print('R2');r2
#not great.... 



#####---- try to plot spatial and pred. eff sizes----#### 
EFB.model.inla <- pc_fixed$summary.random[2:9]
EFB.model.inla <- do.call(rbind,lapply(EFB.model.inla,unlist))

#plot effect sizes 
THEME <- theme(legend.position='none', 
               axis.text.x = element_text(size = 12,colour = "black"),
               axis.text.y = element_text(size = 24, colour = "black"),
               axis.title.x = element_text(vjust = -0.35),
               #legend.key.height= unit(1,"in"),
               axis.title.y = element_text(vjust = 1.2)) #+ theme_bw()

library(ggdist)
#theme_set(theme_ggdist())
library(viridis)

Efxplot <- function(modellist, sig = TRUE, ModelNames = NULL){
  require(tidyverse);require(ggplot2)
  graphlist <- list()
  for(i in 1:length(modellist)){
    model <- modellist[[i]]
    
    if(class(model) == "matrix"){
      graph <- as.data.frame(model)
      colnames(graph)[which(colnames(graph)%in%c("0.025quant","0.5quant","0.975quant"))] <- c("Lower","Middle","Upper") 
      colnames(graph)[which(colnames(graph)%in%c("mean"))] <- c("Estimate")
    }
    
    
    graph$Model <- i
    graph$Factor <- rownames(graph)
    
    graphlist[[i]] <- graph
    
  }
  
  graph <- bind_rows(graphlist)
  
  graph$starloc <- NA
  
  min <- min(graph$Lower,na.rm = T) #,na.rm = T
  max <- max(graph$Upper,na.rm = T)#,na.rm = T
  
  if(sig == TRUE){
    graph$starloc <- max + (max-min)/10
  }
  
  #graph$Sig <- with(graph,ifelse((Lower<0&Upper<0)|(Lower>0&Upper>0),"*",""))
  
  graph$Model <- as.factor(graph$Model)
  
  if(!is.null(ModelNames)){
    levels(graph$Model) <- ModelNames
  }
  
  
  ggplot(as.data.frame(graph),aes(x = Factor,y = Estimate,colour = Model)) + 
    geom_pointinterval(position = position_dodge(w = 0.5), aes(ymin=Lower,ymax=Upper),size=10,width=2) + #change width for potential icons
    #geom_pointinterval(aes(ymin=Lower,ymax=Upper), size=8,width=6, color="black") + 
    # geom_point(position = position_dodge(w = 0.5)) + 
    # geom_errorbar(position = position_dodge(w = 0.5),aes(ymin = Lower,ymax = Upper),size = 0.8,width = 0.6) +
    #geom_errorbar(position = position_dodge(w = 0.5),aes(ymin=Estimate+Middle,ymax=Estimate-Middle),size = 0.9,width = 0.3) +
    #scale_color_viridis(discrete=TRUE, option="viridis") +
    THEME + labs(x = NULL, y = "Effect Size (Odds Ratio)") + geom_hline(aes(yintercept = 0),lty = 2,size=1) +
    coord_flip() #+ ylim(0.5,1.075)
  #geom_text(aes(label = Sig,y = starloc),position = position_dodge(w = 0.5))
}

ModelList<- list(EFB.model.inla)
plot_nT <- Efxplot(ModelList,ModelNames=c("Spatial"))
#works! need to plogis() everything tho 

####----plogis plot----#### 
depth_func <- function(x) {(plogis(x)*sd(efb_data$wtr__))/mean(efb_data$wtr__)} #/sd(efb_data$wtr_dp_)
typha_func <- function(x) {(plogis(x)*sd(efb_data$typ_cover))/mean(efb_data$typ_cover)} #/sd(efb_data$typh_cm)
boat_func <- function(x)  {(plogis(x)*sd(efb_data$NEAR_DIST))/mean(efb_data$NEAR_DIST)} #/sd(efb_data$NEAR_DIST)
fetch_func <- function(x) {(plogis(x)*sd(efb_data$MeanFetch))/mean(efb_data$MeanFetch)} #/sd(efb_data$MeanFetch)

# depth_func <- function(x) {plogis(x)/sd(efb_data$wtr_dp_)} #sd(efb_data$wtr_dp_)
# typha_func <- function(x) {plogis(x)/sd(efb_data$typh_cm)} #sd(efb_data$typh_cm)
# boat_func <- function(x)  {plogis(x)/sd(efb_data$NEAR_DIST)} #/sd(efb_data$NEAR_DIST)
# fetch_func <- function(x) {plogis(x)/sd(efb_data$MeanFetch)}#/sd(efb_data$MeanFetch)

plogs_func <- function (x) {
  df <- depth_func(x[1:2,]) 
  tf <- typha_func(x[3:4,])
  bf <- boat_func(x[5:6,])
  ff <- fetch_func(x[7:8,])
  final <- rbind(df,tf,bf,ff)
  return(final)
}

EFB.model.inlaPL <- plogs_func(EFB.model.inla)

ModelList_T <- list(EFB.model.inlaPL)
plot_T <- Efxplot(ModelList_T,ModelNames=c("Spatial"))

png("Figure3.png",width=1000,height=700,res=600)
plot_T
dev.off()

plot_nT 
ggsave("Figure3.jpg",width=8,height=6,units=c("in"),dpi=600)


require(stats)
x <- matrix(1:10, ncol = 2)
(centered.x <- scale(x, scale = FALSE))
cov(centered.scaled.x <- scale(x)) # all 1

####----misc.----####
shared_spat <- inla.tmarginal(function(x) x, pc_fixed$marginals.hyperpar$`Beta for s.index_mZ`)

shared_spat <- as.data.frame(shared_spat)

mean(shared_spat$x)


# marg.variance <- inla.tmarginal(function(x) 1/x,
#                                 pc_clog$marginals.hyperpar$`Precision for idY`)
# 
# m <- inla.emarginal(function(x) x, marg.variance)
# m
# 
# mm <- inla.emarginal(function(x) x^2, marg.variance)
# sqrt(mm - m^2)
# 
# inla.qmarginal(c(0.025, 0.5, 0.975), marg.variance)
# 
# inla.zmarginal(marg.variance)
# 
# inla.zmarginal(pc_clog$marginals.random$idY$index.1)
# 
# # plot(inla.tmarginal(function(x) plogis(x)/sd(efb_data$EFB_cover), pc_clog$marginals.fixed$b0Y))
# # plot(inla.tmarginal(function(x) plogis(x)/sd(efb_data$hyd_bin), pc_clog$marginals.fixed$b0Y))
# 
# inla.tmarginal(function(x) plogis(x)/sd(efb_data$wtr_dp_), pc_clog$marginals.random$idY$index.1)
# 
# inla.tmarginal(function(x) plogis(pc_clog$marginals.random$idY$index.1[,1]) / sd(efb_data$wtr_dp_)) 
# 
# list_marginals <- as.list(as.data.frame(pc_clog$marginals.random[2:9]))
# 
# depth_y <- plot(inla.tmarginal(function(x) plogis(x)/sd(efb_data$wtr_dp_), pc_clog$marginals.random$idY$index.1), type='l')
# depth_z <- plot(inla.tmarginal(function(x) plogis(x)/sd(efb_data$wtr_dp_), pc_clog$marginals.random$idZ$index.1), type='l')
# 
# typha_y <- plot(inla.tmarginal(function(x) plogis(x)*sd(efb_data$typh_cm), pc_clog$marginals.random$idY2$index.1), type='l')
# typha_z <- plot(inla.tmarginal(function(x) plogis(x)/sd(efb_data$typh_cm), pc_clog$marginals.random$idZ2$index.1), type='l')
# 
# boats_y <- plot(inla.tmarginal(function(x) plogis(x)/sd(efb_data$NEAR_DIST), pc_clog$marginals.random$idY3$index.1), type='l')
# boats_z <- plot(inla.tmarginal(function(x) plogis(x)/sd(efb_data$NEAR_DIST), pc_clog$marginals.random$idZ3$index.1), type='l')
# 
# fetch_y <- plot(inla.tmarginal(function(x) plogis(x)/sd(efb_data$MeanFetch), pc_clog$marginals.random$idY4$index.1), type='l')
# fetch_z <- plot(inla.tmarginal(function(x) plogis(x)/sd(efb_data$MeanFetch), pc_clog$marginals.random$idZ4$index.1), type='l')
# 
# 
# #y is cover, z is occ 
# list_marginals <- as.list(as.data.frame(pc_clog$marginals.random[2:9]))
# 
# depth_func <- function(x) {plogis(x)/sd(efb_data$wtr_dp_)}
# typha_func <- function(x) plogis(x)/sd(efb_data$typh_cm)
# boat_func <- function(x) {plogis(x)/sd(efb_data$NEAR_DIST)}
# fetch_func <- function(x) {plogis(x)/sd(efb_data$MeanFetch)}
# 
# cov_depth_t <- lapply(list_marginals[1], depth_func)
# occ_depth_t <- lapply(list_marginals[3], depth_func)
# cov_typha_t <- lapply(list_marginals[5],typha_func)
# occ_typha_t <- lapply(list_marginals[7],typha_func)
# cov_boat_t <- lapply(list_marginals[9],boat_func)
# occ_boat_t <- lapply(list_marginals[11],boat_func)
# cov_fetch_t <- lapply(list_marginals[13],fetch_func)
# occ_fetch_t <- lapply(list_marginals[15],fetch_func)
# 
# cov_depth_y <- list_marginals[2]
# occ_depth_y <- list_marginals[4]
# cov_typha_y <- list_marginals[6]
# occ_typha_y <- list_marginals[8]
# cov_boat_y <- list_marginals[10]
# occ_boat_y <- list_marginals[12]
# cov_fetch_y <- list_marginals[14]
# occ_fetch_y <- list_marginals[16]
# 
# marg_df <- cbind(cov_depth_t$idY.index.1.x,cov_depth_y$idY.index.1.y,occ_depth_t$idZ.index.1.x,occ_depth_y$idZ.index.1.y,
#                  cov_typha_t$idY2.index.1.x, cov_typha_y$idY2.index.1.y, occ_typha_t$idZ2.index.1.x,occ_typha_y$idZ2.index.1.y,
#                  cov_boat_t$idY3.index.1.x,cov_boat_y$idY3.index.1.y, occ_boat_t$idZ3.index.1.x, occ_boat_y$idZ3.index.1.y, 
#                  cov_fetch_t$idY4.index.1.x, cov_fetch_y$idY4.index.1.y, occ_fetch_t$idZ4.index.1.x, occ_fetch_y$idZ4.index.1.y) 
# 
# 
# marg_df <- as.data.frame(marg_df)
# names(marg_df)[1:16] <- c("cov_depth_x","cov_depth_y","occ_depth_x", "occ_depth_y",
#                           "cov_typha_x","cov_typha_y","occ_typha_x","occ_typha_y",
#                           "cov_boat_x","cov_boat_y","occ_boat_x","occ_typha_y",
#                           "cov_fetch_x","cov_fetch_y","occ_fetch_x","occ_fetch_y")
# #cover intercept = 0.1990273
# 
# cov_depth_plot <- ggplot(marg_df[1:2],aes(x=marg_df$cov_depth_x,y=marg_df$cov_depth_y)) +
#   geom_line() +
#   labs(x="",y="Density") + 
#   geom_vline(xintercept = mean(marg_df$cov_depth_x), linetype="dashed") + 
#   ylim(0,65) + 
#   theme_bw()
# 
# occ_depth_plot <- ggplot(marg_df[3:4],aes(x=marg_df$occ_depth_x,y=marg_df$occ_depth_y)) +
#   geom_line() +
#   labs(x="",y="Density") + 
#   geom_vline(xintercept = mean(marg_df$occ_depth_x), linetype="dashed") + 
#   ylim(0,65) + 
#   theme_bw()
# 
# depth_grid <- grid.arrange(cov_depth_plot,occ_depth_plot,ncol=1)
# 
# cov_typha_plot <- ggplot(marg_df[5:6],aes(x=marg_df$cov_typha_x,y=marg_df$cov_typha_y)) +
#   geom_line() +
#   labs(x="",y="") + 
#   geom_vline(xintercept = mean(marg_df$cov_typha_x), linetype="dashed") +
#   ylim(0,65) + 
#   theme_bw()
# 
# occ_typha_plot <- ggplot(marg_df[7:8],aes(x=marg_df$occ_typha_x,y=marg_df$occ_typha_y)) +
#   geom_line() +
#   labs(x="",y="") + 
#   geom_vline(xintercept = mean(marg_df$occ_typha_x), linetype="dashed") +
#   ylim(0,65) + 
#   theme_bw()
# 
# typha_grid <- grid.arrange(cov_typha_plot,occ_typha_plot,ncol=1)
# 
# cov_boat_plot <- ggplot(marg_df[9:10],aes(x=marg_df$cov_boat_x,y=marg_df$cov_boat_y)) +
#   geom_line() +
#   labs(x="",y="") + 
#   geom_vline(xintercept = mean(marg_df$cov_boat_x), linetype="dashed") +
#   ylim(0,65) + 
#   theme_bw()
# 
# occ_boat_plot <- ggplot(marg_df[11:12],aes(x=marg_df$occ_boat_x,y=marg_df$occ_boat_y)) +
#   geom_line() +
#   labs(x="",y="") + 
#   geom_vline(xintercept = mean(marg_df$occ_boat_x), linetype="dashed") + 
#   ylim(0,65) + 
#   theme_bw()
# 
# boat_grid <- grid.arrange(cov_boat_plot,occ_typha_plot,ncol=1)
# 
# cov_fetch_plot <- ggplot(marg_df[13:14],aes(x=marg_df$cov_fetch_x,y=marg_df$cov_fetch_y)) +
#   geom_line() +
#   labs(x="",y="") + 
#   geom_vline(xintercept = mean(marg_df$cov_fetch_x), linetype="dashed") + 
#   ylim(0,65) + 
#   theme_bw()
# 
# occ_fetch_plot <- ggplot(marg_df[15:16],aes(x=marg_df$occ_fetch_x,y=marg_df$occ_fetch_y)) +
#   geom_line() +
#   labs(x="",y="") + 
#   geom_vline(xintercept = mean(marg_df$occ_fetch_x), linetype="dashed") +
#   ylim(0,65) + 
#   theme_bw()
# 
# fetch_grid <- grid.arrange(cov_fetch_plot,occ_fetch_plot,ncol=1)
# 
# png("Post_Predictors.png", width=600,height=1300)
# post_all <- grid.arrange(depth_grid,boat_grid,fetch_grid,typha_grid,ncol=4)
# dev.off()
# 
# # list_marginals <- pc_clog$marginals.random[2:9]
# # marginals <- data.frame(do.call(rbind, list_marginals))
# # marginals$id <- rep(names(list_marginals),
# #                           times = sapply(list_marginals, nrow))
# # 
# # #calculate exceedance probabiliteis 
# marg <- pc_clog$marginals.random$idY$index.1
# 1 - inla.pmarginal(q = 0, marginal = marg)
# 
# sapply(pc_clog$marginals.random,
#        FUN = function(marg){1-inla.pmarginal(q = 0, marginal = marg)})
# 
# #highest posterior density 
# inla.hpdmarginal(0.95, pc_clog$marginals.hyperpar$`Precision for idY`)