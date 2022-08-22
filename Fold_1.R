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
require('devtools')
require('INLAutils')
require('inlabru')
require('sf')
require('scales')
require('tidyverse')
require('ggplot2')
require('RColorBrewer')
require('blockCV')
require('automap')
library(naniar)

setwd("/home/ljochems/Fold_1/")
efb_fine <- read.csv("FullVegData11_21.csv")
#remove weird NA's at end of spreadsheet. hope that hasn't affected models? 
efb_fine <- drop_na(efb_fine)

plot(efb_fine$hyd_bin~efb_fine$wtr__)
plot(efb_fine$EFB_cover~efb_fine$wtr__)
plot(efb_fine$hyd_bin~efb_fine$typ_cover)
plot(efb_fine$EFB_cover~efb_fine$typ_cover)
#won't do quadratics for now

#make into spdf
coordinates(efb_fine) <- ~UTM_X + UTM_Y
proj4string(efb_fine) <- "+proj=utm +zone=16,17 +datum=WGS84 +units=m +no_defs"

GLbuffer <- shapefile("Great_Lakes_Buffer.shp")
proj4string(GLbuffer) <- "+proj=longlat +datum=WGS84 +no_defs"
GL_utm <- spTransform(GLbuffer, CRS("+proj=utm +zone=16,17 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))


####---determine range and create blocks----#### 
##----just envs by themselves 
BL <- autofitVariogram(NEAR_DIST~1,
                       efb_fine)
plot(BL)

fetch <- autofitVariogram(MeanFetch~1,
                          efb_fine)
plot(fetch)

typ <- autofitVariogram(typ_cover~1,
                        efb_fine)
plot(typ)

depth <- autofitVariogram(wtr__~1,
                          efb_fine)
plot(depth)

# ranges <- rbind(BL$var_model$range,
#                 fetch$var_model$range,
#                 typ$var_model$range,
#                 depth$var_model$range)
# 
# mean_range <- function(x) {
#   avg <- mean(x[,2]) %>% print()
#   med <- median(x[,2]) %>% print()
# }
# 
# mean_range(ranges)
#avg 31839.11 m, better than 220 km? 
#median 10707 m, which to use? 

# #make a spatialblock arrangement 
sb_pred <- spatialBlock(speciesData = efb_fine,
                        species = "hyd_bin",
                        theRange = 10707,
                        k = 5,
                        showBlocks = TRUE,
                        biomod2Format = TRUE)
# #i think this is the way to go 
# 
# #####----set up cv loop for inla------###### 
# fine_folds <- sb_pred$foldID
# 
# efb_fine <- as.data.frame(efb_fine)
# efb_fine <- efb_fine
# efb_fine$folds <- fine_folds
#write.csv(efb_fine,"FullVeg_Folds.csv")


#####--- see if we can get this to work (juanmi's way)-----##### 
cutoff <- 8.048969e-02
#build the mesh
mesh1 <- inla.mesh.2d(boundary=GL_utm,loc=efb_fine[c(43,44)],
                      cutoff=cutoff, max.edge=c(222000,888000),
                      crs = crs(GL_utm))
plot(mesh1);points(efb_fine[c(43,44)])
spde <- inla.spde2.matern(mesh1)

s.index <- inla.spde.make.index("s.index_mY", n.spde = spde$n.spde)

efb_fine <- as.data.frame(efb_fine)
A_matrix <- inla.spde.make.A(mesh1, loc = as.matrix(efb_fine[,c(43,44)]))

A_dim <- dim(A_matrix)


folds <- sb_pred$folds

for (k in seq_along(folds)) {
  assign(paste("IDtrainSet_",k,sep = ""),unlist(folds[[1]][1]))
}

for(k in seq_along(folds)) {
  assign(paste("IDtestSet_",k,sep =""),unlist(folds[[1]][2]))
}

train_list <- list(IDtrainSet_1,IDtrainSet_2,IDtrainSet_3,
                   IDtrainSet_4,IDtrainSet_5)
test_list <- list(IDtestSet_1,IDtestSet_2,IDtestSet_3,
                  IDtestSet_4,IDtestSet_5)

#match and subset to respective training and testing datasets 
train_dfs <- lapply(train_list, function(i) efb_fine[match(i,efb_fine$ï..OID_),]) 
test_dfs <- lapply(test_list, function(i) efb_fine[match(i,efb_fine$ï..OID_),])
#second row of each has all NAs.. hopefully won't affect things down the road 

#set PA to NAs 
test_dfs <- lapply(seq_along(test_dfs),function(i,x) {x[[i]] %>% 
    mutate(hyd_bin = ifelse(as.numeric(hyd_bin) >= 0, "NA",as.numeric(hyd_bin)))},x=test_dfs)
#set abundance to NAs 
test_dfs <- lapply(seq_along(test_dfs),function(i,x) {x[[i]] %>% 
    mutate(EFB_cover = ifelse(as.numeric(EFB_cover) >= 0, "NA",as.numeric(EFB_cover)))},x=test_dfs)

#combine test_dfs to training dfs 
final_dfs <- mapply(rbind,train_dfs,test_dfs,SIMPLIFY = FALSE)

#make sure response vars have desired training and testing NAs 
sapply(final_dfs[[1]],tail)
#everything in quotes? characters? 
str(final_dfs[[1]])
#yes, but seems to change when you make inla stack 

final_dfs <- lapply(final_dfs,function(df) mutate_at(df, .vars = c("hyd_bin","EFB_cover"),as.double))
str(final_dfs[[1]])

z <- final_dfs[[1]]$hyd_bin
y <- ifelse(final_dfs[[1]]$EFB_cover > 0, yes=final_dfs[[1]]$EFB_cover,no = NA)
y_noNA <- y[!is.na(y)]
cens <- 0.001
y_noNA[y_noNA >= 1-cens] <- 1-cens
y[y >= 1-cens] <- 1-cens
    
stack.EFB_y <- inla.stack(data = list(alldata = cbind(y,NA)), 
                          A = list(A_matrix, 1),
                          effects = list(s.index_mY = spde$n.spde,
                                         list(b0Y = rep(1, nrow(final_dfs[[1]])),
                                              data.frame(depth=scale(final_dfs[[1]]$wtr__)[,1]),data.frame(typha=scale(final_dfs[[1]]$typ_cover)[,1]),
                                              data.frame(boats=scale(final_dfs[[1]]$NEAR_DIST)[,1]), data.frame(fetch=scale(final_dfs[[1]]$MeanFetch)[,1]),
                                              idY = rep(1,nrow(final_dfs[[1]])), idY2 = rep(1,nrow(final_dfs[[1]])),idY3 = rep(1,nrow(final_dfs[[1]])),
                                              idY4 = rep(1,nrow(final_dfs[[1]])))), 
                            tag = "Beta (EFB Cover)")
    
stack.EFB_z <- inla.stack(data = list(alldata = cbind(NA,z)), 
                          A = list(A_matrix, 1),
                          effects = list(s.index_mZ = spde$n.spde,
                                         list(b0Z = rep(1, nrow(final_dfs[[1]])),
                                              data.frame(depth=scale(final_dfs[[1]]$wtr__)[,1]),data.frame(typha=scale(final_dfs[[1]]$typ_cover)[,1]),
                                              data.frame(boats=scale(final_dfs[[1]]$NEAR_DIST)[,1]), data.frame(fetch=scale(final_dfs[[1]]$MeanFetch)[,1]),
                                                  idZ = rep(1,nrow(final_dfs[[1]])), idZ2 = rep(1,nrow(final_dfs[[1]])),idZ3 = rep(1,nrow(final_dfs[[1]])),
                                                  idZ4 = rep(1,nrow(final_dfs[[1]])))), 
                              tag = "Bernoulli (EFB Occurence)")
    
stackm <- inla.stack(stack.EFB_y,stack.EFB_z)
    
prec.prior <- list(prec = list(prec = list(initial = 40000, fixed = TRUE)))
    
formula_all <- alldata ~ 0 + b0Y + b0Z +
  f(s.index_mY,model=spde) +
  f(s.index_mZ, copy = "s.index_mY", hyper = list(beta = list(fixed = FALSE))) +
  f(idY, depth, hyper = prec.prior) +
  f(idZ, depth, hyper = prec.prior) +
  f(idY2, typha, hyper = prec.prior) +
  f(idZ2, typha, hyper = prec.prior) +
  f(idY3, boats, hyper = prec.prior) +
  f(idZ3, boats, hyper = prec.prior) +
  f(idY4, fetch, hyper = prec.prior) +
  f(idZ4, fetch, hyper = prec.prior) 
    
    
EFB.hurdlemodel.inla <- inla(formula_all,
                             data = inla.stack.data(stackm),
                             control.predictor = list(A = inla.stack.A(stackm), link = c(rep(1,length(y)), rep(2,length(z))), compute = TRUE),
                             control.compute = list(dic = T, waic = T, config = T,
                                                    hyperpar = T, return.marginals=T), 
                             family = c("beta", "binomial"),
                             control.family = list(list(link = 'logit'), list(link = 'cloglog')), 
                             control.fixed = list(mean = c(0,0), prec = c(list(prior = "pc.prec", param = c(0.2, 0.85)), list(prior = "pc.prec", param = c(0.2, 0.85)))),
                             control.inla(list(strategy="laplace")),
                             verbose = T)
    
EFB.hurdlemodel.inla.complete <- list(summary.fixed = EFB.hurdlemodel.inla$summary.fixed,
                                      summary.hyperpar = EFB.hurdlemodel.inla$summary.hyperpar,
                                      summary.fitted.values = EFB.hurdlemodel.inla$summary.fitted.values,
                                      summary.random = EFB.hurdlemodel.inla$summary.random,
                                      marginals.fixed = EFB.hurdlemodel.inla$marginals.fixed,
                                      marginals.random = EFB.hurdlemodel.inla$marginals.random,
                                      marginals.hyperpar = EFB.hurdlemodel.inla$marginals.hyperpar,
                                      internal.marginals.hyperpar = EFB.hurdlemodel.inla$internal.marginals.hyperpar,
                                      marginals.spde2.blc = EFB.hurdlemodel.inla$marginals.spde2.blc,
                                      marginals.spde3.blc = EFB.hurdlemodel.inla$marginals.spde3.blc,
                                      internal.marginals.hyperpar = EFB.hurdlemodel.inla$internal.marginals.hyperpar)
    
    
save(EFB.hurdlemodel.inla.complete, file = "Fold1.RData")
    
