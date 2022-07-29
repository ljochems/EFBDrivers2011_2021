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
library(blockCV)
library(automap)

efb_data <- read.csv("FullVegData11_21.csv")
#remove weird NA's at end of spreadsheet. hope that hasn't affected models? 
efb_data <- drop_na(efb_data)

#make into spdf
coordinates(efb_data) <- ~UTM_X + UTM_Y
proj4string(efb_data) <- "+proj=utm +zone=16,17 +datum=WGS84 +units=m +no_defs"

GLbuffer <- shapefile("Great_Lakes_Buffer.shp")
proj4string(GLbuffer) <- "+proj=longlat +datum=WGS84 +no_defs"
GL_utm <- spTransform(GLbuffer, CRS("+proj=utm +zone=16,17 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))


####---determine range and create blocks----#### 
##----just envs by themselves 
BL <- autofitVariogram(NEAR_DIST~1, 
                       efb_data)
plot(BL)

fetch <- autofitVariogram(MeanFetch~1, 
                          efb_data)
plot(fetch)

typ <- autofitVariogram(typ_cover~1, 
                        efb_data)
plot(typ)

depth <- autofitVariogram(wtr__~1, 
                          efb_data)
plot(depth)

ranges <- rbind(BL$var_model$range,
                fetch$var_model$range,
                typ$var_model$range,
                depth$var_model$range)

mean_range <- function(x) {
  avg <- mean(x[,2]) %>% print()
  med <- median(x[,2]) %>% print()
}

mean_range(ranges)
#avg 31839.11 m, better than 220 km? 
#median 10706 m, which to use? 

#make a spatialblock arrangement 
sb_pred <- spatialBlock(speciesData = efb_data,
                        species = "hyd_bin",
                        theRange = 10706,
                        k = 5,
                        showBlocks = TRUE,
                        biomod2Format = TRUE)
#i think this is the way to go 

#####----set up cv loop for inla------###### 
fine_folds <- sb_pred$foldID


#attach each fold config. to separate dataframes 
efb_data <- as.data.frame(efb_data)

efb_fine <- efb_data
efb_fine$folds <- fine_folds

# detect coordinates duplicated
dup <- efb_data[duplicated(efb_data[,c(43,44)]),c(43,44)]
dup

#make mesh 
dis <- as.matrix(dist(coordinates(efb_fine[,c(43,44)])))
diag(dis) <- NA
dis2 <- dis[lower.tri(dis, diag = FALSE)] #returns lower triangle of matrix distance values 
range(dis2) #very small min distance (transects ~0.1-0.2m apart)
hist(dis2) #in degrees, convert to km but from conversion factor at equator, need to adjust for great lakes region 
cutoff <- 8.048969e-02

# build the mesh
mesh1 <- inla.mesh.2d(boundary=GL_utm,loc=efb_fine[c(43,44)],cutoff=cutoff, max.edge=c(222000,888000),
                      crs = crs(GL_utm)) # may need to convert to m, max.edge = c(2,8)
# plot the mesh and sampled locations
plot(mesh1);points(efb_fine[,c(43,44)], col = "red", pch = 16)
#looks decent with more or less all points on vertices 
#do I need a outer boundary? 
summary(mesh1)
print("Number of nodes");mesh1$n

# SPDE model, spatial random effect considered as gaussian field  and GMRF to get cov matrices 
#with many zeros based via Stochastic partial differential equation 
spde <- inla.spde2.matern(mesh1)

# to assist with keeping track of which elements relate to what effects, create an index, represents mesh
s.index <- inla.spde.make.index("s.index_mY", n.spde = spde$n.spde)

# to create a projector matrix (links the m mesh vertices to the n responses)
A_matrix <- inla.spde.make.A(mesh1, loc = as.matrix(efb_fine[,c(43,44)]))
# number of observations x number of nodes
A_dim <- dim(A_matrix)
print("dim of A matrix"); A_dim
# how many points are on nodes: 1s because there is only 1 value diferent from zero
node_pts <- table(apply(A_matrix,1,nnzero))
print("Observations on nodes - within triangles"); node_pts
#all obs on nodes! (all but 1 obs for 2011-21)

# how many columns (nodes) have no points (either on or within the triangles)
nodes_nopts <- table(apply(A_matrix, 2, sum) > 0)
print("Nodes not associated with any observation - associated"); nodes_nopts
# FALSE because 0 > 0 is FALSE

pts_assoc <- table(apply(A_matrix, 2, nnzero))
print("How nodes are associated with observations"); pts_assoc

nzero_wt <- which(apply(A_matrix, 2, sum)>0)[1]
print("The first node with a nonzero weight"); nzero_wt

#####-----create inla stacks-----##### 
#create beta (y) amd Bernoulli (z) response vars (from Juanmi's code)
z <- as.numeric(efb_fine$hyd_bin)
y <- ifelse(efb_fine$EFB_cover > 0, yes=efb_fine$EFB_cover,no = NA)

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

#one stack for each 
#need to fit with intercept, need to provide response variable as list
#stack for Beta process (EFB cover)
stack.EFB_y <- inla.stack(data = list(alldata = cbind(y,NA)), 
                          A = list(A_matrix, 1),
                          effects = list(s.index_mY = spde$n.spde,
                                         list(b0Y = rep(1, nrow(efb_fine)),
                                              data.frame(depth=scale(efb_fine$wtr__)[,1]),data.frame(typha=scale(efb_fine$typ_cover)[,1]),
                                              data.frame(boats=scale(efb_fine$NEAR_DIST)[,1]), data.frame(fetch=scale(efb_fine$MeanFetch)[,1]),
                                              idY = rep(1,nrow(efb_fine)), idY2 = rep(1,nrow(efb_fine)),idY3 = rep(1,nrow(efb_fine)),
                                              idY4 = rep(1,nrow(efb_fine)))), 
                          tag = "Beta (EFB Cover)")


#stack for Bernoulli process (EFB occurence)
stack.EFB_z <- inla.stack(data = list(alldata = cbind(NA,z)), 
                          A = list(A_matrix, 1),
                          effects = list(s.index_mZ = spde$n.spde,
                                         list(b0Z = rep(1, nrow(efb_fine)),
                                              data.frame(depth=scale(efb_fine$wtr__)[,1]),data.frame(typha=scale(efb_fine$typ_cover)[,1]),
                                              data.frame(boats=scale(efb_fine$NEAR_DIST)[,1]), data.frame(fetch=scale(efb_fine$MeanFetch)[,1]),
                                              idZ = rep(1,nrow(efb_fine)), idZ2 = rep(1,nrow(efb_fine)),idZ3 = rep(1,nrow(efb_fine)),
                                              idZ4 = rep(1,nrow(efb_fine)))), 
                          tag = "Bernoulli (EFB Occurence)")
#high corr bw fetch and REI, so just use fetch for now 

stackm <- inla.stack(stack.EFB_y,stack.EFB_z)
stackm_form <- inla.stack.data(stackm)

######----model formulation and fitting----#### 
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


trainSet <- which(efb_fine$folds !=1)
stackm$data$data[trainSet,]

for(k in seq_len(5)){
  # extracting the training and testing indices
  # this way only works with foldID
  trainSet <- which(folds != k) # training set indices
  testSet <- which(folds == k) # testing set indices
  EFB.hurdlemodel.inla <- inla(formula_all,
                               data = inla.stack.data(stackm$data$data[trainSet,]),
                               control.predictor = list(A = inla.stack.A(stackm), link = c(rep(1,length(y)), rep(2,length(z))), compute = TRUE),
                               control.compute = list(dic = T, waic = T, config = T,
                                                      hyperpar = T, return.marginals=T), 
                               family = c("beta", "binomial"),
                               control.family = list(list(link = 'logit'), list(link = 'cloglog')), 
                               control.fixed = list(mean = c(0,0), prec = c(list(prior = "pc.prec", param = c(0.2, 0.85)), list(prior = "pc.prec", param = c(0.2, 0.85)))),
                               control.inla(list(strategy="laplace")),
                               verbose = T)
   # fitting a maxent model using linear, quadratic and hinge features
  mx <- maxnet(p = pb[trainSet], 
               data = mydata[trainSet, ], 
               maxnet.formula(p = pb[trainSet], 
                              data = mydata[trainSet, ], 
                              classes = "default"))
  testTable <- pb_data[testSet, ] # a table for testing predictions and reference data
  testTable$pred <- predict(mx, mydata[testSet, ], type = "cloglog") # predict the test set
  # calculate area under the ROC curve
  precrec_obj <- evalmod(scores = testTable$pred, labels = testTable$Species)
  AUCs[k] <- auc(precrec_obj)[1,4] # extract AUC-ROC
}


