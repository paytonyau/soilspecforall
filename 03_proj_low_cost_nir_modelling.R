# Factor: 
# Modelling: (70 cal. and 30 val.)
# 1. Only NIR from ASD 
# 2. Only NIR from NIRScan 
# 3. Only NIR from ASD using all available data and then predict spec from NIRScan
# 4. Adjusted NIR from ASD (using opendata) and then predict using NIRScan spec data
# 5. Only vis-NIR from ASD using all available data and then predict spec from ASD

# Libraries ----
if(!require()){install.packages(''); require()}
if(!require(parallel)){install.packages('parallel'); require(parallel)}
if(!require(doParallel)){install.packages('doParallel'); require(doParallel)}
if(!require(pls)){install.packages('pls'); require(pls)}
if(!require(Cubist)){install.packages('Cubist'); require(Cubist)}
if(!require(caret)){install.packages('caret'); require(caret)}
if(!require(ggplot2)){install.packages('ggplot2'); require(ggplot2)}

# Good of Fitness (GOOF) to get model evaluation
goof <- function(observed,predicted, plot.it = FALSE, type="DSM"){
  # Coefficient of determination
  rLM <- lm(predicted ~ observed)
  MEC <- as.matrix(summary(rLM)$adj.r.squared)
  
  # Standard error of prediction ^2
  SEP2 <- mean((observed - predicted)^2)
  
  # Standard error of prediction
  SEP <- sqrt(SEP2)
  
  #Bias
  bias <- mean(predicted) - mean(observed)
  
  # residual  variance
  SEP2c <- sum(((predicted - bias - observed)^2) / length(observed))
  SEPc <- sqrt(SEP2c)
  
  # ratio of performance to deviation
  RPD <- sd(observed) / SEP
  
  # Ratio of performance to interquartile distance
  IQ <- c(quantile(observed))[4] - c(quantile(observed))[2]
  RPIQ <- IQ / SEP
  
  # Concordance
  mx <- mean(observed)
  my <- mean(predicted)
  s2x <- var(observed)
  s2y <- var(predicted)
  sxy <- mean((observed-mx) * (predicted-my))
  ccc <- 2 * sxy / (s2x + s2y + (mx - my)^2)
  
  if (plot.it==TRUE){eqscplot(observed, predicted)
    abline(a = 0, b = 1, col = "brown4")}
  
  if (type == "DSM"){ gf <- data.frame(MEC=MEC, concordance=ccc, MSE=SEP2, RMSE=SEP, bias=bias, row.names=NULL)}
  else if (type == "spec"){ gf <- data.frame(MEC=MEC, concordance=ccc, MSE=SEP2, RMSE=SEP, bias=bias, 
                                             MSEc=SEP2c,RMSEc=SEPc, RPD=RPD, RPIQ=RPIQ, row.names=NULL)}
  else {stop("ERROR: Revise the type of output you require. Select from either DSM or spec")} 
  
  return(gf)
} 

# Setting the folder path ----
dir.p = 'D:/project_lowcost_nir/database' # homeoffice/

# Importing soil data ----
db.soil = read.csv(paste0(dir.p, '/tables/db_model/db_soil_oscar.csv'))
db.soil
db.soil$id = sprintf("osc_%02d", db.soil$id)
db.soil = db.soil[1:102,]
db.soil = db.soil[!(db.soil$tn.perc == "#N/A" & db.soil$tc.perc == "#N/A"), ]

# 1st model: ASD vis-NIR range (400 - 2450 nm) ----
db.asd.visnir.cal = read.csv(paste0(dir.p, '/tables/db_model/asd_visnir_sg_cal.csv'))[-1]
db.asd.visnir.cal[1:5, c(1:4, ncol(db.asd.visnir.cal))]
colnames(db.asd.visnir.cal) = gsub("X", "", colnames(db.asd.visnir.cal))

### filtering from 400 - 2450 nm
which(colnames(db.asd.visnir.cal) == "400") 
which(colnames(db.asd.visnir.cal) == "2450") 
colnames(db.asd.visnir.cal)[c(47, 2097)]
db.asd.visnir.cal = db.asd.visnir.cal[, c(1, 47:2097)]

db.soil.asd.visnir.cal = merge(db.soil, db.asd.visnir.cal, by="id")
db.soil.asd.visnir.cal[1:5, c(1:10, ncol(db.soil.asd.visnir.cal))]
str(db.soil.asd.visnir.cal[, 1:10])
db.soil.asd.visnir.cal$tn.perc = as.numeric(db.soil.asd.visnir.cal$tn.perc)
db.soil.asd.visnir.cal$tc.perc = as.numeric(db.soil.asd.visnir.cal$tc.perc)

cores = parallel::detectCores()-2

## PLSR
colnames(db.soil.asd.visnir.cal)[1:20]
db.tc.asd.visnir.pls = db.soil.asd.visnir.cal[, c(6, 9:ncol(db.soil.asd.visnir.cal))]
colnames(db.tc.asd.visnir.pls)[1:20]

pls::pls.options(parallel = makeCluster(cores, type = 'PSOCK'))
system.time(tc.asd.visnir.pls.m <- pls::plsr(tc.perc~., 
                                             ncomp = 15, 
                                             data = db.tc.asd.visnir.pls, 
                                             valid = 'CV', 
                                             scale = F)) 

stopCluster(pls.options()$parallel)

summary(tc.asd.visnir.pls.m)
tc.asd.visnir.pls.m$validation
ncomp.sel = as.data.frame(tc.asd.visnir.pls.m$validation$PRESS)
colnames(ncomp.sel) = gsub("comps", "", names(ncomp.sel))
ncomp = as.numeric(colnames(ncomp.sel[which.min(ncomp.sel[1,])]))

jpeg(paste0(dir.p, '/Results/figs_pls/ncomp/soc_asd_visnir_PLSR_Optimal_N_Comp.jpeg'), width = 10, height = 10, units = 'in', res = 300)
par(mar = c(6, 6, 6, 6))
plot(RMSEP(tc.asd.visnir.pls.m),
     main = '',
     xlim = c(0, 15),
     xaxt = 'n',
     frame.plot = FALSE,
     type = 'l',
     lty = 1:2,
     col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5),
     xlab =  'Number of components', # expression(paste('Wavelength(cm'^'-1',')')),
     ylab = 'RMSEP',
     las=1, font=2)
axis(1, at = seq(0, 15, by = 1), labels = T)
mtext('SOC SG-0D')
legend("topright",legend=c("Cross-Validation (CV)","Adjusted CV"),
       lty = 1:2,
       col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5),
       bty = 'n',
       x.intersp = 0.25,
       seg.len = 1)
points(ncomp, 0.2076, col = 'red', pch =20)
dev.off()

## Predicting
db.asd.visnir.val = read.csv(paste0(dir.p, '/tables/db_model/asd_visnir_sg_val.csv'))[-1]
db.asd.visnir.val[1:5, c(1:4, ncol(db.asd.visnir.val))]
colnames(db.asd.visnir.val) = gsub("X", "", colnames(db.asd.visnir.val))

### filtering from 400 - 2450 nm
which(colnames(db.asd.visnir.val) == "400") 
which(colnames(db.asd.visnir.val) == "2450") 
colnames(db.asd.visnir.val)[c(47, 2097)]
db.asd.visnir.val = db.asd.visnir.val[, c(1, 47:2097)]

db.soil.asd.visnir.val = merge(db.soil, db.asd.visnir.val, by="id")
db.soil.asd.visnir.val[1:5, c(1:10, ncol(db.soil.asd.visnir.val))]

### Creating a validation table for all results
db.soil.val.all = as.data.frame(db.soil.asd.visnir.val[, 1:6])
str(db.soil.val.all)
#write.csv(db.soil.val.all[,c(1, 5:6)], paste0(dir.p, '/tables/db_soil_asd_nir_val.csv'))
db.soil.val.all$tn.perc = as.numeric(db.soil.val.all$tn.perc)
db.soil.val.all$tc.perc = as.numeric(db.soil.val.all$tc.perc)

db.soil.val.all$tc.asd.visnir.pls = stats::predict(tc.asd.visnir.pls.m, 
                                                   newdata = db.asd.visnir.val[, 2:ncol(db.asd.visnir.val)],
                                                   ncomp = ncomp)

db.soil.val.all$tc.asd.visnir.pls = as.numeric(db.soil.val.all$tc.asd.visnir.pls[, 1, 1])

metric.tc.asd.visnir.pls = goof(observed = db.soil.val.all$tc.perc, 
                           predicted = db.soil.val.all$tc.asd.visnir.pls,
                           type = 'spec')

jpeg(paste0(dir.p, '/Results/figs_pls/soc_sg_asd_visnir_PLSR_Obs_Pred.jpeg'), width = 10, height = 10, units = 'in', res = 300)
par(mar = c(6, 6, 6, 6))
summary(db.soil.val.all[, c('tc.asd.visnir.pls', 'tc.perc')])
plot(db.soil.val.all$tc.asd.visnir.pls, db.soil.val.all$tc.perc,
     xlab = 'Predicted (%)', #expression(paste('Predicted ', '(', 'mg', ' kg'^'-1',')')),
     ylab = 'Observed (%)', #expression(paste('Observed ', '(', 'mg', ' kg'^'-1',')')),
     xlim = c(0, 2),
     ylim = c(0, 2),
     type = 'p',
     pch =  '🔅', # 16,
     bty = 'L',
     cex = 0.9,
     cex.lab= 1.1,
     cex.axis=1.2,
     col = rgb(red=0.0 , green = 0.0, blue = 0.0, alpha = .8),
     family = 'A')
abline(0, 1, col='black', lty=1, lwd=1)
abline(lm(db.soil.val.all$tc.asd.visnir.pls ~ db.soil.val.all$tc.perc), col = 'red', lty=2, lwd=1)
text(1.8, 0.5, bquote(RMSE == .(round(metric.tc.asd.visnir.pls$RMSE, 2))))
text(1.8, 0.4, bquote(MEC == .(round(metric.tc.asd.visnir.pls$MEC, 2))))
text(1.8, 0.3, bquote(CCC == .(round(metric.tc.asd.visnir.pls$concordance, 2))))
text(1.8, 0.2, bquote(bias == .(round(metric.tc.asd.visnir.pls$bias, 2))))
text(1.8, 0.1, bquote(RPIQ == .(round(metric.tc.asd.visnir.pls$RPIQ, 2))))
dev.off()

## Cubist
ctrl = caret::trainControl(method = 'repeatedcv', number = 10, repeats = 10, search = 'random')

cores # it has the number of cores from parallel::detectCores()-2
cl = parallel::makeCluster(cores) # Create a cluster using all cores
doParallel::registerDoParallel(cl)

#colnames(db.soil.asd.visnir.cal)[1:20]
system.time(tc.asd.visnir.cub.m <- caret::train(x=db.soil.asd.visnir.cal[ ,9:ncol(db.soil.asd.visnir.cal)], 
                                                y=db.soil.asd.visnir.cal$tc.perc,
                                                method = "cubist", 
                                                trControl = ctrl,
                                                importance = T))  
stopCluster(cl)

print(tc.asd.visnir.cub.m)

tc.asd.visnir.cub.m.var = varImp(tc.asd.visnir.cub.m, scale = T)

jpeg(paste0(dir.p, '/results/figs_cub/varimp/soc_sg_asd_visnir_Var_Importance_Cub.jpeg'), width = 10, height = 10, units = 'in', res = 300)
plot(tc.asd.visnir.cub.m.var, top = 20, main = "Original data")
dev.off()


## Predicting
db.soil.asd.visnir.val = merge(db.soil, db.asd.visnir.val, by="id")
db.soil.asd.visnir.val[1:5, c(1:10, ncol(db.soil.asd.visnir.val))]
db.soil.asd.visnir.val[1:5,1:10]

db.soil.val.all$tc.asd.visnir.cub = stats::predict(tc.asd.visnir.cub.m, 
                                                   db.soil.asd.visnir.val[, 9:ncol(db.soil.asd.visnir.val)]) # Predicting using external data

metric.tc.asd.visnir.cub = goof(observed = db.soil.val.all$tc.perc, 
                                predicted = db.soil.val.all$tc.asd.visnir.cub,
                                type = 'spec')


jpeg(paste0(dir.p, '/Results/figs_cub/soc_sg_asd_visnir_Cub_Obs_Pred.jpeg'), width = 10, height = 10, units = 'in', res = 300)
par(mar = c(6, 6, 6, 6))
summary(db.soil.val.all[, c('tc.asd.visnir.cub', 'tc.perc')])
plot(db.soil.val.all$tc.asd.visnir.cub, db.soil.val.all$tc.perc,
     xlab = 'Predicted (%)', #expression(paste('Predicted ', '(', 'mg', ' kg'^'-1',')')),
     ylab = 'Observed (%)', #expression(paste('Observed ', '(', 'mg', ' kg'^'-1',')')),
     xlim = c(0, 2),
     ylim = c(0, 2),
     type = 'p',
     pch =  '🔅', # 16,
     bty = 'L',
     cex = 0.9,
     cex.lab= 1.1,
     cex.axis=1.2,
     col = rgb(red=0.0 , green = 0.0, blue = 0.0, alpha = .8),
     family = 'A')
abline(0, 1, col='black', lty=1, lwd=1)
abline(lm(db.soil.val.all$tc.asd.visnir.cub ~ db.soil.val.all$tc.perc), col = 'red', lty=2, lwd=1)
text(1.8, 0.5, bquote(RMSE == .(round(metric.tc.asd.visnir.cub$RMSE, 2))))
text(1.8, 0.4, bquote(MEC == .(round(metric.tc.asd.visnir.cub$MEC, 2))))
text(1.8, 0.3, bquote(CCC == .(round(metric.tc.asd.visnir.cub$concordance, 2))))
text(1.8, 0.2, bquote(bias == .(round(metric.tc.asd.visnir.cub$bias, 2))))
text(1.8, 0.1, bquote(RPIQ == .(round(metric.tc.asd.visnir.cub$RPIQ, 2))))
dev.off()


# 2nd model: ASD NIR range (955 - 1645 nm) ----
db.asd.nir.cal = read.csv(paste0(dir.p, '/tables/db_model/db_target_asd_spec_cal.csv'))[-1]
db.asd.nir.cal[1:5, c(1:4, ncol(db.asd.nir.cal))]
colnames(db.asd.nir.cal) = gsub("NIR.", "", colnames(db.asd.nir.cal))
which(colnames(db.asd.nir.cal)=="1645")
db.asd.nir.cal = db.asd.nir.cal[, c(1, 3:141)]

db.soil.asd.nir.cal = merge(db.soil, db.asd.nir.cal, by="id")
db.soil.asd.nir.cal[1:5, c(1:10, ncol(db.soil.asd.nir.cal))]
str(db.soil.asd.nir.cal[, 1:10])
db.soil.asd.nir.cal$tn.perc = as.numeric(db.soil.asd.nir.cal$tn.perc)
db.soil.asd.nir.cal$tc.perc = as.numeric(db.soil.asd.nir.cal$tc.perc)

## PLSR
colnames(db.soil.asd.nir.cal)[1:20]
db.tc.asd.nir.pls = db.soil.asd.nir.cal[, c(6, 9:ncol(db.soil.asd.nir.cal))]
colnames(db.tc.asd.nir.pls)[1:20]

pls::pls.options(parallel = makeCluster(cores, type = 'PSOCK'))
system.time(tc.asd.nir.pls.m <- pls::plsr(tc.perc~., 
                                             ncomp = 15, 
                                             data = db.tc.asd.nir.pls, 
                                             valid = 'CV', 
                                             scale = F)) 

stopCluster(pls.options()$parallel)

summary(tc.asd.nir.pls.m)
tc.asd.nir.pls.m$validation
ncomp.sel = as.data.frame(tc.asd.nir.pls.m$validation$PRESS)
colnames(ncomp.sel) = gsub("comps", "", names(ncomp.sel))
ncomp = as.numeric(colnames(ncomp.sel[which.min(ncomp.sel[1,])]))

jpeg(paste0(dir.p, '/Results/figs_pls/ncomp/soc_asd_nir_PLSR_Optimal_N_Comp.jpeg'), width = 10, height = 10, units = 'in', res = 300)
par(mar = c(6, 6, 6, 6))
plot(RMSEP(tc.asd.nir.pls.m),
     main = '',
     xlim = c(0, 15),
     xaxt = 'n',
     frame.plot = FALSE,
     type = 'l',
     lty = 1:2,
     col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5),
     xlab =  'Number of components', # expression(paste('Wavelength(cm'^'-1',')')),
     ylab = 'RMSEP',
     las=1, font=2)
axis(1, at = seq(0, 15, by = 1), labels = T)
mtext('SOC SG-0D: ASD NIR')
legend("topright",legend=c("Cross-Validation (CV)","Adjusted CV"),
       lty = 1:2,
       col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5),
       bty = 'n',
       x.intersp = 0.25,
       seg.len = 1)
points(ncomp, 0.1718, col = 'red', pch =20)
dev.off()

## Predicting
db.asd.nir.val = read.csv(paste0(dir.p, '/tables/db_model/db_target_asd_spec_val.csv'))[-1]
db.asd.nir.val[1:5, c(1:4, ncol(db.asd.nir.val))]
colnames(db.asd.nir.val) = gsub("NIR.", "", colnames(db.asd.nir.val))
which(colnames(db.asd.nir.val)=="1645")
db.asd.nir.val = db.asd.nir.val[, c(1, 3:141)]

db.soil.val.all$tc.asd.nir.pls = stats::predict(tc.asd.nir.pls.m, 
                                                   newdata = db.asd.nir.val[, 2:ncol(db.asd.nir.val)],
                                                   ncomp = ncomp)

db.soil.val.all$tc.asd.nir.pls = as.numeric(db.soil.val.all$tc.asd.nir.pls[, 1, 1])

metric.tc.asd.nir.pls = goof(observed = db.soil.val.all$tc.perc, 
                                predicted = db.soil.val.all$tc.asd.nir.pls,
                                type = 'spec')

jpeg(paste0(dir.p, '/Results/figs_pls/soc_sg_asd_nir_PLSR_Obs_Pred.jpeg'), width = 10, height = 10, units = 'in', res = 300)
par(mar = c(6, 6, 6, 6))
summary(db.soil.val.all[, c('tc.asd.nir.pls', 'tc.perc')])
plot(db.soil.val.all$tc.asd.nir.pls, db.soil.val.all$tc.perc,
     xlab = 'Predicted (%)', #expression(paste('Predicted ', '(', 'mg', ' kg'^'-1',')')),
     ylab = 'Observed (%)', #expression(paste('Observed ', '(', 'mg', ' kg'^'-1',')')),
     xlim = c(0, 2),
     ylim = c(0, 2),
     type = 'p',
     pch =  '🔅', # 16,
     bty = 'L',
     cex = 0.9,
     cex.lab= 1.1,
     cex.axis=1.2,
     col = rgb(red=0.0 , green = 0.0, blue = 0.0, alpha = .8),
     family = 'A')
abline(0, 1, col='black', lty=1, lwd=1)
abline(lm(db.soil.val.all$tc.asd.nir.pls ~ db.soil.val.all$tc.perc), col = 'red', lty=2, lwd=1)
text(1.8, 0.5, bquote(RMSE == .(round(metric.tc.asd.nir.pls$RMSE, 2))))
text(1.8, 0.4, bquote(MEC == .(round(metric.tc.asd.nir.pls$MEC, 2))))
text(1.8, 0.3, bquote(CCC == .(round(metric.tc.asd.nir.pls$concordance, 2))))
text(1.8, 0.2, bquote(bias == .(round(metric.tc.asd.nir.pls$bias, 2))))
text(1.8, 0.1, bquote(RPIQ == .(round(metric.tc.asd.nir.pls$RPIQ, 2))))
dev.off()

## Cubist
cores # it has the number of cores from parallel::detectCores()-2
cl = parallel::makeCluster(cores) # Create a cluster using all cores
doParallel::registerDoParallel(cl)

system.time(tc.asd.nir.cub.m <- caret::train(x=db.soil.asd.nir.cal[ ,9:ncol(db.soil.asd.nir.cal)], 
                                                y=db.soil.asd.nir.cal$tc.perc,
                                                method = "cubist", 
                                                trControl = ctrl,
                                                importance = T))  
stopCluster(cl)

print(tc.asd.nir.cub.m)

tc.asd.nir.cub.m.var = varImp(tc.asd.nir.cub.m, scale = T)

jpeg(paste0(dir.p, '/results/figs_cub/varimp/soc_sg_asd_nir_Var_Importance_Cub.jpeg'), width = 10, height = 10, units = 'in', res = 300)
plot(tc.asd.nir.cub.m.var, top = 20, main = "Original data")
dev.off()


## Predicting
db.soil.val.all$tc.asd.nir.cub = stats::predict(tc.asd.nir.cub.m, 
                                                   db.asd.nir.val[, 2:ncol(db.asd.nir.val)]) # Predicting using external data

metric.tc.asd.nir.cub = goof(observed = db.soil.val.all$tc.perc, 
                                predicted = db.soil.val.all$tc.asd.nir.cub,
                                type = 'spec')


jpeg(paste0(dir.p, '/Results/figs_cub/soc_sg_asd_nir_Cub_Obs_Pred.jpeg'), width = 10, height = 10, units = 'in', res = 300)
par(mar = c(6, 6, 6, 6))
summary(db.soil.val.all[, c('tc.asd.nir.cub', 'tc.perc')])
plot(db.soil.val.all$tc.asd.nir.cub, db.soil.val.all$tc.perc,
     xlab = 'Predicted (%)', #expression(paste('Predicted ', '(', 'mg', ' kg'^'-1',')')),
     ylab = 'Observed (%)', #expression(paste('Observed ', '(', 'mg', ' kg'^'-1',')')),
     xlim = c(0, 2),
     ylim = c(0, 2),
     type = 'p',
     pch =  '🔅', # 16,
     bty = 'L',
     cex = 0.9,
     cex.lab= 1.1,
     cex.axis=1.2,
     col = rgb(red=0.0 , green = 0.0, blue = 0.0, alpha = .8),
     family = 'A')
abline(0, 1, col='black', lty=1, lwd=1)
abline(lm(db.soil.val.all$tc.asd.nir.cub ~ db.soil.val.all$tc.perc), col = 'red', lty=2, lwd=1)
text(1.8, 0.5, bquote(RMSE == .(round(metric.tc.asd.nir.cub$RMSE, 2))))
text(1.8, 0.4, bquote(MEC == .(round(metric.tc.asd.nir.cub$MEC, 2))))
text(1.8, 0.3, bquote(CCC == .(round(metric.tc.asd.nir.cub$concordance, 2))))
text(1.8, 0.2, bquote(bias == .(round(metric.tc.asd.nir.cub$bias, 2))))
text(1.8, 0.1, bquote(RPIQ == .(round(metric.tc.asd.nir.cub$RPIQ, 2))))
dev.off()


# 3rd model: ISC NIRScan NIR range (955 - 1645 nm) ----
db.isc.nir.cal = read.csv(paste0(dir.p, '/tables/db_model/db_source_isc_spec_cal.csv'))[-1]
db.isc.nir.cal[1:5, c(1:4, ncol(db.isc.nir.cal))]
colnames(db.isc.nir.cal) = gsub("NIR.", "", colnames(db.isc.nir.cal))
which(colnames(db.isc.nir.cal) == "1645")
db.isc.nir.cal = db.isc.nir.cal[, c(1, 3:141)]

db.soil.isc.nir.cal = merge(db.soil, db.isc.nir.cal, by="id")
db.soil.isc.nir.cal[1:5, c(1:10, ncol(db.soil.isc.nir.cal))]
str(db.soil.isc.nir.cal[, 1:10])
db.soil.isc.nir.cal$tn.perc = as.numeric(db.soil.isc.nir.cal$tn.perc)
db.soil.isc.nir.cal$tc.perc = as.numeric(db.soil.isc.nir.cal$tc.perc)

## PLSR
colnames(db.soil.isc.nir.cal)[1:20]
db.tc.isc.nir.pls = db.soil.isc.nir.cal[, c(6, 9:ncol(db.soil.isc.nir.cal))]
colnames(db.tc.isc.nir.pls)[1:20]

pls::pls.options(parallel = makeCluster(cores, type = 'PSOCK'))
system.time(tc.isc.nir.pls.m <- pls::plsr(tc.perc~., 
                                          ncomp = 15, 
                                          data = db.tc.isc.nir.pls, 
                                          valid = 'CV', 
                                          scale = F)) 

stopCluster(pls.options()$parallel)

summary(tc.isc.nir.pls.m)
tc.isc.nir.pls.m$validation
ncomp.sel = as.data.frame(tc.isc.nir.pls.m$validation$PRESS)
colnames(ncomp.sel) = gsub("comps", "", names(ncomp.sel))
ncomp = as.numeric(colnames(ncomp.sel[which.min(ncomp.sel[1,])]))

jpeg(paste0(dir.p, '/Results/figs_pls/ncomp/soc_isc_nir_PLSR_Optimal_N_Comp.jpeg'), width = 10, height = 10, units = 'in', res = 300)
par(mar = c(6, 6, 6, 6))
plot(RMSEP(tc.isc.nir.pls.m),
     main = '',
     xlim = c(0, 15),
     xaxt = 'n',
     frame.plot = FALSE,
     type = 'l',
     lty = 1:2,
     col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5),
     xlab =  'Number of components', # expression(paste('Wavelength(cm'^'-1',')')),
     ylab = 'RMSEP',
     las=1, font=2)
axis(1, at = seq(0, 15, by = 1), labels = T)
mtext('SOC SG-0D: isc NIR')
legend("topright",legend=c("Cross-Validation (CV)","Adjusted CV"),
       lty = 1:2,
       col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5),
       bty = 'n',
       x.intersp = 0.25,
       seg.len = 1)
points(ncomp, 0.3211, col = 'red', pch =20)
dev.off()

## Predicting
db.isc.nir.val = read.csv(paste0(dir.p, '/tables/db_model/db_source_isc_spec_val.csv'))[-1]
db.isc.nir.val[1:5, c(1:4, ncol(db.isc.nir.val))]
colnames(db.isc.nir.val) = gsub("NIR.", "", colnames(db.isc.nir.val))
which(colnames(db.isc.nir.val)== "1645")
db.isc.nir.val = db.isc.nir.val[, c(1, 3:141)]

db.soil.val.all$tc.isc.nir.pls = stats::predict(tc.isc.nir.pls.m, 
                                                newdata = db.isc.nir.val[, 2:ncol(db.isc.nir.val)],
                                                ncomp = ncomp)

db.soil.val.all$tc.isc.nir.pls = as.numeric(db.soil.val.all$tc.isc.nir.pls[, 1, 1])

metric.tc.isc.nir.pls = goof(observed = db.soil.val.all$tc.perc, 
                             predicted = db.soil.val.all$tc.isc.nir.pls,
                             type = 'spec')

jpeg(paste0(dir.p, '/Results/figs_pls/soc_sg_isc_nir_PLSR_Obs_Pred.jpeg'), width = 10, height = 10, units = 'in', res = 300)
par(mar = c(6, 6, 6, 6))
summary(db.soil.val.all[, c('tc.isc.nir.pls', 'tc.perc')])
plot(db.soil.val.all$tc.isc.nir.pls, db.soil.val.all$tc.perc,
     xlab = 'Predicted (%)', #expression(paste('Predicted ', '(', 'mg', ' kg'^'-1',')')),
     ylab = 'Observed (%)', #expression(paste('Observed ', '(', 'mg', ' kg'^'-1',')')),
     xlim = c(0, 2),
     ylim = c(0, 2),
     type = 'p',
     pch =  '🔅', # 16,
     bty = 'L',
     cex = 0.9,
     cex.lab= 1.1,
     cex.axis=1.2,
     col = rgb(red=0.0 , green = 0.0, blue = 0.0, alpha = .8),
     family = 'A')
abline(0, 1, col='black', lty=1, lwd=1)
abline(lm(db.soil.val.all$tc.isc.nir.pls ~ db.soil.val.all$tc.perc), col = 'red', lty=2, lwd=1)
text(1.8, 0.5, bquote(RMSE == .(round(metric.tc.isc.nir.pls$RMSE, 2))))
text(1.8, 0.4, bquote(MEC == .(round(metric.tc.isc.nir.pls$MEC, 2))))
text(1.8, 0.3, bquote(CCC == .(round(metric.tc.isc.nir.pls$concordance, 2))))
text(1.8, 0.2, bquote(bias == .(round(metric.tc.isc.nir.pls$bias, 2))))
text(1.8, 0.1, bquote(RPIQ == .(round(metric.tc.isc.nir.pls$RPIQ, 2))))
dev.off()

## Cubist
cores # it has the number of cores from parallel::detectCores()-2
cl = parallel::makeCluster(cores) # Create a cluster using all cores
doParallel::registerDoParallel(cl)

system.time(tc.isc.nir.cub.m <- caret::train(x=db.soil.isc.nir.cal[ ,9:ncol(db.soil.isc.nir.cal)], 
                                             y=db.soil.isc.nir.cal$tc.perc,
                                             method = "cubist", 
                                             trControl = ctrl,
                                             importance = T))  
stopCluster(cl)

print(tc.isc.nir.cub.m)

tc.isc.nir.cub.m.var = varImp(tc.isc.nir.cub.m, scale = T)

jpeg(paste0(dir.p, '/results/figs_cub/varimp/soc_sg_isc_nir_Var_Importance_Cub.jpeg'), width = 10, height = 10, units = 'in', res = 300)
plot(tc.isc.nir.cub.m.var, top = 20, main = "Original data")
dev.off()


## Predicting
db.soil.val.all$tc.isc.nir.cub = stats::predict(tc.isc.nir.cub.m, 
                                                db.isc.nir.val[, 2:ncol(db.isc.nir.val)]) # Predicting using external data

metric.tc.isc.nir.cub = goof(observed = db.soil.val.all$tc.perc, 
                             predicted = db.soil.val.all$tc.isc.nir.cub,
                             type = 'spec')


jpeg(paste0(dir.p, '/Results/figs_cub/soc_sg_isc_nir_Cub_Obs_Pred.jpeg'), width = 10, height = 10, units = 'in', res = 300)
par(mar = c(6, 6, 6, 6))
summary(db.soil.val.all[, c('tc.isc.nir.cub', 'tc.perc')])
plot(db.soil.val.all$tc.isc.nir.cub, db.soil.val.all$tc.perc,
     xlab = 'Predicted (%)', #expression(paste('Predicted ', '(', 'mg', ' kg'^'-1',')')),
     ylab = 'Observed (%)', #expression(paste('Observed ', '(', 'mg', ' kg'^'-1',')')),
     xlim = c(0, 2),
     ylim = c(0, 2),
     type = 'p',
     pch =  '🔅', # 16,
     bty = 'L',
     cex = 0.9,
     cex.lab= 1.1,
     cex.axis=1.2,
     col = rgb(red=0.0 , green = 0.0, blue = 0.0, alpha = .8),
     family = 'A')
abline(0, 1, col='black', lty=1, lwd=1)
abline(lm(db.soil.val.all$tc.isc.nir.cub ~ db.soil.val.all$tc.perc), col = 'red', lty=2, lwd=1)
text(1.8, 0.5, bquote(RMSE == .(round(metric.tc.isc.nir.cub$RMSE, 2))))
text(1.8, 0.4, bquote(MEC == .(round(metric.tc.isc.nir.cub$MEC, 2))))
text(1.8, 0.3, bquote(CCC == .(round(metric.tc.isc.nir.cub$concordance, 2))))
text(1.8, 0.2, bquote(bias == .(round(metric.tc.isc.nir.cub$bias, 2))))
text(1.8, 0.1, bquote(RPIQ == .(round(metric.tc.isc.nir.cub$RPIQ, 2))))
dev.off()


# 4th model: ASD SNV NIR range (955 - 1645 nm) ----
db.asd.snv.nir.cal = read.csv(paste0(dir.p, '/tables/db_model/db_target_asd_spec_snv_cal.csv'))[-1]
db.asd.snv.nir.cal[1:5, c(1:4, ncol(db.asd.snv.nir.cal))]
colnames(db.asd.snv.nir.cal) = gsub("NIR.", "", colnames(db.asd.snv.nir.cal))
which(colnames(db.asd.snv.nir.cal)=="1645")
db.asd.snv.nir.cal = db.asd.snv.nir.cal[, c(1, 3:141)]

db.soil.asd.snv.nir.cal = merge(db.soil, db.asd.snv.nir.cal, by="id")
db.soil.asd.snv.nir.cal[1:5, c(1:10, ncol(db.soil.asd.snv.nir.cal))]
str(db.soil.asd.snv.nir.cal[, 1:10])
db.soil.asd.snv.nir.cal$tn.perc = as.numeric(db.soil.asd.snv.nir.cal$tn.perc)
db.soil.asd.snv.nir.cal$tc.perc = as.numeric(db.soil.asd.snv.nir.cal$tc.perc)

## PLSR
colnames(db.soil.asd.snv.nir.cal)[1:20]
db.tc.asd.snv.nir.pls = db.soil.asd.snv.nir.cal[, c(6, 9:ncol(db.soil.asd.snv.nir.cal))]
colnames(db.tc.asd.snv.nir.pls)[1:20]

pls::pls.options(parallel = makeCluster(cores, type = 'PSOCK'))
system.time(tc.asd.snv.nir.pls.m <- pls::plsr(tc.perc~., 
                                          ncomp = 15, 
                                          data = db.tc.asd.snv.nir.pls, 
                                          valid = 'CV', 
                                          scale = F)) 

stopCluster(pls.options()$parallel)

summary(tc.asd.snv.nir.pls.m)
tc.asd.snv.nir.pls.m$validation
ncomp.sel = as.data.frame(tc.asd.snv.nir.pls.m$validation$PRESS)
colnames(ncomp.sel) = gsub("comps", "", names(ncomp.sel))
ncomp = as.numeric(colnames(ncomp.sel[which.min(ncomp.sel[1,])]))

jpeg(paste0(dir.p, '/Results/figs_pls/ncomp/soc_asd_SNV_nir_PLSR_Optimal_N_Comp.jpeg'), width = 10, height = 10, units = 'in', res = 300)
par(mar = c(6, 6, 6, 6))
plot(RMSEP(tc.asd.snv.nir.pls.m),
     main = '',
     xlim = c(0, 15),
     xaxt = 'n',
     frame.plot = FALSE,
     type = 'l',
     lty = 1:2,
     col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5),
     xlab =  'Number of components', # expression(paste('Wavelength(cm'^'-1',')')),
     ylab = 'RMSEP',
     las=1, font=2)
axis(1, at = seq(0, 15, by = 1), labels = T)
mtext('SOC SG-0D: ASD NIR')
legend("topright",legend=c("Cross-Validation (CV)","Adjusted CV"),
       lty = 1:2,
       col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5),
       bty = 'n',
       x.intersp = 0.25,
       seg.len = 1)
points(ncomp, 0.1616, col = 'red', pch =20)
dev.off()

## Predicting
db.asd.snv.nir.val = read.csv(paste0(dir.p, '/tables/db_model/db_target_asd_spec_snv_val.csv'))[-1]
db.asd.snv.nir.val[1:5, c(1:4, ncol(db.asd.snv.nir.val))]
colnames(db.asd.snv.nir.val) = gsub("NIR.", "", colnames(db.asd.snv.nir.val))
which(colnames(db.asd.snv.nir.val)=="1645")
db.asd.snv.nir.val = db.asd.snv.nir.val[, c(1, 3:141)]

db.soil.val.all$tc.asd.snv.nir.pls = stats::predict(tc.asd.snv.nir.pls.m, 
                                                newdata = db.asd.snv.nir.val[, 2:ncol(db.asd.snv.nir.val)],
                                                ncomp = ncomp)

db.soil.val.all$tc.asd.snv.nir.pls = as.numeric(db.soil.val.all$tc.asd.snv.nir.pls[, 1, 1])

metric.tc.asd.snv.nir.pls = goof(observed = db.soil.val.all$tc.perc, 
                             predicted = db.soil.val.all$tc.asd.snv.nir.pls,
                             type = 'spec')

jpeg(paste0(dir.p, '/Results/figs_pls/soc_SNV_asd_nir_PLSR_Obs_Pred.jpeg'), width = 10, height = 10, units = 'in', res = 300)
par(mar = c(6, 6, 6, 6))
summary(db.soil.val.all[, c('tc.asd.snv.nir.pls', 'tc.perc')])
plot(db.soil.val.all$tc.asd.snv.nir.pls, db.soil.val.all$tc.perc,
     xlab = 'Predicted (%)', #expression(paste('Predicted ', '(', 'mg', ' kg'^'-1',')')),
     ylab = 'Observed (%)', #expression(paste('Observed ', '(', 'mg', ' kg'^'-1',')')),
     xlim = c(0, 2),
     ylim = c(0, 2),
     type = 'p',
     pch =  '🔅', # 16,
     bty = 'L',
     cex = 0.9,
     cex.lab= 1.1,
     cex.axis=1.2,
     col = rgb(red=0.0 , green = 0.0, blue = 0.0, alpha = .8),
     family = 'A')
abline(0, 1, col='black', lty=1, lwd=1)
abline(lm(db.soil.val.all$tc.asd.snv.nir.pls ~ db.soil.val.all$tc.perc), col = 'red', lty=2, lwd=1)
text(1.8, 0.5, bquote(RMSE == .(round(metric.tc.asd.snv.nir.pls$RMSE, 2))))
text(1.8, 0.4, bquote(MEC == .(round(metric.tc.asd.snv.nir.pls$MEC, 2))))
text(1.8, 0.3, bquote(CCC == .(round(metric.tc.asd.snv.nir.pls$concordance, 2))))
text(1.8, 0.2, bquote(bias == .(round(metric.tc.asd.snv.nir.pls$bias, 2))))
text(1.8, 0.1, bquote(RPIQ == .(round(metric.tc.asd.snv.nir.pls$RPIQ, 2))))
dev.off()

## Cubist
cores # it has the number of cores from parallel::detectCores()-2
cl = parallel::makeCluster(cores) # Create a cluster using all cores
doParallel::registerDoParallel(cl)

system.time(tc.asd.snv.nir.cub.m <- caret::train(x=db.soil.asd.snv.nir.cal[ ,9:ncol(db.soil.asd.snv.nir.cal)], 
                                             y=db.soil.asd.snv.nir.cal$tc.perc,
                                             method = "cubist", 
                                             trControl = ctrl,
                                             importance = T))  
stopCluster(cl)

print(tc.asd.snv.nir.cub.m)

tc.asd.snv.nir.cub.m.var = varImp(tc.asd.snv.nir.cub.m, scale = T)

jpeg(paste0(dir.p, '/results/figs_cub/varimp/soc_SNV_asd_nir_Var_Importance_Cub.jpeg'), width = 10, height = 10, units = 'in', res = 300)
plot(tc.asd.snv.nir.cub.m.var, top = 20, main = "Original data")
dev.off()


## Predicting
db.soil.val.all$tc.asd.snv.nir.cub = stats::predict(tc.asd.snv.nir.cub.m, 
                                                db.asd.snv.nir.val[, 2:ncol(db.asd.snv.nir.val)]) # Predicting using external data

metric.tc.asd.snv.nir.cub = goof(observed = db.soil.val.all$tc.perc, 
                             predicted = db.soil.val.all$tc.asd.snv.nir.cub,
                             type = 'spec')


jpeg(paste0(dir.p, '/Results/figs_cub/soc_SNV_asd_nir_Cub_Obs_Pred.jpeg'), width = 10, height = 10, units = 'in', res = 300)
par(mar = c(6, 6, 6, 6))
summary(db.soil.val.all[, c('tc.asd.snv.nir.cub', 'tc.perc')])
plot(db.soil.val.all$tc.asd.snv.nir.cub, db.soil.val.all$tc.perc,
     xlab = 'Predicted (%)', #expression(paste('Predicted ', '(', 'mg', ' kg'^'-1',')')),
     ylab = 'Observed (%)', #expression(paste('Observed ', '(', 'mg', ' kg'^'-1',')')),
     xlim = c(0, 2),
     ylim = c(0, 2),
     type = 'p',
     pch =  '🔅', # 16,
     bty = 'L',
     cex = 0.9,
     cex.lab= 1.1,
     cex.axis=1.2,
     col = rgb(red=0.0 , green = 0.0, blue = 0.0, alpha = .8),
     family = 'A')
abline(0, 1, col='black', lty=1, lwd=1)
abline(lm(db.soil.val.all$tc.asd.snv.nir.cub ~ db.soil.val.all$tc.perc), col = 'red', lty=2, lwd=1)
text(1.8, 0.5, bquote(RMSE == .(round(metric.tc.asd.snv.nir.cub$RMSE, 2))))
text(1.8, 0.4, bquote(MEC == .(round(metric.tc.asd.snv.nir.cub$MEC, 2))))
text(1.8, 0.3, bquote(CCC == .(round(metric.tc.asd.snv.nir.cub$concordance, 2))))
text(1.8, 0.2, bquote(bias == .(round(metric.tc.asd.snv.nir.cub$bias, 2))))
text(1.8, 0.1, bquote(RPIQ == .(round(metric.tc.asd.snv.nir.cub$RPIQ, 2))))
dev.off()


# 5th model: ISC NIRScan SNV NIR range (955 - 1645 nm) ----
db.isc.snv.nir.cal = read.csv(paste0(dir.p, '/tables/db_model/db_source_isc_spec_snv_cal.csv'))[-1]
db.isc.snv.nir.cal[1:5, c(1:4, ncol(db.isc.snv.nir.cal))]
colnames(db.isc.snv.nir.cal) = gsub("NIR.", "", colnames(db.isc.snv.nir.cal))
which(colnames(db.isc.snv.nir.cal)=="1645")
db.isc.snv.nir.cal = db.isc.snv.nir.cal[, c(1, 3:141)]

db.soil.isc.snv.nir.cal = merge(db.soil, db.isc.snv.nir.cal, by="id")
db.soil.isc.snv.nir.cal[1:5, c(1:10, ncol(db.soil.isc.snv.nir.cal))]
str(db.soil.isc.snv.nir.cal[, 1:10])
db.soil.isc.snv.nir.cal$tn.perc = as.numeric(db.soil.isc.snv.nir.cal$tn.perc)
db.soil.isc.snv.nir.cal$tc.perc = as.numeric(db.soil.isc.snv.nir.cal$tc.perc)

## PLSR
colnames(db.soil.isc.snv.nir.cal)[1:20]
db.tc.isc.snv.nir.pls = db.soil.isc.snv.nir.cal[, c(6, 9:ncol(db.soil.isc.snv.nir.cal))]
colnames(db.tc.isc.snv.nir.pls)[1:20]

pls::pls.options(parallel = makeCluster(cores, type = 'PSOCK'))
system.time(tc.isc.snv.nir.pls.m <- pls::plsr(tc.perc~., 
                                              ncomp = 15, 
                                              data = db.tc.isc.snv.nir.pls, 
                                              valid = 'CV', 
                                              scale = F)) 

stopCluster(pls.options()$parallel)

summary(tc.isc.snv.nir.pls.m)
tc.isc.snv.nir.pls.m$validation
ncomp.sel = as.data.frame(tc.isc.snv.nir.pls.m$validation$PRESS)
colnames(ncomp.sel) = gsub("comps", "", names(ncomp.sel))
ncomp = as.numeric(colnames(ncomp.sel[which.min(ncomp.sel[1,])]))

jpeg(paste0(dir.p, '/Results/figs_pls/ncomp/soc_isc_SNV_nir_PLSR_Optimal_N_Comp.jpeg'), width = 10, height = 10, units = 'in', res = 300)
par(mar = c(6, 6, 6, 6))
plot(RMSEP(tc.isc.snv.nir.pls.m),
     main = '',
     xlim = c(0, 15),
     xaxt = 'n',
     frame.plot = FALSE,
     type = 'l',
     lty = 1:2,
     col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5),
     xlab =  'Number of components', # expression(paste('Wavelength(cm'^'-1',')')),
     ylab = 'RMSEP',
     las=1, font=2)
axis(1, at = seq(0, 15, by = 1), labels = T)
mtext('SOC SG-0D: isc NIR')
legend("topright",legend=c("Cross-Validation (CV)","Adjusted CV"),
       lty = 1:2,
       col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5),
       bty = 'n',
       x.intersp = 0.25,
       seg.len = 1)
points(ncomp, 0.1616, col = 'red', pch =20)
dev.off()

## Predicting
db.isc.snv.nir.val = read.csv(paste0(dir.p, '/tables/db_model/db_source_isc_spec_snv_val.csv'))[-1]
db.isc.snv.nir.val[1:5, c(1:4, ncol(db.isc.snv.nir.val))]
colnames(db.isc.snv.nir.val) = gsub("NIR.", "", colnames(db.isc.snv.nir.val))
which(colnames(db.isc.snv.nir.val)=="1645")
db.isc.snv.nir.val = db.isc.snv.nir.val[, c(1, 3:141)]

db.soil.val.all$tc.isc.snv.nir.pls = stats::predict(tc.isc.snv.nir.pls.m, 
                                                    newdata = db.isc.snv.nir.val[, 2:ncol(db.isc.snv.nir.val)],
                                                    ncomp = ncomp)

db.soil.val.all$tc.isc.snv.nir.pls = as.numeric(db.soil.val.all$tc.isc.snv.nir.pls[, 1, 1])

metric.tc.isc.snv.nir.pls = goof(observed = db.soil.val.all$tc.perc, 
                                 predicted = db.soil.val.all$tc.isc.snv.nir.pls,
                                 type = 'spec')

jpeg(paste0(dir.p, '/Results/figs_pls/soc_SNV_isc_nir_PLSR_Obs_Pred.jpeg'), width = 10, height = 10, units = 'in', res = 300)
par(mar = c(6, 6, 6, 6))
summary(db.soil.val.all[, c('tc.isc.snv.nir.pls', 'tc.perc')])
plot(db.soil.val.all$tc.isc.snv.nir.pls, db.soil.val.all$tc.perc,
     xlab = 'Predicted (%)', #expression(paste('Predicted ', '(', 'mg', ' kg'^'-1',')')),
     ylab = 'Observed (%)', #expression(paste('Observed ', '(', 'mg', ' kg'^'-1',')')),
     xlim = c(0, 2),
     ylim = c(0, 2),
     type = 'p',
     pch =  '🔅', # 16,
     bty = 'L',
     cex = 0.9,
     cex.lab= 1.1,
     cex.axis=1.2,
     col = rgb(red=0.0 , green = 0.0, blue = 0.0, alpha = .8),
     family = 'A')
abline(0, 1, col='black', lty=1, lwd=1)
abline(lm(db.soil.val.all$tc.isc.snv.nir.pls ~ db.soil.val.all$tc.perc), col = 'red', lty=2, lwd=1)
text(1.8, 0.5, bquote(RMSE == .(round(metric.tc.isc.snv.nir.pls$RMSE, 2))))
text(1.8, 0.4, bquote(MEC == .(round(metric.tc.isc.snv.nir.pls$MEC, 2))))
text(1.8, 0.3, bquote(CCC == .(round(metric.tc.isc.snv.nir.pls$concordance, 2))))
text(1.8, 0.2, bquote(bias == .(round(metric.tc.isc.snv.nir.pls$bias, 2))))
text(1.8, 0.1, bquote(RPIQ == .(round(metric.tc.isc.snv.nir.pls$RPIQ, 2))))
dev.off()

## Cubist
cores # it has the number of cores from parallel::detectCores()-2
cl = parallel::makeCluster(cores) # Create a cluster using all cores
doParallel::registerDoParallel(cl)

system.time(tc.isc.snv.nir.cub.m <- caret::train(x=db.soil.isc.snv.nir.cal[ ,9:ncol(db.soil.isc.snv.nir.cal)], 
                                                 y=db.soil.isc.snv.nir.cal$tc.perc,
                                                 method = "cubist", 
                                                 trControl = ctrl,
                                                 importance = T))  
stopCluster(cl)

print(tc.isc.snv.nir.cub.m)

tc.isc.snv.nir.cub.m.var = varImp(tc.isc.snv.nir.cub.m, scale = T)

jpeg(paste0(dir.p, '/results/figs_cub/varimp/soc_SNV_isc_nir_Var_Importance_Cub.jpeg'), width = 10, height = 10, units = 'in', res = 300)
plot(tc.isc.snv.nir.cub.m.var, top = 20, main = "Original data")
dev.off()


## Predicting
db.soil.val.all$tc.isc.snv.nir.cub = stats::predict(tc.isc.snv.nir.cub.m, 
                                                    db.isc.snv.nir.val[, 2:ncol(db.isc.snv.nir.val)]) # Predicting using external data

metric.tc.isc.snv.nir.cub = goof(observed = db.soil.val.all$tc.perc, 
                                 predicted = db.soil.val.all$tc.isc.snv.nir.cub,
                                 type = 'spec')


jpeg(paste0(dir.p, '/Results/figs_cub/soc_SNV_isc_nir_Cub_Obs_Pred.jpeg'), width = 10, height = 10, units = 'in', res = 300)
par(mar = c(6, 6, 6, 6))
summary(db.soil.val.all[, c('tc.isc.snv.nir.cub', 'tc.perc')])
plot(db.soil.val.all$tc.isc.snv.nir.cub, db.soil.val.all$tc.perc,
     xlab = 'Predicted (%)', #expression(paste('Predicted ', '(', 'mg', ' kg'^'-1',')')),
     ylab = 'Observed (%)', #expression(paste('Observed ', '(', 'mg', ' kg'^'-1',')')),
     xlim = c(0, 2),
     ylim = c(0, 2),
     type = 'p',
     pch =  '🔅', # 16,
     bty = 'L',
     cex = 0.9,
     cex.lab= 1.1,
     cex.axis=1.2,
     col = rgb(red=0.0 , green = 0.0, blue = 0.0, alpha = .8),
     family = 'A')
abline(0, 1, col='black', lty=1, lwd=1)
abline(lm(db.soil.val.all$tc.isc.snv.nir.cub ~ db.soil.val.all$tc.perc), col = 'red', lty=2, lwd=1)
text(1.8, 0.5, bquote(RMSE == .(round(metric.tc.isc.snv.nir.cub$RMSE, 2))))
text(1.8, 0.4, bquote(MEC == .(round(metric.tc.isc.snv.nir.cub$MEC, 2))))
text(1.8, 0.3, bquote(CCC == .(round(metric.tc.isc.snv.nir.cub$concordance, 2))))
text(1.8, 0.2, bquote(bias == .(round(metric.tc.isc.snv.nir.cub$bias, 2))))
text(1.8, 0.1, bquote(RPIQ == .(round(metric.tc.isc.snv.nir.cub$RPIQ, 2))))
dev.off()


# Predicting Transformed validation spectra within the best models ----
## Substracted correction SG
db.isc.nir.subcorr.val = read.csv(paste0(dir.p, '/tables/db_model/db_source_spec_subcorr_val.csv'))[-1]
db.isc.nir.subcorr.val[1:5, c(1:4, ncol(db.isc.nir.subcorr.val))]
colnames(db.isc.nir.subcorr.val) = gsub("NIR.", "", colnames(db.isc.nir.subcorr.val))
which(colnames(db.isc.nir.subcorr.val)=="1645")
db.isc.nir.subcorr.val = db.isc.nir.subcorr.val[, c(1, 3:141)]

### 2nd model PLSR was the best fit with 15 pc.
db.soil.val.all$tc.isc.nir.subcorr.pls = stats::predict(tc.asd.nir.pls.m,
                                                        newdata = db.isc.nir.subcorr.val[, 2:ncol(db.isc.nir.subcorr.val)],
                                                        ncomp = 15)

db.soil.val.all$tc.isc.nir.subcorr.pls = as.numeric(db.soil.val.all$tc.isc.nir.subcorr.pls[, 1, 1])

metric.tc.isc.nir.subcorr.pls = goof(observed = db.soil.val.all$tc.perc, 
                                 predicted = db.soil.val.all$tc.isc.nir.subcorr.pls,
                                 type = 'spec')

jpeg(paste0(dir.p, '/Results/figs_pls/soc_isc_nir_SUBCORRECTED_PLSR_Obs_Pred.jpeg'), width = 10, height = 10, units = 'in', res = 300)
par(mar = c(6, 6, 6, 6))
summary(db.soil.val.all[, c('tc.isc.nir.subcorr.pls', 'tc.perc')])
plot(db.soil.val.all$tc.isc.nir.subcorr.pls, db.soil.val.all$tc.perc,
     xlab = 'Predicted (%)', #expression(paste('Predicted ', '(', 'mg', ' kg'^'-1',')')),
     ylab = 'Observed (%)', #expression(paste('Observed ', '(', 'mg', ' kg'^'-1',')')),
     xlim = c(0, 3),
     ylim = c(0, 3),
     type = 'p',
     pch =  '🔅', # 16,
     bty = 'L',
     cex = 0.9,
     cex.lab= 1.1,
     cex.axis=1.2,
     col = rgb(red=0.0 , green = 0.0, blue = 0.0, alpha = .8),
     family = 'A')
abline(0, 1, col='black', lty=1, lwd=1)
abline(lm(db.soil.val.all$tc.isc.nir.subcorr.pls ~ db.soil.val.all$tc.perc), col = 'red', lty=2, lwd=1)
text(0.4, 2.9, bquote(RMSE == .(round(metric.tc.isc.nir.subcorr.pls$RMSE, 2))))
text(0.4, 2.8, bquote(MEC == .(round(metric.tc.isc.nir.subcorr.pls$MEC, 2))))
text(0.4, 2.7, bquote(CCC == .(round(metric.tc.isc.nir.subcorr.pls$concordance, 2))))
text(0.4, 2.6, bquote(bias == .(round(metric.tc.isc.nir.subcorr.pls$bias, 2))))
text(0.4, 2.5, bquote(RPIQ == .(round(metric.tc.isc.nir.subcorr.pls$RPIQ, 2))))
dev.off()

## Substracted correction SG-SNV
db.isc.nir.snv.subcorr.val = read.csv(paste0(dir.p, '/tables/db_model/db_source_spec_snv_subcorr_val.csv'))[-1]
db.isc.nir.snv.subcorr.val[1:5, c(1:4, ncol(db.isc.nir.snv.subcorr.val))]
colnames(db.isc.nir.snv.subcorr.val) = gsub("NIR.", "", colnames(db.isc.nir.snv.subcorr.val))
which(colnames(db.isc.nir.snv.subcorr.val)=="1645")
db.isc.nir.snv.subcorr.val = db.isc.nir.snv.subcorr.val[, c(1, 3:141)]

### 4th model Cubist was the best fit
db.soil.val.all$tc.isc.nir.snv.subcorr.cub = stats::predict(tc.asd.snv.nir.cub.m,
                                                        newdata = db.isc.nir.snv.subcorr.val[, 2:ncol(db.isc.nir.snv.subcorr.val)])

metric.tc.isc.nir.snv.subcorr.cub = goof(observed = db.soil.val.all$tc.perc, 
                                     predicted = db.soil.val.all$tc.isc.nir.snv.subcorr.cub,
                                     type = 'spec')

jpeg(paste0(dir.p, '/Results/figs_cub/soc_isc_nir_snv_SUBCORRECTED_CUB_Obs_Pred.jpeg'), width = 10, height = 10, units = 'in', res = 300)
par(mar = c(6, 6, 6, 6))
summary(db.soil.val.all[, c('tc.isc.nir.snv.subcorr.cub', 'tc.perc')])
plot(db.soil.val.all$tc.isc.nir.snv.subcorr.cub, db.soil.val.all$tc.perc,
     xlab = 'Predicted (%)', #expression(paste('Predicted ', '(', 'mg', ' kg'^'-1',')')),
     ylab = 'Observed (%)', #expression(paste('Observed ', '(', 'mg', ' kg'^'-1',')')),
     xlim = c(0, 3),
     ylim = c(0, 3),
     type = 'p',
     pch =  '🔅', # 16,
     bty = 'L',
     cex = 0.9,
     cex.lab= 1.1,
     cex.axis=1.2,
     col = rgb(red=0.0 , green = 0.0, blue = 0.0, alpha = .8),
     family = 'A')
abline(0, 1, col='black', lty=1, lwd=1)
abline(lm(db.soil.val.all$tc.isc.nir.snv.subcorr.cub ~ db.soil.val.all$tc.perc), col = 'red', lty=2, lwd=1)
text(0.4, 2.9, bquote(RMSE == .(round(metric.tc.isc.nir.snv.subcorr.cub$RMSE, 2))))
text(0.4, 2.8, bquote(MEC == .(round(metric.tc.isc.nir.snv.subcorr.cub$MEC, 2))))
text(0.4, 2.7, bquote(CCC == .(round(metric.tc.isc.nir.snv.subcorr.cub$concordance, 2))))
text(0.4, 2.6, bquote(bias == .(round(metric.tc.isc.nir.snv.subcorr.cub$bias, 2))))
text(0.4, 2.5, bquote(RPIQ == .(round(metric.tc.isc.nir.snv.subcorr.cub$RPIQ, 2))))
dev.off()

## PDS correction SG-SNV
db.isc.nir.snv.pdscorr.val = read.csv(paste0(dir.p, '/tables/db_model/db_source_spec_snv_pdscorr_val.csv'))[-1]
db.isc.nir.snv.pdscorr.val[1:5, c(1:4, ncol(db.isc.nir.snv.pdscorr.val))]
colnames(db.isc.nir.snv.pdscorr.val) = gsub("NIR.", "", colnames(db.isc.nir.snv.pdscorr.val))

### 4th model Cubist was the best fit
db.soil.val.all$tc.isc.nir.snv.pdscorr.cub = stats::predict(tc.asd.snv.nir.cub.m,
                                                            newdata = db.isc.nir.snv.pdscorr.val[, 2:ncol(db.isc.nir.snv.pdscorr.val)])

metric.tc.isc.nir.snv.pdscorr.cub = goof(observed = db.soil.val.all$tc.perc, 
                                         predicted = db.soil.val.all$tc.isc.nir.snv.pdscorr.cub,
                                         type = 'spec')

jpeg(paste0(dir.p, '/Results/figs_cub/soc_isc_nir_snv_PDSCORRECTED_CUB_Obs_Pred.jpeg'), width = 10, height = 10, units = 'in', res = 300)
par(mar = c(6, 6, 6, 6))
summary(db.soil.val.all[, c('tc.isc.nir.snv.pdscorr.cub', 'tc.perc')])
plot(db.soil.val.all$tc.isc.nir.snv.pdscorr.cub, db.soil.val.all$tc.perc,
     xlab = 'Predicted (%)', #expression(paste('Predicted ', '(', 'mg', ' kg'^'-1',')')),
     ylab = 'Observed (%)', #expression(paste('Observed ', '(', 'mg', ' kg'^'-1',')')),
     xlim = c(0, 3),
     ylim = c(0, 3),
     type = 'p',
     pch =  '🔅', # 16,
     bty = 'L',
     cex = 0.9,
     cex.lab= 1.1,
     cex.axis=1.2,
     col = rgb(red=0.0 , green = 0.0, blue = 0.0, alpha = .8),
     family = 'A')
abline(0, 1, col='black', lty=1, lwd=1)
abline(lm(db.soil.val.all$tc.isc.nir.snv.pdscorr.cub ~ db.soil.val.all$tc.perc), col = 'red', lty=2, lwd=1)
text(0.4, 2.9, bquote(RMSE == .(round(metric.tc.isc.nir.snv.pdscorr.cub$RMSE, 2))))
text(0.4, 2.8, bquote(MEC == .(round(metric.tc.isc.nir.snv.pdscorr.cub$MEC, 2))))
text(0.4, 2.7, bquote(CCC == .(round(metric.tc.isc.nir.snv.pdscorr.cub$concordance, 2))))
text(0.4, 2.6, bquote(bias == .(round(metric.tc.isc.nir.snv.pdscorr.cub$bias, 2))))
text(0.4, 2.5, bquote(RPIQ == .(round(metric.tc.isc.nir.snv.pdscorr.cub$RPIQ, 2))))
dev.off()

# Rounding table -----
# Round the specified columns in the dataframe
colnames(db.soil.val.all)
db.soil.val.all[, 5:ncol(db.soil.val.all)] = round(db.soil.val.all[, 5:ncol(db.soil.val.all)], 2)
write.csv(db.soil.val.all, paste0(dir.p, '/tables/db_carbon_nitrogen_val.csv'))
# Saving R workspace ----
save.image(file = paste0(dir.p, '/r_code/proj_low_cost_nir_modelling.RData'))
