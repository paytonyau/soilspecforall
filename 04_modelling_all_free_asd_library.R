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

# Importing more ASD data and testing without the transferred data -----
## Paulinenaue / Peatlands
asd.paul = read.csv(paste0(dir.p, "/raw_data/paulinenaue/db_soil_spec_epsg32633.csv"), header = T, sep = ",")[-1]
head(asd.paul[, 1:10], 5) 

colnames(asd.paul)[1:42]
asd.paul = asd.paul[, c(1, 11:12, 42:ncol(asd.paul))]
head(asd.paul[, c(1:10, ncol(asd.paul))], 5) # how it should look like

asd.paul = asd.paul[complete.cases(asd.paul), ] # 262 obs
colnames(asd.paul) = gsub("X", "", names(asd.paul))

matplot(as.numeric(colnames(asd.paul)[4:ncol(asd.paul)]), 
        t(asd.paul[,4:ncol(asd.paul)]),
        xlim = c(350, 2500),
        xaxt = 'n',
        frame.plot = FALSE,
        type = 'l',
        lty = 1,
        col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5),
        xlab =  'nm', #expression(paste('Wavelength(cm'^'-1',')'))
        ylab = 'Reflectance',
        las=1, font=2)
axis(1, at = seq(350, 2500, by = 50), labels = T, las = 2)
mtext('ASD Spectra: Peatlands')
grid(lwd=1, nx=50, ny=50)

matplot(as.numeric(colnames(asd.vnir.osc)[4:ncol(asd.vnir.osc)]), 
        t(asd.vnir.osc[,4:ncol(asd.vnir.osc)]),
        xlim = c(350, 2500),
        xaxt = 'n',
        frame.plot = FALSE,
        type = 'l',
        lty = 1,
        col = rgb(red = 1.0, green = 0.5, blue = 0.5, alpha = 0.5),
        xlab =  'nm', #expression(paste('Wavelength(cm'^'-1',')'))
        ylab = 'Reflectance',
        las=1, font=2,
        add = T)

## Merging datasets for preprocessing ----
asd.paul[1:5, c(1:10, ncol(asd.paul))]
colnames(asd.paul)[1] = "id"
asd.paul = asd.paul[, c(1, 3, 2, 4:ncol(asd.paul))]
colnames(asd.paul)[2:3] = c("tn.perc", "tc.perc")
asd.paul$tn.perc = asd.paul$tn.perc/10
nrow(asd.paul)

asd.osc[1:5, c(1:10, ncol(asd.osc))]
nrow(asd.osc)

asd.chem.de = rbind(asd.paul, asd.osc)
asd.chem.de[1:5, c(1:10, ncol(asd.chem.de))]

## Splice Correction applicable for ASD ----
asd.chem.de[1:5, c(1:10, ncol(asd.chem.de))]

db.asd.paul.spl = asd.paul
colnames(db.asd.paul.spl)[c(1:4, ncol(db.asd.paul.spl))]
head(db.asd.paul.spl[, c(1:4, ncol(db.asd.paul.spl))])
rownames(db.asd.paul.spl) = db.asd.paul.spl$id

class(db.asd.paul.spl)
colnames(db.asd.paul.spl)[1:10]
wav.asd = as.numeric(colnames(db.asd.paul.spl[, 4:ncol(db.asd.paul.spl)]))
db.asd.paul.spl.m = as.matrix(db.asd.paul.spl[, 4:ncol(db.asd.paul.spl)]) # assembling matrix
class(db.asd.paul.spl.m) # checking class

db.asd.paul.spl.c = prospectr::spliceCorrection(db.asd.paul.spl.m, wav.asd, 
                                               splice = c(1000, 1800), interpol.bands = 10)
class(db.asd.paul.spl.c)
colnames(db.asd.paul.spl.c)
head(db.asd.paul.spl.c[1:5, 1:5])

## Resampling from 1nm to 5 nm 
db.asd.paul.spl.c[1:5, c(1, ncol(db.asd.paul.spl.c))]
db.asd.paul.r = prospectr::resample(db.asd.paul.spl.c, 
                                    wav.asd, 
                                    seq(435, 2415, 5), 
                                    interpol = "spline") # window size of 5 nm


class(db.asd.paul.r)

## Savitzky-Golay
db.asd.paul.sg = prospectr::savitzkyGolay(db.asd.paul.r[, 1:ncol(db.asd.paul.r)],
                                         m = 0, # differentiation order (derivative)
                                         p = 2, # polynomial order
                                         w = 11 # window size
) 
class(db.asd.paul.sg)
head(db.asd.paul.sg[1:5, 1:5])
which(colnames(db.asd.paul.sg) == "950") 
which(colnames(db.asd.paul.sg) == "1650") 
db.asd.paul.sg = db.asd.paul.sg[, 99:239]
db.asd.paul.sg[1:5, c(1, ncol(db.asd.paul.sg))]

## SNV - Select the range before perform it!!!
db.asd.paul.snv = prospectr::standardNormalVariate(db.asd.paul.sg)

### Transforming to DF 
db.asd.paul.d = as.data.frame(db.asd.paul.snv) # converting back to DF
db.asd.paul.d[1:20, c(1:5, ncol(db.asd.paul.d))]
db.asd.paul.d$id = rownames(db.asd.paul.d)
db.asd.paul.d = db.asd.paul.d %>%
  dplyr::select(c('id'), everything())
head(db.asd.paul.d[1:5, c(1:2, ncol(db.asd.paul.d))])

### filtering from 950 - 1650 nm
#which(colnames(db.asd.paul.d) == "950") 
#which(colnames(db.asd.paul.d) == "1650") 
#colnames(db.asd.paul.d)[c(1,100, 240)]
#db.asd.paul.d = db.asd.paul.d[, c(1, 100:240)]
head(db.asd.paul.d[1:5, c(1:2, ncol(db.asd.paul.d))])

### Merging with soil data
#asd.chem.de
asd.paul[1:5, c(1:5, ncol(asd.paul))]
nrow(asd.paul) #262 obs
#head(asd.chem.de[1:5, c(1:4, ncol(asd.chem.de))])
head(db.asd.paul.d[1:5, c(1:2, ncol(db.asd.paul.d))])
db.soil.asd.nir = merge(asd.paul[, 1:3], db.asd.paul.d, by="id")
head(db.soil.asd.nir[1:5, c(1:5, ncol(db.soil.asd.nir))])

#jpeg(paste0(dir.p, '/results/figs/splice_asd_osc_vnir.jpeg'), width = 10, height = 10, units = 'in', res = 300)
matplot(as.numeric(colnames(db.soil.asd.nir)[4:ncol(db.soil.asd.nir)]), 
        t(db.soil.asd.nir[,4:ncol(db.soil.asd.nir)]),
        xlim = c(950, 1650),
        ylim = c(0,1),
        xaxt = 'n',
        frame.plot = FALSE,
        type = 'l',
        lty = 1,
        col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5),
        xlab =  'nm', #expression(paste('Wavelength(cm'^'-1',')'))
        ylab = 'Reflectance',
        las=1, font=2)
axis(1, at = seq(950, 1650, by = 10), labels = T, las = 2)
mtext('ASD Spectra: SG-SNV (955 - 1645 nm)')
grid(lwd=1, nx=50, ny=50)
#dev.off()

## seeing where the target fits in 
asd.target = read.csv(paste0(dir.p, "/tables/db_model/db_target_asd_spec_cal.csv"), header = T, sep = ",")[-1]
head(asd.target[1:5, c(1:3, ncol(asd.target))])
colnames(asd.target) = gsub("NIR.", "", names(asd.target))
#asd.target = asd.target[, c(1, 3:141)]

matplot(as.numeric(colnames(asd.target)[2:ncol(asd.target)]), 
        t(asd.target[,2:ncol(asd.target)]),
        #xlim = c(400, 2450),
        #ylim = c(0,0.6),
        xaxt = 'n',
        frame.plot = FALSE,
        type = 'l',
        lty = 1,
        lwd = 1,
        col = rgb(red = 0.13, green = 0.55, blue = 0.13, alpha = 0.5),
        las=1, font=2,
        add = T)


# Adding local data to regional model -----
## Oscar 
db.soil = read.csv(paste0(dir.p, '/tables/db_model/db_soil_oscar.csv'))
db.soil
db.soil$id = sprintf("osc_%02d", db.soil$id)
db.soil = db.soil[1:102,]
db.soil = db.soil[!(db.soil$tn.perc == "#N/A" & db.soil$tc.perc == "#N/A"), ]
db.soil[1:5, 1:8]

head(asd.target[1:5, c(1:3, ncol(asd.target))])

asd.soil.nir.osc = merge(db.soil[, c(1, 5:6)], asd.target, by="id")

#which(colnames(asd.vnir.osc) == "435") 
#which(colnames(asd.vnir.osc) == "2415") 
asd.soil.nir.osc[1:5, c(1:10, ncol(asd.soil.nir.osc))]
nrow(asd.soil.nir.osc)
ncol(asd.soil.nir.osc)
db.soil.asd.nir[1:5, c(1:10, ncol(db.soil.asd.nir))]
nrow(db.soil.asd.nir)
ncol(db.soil.asd.nir)

asd.chem.de = rbind(db.soil.asd.nir, asd.soil.nir.osc)

write.csv(asd.chem.de, paste0(dir.p, "/tables/db_model/db_soil_germany_asd_nir_snv.csv"))

### Exporting the table ----
write.csv(db.soil.asd.nir, paste0(dir.p, '/tables/db_soil_asd_peat_nir.csv'))

# OSSSL model: Testing ASD large library SNV NIR range (955 - 1645 nm) ----
str(db.soil.asd.nir[, 1:10])
#colnames(db.soil.asd.nir)[2] = "tc.perc"

ctrl = caret::trainControl(method = 'repeatedcv', number = 10, repeats = 10, search = 'random')
cores = detectCores()-2
#asd.chem.de
#asd.chem.de[1:5, c(1:10, ncol(asd.chem.de))]
class(asd.chem.de)
str(asd.chem.de)
asd.chem.de$tc.perc = as.numeric(asd.chem.de$tc.perc)
## Cubist
cores # it has the number of cores from parallel::detectCores()-2
cl = parallel::makeCluster(cores) # Create a cluster using all cores
doParallel::registerDoParallel(cl)

system.time(ossl.tc.m <- caret::train(x=asd.chem.de[ ,4:ncol(asd.chem.de)], #ncol(db.soil.asd.nir)
                                                 y=asd.chem.de$tc.perc,
                                                 method = "cubist", 
                                                 trControl = ctrl,
                                                 importance = T))  
stopCluster(cl)

print(ossl.tc.m)

ossl.tc.m.var = varImp(ossl.tc.m, scale = T)

#jpeg(paste0(dir.p, '/results/figs_cub/varimp/soc_SNV_isc_nir_Var_Importance_Cub.jpeg'), width = 10, height = 10, units = 'in', res = 300)
plot(ossl.tc.m.var, top = 20, main = "Original data")
#dev.off()

## Predicting
### ASD data
db.soil.val.all = read.csv(paste0(dir.p, "/tables/db_soil_asd_nir_val.csv"), header = T, sep = ",")[-1] 
str(db.soil.val.all)

db.asd.nir.val = read.csv(paste0(dir.p, "/tables/db_model/db_target_asd_spec_val.csv"), header = T, sep = ",")[-1] 
str(db.asd.nir.val)
colnames(db.asd.nir.val) = gsub("NIR.", "", names(db.asd.nir.val))
db.asd.nir.val[1:5, c(1:5, ncol(db.asd.nir.val))]
db.soil.val.all$tc.asd = stats::predict(ossl.tc.m,
                                       db.asd.nir.val[, 2:ncol(db.asd.nir.val)]) # Predicting using external data / ncol(db.asd.nir.val)

db.soil.val.all$tc.asd = round(db.soil.val.all$tc.asd, 2)

asd.nir.val.metrics = goof(observed = db.soil.val.all$tc.perc, 
                           predicted = db.soil.val.all$tc.asd,
                           type = 'spec')


#jpeg(paste0(dir.p, '/Results/figs_cub/soc_SNV_isc_nir_Cub_Obs_Pred.jpeg'), width = 10, height = 10, units = 'in', res = 300)
par(mar = c(6, 6, 6, 6))
summary(db.soil.val.all[, c('tc.asd', 'tc.perc')])
plot(db.soil.val.all$tc.asd, db.soil.val.all$tc.perc,
     xlab = 'Predicted (%)', #expression(paste('Predicted ', '(', 'mg', ' kg'^'-1',')')),
     ylab = 'Observed (%)', #expression(paste('Observed ', '(', 'mg', ' kg'^'-1',')')),
     xlim = c(0, 30),
     ylim = c(0, 30),
     type = 'p',
     pch =  '🔅', # 16,
     bty = 'L',
     cex = 0.9,
     cex.lab= 1.1,
     cex.axis=1.2,
     col = rgb(red=0.0 , green = 0.0, blue = 0.0, alpha = .8),
     family = 'A')
abline(0, 1, col='black', lty=1, lwd=1)
abline(lm(db.soil.val.all$tc.asd ~ db.soil.val.all$tc.perc), col = 'red', lty=2, lwd=1)
#text(1.8, 0.5, bquote(RMSE == .(round(asd.nir.val.metrics$RMSE, 2))))
#text(1.8, 0.4, bquote(MEC == .(round(asd.nir.val.metrics$MEC, 2))))
#text(1.8, 0.3, bquote(CCC == .(round(asd.nir.val.metrics$concordance, 2))))
#text(1.8, 0.2, bquote(bias == .(round(asd.nir.val.metrics$bias, 2))))
#text(1.8, 0.1, bquote(RPIQ == .(round(asd.nir.val.metrics$RPIQ, 2))))
#dev.off()

### ISC data
db.isc.nir.val = read.csv(paste0(dir.p, "/tables/db_model/db_source_spec_subcorr_val.csv"), header = T, sep = ",")[-1] 
str(db.isc.nir.val)
colnames(db.isc.nir.val) = gsub("NIR.", "", names(db.isc.nir.val))
db.isc.nir.val[1:5, c(1:5, ncol(db.isc.nir.val))]

db.soil.val.all$tc.isc = stats::predict(ossl.tc.m,
                                        db.isc.nir.val[, 2:ncol(db.isc.nir.val)]) # Predicting using external data

db.soil.val.all$tc.isc = round(db.soil.val.all$tc.isc, 2)

isc.nir.val.metrics = goof(observed = db.soil.val.all$tc.perc, 
                           predicted = db.soil.val.all$tc.isc,
                           type = 'spec')


#jpeg(paste0(dir.p, '/Results/figs_cub/soc_SNV_isc_nir_Cub_Obs_Pred.jpeg'), width = 10, height = 10, units = 'in', res = 300)
par(mar = c(6, 6, 6, 6))
summary(db.soil.val.all[, c('tc.isc', 'tc.perc')])
plot(db.soil.val.all$tc.asd, db.soil.val.all$tc.perc,
     xlab = 'Predicted (%)', #expression(paste('Predicted ', '(', 'mg', ' kg'^'-1',')')),
     ylab = 'Observed (%)', #expression(paste('Observed ', '(', 'mg', ' kg'^'-1',')')),
     xlim = c(0, 10),
     ylim = c(0, 10),
     type = 'p',
     pch =  '🔅', # 16,
     bty = 'L',
     cex = 0.9,
     cex.lab= 1.1,
     cex.axis=1.2,
     col = rgb(red=0.0 , green = 0.0, blue = 0.0, alpha = .8),
     family = 'A')
abline(0, 1, col='black', lty=1, lwd=1)
abline(lm(db.soil.val.all$tc.asd ~ db.soil.val.all$tc.perc), col = 'red', lty=2, lwd=1)
#text(1.8, 0.5, bquote(RMSE == .(round(asd.nir.val.metrics$RMSE, 2))))
#text(1.8, 0.4, bquote(MEC == .(round(asd.nir.val.metrics$MEC, 2))))
#text(1.8, 0.3, bquote(CCC == .(round(asd.nir.val.metrics$concordance, 2))))
#text(1.8, 0.2, bquote(bias == .(round(asd.nir.val.metrics$bias, 2))))
#text(1.8, 0.1, bquote(RPIQ == .(round(asd.nir.val.metrics$RPIQ, 2))))
#dev.off()

matplot(as.numeric(colnames(asd.chem.de)[4:ncol(asd.chem.de)]), 
        t(asd.chem.de[,4:ncol(asd.chem.de)]),
        xlim = c(950, 1650),
        ylim = c(0,1),
        xaxt = 'n',
        frame.plot = FALSE,
        type = 'l',
        lty = 1,
        col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5),
        xlab =  'nm', #expression(paste('Wavelength(cm'^'-1',')'))
        ylab = 'Reflectance',
        las=1, font=2)
axis(1, at = seq(950, 1650, by = 10), labels = T, las = 2)
mtext('ASD Spectra: SG-SNV (955 - 1645 nm)')
grid(lwd=1, nx=50, ny=50)
#dev.off()

## seeing where the source fits in 

matplot(as.numeric(colnames(db.isc.nir.val)[2:ncol(db.isc.nir.val)]), 
        t(db.isc.nir.val[,2:ncol(db.isc.nir.val)]),
        #xlim = c(400, 2450),
        #ylim = c(0,0.6),
        xaxt = 'n',
        frame.plot = FALSE,
        type = 'l',
        lty = 1,
        lwd = 1,
        col = rgb(red = 0.13, green = 0.55, blue = 0.13, alpha = 0.5),
        las=1, font=2,
        add = T)


# Saving R workspace ----
save.image(file = paste0(dir.p, '/r_code/modelling_all_free_asd_library.RData'))
