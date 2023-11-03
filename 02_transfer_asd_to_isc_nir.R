# Libraries ----
if(!require(dplyr)){install.packages('dplyr'); require(dplyr)}
if(!require(tidyr)){install.packages('tidyr'); require(tidyr)}
if(!require(prospectr)){install.packages('prospectr'); require(prospectr)}
if(!require(pls)){install.packages('pls'); require(pls)}
if(!require(resemble)){install.packages('resemble'); require(resemble)}
if(!require(clhs)){install.packages('clhs'); require(clhs)}
if(!require(ggplot2)){install.packages('ggplot2'); require(ggplot2)}

# Setting the folder path ----
dir.p = 'D:/project_lowcost_nir/database' # homeoffice/

# Importing data from ASD spectrometer ----
asd.osc = read.delim(paste0(dir.p, "/raw_data/oscar/asd/osc_asd.txt"), header = T, sep = ";")
head(asd.osc[, 1:10], 5)

# Reorganise the dataframe
asd.osc.df = as.data.frame(t(asd.osc))
colnames(asd.osc.df) = as.integer(asd.osc.df[1, ])
asd.osc.df = asd.osc.df[-1,]
asd.osc.df$id = rownames(asd.osc.df)

asd.osc.df = asd.osc.df %>%
  dplyr::select(c('id'), everything())

asd.osc.df$id = stringr::str_extract(asd.osc.df$id, "(?<=_)(\\d+)(?=\\_)")

asd.osc.df$id = paste("osc", asd.osc.df$id, sep = "_")

asd.osc.df[1:5, c(1:5, ncol(asd.osc.df))]

asd.osc.vnir = aggregate(.~id, data = asd.osc.df, FUN = mean) # averaging the 3 repetitions

asd.osc.vnir[1:5, c(1:5, ncol(asd.osc.vnir))]

write.csv(asd.osc.vnir, paste0(dir.p, '/tables/asd_vnir_oscar_100pts.csv'))

# Importing data from ISC NIRScan spectrometer ----
isc.osc = read.csv(paste0(dir.p, "/raw_data/oscar/isc/combinedSpec.csv"), header = T, sep = ",")[-1]
head(isc.osc[, 1:10], 5)

# Reorganise the dataframe
isc.osc.df = as.data.frame(t(isc.osc))
colnames(isc.osc.df) = as.integer(isc.osc.df[1, ])
isc.osc.df = isc.osc.df[-1,]

isc.osc.df$id = rownames(isc.osc.df)

isc.osc.df = isc.osc.df %>%
  dplyr::select(c('id'), everything())

isc.osc.df[1:5, c(1:5, ncol(isc.osc.df))]

isc.osc.nir = isc.osc.df

rownames(isc.osc.nir) = NULL

isc.osc.nir[1:5, c(1:5, ncol(isc.osc.nir))]

write.csv(isc.osc.nir, paste0(dir.p, '/tables/isc_nir_oscar_100pts.csv'))

# Plotting both spectra ----
## ASD
asd.osc.vnir[1:5, c(1:5, ncol(asd.osc.vnir))]

jpeg(paste0(dir.p, '/results/figs/original_asd_osc_vnir.jpeg'), width = 10, height = 10, units = 'in', res = 300)
matplot(as.numeric(colnames(asd.osc.vnir)[2:ncol(asd.osc.vnir)]), 
        t(asd.osc.vnir[,2:ncol(asd.osc.vnir)]),
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
mtext('Original ASD Spectra')
grid(lwd=1, nx=50, ny=50)
dev.off()

isc.osc.nir[1:5, c(1:5, ncol(isc.osc.nir))]

jpeg(paste0(dir.p, '/results/figs/original_isc_osc_nir.jpeg'), width = 10, height = 10, units = 'in', res = 300)
matplot(as.numeric(colnames(isc.osc.nir)[2:ncol(isc.osc.nir)]), 
        t(isc.osc.nir[,2:ncol(isc.osc.nir)]),
        xlim = c(350, 2500),
        xaxt = 'n',
        frame.plot = FALSE,
        type = 'l',
        lty = 1,
        col = rgb(red = 0.0, green = 0.4, blue = 0.6, alpha = 0.4),
        xlab =  'nm', #expression(paste('Wavelength(cm'^'-1',')'))
        ylab = 'Reflectance',
        las=1, font=2,
        add = T)
axis(1, at = seq(350, 2500, by = 50), labels = T, las = 2)
mtext('Original NIRScan Spectra')
grid(lwd=1, nx=50, ny=50)
dev.off()

## All together
jpeg(paste0(dir.p, '/results/figs/original_asd_vnir_isc_nir_osc.jpeg'), width = 10, height = 10, units = 'in', res = 300)
matplot(as.numeric(colnames(asd.osc.vnir)[2:ncol(asd.osc.vnir)]), 
        t(asd.osc.vnir[1:3,2:ncol(asd.osc.vnir)]),
        #xlim = c(400, 2450),
        ylim = c(0,1),
        xaxt = 'n',
        frame.plot = FALSE,
        type = 'l',
        lty = 1,
        col = c("red", "blue", "grey"), #rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3),
        xlab =  'Wavelength (nm)', #expression(paste('Wavenumber(cm'^'-1',')'))
        ylab = 'Reflectance Factor (x100)',
        las=1, font=2)
axis(1, at = seq(400, 2450, by = 50), labels = T, las = 2)

matplot(as.numeric(colnames(isc.osc.nir)[2:ncol(isc.osc.nir)]), 
        t(isc.osc.nir[1:3,2:ncol(isc.osc.nir)]),
        #xlim = c(400, 2450),
        #ylim = c(0,0.6),
        xaxt = 'n',
        frame.plot = FALSE,
        type = 'l',
        lty = 1,
        lwd = 2,
        col = c("red", "blue", "grey"), #rgb(red = 0.0, green = 0.4, blue = 0.6, alpha = 0.4),
        las=1, font=2,
        add = T)
legend("topright", legend = c("Sample 1", "Sample 2", "Sample 3"), col = c("red", "blue", "grey"), lty = 1, bty = "n")
dev.off()

# Pre-processing the spectra ----
## Splice Correction applicable for ASD only
asd.osc.vnir[1:5, c(1:5, ncol(asd.osc.vnir))]

db.spl.asd.osc = asd.osc.vnir
colnames(db.spl.asd.osc)[c(1:4, ncol(db.spl.asd.osc))]
rownames(db.spl.asd.osc) = db.spl.asd.osc$id

head(db.spl.asd.osc[, c(1:4, ncol(db.spl.asd.osc))])

class(db.spl.asd.osc)
colnames(db.spl.asd.osc)[1:10]
wav.spl.asd = as.numeric(colnames(db.spl.asd.osc[, 2:ncol(db.spl.asd.osc)]))
db.spl.asd.osc.m = as.matrix(db.spl.asd.osc[, 2:ncol(db.spl.asd.osc)]) # assembling matrix
class(db.spl.asd.osc.m) # checking class

db.spl.asd.osc.c = prospectr::spliceCorrection(db.spl.asd.osc.m, wav.spl.asd, 
                                               splice = c(1000, 1800), interpol.bands = 10)
class(db.spl.asd.osc.c)

db.spl.asd.osc.d = as.data.frame(db.spl.asd.osc.c) # converting back to DF
db.spl.asd.osc.d[1:20, 1:5]
db.spl.asd.osc.d$id = rownames(db.spl.asd.osc.d)
db.spl.asd.osc.d = db.spl.asd.osc.d %>%
  dplyr::select(c('id'), everything())

colnames(db.spl.asd.osc.d)[52:53]

db.spl.asd.osc.d[1:5, c(1:3, 52:53, ncol(db.spl.asd.osc.d))]

jpeg(paste0(dir.p, '/results/figs/splice_asd_osc_vnir.jpeg'), width = 10, height = 10, units = 'in', res = 300)
matplot(as.numeric(colnames(db.spl.asd.osc.d)[2:ncol(db.spl.asd.osc.d)]), 
        t(db.spl.asd.osc.d[,2:ncol(db.spl.asd.osc.d)]),
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
mtext('Splice Correction ASD Spectra')
grid(lwd=1, nx=50, ny=50)
dev.off()

write.csv(db.spl.asd.osc.d, paste0(dir.p, '/tables/asd_vnir_splicecorrection_oscar_100pts.csv'))

## Savitzky-Golay
### ASD
db.spl.asd.osc.d[1:5, c(1:5, ncol(db.spl.asd.osc.d))]
colnames(db.spl.asd.osc.d[, 1:20])

db.asd.osc.sg = prospectr::savitzkyGolay(db.spl.asd.osc.d[, 2:ncol(db.spl.asd.osc.d)],
                                         m = 0, # differentiation order (derivative)
                                         p = 2, # polynomial order
                                         w = 11 # window size
                                         ) 
class(db.asd.osc.sg)
db.asd.osc.sg = as.data.frame(db.asd.osc.sg)
db.asd.osc.sg[1:5, c(1:5, ncol(db.asd.osc.sg))]
db.asd.osc.sg$id = rownames(db.asd.osc.sg)
head(db.asd.osc.sg[, c(1:5, ncol(db.asd.osc.sg))])
#colnames(db.train.sg) = gsub("X", "", names(db.train.sg))

db.asd.osc.sg = db.asd.osc.sg %>%
  select(c("id"), everything())

head(db.asd.osc.sg[, c(1:6, ncol(db.asd.osc.sg))])

jpeg(paste0(dir.p, '/results/figs/sg_2ord_11pts_0der_asd_osc_vnir.jpeg'), width = 10, height = 10, units = 'in', res = 300)
matplot(as.numeric(colnames(db.asd.osc.sg)[2:ncol(db.asd.osc.sg)]), 
        t(db.asd.osc.sg[,2:ncol(db.asd.osc.sg)]),
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
mtext('ASD Spectra: Savitzky-Golay 2nd order polynomial and window size of 11 smoothing points')
grid(lwd=1, nx=50, ny=50)
dev.off()

write.csv(db.asd.osc.sg, paste0(dir.p, '/tables/asd_vnir_sg_2p_w11_0d_oscar_100pts.csv'))

### NIRScan
isc.osc.nir[1:5, c(1:5, ncol(isc.osc.nir))]
colnames(isc.osc.nir[, 1:20])

rownames(isc.osc.nir) = isc.osc.nir$id

db.isc.osc.sg = prospectr::savitzkyGolay(isc.osc.nir[, 2:ncol(isc.osc.nir)],
                                         m = 0, # differentiation order (derivative)
                                         p = 2, # polynomial order
                                         w = 11 # window size
) 
class(db.isc.osc.sg)
db.isc.osc.sg = as.data.frame(db.isc.osc.sg)
db.isc.osc.sg[1:5, c(1:5, ncol(db.isc.osc.sg))]
db.isc.osc.sg$id = rownames(db.isc.osc.sg)
head(db.isc.osc.sg[, c(1:5, ncol(db.isc.osc.sg))])
#colnames(db.train.sg) = gsub("X", "", names(db.train.sg))

db.isc.osc.sg = db.isc.osc.sg %>%
  select(c("id"), everything())

head(db.isc.osc.sg[, c(1:6, ncol(db.isc.osc.sg))])

jpeg(paste0(dir.p, '/results/figs/sg_2ord_11pts_0der_isc_osc_vnir.jpeg'), width = 10, height = 10, units = 'in', res = 300)
matplot(as.numeric(colnames(db.isc.osc.sg)[2:ncol(db.isc.osc.sg)]), 
        t(db.isc.osc.sg[,2:ncol(db.isc.osc.sg)]),
        xlim = c(350, 2500),
        #ylim = c(0.0, 0.7),
        xaxt = 'n',
        frame.plot = FALSE,
        type = 'l',
        lty = 1,
        col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5),
        xlab =  'nm', #expression(paste('Wavelength(cm'^'-1',')'))
        ylab = 'Reflectance',
        las=1, font=2)
axis(1, at = seq(350, 2500, by = 50), labels = T, las = 2)
mtext('NIRScan Spectra: Savitzky-Golay 2nd order polynomial and window size of 11 smoothing points')
grid(lwd=1, nx=50, ny=50)
dev.off()

write.csv(db.isc.osc.sg, paste0(dir.p, '/tables/isc_vnir_sg_2p_w11_0d_oscar_100pts.csv'))

## All in one
jpeg(paste0(dir.p, '/results/figs/sg_asd_vnir_isc_nir_osc.jpeg'), width = 10, height = 10, units = 'in', res = 300)
matplot(as.numeric(colnames(db.asd.osc.sg)[2:ncol(db.asd.osc.sg)]), 
        t(db.asd.osc.sg[1:3,2:ncol(db.asd.osc.sg)]),
        #xlim = c(400, 2450),
        ylim = c(0,1),
        xaxt = 'n',
        frame.plot = FALSE,
        type = 'l',
        lty = 1,
        col = c("red", "blue", "grey"), #rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3), 
        xlab =  'Wavelength (nm)', #expression(paste('Wavenumber(cm'^'-1',')'))
        ylab = 'Reflectance Factor (x100)',
        las=1, font=2)
axis(1, at = seq(350, 2500, by = 50), labels = T, las = 2)

matplot(as.numeric(colnames(db.isc.osc.sg)[2:ncol(db.isc.osc.sg)]), 
        t(db.isc.osc.sg[1:3,2:ncol(db.isc.osc.sg)]),
        #xlim = c(400, 2450),
        #ylim = c(0,0.6),
        xaxt = 'n',
        frame.plot = FALSE,
        type = 'l',
        lty = 1,
        lwd = 2,
        col = c("red", "blue", "grey"), #rgb(red = 0.0, green = 0.4, blue = 0.6, alpha = 0.4),
        las=1, font=2,
        add = T)
legend("topright", legend = c("Sample 1", "Sample 2", "Sample 3"), col = c("red", "blue", "grey"), lty = 1, bty = "n")
dev.off()

# Selecting the same spectral range (950 - 1650 nm) ----
## ASD
db.asd.osc.sg[1:3, c(596:597, 1296:1297)]
db.asd.osc.ok = db.asd.osc.sg[, 597:1297]
db.asd.osc.ok[1:3, c(1:4, ncol(db.asd.osc.ok))]

## NIRScan
db.isc.osc.sg[1:3, c(9:10, 209:211)]
db.isc.osc.ok = db.isc.osc.sg[, 9:209]
db.isc.osc.ok[1:3, c(1:4, ncol(db.isc.osc.ok))]


# Resampling them to the same window size of 5 nm ----
## ASD from 1 nm
wav.asd = as.numeric(colnames(db.asd.osc.ok)) # retrieving wavelengths

db.asd.osc.r = as.matrix(db.asd.osc.ok)

#rownames(db.asd.osc.r) = db.asd.osc.ok$id

db.asd.osc.r = prospectr::resample(db.asd.osc.r, 
                                   wav.asd, 
                                   seq(950, 1650, 5), 
                                   interpol = "spline") # window size of 5 nm
class(db.asd.osc.r)
db.asd.osc.r = as.data.frame(db.asd.osc.r)

db.asd.osc.r$id = rownames(db.asd.osc.r)
db.asd.osc.r = db.asd.osc.r %>%
  select(c("id"), everything())

db.asd.osc.r[1:5, c(1:5, ncol(db.asd.osc.r))]

## NIRScan from 3 nm
wav.isc = as.numeric(colnames(db.isc.osc.ok)) # retrieving wavelengths

db.isc.osc.r = as.matrix(db.isc.osc.ok)

#rownames(db.isc.osc.r) = db.isc.osc.ok$id

db.isc.osc.r = prospectr::resample(db.isc.osc.r, 
                                   wav.isc, 
                                   seq(950, 1650, 5), 
                                   interpol = "spline") # window size of 5 nm
class(db.isc.osc.r)
db.isc.osc.r = as.data.frame(db.isc.osc.r)

db.isc.osc.r$id = rownames(db.isc.osc.r)
db.isc.osc.r = db.isc.osc.r %>%
  select(c("id"), everything())

db.isc.osc.r[1:5, c(1:5, ncol(db.isc.osc.r))]

## Plot all together
jpeg(paste0(dir.p, '/results/figs/resampled_asd_isc_nir_osc.jpeg'), width = 10, height = 10, units = 'in', res = 300)
matplot(as.numeric(colnames(db.asd.osc.r)[2:ncol(db.asd.osc.r)]), 
        t(db.asd.osc.r[,2:ncol(db.asd.osc.r)]),
        #xlim = c(400, 2450),
        ylim = c(0,1),
        xaxt = 'n',
        frame.plot = FALSE,
        type = 'l',
        lty = 1,
        col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3), #c("red", "blue", "grey"),
        xlab =  'Wavelength (nm)', #expression(paste('Wavenumber(cm'^'-1',')'))
        ylab = 'Reflectance Factor (x100)',
        las=1, font=2)
axis(1, at = seq(950, 1650, by = 10), labels = T, las = 2)

matplot(as.numeric(colnames(db.isc.osc.r)[2:ncol(db.isc.osc.r)]), 
        t(db.isc.osc.r[,2:ncol(db.isc.osc.r)]),
        #xlim = c(400, 2450),
        #ylim = c(0,0.6),
        xaxt = 'n',
        frame.plot = FALSE,
        type = 'l',
        lty = 1,
        lwd = 2,
        col = rgb(red = 0.6, green = 0.4, blue = 0.2, alpha = 0.4), #c("red", "blue", "grey"), 
        las=1, font=2,
        add = F)
legend("topright", legend = c("ASD", "NIRScan"), col = c(rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3),
                                                         rgb(red = 0.6, green = 0.4, blue = 0.2, alpha = 0.4)), 
       lty = 1, bty = "n")
dev.off()

# Splitting the dataset using cLHS ----
db.asd.osc.r
colnames(db.asd.osc.r)[1:20]
db.asd.osc.r[1:5, 1:5]

#### Performing PC analysis before splitting the data into cal and val using ASD as standard
maxexplvar = 0.999 # indicating the max amount of cumulative variance explained
asd.osc.pc.spec = resemble::pc_projection(Xr = as.matrix(db.asd.osc.r[, -1]),
                                                 pc_selection = list("cumvar", maxexplvar),
                                                 center = T,
                                                 scale = F)


### Conditioned Latin Hypercube Sampling (cLHS)
nrow(db.asd.osc.r) # 100 obs.
class(asd.osc.pc.spec$scores)

# Soil attributes with 262 obs.
asd.osc.clhs = clhs::clhs(x = as.data.frame(asd.osc.pc.spec$scores),
                          size = round(0.70*nrow(db.asd.osc.r)), 
                          iter = 1000, simple = FALSE)

jpeg(paste0(dir.p, '/results/figs/asd_osc_clhs.jpeg'),  
     width = 10, height = 10, units = 'in' , res = 300)

## plot the 2 PCs
# Create a data frame for the first two PCs and calibration samples
pc_scores = as.data.frame(asd.osc.pc.spec$scores[, 1:2])
calibration_samples = as.data.frame(asd.osc.pc.spec$scores[, 1:2][asd.osc.clhs$index_samples, ])

# Add explained variances and cumulative explained variances to the plot
explained_var = as.data.frame(asd.osc.pc.spec$variance$x_var[2, ]) #explained_var
cumulative_explained_var = as.data.frame(asd.osc.pc.spec$variance$x_var[3, ]) #cumulative_explained_var

# Plot the first two PCs
ggplot(pc_scores, aes(pc_1, pc_2)) +
  geom_point(color = "gray", alpha = 0.5, size = 3) +
  geom_point(data = calibration_samples, color = "blue", size = 3) +
  theme_minimal() +
  labs(x = "PC 1 (99.5%)", y = "PC 2 (0.40%)") +
  ggtitle("Scatter Plot of PC 1 and PC 2") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") #+
  #geom_text(aes(label = paste0("Explained variance: ", round(explained_var, 3))), 
  #          x = max(pc_scores$pc_1), y = max(pc_scores$pc_2), hjust = 1, vjust = 1) +
  #geom_text(aes(label = paste0("Cumulative explained variance: ", round(cumulative_explained_var, 3))), 
  #          x = max(pc_scores$pc_1), y = max(pc_scores$pc_2) - 0.01, hjust = 1, vjust = 1)

dev.off()

# ASD calibration
db.asd.osc.r.cal = db.asd.osc.r[asd.osc.clhs$index_samples, ] 
nrow(db.asd.osc.r.cal) # 70% / 70
write.csv(db.asd.osc.r.cal, paste0(dir.p, '/tables/asd_nir_sg_cal_oscar_70pts.csv'))

db.asd.osc.sg.cal = db.asd.osc.sg[db.asd.osc.sg$id %in% db.asd.osc.r.cal$id,]
write.csv(db.asd.osc.sg.cal, paste0(dir.p, '/tables/db_model/asd_visnir_sg_cal_oscar_70pts.csv'))

db.asd.osc.r.val = db.asd.osc.r[-asd.osc.clhs$index_samples, ] 
nrow(db.asd.osc.r.val) # 30% / 30
write.csv(db.asd.osc.r.val, paste0(dir.p, '/tables/asd_nir_sg_cal_oscar_30pts.csv'))

db.asd.osc.sg.val = db.asd.osc.sg[db.asd.osc.sg$id %in% db.asd.osc.r.val$id,]
write.csv(db.asd.osc.sg.val, paste0(dir.p, '/tables/db_model/asd_visnir_sg_val_oscar_70pts.csv'))


# Selecting the same samples from ISC NIRScan
db.isc.osc.r
colnames(db.isc.osc.r)
colnames(db.asd.osc.r.cal)
db.isc.osc.r.cal = db.isc.osc.r[db.isc.osc.r$id %in% db.asd.osc.r.cal$id,]
write.csv(db.isc.osc.r.cal, paste0(dir.p, '/tables/isc_nir_sg_cal_oscar_70pts.csv'))

db.isc.osc.r.val = db.isc.osc.r[db.isc.osc.r$id %in% db.asd.osc.r.val$id,]
write.csv(db.isc.osc.r.val, paste0(dir.p, '/tables/isc_nir_sg_val_oscar_30pts.csv'))

# Organising the data into source and target calibration and validation and performing SNV ----
## CALIBRATION
### SOURCE
source.spectra = data.frame(db.isc.osc.r.cal[, 1], NIR = I(db.isc.osc.r.cal[,-1]))
colnames(source.spectra)[1] = "id"
write.csv(source.spectra, paste0(dir.p, '/tables/db_model/db_source_isc_spec_cal.csv'))
#### SNV
source.snv = source.spectra
source.snv$NIR = prospectr::standardNormalVariate(source.snv$NIR)
colnames(source.snv)
write.csv(source.snv, paste0(dir.p, '/tables/db_model/db_source_isc_spec_snv_cal.csv'))

### TARGET
target.spectra = data.frame(db.asd.osc.r.cal[, 1], NIR = I(db.asd.osc.r.cal[,-1]))
colnames(target.spectra)[1] = "id"
write.csv(target.spectra, paste0(dir.p, '/tables/db_model/db_target_asd_spec_cal.csv'))
#### SNV
target.snv = target.spectra
target.snv$NIR = prospectr::standardNormalVariate(target.snv$NIR)
colnames(target.snv)
write.csv(target.snv, paste0(dir.p, '/tables/db_model/db_target_asd_spec_snv_cal.csv'))

wavelength = as.numeric(substring(colnames(source.snv$NIR), 1, 7))
str(source.snv)
class(source.snv)

## VALIDATION
### SOURCE
source.spectra.val = data.frame(db.isc.osc.r.val[, 1], NIR = I(db.isc.osc.r.val[,-1]))
colnames(source.spectra.val)[1] = "id"
write.csv(source.spectra.val, paste0(dir.p, '/tables/db_model/db_source_isc_spec_val.csv'))
### SNV
source.snv.val = source.spectra.val
source.snv.val$NIR = prospectr::standardNormalVariate(source.snv.val$NIR)
colnames(source.snv.val)
write.csv(source.snv.val, paste0(dir.p, '/tables/db_model/db_source_isc_spec_snv_val.csv'))

### TARGET
target.spectra.val = data.frame(db.asd.osc.r.val[, 1], NIR = I(db.asd.osc.r.val[,-1]))
colnames(target.spectra.val)[1] = "id"
write.csv(target.spectra.val, paste0(dir.p, '/tables/db_model/db_target_asd_spec_val.csv'))

#### SNV
target.snv.val = target.spectra.val
target.snv.val$NIR = prospectr::standardNormalVariate(target.snv.val$NIR)
colnames(target.snv.val)
write.csv(target.snv.val, paste0(dir.p, '/tables/db_model/db_target_asd_spec_snv_val.csv'))

# Substraction transfer, spectral outliers removed by IQR -----
common.ids = intersect(db.asd.osc.r.cal$id, db.isc.osc.r.cal$id)

## Through substraction
subtracted.spectra.asd.isc = data.frame(id = common.ids)

for (col in names(target.spectra$NIR)) {  # Exclude the ID column
  subtracted.spectra.asd.isc[[col]] = round(target.spectra$NIR[target.spectra$id %in% common.ids, col] - 
                                      source.spectra$NIR[source.spectra$id %in% common.ids, col],3)
}

# Print the subtracted dataframe
print(subtracted.spectra.asd.isc)

summary(subtracted.spectra.asd.isc[, c("1400", "1415")])
boxplot(subtracted.spectra.asd.isc[, "1415"])
# Removing outliers using as reference water band 1415 nm
quartiles = quantile(subtracted.spectra.asd.isc[, "1415"], probs=c(.25, .75), rm.na=T)
iqr = IQR(subtracted.spectra.asd.isc[, "1415"])
lv = quartiles[1]-1.5*iqr
uv = quartiles[2]+1.5*iqr

subtracted.spectra.asd.isc.noout = subset(subtracted.spectra.asd.isc, subtracted.spectra.asd.isc[, "1415"] > lv & subtracted.spectra.asd.isc[, "1415"] < uv)
boxplot(subtracted.spectra.asd.isc.noout[, "1415"])
hist(subtracted.spectra.asd.isc.noout[, "1415"])
summary(subtracted.spectra.asd.isc.noout[, c("1400", "1415")])
isc.subtracted.median = apply(subtracted.spectra.asd.isc.noout[-1], 2, median)
isc.factor = isc.subtracted.median['1415'] # Water band 1400 nm or *1415nm

source.spectra.factor.cor = source.spectra

source.spectra.factor.cor$NIR = source.spectra.factor.cor$NIR + isc.factor

## All in one Calibration
jpeg(paste0(dir.p, '/results/figs/MEAN_SPEC_SIGN_factor_correction_nir_asd_isc_osc_calibration.jpeg'), width = 10, height = 10, units = 'in', res = 300)
matplot(wavelength, 
        colMeans(target.spectra$NIR),
        #xlim = c(400, 2450),
        ylim = c(0,1),
        xaxt = 'n',
        frame.plot = FALSE,
        type = 'l',
        lty = 3,
        lwd = 2,
        col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.8),
        xlab =  'Wavelength (nm)', #expression(paste('Wavenumber(cm'^'-1',')'))
        ylab = 'Reflectance Factor (x100)',
        las=1, font=2)
axis(1, at = seq(950, 1650, by = 10), labels = T, las = 2)

matplot(wavelength, 
        colMeans(source.spectra$NIR),
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

matplot(wavelength, 
        colMeans(source.spectra.factor.cor$NIR),
        #xlim = c(400, 2450),
        #ylim = c(0,0.6),
        xaxt = 'n',
        frame.plot = FALSE,
        type = 'l',
        lty = 1,
        lwd = 1,
        col = rgb(red = 1.0, green = 0.0, blue = 0.0, alpha = 0.5),
        las=1, font=2,
        add = T)

legend("topleft",bty="n", legend=c("Target: ASD data","Source: ISC NIRScan data", "Index: Fitting source in target data"),
       col=c(rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3),
             "forestgreen",
             "red"),
       lty=c(1),cex=1,lwd=3)

dev.off()

## VALIDATION
source.spectra.factor.cor.val = source.spectra.val

source.spectra.factor.cor.val$NIR = source.spectra.factor.cor.val$NIR + isc.factor

## All in one Validation
jpeg(paste0(dir.p, '/results/figs/MEAN_SPEC_SIGN_factor_correction_nir_asd_isc_osc_validation.jpeg'), width = 10, height = 10, units = 'in', res = 300)
matplot(wavelength, 
        colMeans(target.spectra.val$NIR),
        #xlim = c(400, 2450),
        ylim = c(0,1),
        xaxt = 'n',
        frame.plot = FALSE,
        type = 'l',
        lty = 3,
        lwd = 2,
        col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.8),
        xlab =  'Wavelength (nm)', #expression(paste('Wavenumber(cm'^'-1',')'))
        ylab = 'Reflectance Factor (x100)',
        las=1, font=2)
axis(1, at = seq(950, 1650, by = 10), labels = T, las = 2)

matplot(wavelength, 
        colMeans(source.spectra.val$NIR),
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

matplot(wavelength, 
        colMeans(source.spectra.factor.cor.val$NIR),
        #xlim = c(400, 2450),
        #ylim = c(0,0.6),
        xaxt = 'n',
        frame.plot = FALSE,
        type = 'l',
        lty = 1,
        lwd = 1,
        col = rgb(red = 1.0, green = 0.0, blue = 0.0, alpha = 0.5),
        las=1, font=2,
        add = T)

legend("topleft",bty="n", legend=c("Target: ASD data","Source: ISC NIRScan data", "Index: Fitting source in target data"),
       col=c(rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3),
             "forestgreen",
             "red"),
       lty=c(1),cex=1,lwd=3)

dev.off()

### Exporting data.frame
source.spectra.factor.cor
write.csv(source.spectra.factor.cor, paste0(dir.p, '/tables/db_model/db_source_spec_subcorr_cal.csv'))
source.spectra.factor.cor.val
write.csv(source.spectra.factor.cor.val, paste0(dir.p, '/tables/db_model/db_source_spec_subcorr_val.csv'))

## SG-SNV CALIBRATION
subtracted.snv.asd.isc = data.frame(id = common.ids)
#colnames(target.snv$NIR)
for (col in colnames(target.snv$NIR)) {  # Exclude the ID column
  subtracted.snv.asd.isc[[col]] = round(target.snv$NIR[target.snv$id %in% common.ids, col] - 
                                              source.snv$NIR[source.snv$id %in% common.ids, col],3)
}

# Print the subtracted dataframe
print(subtracted.snv.asd.isc)

summary(subtracted.snv.asd.isc[, c("1400", "1415")])
boxplot(subtracted.snv.asd.isc[, "1415"])
# Removing outliers using as reference water band 1415 nm
quartiles.snv = quantile(subtracted.snv.asd.isc[, "1415"], probs=c(.25, .75), rm.na=T)
iqr.snv = IQR(subtracted.snv.asd.isc[, "1415"])
lv.snv = quartiles[1]-1.5*iqr
uv.snv = quartiles[2]+1.5*iqr

subtracted.snv.asd.isc.noout = subset(subtracted.snv.asd.isc, subtracted.snv.asd.isc[, "1415"] > lv.snv & subtracted.snv.asd.isc[, "1415"] < uv.snv)
boxplot(subtracted.snv.asd.isc.noout[, "1415"])
hist(subtracted.snv.asd.isc.noout[, "1415"])
summary(subtracted.snv.asd.isc.noout[, c("1400", "1415")])
isc.subtracted.snv.median = apply(subtracted.snv.asd.isc.noout[-1], 2, median)
isc.factor.snv = isc.subtracted.snv.median['1415'] # Water band 1400 nm or *1415nm

source.snv.factor.cor = source.snv

source.snv.factor.cor$NIR = source.snv.factor.cor$NIR + isc.factor.snv

## All in one Calibration
jpeg(paste0(dir.p, '/results/figs/factor_correction_nir_asd_isc_SNV_osc_calibration.jpeg'), width = 10, height = 10, units = 'in', res = 300)
matplot(wavelength, 
        colMeans(target.snv$NIR),
        #xlim = c(400, 2450),
        ylim = c(-3,3),
        xaxt = 'n',
        frame.plot = FALSE,
        type = 'l',
        lty = 3,
        lwd = 2,
        col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.8),
        xlab =  'Wavelength (nm)', #expression(paste('Wavenumber(cm'^'-1',')'))
        ylab = 'SG-SNV',
        las=1, font=2)
axis(1, at = seq(950, 1650, by = 10), labels = T, las = 2)

matplot(wavelength, 
        colMeans(source.snv$NIR),
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

matplot(wavelength, 
        colMeans(source.snv.factor.cor$NIR),
        #xlim = c(400, 2450),
        #ylim = c(0,0.6),
        xaxt = 'n',
        frame.plot = FALSE,
        type = 'l',
        lty = 1,
        lwd = 1,
        col = rgb(red = 1.0, green = 0.0, blue = 0.0, alpha = 0.5),
        las=1, font=2,
        add = T)

legend("topleft",bty="n", legend=c("Target: ASD data","Source: ISC NIRScan data", "Index: Fitting source in target data"),
       col=c(rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3),
             "forestgreen",
             "red"),
       lty=c(1),cex=1,lwd=3)

dev.off()

## VALIDATION
source.snv.factor.cor.val = source.snv.val

source.snv.factor.cor.val$NIR = source.snv.factor.cor.val$NIR + isc.factor.snv

## All in one Validation
jpeg(paste0(dir.p, '/results/figs/factor_correction_nir_asd_isc_SNV_osc_validation.jpeg'), width = 10, height = 10, units = 'in', res = 300)
matplot(wavelength, 
        colMeans(target.snv.val$NIR),
        #xlim = c(400, 2450),
        ylim = c(-3,3),
        xaxt = 'n',
        frame.plot = FALSE,
        type = 'l',
        lty = 3,
        lwd = 2,
        col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.8),
        xlab =  'Wavelength (nm)', #expression(paste('Wavenumber(cm'^'-1',')'))
        ylab = 'Reflectance Factor (x100)',
        las=1, font=2)
axis(1, at = seq(950, 1650, by = 10), labels = T, las = 2)

matplot(wavelength, 
        colMeans(source.snv.val$NIR),
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

matplot(wavelength, 
        colMeans(source.snv.factor.cor.val$NIR),
        #xlim = c(400, 2450),
        #ylim = c(0,0.6),
        xaxt = 'n',
        frame.plot = FALSE,
        type = 'l',
        lty = 1,
        lwd = 1,
        col = rgb(red = 1.0, green = 0.0, blue = 0.0, alpha = 0.5),
        las=1, font=2,
        add = T)

legend("topleft",bty="n", legend=c("Target: ASD data","Source: ISC NIRScan data", "Index: Fitting source in target data"),
       col=c(rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3),
             "forestgreen",
             "red"),
       lty=c(1),cex=1,lwd=3)
 
dev.off()

### Exporting data.frame
source.snv.factor.cor
write.csv(source.snv.factor.cor, paste0(dir.p, '/tables/db_model/db_source_spec_snv_subcorr_cal.csv'))
source.snv.factor.cor.val
write.csv(source.snv.factor.cor.val, paste0(dir.p, '/tables/db_model/db_source_spec_snv_subcorr_val.csv'))

# Piecewise Direct Standardisation (PDS) ----
# Define the Piecewise Direct Standardization (PDS) function with PLS regression
perform_pds <- function(target, source, MWsize, Ncomp, wavelength){
  
  require(pls)
  
  #Loop Initialization:
  i<-MWsize
  k<-i-1
  #Creation of an empty P matrix:
  P<-matrix(0,nrow=ncol(target),ncol=ncol(target)-(2*i)+2)
  InterceptReg<-c()
  
  while(i<=(ncol(target)-k)){
    
    #PLS regression:
    fit<- plsr(target[,i] ~ as.matrix(source[,(i-k):(i+k)]),
               ncomp=Ncomp, scale=F, method="oscorespls")
    
    #Extraction of the regression coefficients:
    coefReg<-as.numeric(coef(fit, ncomp=Ncomp, intercept = TRUE))
    InterceptReg<-c(InterceptReg,coefReg[1])
    coefReg<-coefReg[2:length(coefReg)]
    
    #Add coefficients to the transfer matrix:
    P[(i-k):(i+k),i-k]<-t(coefReg)
    
    rm(coefReg,fit)
    i<-i+1
    
    #Diplay progression:
    cat("\r",paste(round(i/ncol(target)*100)," %",sep=""))}
  
  P<-data.frame(matrix(0,nrow=ncol(target),ncol=k), P,
                matrix(0,nrow=ncol(target),ncol=k))
  InterceptReg<-c(rep(0,k),InterceptReg,rep(0,k)) 
  
  Output<-list(P = P , Intercept = InterceptReg)
  
  return(Output)}

# Apply PDS to transfer spectra from source to target instrument
Ncomp <- 2 # 3
MWsize <- 2 # 8
#wavelength = as.numeric(substring(colnames(target.snv$NIR), 1, 7))
transformed.pds.target <- perform_pds(target.snv$NIR, source.snv$NIR, MWsize, Ncomp, wavelength)

## Matrix that will help to transfer spec from source to target using SNV data
source.snv.pds.cal <- as.matrix(source.snv$NIR)%*%as.matrix(transformed.pds.target$P)
source.snv.pds.cal <- sweep(source.snv.pds.cal, 2, as.numeric(t(transformed.pds.target$Intercept)), "+")
source.snv.pds.cal[1:5, 1:141]
colnames(source.snv.pds.cal[, 1:141])
source.snv.pds.cal <- data.frame(source.snv$id, NIR = I(source.snv.pds.cal[, 2:140]))
colnames(source.snv.pds.cal)[1] = "id"


## All in one graph - calibration
jpeg(paste0(dir.p, '/results/figs/MEAN_SPEC_SIGN_pds_correction_nir_asd_isc_osc_cal.jpeg'), width = 10, height = 10, units = 'in', res = 300)
matplot(wavelength[2:140],
        colMeans(target.snv$NIR[,2:140]),
        #xlim = c(400, 2450),
        ylim = c(-3,3),
        xaxt = 'n',
        frame.plot = FALSE,
        type = 'l',
        lty = 3,
        lwd = 2,
        col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.8), 
        xlab =  'Wavelength (nm)', #expression(paste('Wavenumber(cm'^'-1',')'))
        ylab = 'SG-SNV',
        las=1, font=2)
axis(1, at = seq(950, 1650, by = 10), labels = T, las = 2)

matplot(wavelength[2:140],
        colMeans(source.snv$NIR[, 2:140]),
        #xlim = c(400, 2450),
        #ylim = c(0,0.6),
        xaxt = 'n',
        frame.plot = FALSE,
        type = 'l',
        lty = 1,
        lwd = 1,
        col = rgb(red = 0.13, green = 0.55, blue = 0.13, alpha = 0.8),
        las=1, font=2,
        add = T)

matplot(wavelength[2:140],
        colMeans(source.snv.pds.cal$NIR),
        #xlim = c(400, 2450),
        #ylim = c(0,0.6),
        xaxt = 'n',
        frame.plot = FALSE,
        type = 'l',
        lty = 1,
        lwd = 1,
        col = rgb(red = 1.0, green = 0.0, blue = 0.0, alpha = 0.8),
        las=1, font=2,
        add = T)
legend("topleft",bty="n", legend=c("Target: ASD data","Source: ISC NIRScan data", "PDS: Fitting source in target data"),
       col=c(rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3),
             "forestgreen",
             "red"),
       lty=c(1),cex=1,lwd=3)

dev.off()

# VALIDATION
## Matrix that will help to transfer spec from source to target using SNV data
source.snv.pds.val <- as.matrix(source.snv.val$NIR)%*%as.matrix(transformed.pds.target$P)
source.snv.pds.val <- sweep(source.snv.pds.val, 2, as.numeric(t(transformed.pds.target$Intercept)), "+")
source.snv.pds.val[1:5, 1:141]
colnames(source.snv.pds.val[, 1:141])
source.snv.pds.val <- data.frame(source.snv.val$id, NIR = I(source.snv.pds.val[, 2:140]))
colnames(source.snv.pds.val)[1] = "id"


## All in one graph - calibration
jpeg(paste0(dir.p, '/results/figs/MEAN_SPEC_SIGN_pds_correction_nir_asd_isc_osc_val.jpeg'), width = 10, height = 10, units = 'in', res = 300)
matplot(wavelength[2:140],
        colMeans(target.snv.val$NIR[,2:140]),
        #xlim = c(400, 2450),
        ylim = c(-3,3),
        xaxt = 'n',
        frame.plot = FALSE,
        type = 'l',
        lty = 3,
        lwd = 2,
        col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.8), 
        xlab =  'Wavelength (nm)', #expression(paste('Wavenumber(cm'^'-1',')'))
        ylab = 'SG-SNV',
        las=1, font=2)
axis(1, at = seq(950, 1650, by = 10), labels = T, las = 2)

matplot(wavelength[2:140],
        colMeans(source.snv.val$NIR[, 2:140]),
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

matplot(wavelength[2:140],
        colMeans(source.snv.pds.val$NIR),
        #xlim = c(400, 2450),
        #ylim = c(0,0.6),
        xaxt = 'n',
        frame.plot = FALSE,
        type = 'l',
        lty = 1,
        lwd = 1,
        col = rgb(red = 1.0, green = 0.0, blue = 0.0, alpha = 0.8),
        las=1, font=2,
        add = T)
legend("topleft",bty="n", legend=c("Target: ASD data","Source: ISC NIRScan data", "PDS: Fitting source in target data"),
       col=c(rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3),
             "forestgreen",
             "red"),
       lty=c(1),cex=1,lwd=3)

dev.off()

### Exporting data.frame
colnames(source.snv.pds.cal$NIR) = wavelength[2:140]
write.csv(source.snv.pds.cal, paste0(dir.p, '/tables/db_model/db_source_spec_snv_pdscorr_cal.csv'))
source.snv.pds.val
colnames(source.snv.pds.val$NIR) = wavelength[2:140]
write.csv(source.snv.pds.val, paste0(dir.p, '/tables/db_model/db_source_spec_snv_pdscorr_val.csv'))


# Other method ----

# Saving R workspace ----
save.image(file = paste0(dir.p, '/r_code/proj_low_cost_nir.RData'))
