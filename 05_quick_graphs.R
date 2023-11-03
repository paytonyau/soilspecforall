 
# Setting up local folder ----
dir.l = 'C:/BSSS_project_lowcost_nir/database/'

# R packages / Libraries ----
if(!require(dplyr)){install.packages('dplyr'); require(dplyr)}
if(!require(tidyr)){install.packages('tidyr'); require(tidyr)}
if(!require(prospectr)){install.packages('prospectr'); require(prospectr)}
if(!require(pls)){install.packages('pls'); require(pls)}
if(!require(resemble)){install.packages('resemble'); require(resemble)}
if(!require(clhs)){install.packages('clhs'); require(clhs)}
if(!require(ggplot2)){install.packages('ggplot2'); require(ggplot2)}
if(!require(ggridges)){install.packages('ggridges'); require(ggridges)}
if(!require(patchwork)){install.packages('patchwork'); require(patchwork)}


# soil data
db.soil = read.csv(paste0(dir.l, "/raw_data/db_soil_analysis_all.csv"))
head(db.soil, 5)
summary(db.soil)

# Create a blind-friendly color palette
blind_friendly_colors <- c("#56B4E9", "#E69F00", "#999999", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


carbon.plot <- ggplot(db.soil, aes(x = tc.perc, y = country, fill = landuse)) +
  geom_density_ridges(alpha = 0.7, color = "white") +  # Adjust alpha and set a white outline
  labs(x = "%", y = "", fill = "Land Use") +
  scale_fill_manual(values = blind_friendly_colors) +  # Use the blind-friendly color palette
  ggtitle("Total Carbon") +
  theme_ridges() +
  theme(
    #legend.position="top",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )

carbon.plot

nitrogen.plot <- ggplot(db.soil, aes(x = tn.perc, y = country, fill = landuse)) +
  geom_density_ridges(alpha = 0.7, color = "white") +  # Adjust alpha and set a white outline
  labs(x = "%", y = "", fill = "Land Use") +
  scale_fill_manual(values = blind_friendly_colors) +  # Use the blind-friendly color palette
  ggtitle("Total Nitrogen") +
  theme_ridges() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )

pH.plot <- ggplot(db.soil, aes(x = ph.water, y = country, fill = landuse)) +
  geom_density_ridges(alpha = 0.7, color = "white") +  # Adjust alpha and set a white outline
  labs(x = "adimensional", y = "", fill = "Land Use") +
  scale_fill_manual(values = blind_friendly_colors) +  # Use the blind-friendly color palette
  ggtitle("pH in water") +
  theme_ridges() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )

p.plot <- ggplot(db.soil, aes(x = p.dl.mg_100g_ok, y = country, fill = landuse)) +
  geom_density_ridges(alpha = 0.7, color = "white") +  # Adjust alpha and set a white outline
  labs(x = "mg/100g", y = "", fill = "Land Use") +
  scale_fill_manual(values = blind_friendly_colors) +  # Use the blind-friendly color palette
  ggtitle("Phosphorus") +
  theme_ridges() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )


combined_plots <- (carbon.plot / nitrogen.plot / pH.plot / p.plot) +
  plot_layout(ncol = 2, nrow = 2)

jpeg(paste0(dir.l, '/results/soil_data.jpeg'), width = 10, height = 10, units = 'in', res = 300)

combined_plots

dev.off()

# Plot the specs
db.isc.nir = read.csv(paste0(dir.l, "raw_data/spectra_usa_ge/isc/combinedSpec.csv"))[-1]
head(db.isc.nir, 5)
colnames(db.isc.nir) = gsub("X", "", names(db.isc.nir))

db.asd.vnir = read.csv(paste0(dir.l, "raw_data/spectra_usa_ge/asd/asd_all_usa_ge.csv"))
head(db.asd.vnir[, c(1:5, ncol(db.asd.vnir))], 5)
colnames(db.asd.vnir) = gsub("X", "", names(db.asd.vnir))


#############################
jpeg(paste0(dir.l, '/results/original_isc_nir.jpeg'), width = 10, height = 10, units = 'in', res = 300)
matplot(as.numeric(colnames(db.isc.nir)[2:ncol(db.isc.nir)]), 
        t(db.isc.nir[,2:ncol(db.isc.nir)]),
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
mtext('ISC NIR Spectra')
grid(lwd=1, nx=50, ny=50)
dev.off()

db.isc.nir[1:5, c(1:5, ncol(db.isc.nir))]

jpeg(paste0(dir.l, '/results/original_asd_vnir.jpeg'), width = 10, height = 10, units = 'in', res = 300)
matplot(as.numeric(colnames(db.asd.vnir)[2:ncol(db.asd.vnir)]), 
        t(db.asd.vnir[,2:ncol(db.asd.vnir)]),
        xlim = c(350, 2500),
        xaxt = 'n',
        frame.plot = FALSE,
        type = 'l',
        lty = 1,
        col = rgb(red = 0.0, green = 0.4, blue = 0.6, alpha = 0.4),
        xlab =  'nm', #expression(paste('Wavelength(cm'^'-1',')'))
        ylab = 'Reflectance',
        las=1, font=2,
        add = F)
axis(1, at = seq(350, 2500, by = 50), labels = T, las = 2)
mtext('Original ASD vis-NIR Spectra')
grid(lwd=1, nx=50, ny=50)
dev.off()


db.asd.vnir[1:5, c(1:5, 552, 1352, ncol(db.asd.vnir))]

jpeg(paste0(dir.l, '/results/original_asd_nir.jpeg'), width = 10, height = 10, units = 'in', res = 300)
matplot(as.numeric(colnames(db.asd.vnir)[552:1352]), 
        t(db.asd.vnir[,552:1352]),
        xlim = c(350, 2500),
        xaxt = 'n',
        frame.plot = FALSE,
        type = 'l',
        lty = 1,
        col = rgb(red = 0.0, green = 0.4, blue = 0.6, alpha = 0.4),
        xlab =  'nm', #expression(paste('Wavelength(cm'^'-1',')'))
        ylab = 'Reflectance',
        las=1, font=2,
        add = F)
axis(1, at = seq(350, 2500, by = 50), labels = T, las = 2)
mtext('Original ASD NIR Spectra')
grid(lwd=1, nx=50, ny=50)

dev.off()

# Save R envir. -----
save.image(paste0(dir.l, '/r_code/quick_graphs.RData'))
