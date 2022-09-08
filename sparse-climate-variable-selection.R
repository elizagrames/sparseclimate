# This script is for climate variable selection with Bayesian sparse modeling.

# Data files are curated versions of annual population-level responses and are
# derived from publicly available datasets archived at the DOIs below. 

# Climate data was generated using an alpha version of the climetric package
# (github.com/elizagrames/climetric) using TerraClimate data layers. 

# For questions, see the manuscript or contact Eliza Grames at egrames@unr.edu.

# McNulty 2018: 10.6073/pasta/39f495a15fa0667af635991d749741c4
# Ma et al. 2020: 10.5061/dryad.jdfn2z391
# Ramakers et al. 2018: 10.5061/dryad.35k1n3m
# Wadgymar et al. 2018: 10.5061/dryad.qr5vd
# O'Keefe 2021: 10.6073/pasta/91e3b7c2548a0f2e251729eeacbce312
# Poysa et al. 2019: 10.5061/dryad.b5mkkwh89
# Hinks et al. 2015: 10.5061/dryad.n6g3q
# Visser et al. 2021: 10.5061/dryad.f1vhhmgx6
# Werner et al. 2016: 10.6073/pasta/c174404b0bb5d9a65bc8eccb40db825c
# Huenneke and Browning 2016: 10.6073/pasta/bf98e185c512a8dff9677a611dc05455
# Donoso et al. 2015: 10.5061/dryad.72551
# Wiebe 2020: 10.5061/dryad.bk3j9kd71
# McLean et al. 2020: 10.5061/dryad.zs7h44j56
# Valtonen et al. 2017: 10.5061/dryad.9m6vp
# DeMay and Walters 2019: 10.5061/dryad.0qv86b5
# Cole et al. 2015: 10.5061/dryad.7v1qg
# Lightfoot 2021: 10.6073/pasta/8da81e91177dc6029d147990d1181be3
# Frigerio et al. 2021: 10.5061/dryad.np5hqbztd
# Weed et al. 2010: 10.5061/dryad.ms739

# Workspace set up -------------------------------------------------------------

library(tidyr)
library(ggplot2)
library(ggridges)

# Read in data -----------------------------------------------------------------

# First, generate a list of the data files in the directory
files <- list.files("./", pattern = "csv")

# Run the models to select variables for all datasets --------------------------

# Because we have 19 case studies, but we want summary metrics, we wrap
# everything up in a for loop so we can do the models for each case study and
# save the output to an array, so first we set up placeholder arrays

gammas <- betas <- array(dim = c(length(files), 80))
checks <- c()

for (i in (1:length(files))) {
  dat <- read.csv(files[i])
  covars <- apply(dat[, -c(1:4)], 2, scale)
  # Set up an object to hold all our data that we need for the model
  jags.data <- list(
    y = scale(dat$response)[, 1],
    years = scale(dat$year)[, 1],
    lambda = 0.1,
    covar = covars,
    n.covars = ncol(covars),
    n.obs = nrow(dat)
  )
  
  # Run the model in JAGS
  mod1 <- R2jags::jags(
    data = jags.data,
    parameters.to.save = c("beta0",
                           "beta",
                           "beta.year", 
                           "gamma", "y.sim", "rsq", "disc"),
    model.file = "./jags-horseshoe-ts.R",
    n.chains = 3,
    n.iter = 5000,
    n.burnin = 2000,
    n.thin = 3
  )
  
  summary <- mod1$BUGSoutput$summary
  checkR <- (grepl("beta0", rownames(summary)) | grepl("beta.year", rownames(summary)) | grepl("gamma", rownames(summary)))
  checks[i] <- (all(summary[checkR,'Rhat']<=1.1))
  
  # Save some output for our summary stats
  gammas[i, ] <- mod1$BUGSoutput$mean$gamma
  betas[i, ] <- mod1$BUGSoutput$mean$beta
  
  # We use the McNulty case study as an example plot, so let's save that
  if(grepl("mcnulty", files[i])){
    save(mod1, file="./mcnulty_for_plot.rda")
  }
    
}


# Summaries --------------------------------------------------------------------

# These aren't really of interest in the manuscript since the focus is the method
# But here are some general overviews of the betas 
mean(betas, na.rm=T)
min(betas, na.rm=T)
max(betas, na.rm=T)

# Summary plots ----------------------------------------------------------------

# Note: the following code is rather messy, but it gets the job done and gets 
# all the variables in the right order to be read sensibly in the plot. For the
# figure in the publication, we have done some post-hoc editing in Inkscape to
# make the figure easier to read (e.g. separating the heatmap into blocks)

colnames(gammas) <- colnames(betas) <- colnames(covars)
tmp <- tidyr::gather(as.data.frame(gammas))
tmp$study <- rep(gsub("_cleaned.csv", "", files), ncol(covars))
tmp$value <- tmp$value - .5 # To set 0 as the middle, since 0.5 is expected

# For the main plot, we only want proximate measures, not departures
g1 <- gammas[, grepl("prox", colnames(gammas))]
g1 <- t(g1)

# The following lines are just a kludgy way to add numbers to make sure we sort 
# the rows into the exact order we want for the plot
rownames(g1)[grepl("fall", rownames(g1))] <- paste("9", rownames(g1)[grepl("fall", rownames(g1))], sep = "_")
rownames(g1)[grepl("winter", rownames(g1))] <- paste("7", rownames(g1)[grepl("winter", rownames(g1))], sep = "_")
rownames(g1)[grepl("spring", rownames(g1))] <- paste("6", rownames(g1)[grepl("spring", rownames(g1))], sep = "_")
rownames(g1)[grepl("summer", rownames(g1))] <- paste("5", rownames(g1)[grepl("summer", rownames(g1))], sep = "_")
rownames(g1)[grepl("annual", rownames(g1))] <- paste("3", rownames(g1)[grepl("annual", rownames(g1))], sep = "_")

rownames(g1)[grepl("min.prox", rownames(g1))] <- paste("190", rownames(g1)[grepl("min.prox", rownames(g1))], sep = "")
rownames(g1)[grepl("tmin", rownames(g1)) & !grepl("min.prox", rownames(g1))] <- paste("180", rownames(g1)[grepl("tmin", rownames(g1)) & !grepl("min.prox", rownames(g1))], sep = "")
rownames(g1)[grepl("max.prox", rownames(g1))] <- paste("170", rownames(g1)[grepl("max.prox", rownames(g1))], sep = "")
rownames(g1)[grepl("tmax", rownames(g1)) & !grepl("max.prox", rownames(g1))] <- paste("160", rownames(g1)[grepl("tmax", rownames(g1)) & !grepl("max.prox", rownames(g1))], sep = "")

rownames(g1)[grepl("ppt", rownames(g1))] <- paste("150", rownames(g1)[grepl("ppt", rownames(g1))], sep = "")
rownames(g1)[grepl("soil", rownames(g1))] <- paste("140", rownames(g1)[grepl("soil", rownames(g1))], sep = "")
rownames(g1)[grepl("PDSI", rownames(g1))] <- paste("130", rownames(g1)[grepl("PDSI", rownames(g1))], sep = "")
rownames(g1)[grepl("vpd", rownames(g1))] <- paste("120", rownames(g1)[grepl("vpd", rownames(g1))], sep = "")

# Sort on our new row order
g1 <- g1[order(rownames(g1)), ]

pal <-  c( "#ffffff", "#F2FCFE", "#E5F8FC", "#D8F4FA", "#caf0f8", 
           "#90e0ef", "#00b4d8", "#0077b6", "#03045e")

heatmap(g1,
        scale = "none",
        col = colorRampPalette(pal)(100),
        breaks = seq(0, 1, length.out = 101),
        keep.dendro = F,
        Rowv = NA,
        labCol = gsub("_cleaned.csv", "", files),
        cexCol = 0.5,
        cexRow = 0.5
)

legend(
  "bottomright",
  legend = format(seq(0.1, 1, .1), nsmall = 2),
  col = colorRampPalette(pal)(100)[seq(10, 100, 10)],
  pch = 15,
  bty = "n",
  cex = 2,
  pt.cex = 4
)

# Example plot -----------------------------------------------------------------

# Load the JAGS object we saved when running the analysis to make the figure
load("./mcnulty_for_plot.rda")

# The following code is rather messy, but, it gets things in the same order
# as the main plot

# Gather the gamma and beta data into a nice format for plotting
tmp <- as.data.frame(mod1$BUGSoutput$sims.list$beta[,mod1$BUGSoutput$mean$gamma>0])
colnames(tmp) <- colnames(covars)[gammas[11,]>0]
tmp <- gather(tmp)
tmp$gamma <- rep(gammas[11,], each=3000)-.5

colnames(tmp) <- c("Measurement", "Estimate", "gamma")
tmp <- tmp[grepl("prox", tmp$Measurement),]

# We want things to be in a specific order to be consistent across studies
tmp$Measurement[grepl("fall", tmp$Measurement)] <- paste("9",tmp$Measurement[grepl("fall", tmp$Measurement)], sep="_")
tmp$Measurement[grepl("winter", tmp$Measurement)] <- paste("7",tmp$Measurement[grepl("winter", tmp$Measurement)], sep="_")
tmp$Measurement[grepl("spring", tmp$Measurement)] <- paste("6",tmp$Measurement[grepl("spring", tmp$Measurement)], sep="_")
tmp$Measurement[grepl("summer", tmp$Measurement)] <- paste("5",tmp$Measurement[grepl("summer", tmp$Measurement)], sep="_")
tmp$Measurement[grepl("annual", tmp$Measurement)] <- paste("3",tmp$Measurement[grepl("annual", tmp$Measurement)], sep="_")

tmp$Measurement[grepl("min.prox", tmp$Measurement)] <- paste("190",  tmp$Measurement[grepl("min.prox", tmp$Measurement)], sep="")
tmp$Measurement[grepl("tmin", tmp$Measurement) & !grepl("min.prox", tmp$Measurement)] <- paste("180",  tmp$Measurement[grepl("tmin", tmp$Measurement) & !grepl("min.prox", tmp$Measurement)], sep="")
tmp$Measurement[grepl("max.prox", tmp$Measurement)] <- paste("170", tmp$Measurement[grepl("max.prox", tmp$Measurement)], sep="")
tmp$Measurement[grepl("tmax", tmp$Measurement) & !grepl("max.prox", tmp$Measurement)] <- paste("160", tmp$Measurement[grepl("tmax", tmp$Measurement) & !grepl("max.prox", tmp$Measurement)], sep="")

tmp$Measurement[grepl("ppt", tmp$Measurement)] <- paste("150",  tmp$Measurement[grepl("ppt", tmp$Measurement)], sep="")
tmp$Measurement[grepl("soil", tmp$Measurement)] <- paste("140",  tmp$Measurement[grepl("soil", tmp$Measurement)], sep="")
tmp$Measurement[grepl("PDSI", tmp$Measurement)] <- paste("130",  tmp$Measurement[grepl("PDSI", tmp$Measurement)], sep="")
tmp$Measurement[grepl("vpd", tmp$Measurement)] <- paste("120", tmp$Measurement[grepl("vpd", tmp$Measurement)], sep="")

p1 <- ggplot(tmp, aes(x=Estimate, y=Measurement, fill=gamma)) +
  stat_density_ridges(rel_min_height=0.02,
                      color="black", lwd=0.1) +
  scale_fill_gradient2(low="#FF7F11",
                       high="#0077B6",
                       mid="white",
                       na.value = "grey96",
                       limits=c(-.5, .5)) + 
  xlim(-1,1) + theme(axis.title.y=element_blank(),
                     axis.ticks.y = element_blank(),
                     legend.background =element_blank(),
                     panel.background = element_blank(),
                     panel.grid.major.y = element_line(size=.1, colour="black"),
                     panel.grid.major.x=element_blank())
print(p1)

