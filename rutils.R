library(terra)
library(RColorBrewer)
library(SpaDES.tools)
library(mapview)
library(leaflet)
library(leafem)
library(geojsonsf)
library(sf)

###################################################################
# getRandomPalette
###################################################################
getRandomPalette <- function(nbcol = 1, hex = FALSE, dark = FALSE, pale = FALSE){
  # hues <- c(runif(1, 0, 180), runif(1, 180, 360)) # random hue for the whole set of polygons
  base_hue <- runif(1, 0, 360)
  # Each polygon gets a hue close to the base hue
  hues <- (base_hue + runif(nbcol, -30, 30)) %% 360
  chromas <- runif(nbcol, 0, 100)
  lmin <- 0
  lmax <- 100
  if (dark && !pale) {
    lmin <- 0
    lmax <- 20
  }
  if (!dark && pale) {
    lmin <- 80
    lmax <- 100
  }
  luminances <- runif(nbcol, lmin, lmax)
  cols <- hcl(h = hues, c = chromas, l = luminances)
  if (hex) return(cols)
  return(t(col2rgb(cols)))
}

# debug(getRandomPalette)
# getRandomPalette(2)

###################################################################
# myPlot
###################################################################
myPlot <- function(spatobj, names = NULL) {
  if (!is.list(spatobj)){
    spatobj <- list(spatobj)
  }
  
  if (is.list(names)){
    names <- unlist(names)
  }
  
  # Create base leaflet map
  m <- leaflet()
  m <- addTiles(m)
  
  # Add layer control
  overlay_groups <- if (!is.null(names)) names else paste0("Layer_", seq_along(spatobj))
  m <- addLayersControl(m, 
                        overlayGroups = overlay_groups,
                        options = layersControlOptions(collapsed = FALSE))
  # Loop over spatObj
  for (i in seq_along(spatobj)) {
    sObj <- spatobj[[i]]
    # Set layer name
    layer_name <- if (!is.null(names) && length(names) >= i) names[i] else paste0("Layer_", i)
    
    if ("SpatRaster" %in% class(sObj)){
      # resample the raster if it is too big for leafem/leaflet to handle
      ctab <- coltab(sObj)[[1]]
      aggMethod <- terra::mean
      bytePerPixel <- 8
      if (is.factor(sObj)){
        aggMethod <- terra::modal
        bytePerPixel <- 4
      }
      rastSize <- ncell(sObj) * bytePerPixel
      while (rastSize > 4194304){
        message("Downscaling ", names(sObj)[1], " by a factor of 2 so leaflet can handle it...")
        sObj <- aggregate(sObj, fact = 2, fun = aggMethod, na.rm=TRUE)
        rastSize <- ncell(sObj) * bytePerPixel
      }

      # Convert raster values to colors
      # if the raster has no palette build a random one
      vals <- sort(unique(values(sObj))[, 1])
      browser()
      if (is.null(ctab)){
        ctab <- data.frame(vals, getRandomPalette(length(vals)))
      }
      if (length(vals) > 20){
        pal <- c("#FFFFFF", getRandomPalette(hex = TRUE, dark = TRUE))
        pal <- colorNumeric(palette = pal, domain = ctab$vals, na.color = "transparent")
      }
      else {
        pal <- colorFactor(palette = rgb(ctab$red, ctab$green, ctab$blue, maxColorValue = 255), domain = ctab$vals, na.color = "transparent")
      }
      
      # Add spatObj as overlay
      m <- addRasterImage(m, sObj, colors = pal, opacity = 0.7, group = layer_name, layerId = layer_name)
      m <- addLegend(m, pal = pal, values = values(sObj), title = layer_name, group = layer_name)
      m <- addImageQuery(m, sObj, band = 1, group = layer_name, project = TRUE, type = c("mousemove", "click"), position = "topright", prefix = "Layer"
      )  
    }
    else if ("SpatVector" %in% class(sObj)){
      sObj <- project(sObj, "EPSG:4326")
      m <- addGeoJSON(m, sf_geojson(st_as_sf(sObj)), color = rgb(sample(0:255, 1), sample(0:255, 1), sample(0:255, 1), maxColorValue = 255), weight = 2, fill = FALSE, group = layer_name)
    }
  }

  return(m)
}

###################################################################
# mapViewCol
###################################################################
mapViewCol <-  function(raster, ...) {
  ct <- coltab(raster)
  
  if (is.null(ct)) {
    warning("Raster has no coltab. Using default mapview colors.")
    return(mapview(raster, ...))
  }
  
  if (is.list(ct)) {
    ct <- ct[[1]]  # take first layer
  }
  
  if (inherits(ct, "data.frame")) {
    ct_vals <- ct[, 1]
    rgb_cols <- apply(ct[, 2:4], 1, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
  } else {
    ct_vals <- ct[, 1]
    rgb_cols <- apply(ct, 1, function(x) rgb(x[2], x[3], x[4], maxColorValue = 255))
  }
  
  # Match raster values to colors
  vals <- sort(unique(values(raster)))
  cols <- rgb_cols[match(vals, ct_vals)]
  
  # Plot
  mv <- mapView(raster, col.regions = cols, at = vals, ...)
  # if (!inherits(mv, "Map")) {
  #   class(mv) <- c("leaflet", class(mv))
  # }

  return(mv)
  # mapView(raster, col.regions = cols, at = vals, ...)
}

###################################################################
# getRandomPixelGroupMap
###################################################################
getRandomCategoricalMap<- function(origin = c(1541912, 1072021),
                       width = 1000,
                       pixelsize = 250,
                       crs = "ESRI:102002",
                       nbregion = 200,
                       valuevect = NULL,
                       seed = NULL){
  if (!is.null(seed)){
    set.seed(seed)
  }
  tempRast <- terra::rast(
    crs = crs,
    nrows = width, 
    ncols = width, 
    xmin = origin[1], 
    xmax = origin[1] + width * pixelsize, 
    ymin = origin[2], 
    ymax = origin[2] + width * pixelsize,
    vals = 0
  )
  # mapView(tempRast)
  
  rast <- randomPolygons(ras = tempRast, numTypes = nbregion)
  
  if (!is.null(valuevect)){
    #browser()
    valuevectCnt <- length(unique(valuevect))
    nbCol <- valuevectCnt
    if (nbregion < valuevectCnt){
      message("getRandomPixelGroupMap: Not enough different polygons to assign every values in valueVect")
    }
    else {
      breaks <- seq(0, nbregion, length.out = valuevectCnt + 1)
      reclass_matrix <- cbind(
        from = breaks[-length(breaks)],   # all except last
        to   = breaks[-1],                # all except first
        new  = sort(unique(valuevect))    # new class values
      )
      rast <- terra::classify(rast, rcl = reclass_matrix, include.lowest = TRUE)
    }
  }
  
  vals <- sort(unique(values(rast)))
  vals <- vals[!is.na(vals)]
  
  cols <- getRandomPalette(length(vals))
  coltab(rast) <- cbind(vals, cols)
  
  return(rast)
}
###################################################################
# getRandomCohortData
#  Return a randomly generated cohortdata table
###################################################################
getRandomCohortData <- function(nbPixelGroup, pixelSize, seed = NULL){
  set.seed(seed)
  # Restrict the number of possible species
  maxNbSpecies <- 7
  species <- as.factor(c("Pice_mar", "Pinu_ban", "Popu_tre", "Betu_pap", 
                         "Pice_gla", "Abie_bal", "Acer_rub"))
  speciesProb <- c(0.25, 0.18, 0.17, 0.12, 0.1, 0.1, 0.08)
  speciesMaxAge <- c(275, 175, 90, 110, 275, 225, 175)
  nbSpecies <- min(maxNbSpecies, length(species))
  species <- head(species, nbSpecies)
  speciesProb <- head(speciesProb, nbSpecies)
  speciesMaxAge <- head(speciesMaxAge, nbSpecies)
  
  # Set
  nbSpPerGroupWeights <- c(0.2, 0.3, 0.3, 0.15, 0.05)
  minNbSpPerGroup <- 1
  maxNbSpPerGroup <- min(length(nbSpPerGroupWeights), nbSpecies)
  
  # sum the extra weight together if the number of species is smaller than the 
  # maximum number of species per group
  while (length(nbSpPerGroupWeights) > maxNbSpPerGroup) {
    l <- length(nbSpPerGroupWeights)
    nbSpPerGroupWeights[l - 1] <- nbSpPerGroupWeights[l - 1] + nbSpPerGroupWeights[l]
    nbSpPerGroupWeights <- nbSpPerGroupWeights[-l]
  }
  
  # number of species per group
  repeats <- sample(minNbSpPerGroup:maxNbSpPerGroup, 
                    nbPixelGroup, 
                    replace = TRUE, 
                    prob=nbSpPerGroupWeights)
  #hist(repeats)
  
  # define constants needed to generate random age and biomass
  sdlog <- 0.6
  maxBiomassConst <- 100000 # g/m2 or 100 kg/m2
  growthRate <- 0.05
  inflectionAgeConst <- 0.25
  
  # This function generates a random biomass estimate using a logistic growth model.
  # The estimate depends on the age of the tree and the specified maximum age.
  randomBiomass <- function(ageV, ageMaxV){
    maxBiomassConst / (1 + exp(-growthRate * (ageV - inflectionAgeConst * ageMaxV)))
  }
  
  # # test
  # ageMax <- 300
  # plot(
  #   1:ageMax,
  #   randomBiomass(1:ageMax, rep(ageMax, ageMax))
  # )

  #debug(bio)
  
  # Define a function that generates random ages following a log-normal distribution.
  randomAge <- function(spCode){
    mean_target <- speciesMaxAge[spCode] / 4
    meanlog <- log(mean_target) - (sdlog ^ 2) / 2
    min(round(rlnorm(1, meanlog, sdlog)), 1.25 * speciesMaxAge[spCode])
  }
  
  # test
  # hist(replicate(10000, randomAge(1)))
  
  # build the data table
  pixelGroup <- c()
  speciesCode <- c()
  age <- c()
  ageMax <- c()
  spProb <- c()
  bPerPixel <- c()
  b <-  c()
  
  # for each group
  for (i in seq_along(1:nbPixelGroup)) {
    n <- repeats[i]
    
    pixelGroup <- c(pixelGroup, rep(i, n))
  
    # pick n unique numbers between 1 and 5
    newSpCodes <- sort(sample(1:nbSpecies, size=n, replace=FALSE, prob=speciesProb))
    speciesCode <- c(speciesCode, newSpCodes)
  
    # random age 
    newAges <- sapply(newSpCodes, randomAge)
    age <- c(age, newAges)
    
    # add ageMax to the table
    newAgeMax <- speciesMaxAge[newSpCodes]
    ageMax <-  c(ageMax, newAgeMax)
    
    # add species probability
    newSpProb <- speciesProb[newSpCodes]
    spProb <- c(spProb, newSpProb)
    
    # estimated biomass per pixel
    newBPerPixel <- randomBiomass(newAges, newAgeMax) * pixelSize * pixelSize
    bPerPixel <- c(bPerPixel, newBPerPixel)
    
    # multiply biomass per pixel by the prob of finding the species in the cohort
    newB <- round(newBPerPixel * newSpProb)
    b <- c(b, newB)
  }
  
  # Combine into a data frame
  dt <- data.table(pixelGroup = pixelGroup, 
                   speciesCode = factor(speciesCode, labels = species),
                   age = age,
                   ageMax = ageMax,
                   spProb = spProb,
                   bPerPixel = bPerPixel,
                   B = b)
  # plot(
  #   dt[as.integer(speciesCode) == 1]$age,
  #   bio(dt[as.integer(speciesCode) == 1]$age)
  # )
  # 
  
  return(dt)
}

#getRandomCohortData(100, 250, 20)

# old methods
# set.seed(42)
# species_levels <- c("Pinu_ban", "Pice_mar", "Betu_pap", "Popu_tre", "Acer_rub")
# cohortData <- data.table(do.call(rbind, lapply(pixel_groups, function(pg) {
#   n_cohorts <- sample(1:3, 1)  # 1–3 cohorts per pixel group
#   data.frame(
#     pixelGroup = as.integer(pg),
#     speciesCode = factor(
#       sample(species_levels, n_cohorts, replace = FALSE),
#       levels = species_levels
#     ),
#     age = as.integer(sample(seq(5, 150, by = 5), n_cohorts, replace = TRUE)),
#     B = as.integer((round(runif(n_cohorts, 0.5, 10), 2))),
#     ecoregionGroup = as.factor("1_09")
#     #mortality = round(runif(n_cohorts, 0, 0.5), 2)
#   )
# })))
# 
# 
# set.seed(42)
# cohortData <- data.table(do.call(rbind, lapply(pixel_groups, function(pg) {
#   n_cohorts <- sample(2:3, 1)  # 1–3 cohorts per pixel group
#   data.frame(
#     pixelGroup = as.integer(pg),
#     speciesCode = factor(
#       sample(species_levels, n_cohorts, replace = FALSE),
#       levels = species_levels
#     ),
#     age = as.integer(sample(seq(5, 150, by = 5), n_cohorts, replace = TRUE)),
#     B = as.integer((round(runif(n_cohorts, 0.5, 10), 2))),
#     ecoregionGroup = as.factor("1_09")
#     #mortality = round(runif(n_cohorts, 0, 0.5), 2)
#   )
# })))
