library(SpaDES)
library(terra)
library(RColorBrewer)
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
myPlot <- function(spatobj, names = NULL, labelCols = NULL) {
  if (!is.list(spatobj)){
    spatobj <- list(spatobj)
  }
  
  # Build the list of layer names
  namesV <- c()
  if (is.null(names)){
    for (i in seq_along(spatobj)) {
      if ("SpatVector" %in% class(spatobj[[i]]) || 
          is.null(names(spatobj[[i]])) || 
          all(names(spatobj[[i]]) == "") || 
          length(names(spatobj[[i]])) > 1
         ){
        namesV <- c(namesV, paste("Layer", i))
      }
      else {
        namesV <- c(namesV, names(spatobj[[i]]))
      }
    }
  } 
  else if (is.list(names)){
    namesV <- unlist(names)
  }
  
  # Create base leaflet map
  m <- leaflet()
  m <- addTiles(m)
  
  # Add layer control
  overlay_groups <- if (!is.null(namesV)) namesV else paste0("Layer_", seq_along(spatobj))
  m <- addLayersControl(m, 
                        overlayGroups = overlay_groups,
                        options = layersControlOptions(collapsed = FALSE))
  # Loop over spatobj
  spatVectorIdx <- 0
  for (i in seq_along(spatobj)) {
    sObj <- spatobj[[i]]
    # Set layer name
    layer_name <- if (!is.null(namesV) && length(namesV) >= i) namesV[i] else paste0("Layer_", i)
    
    if ("SpatRaster" %in% class(sObj)){
      # resample the raster if it is too big for leafem/leaflet to handle
      ctab <- coltab(sObj)[[1]]
      aggMethod <- terra::mean
      bytePerPixel <- 8
      # if it's a raster of integer change the aggregate method to modal
      if (is.factor(sObj) || length(unique(values(all(floor(sObj) == sObj), na.rm = TRUE))) == 1){
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
      if (is.factor(sObj)){
        vals <- levels(sObj)[[1]][[2]]
      }
      
      if (is.null(ctab)){
        ctab <- data.frame(vals, getRandomPalette(length(vals)))
      }
      if (length(vals) > 20){
        pal <- c("#FFFFFF", getRandomPalette(hex = TRUE, dark = TRUE))
        pal <- colorNumeric(palette = pal, domain = ctab$vals, na.color = "transparent")
        # Add spatObj as overlay
        m <- addRasterImage(m, sObj, colors = pal, opacity = 0.7, group = layer_name, layerId = layer_name)
        m <- addLegend(m, pal = pal, values = values(sObj), title = layer_name, group = layer_name)
      }
      else {
        # Add spatObj as overlay
        colors <- rgb(ctab$red, ctab$green, ctab$blue, maxColorValue = 255)
        labels <- ctab$vals
        # if (is.factor(sObj)){
        #   labels <- levels(sObj)[[1]][[2]]
        # }
        m <- addRasterImage(m, sObj, colors = colors, opacity = 0.7, group = layer_name, layerId = layer_name)
        m <- addLegend(m, colors = colors, values = values(sObj), labels = labels, title = layer_name, group = layer_name)
      }
      
      m <- addImageQuery(m, sObj, band = 1, group = layer_name, project = TRUE, type = c("mousemove", "click"), position = "topright", prefix = "Layer"
      )  
    }
    else if ("SpatVector" %in% class(sObj)){
      spatVectorIdx <- spatVectorIdx + 1
      sObj <- project(sObj, "EPSG:4326")
      m <- addGeoJSON(
        m, 
        sf_geojson(st_as_sf(sObj)),
        color = rgb(sample(0:255, 1),
                    sample(0:255, 1),
                    sample(0:255, 1),
                    maxColorValue = 255), 
        weight = 2, 
        fill = FALSE, 
        group = layer_name
      )
      
      # add labels if requested
      if (!is.null(labelCols[spatVectorIdx]) && labelCols[spatVectorIdx] %in% namesV(sObj)){
        # we will put each label on the centroid of each polygon
        labelPts <- st_centroid(sf::st_as_sf(sObj))
        # we create a temp column to be used as labels
        labelPts$temp_name <- labelPts[[labelCols[spatVectorIdx]]]
        m <- addLabelOnlyMarkers(
          m,
          data = labelPts,
          #label = ~ .[[labelCols[spatVectorIdx]]], # doesn't work
          label = ~temp_name,
          labelOptions = labelOptions(
            noHide = TRUE, 
            direction = "center",
            style = list(
              "padding" = "1px",
              "margin" = "0px",
              "box-shadow" = "none",
              "border" = "none",
              "background-color" = "transparent",
              "font-weight" = "bold",
              "color" = "black",
              "text-shadow" = paste0(
                "1px 1px 0 rgba(255,255,255,0.6), ",
                "-1px 1px 0 rgba(255,255,255,0.6), ",
                "1px -1px 0 rgba(255,255,255,0.6), ",
                "-1px -1px 0 rgba(255,255,255,0.6)"
              )
            )
          )
        )
      }
    }
  }

  return(m)
}

###################################################################
# mapViewCol
# Build a random color map before passing the raster to mapView
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

  return(mv)
}

###################################################################
# getRandomPixelGroupMap
###################################################################
getRandomCategoricalMap<- function(origin = c(1541912, 1072021),
                       ncol = 1000,
                       nrow = 1000,
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
    ncols = ncol, 
    nrows = nrow, 
    xmin = origin[1], 
    xmax = origin[1] + ncol * pixelsize, 
    ymin = origin[2], 
    ymax = origin[2] + nrow * pixelsize,
    vals = 0
  )
  # mapView(tempRast)
  
  rast <- randomPolygons(ras = tempRast, numTypes = nbregion)
  
  if (!is.null(valuevect)){
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

######################################
# Packages management
######################################
installedPackages <- function(){
  if (!requireNamespace("stringr", quietly = TRUE)) {
    install.packages("stringr")
  }
  library("stringr")
  pkgs <- installed.packages()
  pkg_names <- pkgs[, "Package"]
  
  get_source <- function(pkg) {
    desc <- tryCatch(utils::packageDescription(pkg), error = function(e) return(NULL))
    if (is.null(desc)) return(NA)
    
    # Check common fields
    if (!is.null(desc$Repository)) return(desc$Repository)
    if (!is.null(desc$RemoteType)) return(desc$RemoteType)
    if (!is.null(desc$biocViews)) return("Bioconductor")
    if (!is.null(desc$RemoteUrl)) return(desc$RemoteUrl)
    
    return("Unknown")
  }
  
  sources <- sapply(pkg_names, get_source)
  
  str_trunk_both_ends <- function(str, start_length=10, end_length=10){
    ret = stringr::str_trunc(str, start_length)
    if (nchar(str)[1] > start_length) {
      ret = paste(ret, substring(stringr::str_trunc(str, end_length, "left"), 4), sep="")
    }
    return(ret)
  }
  
  # Build data frame
  df <- data.frame(
    Package = pkg_names,
    Version = pkgs[, "Version"],
    LibPath = str_trunk_both_ends(pkgs[, "LibPath"], 10, 20),
    Source = str_trunk_both_ends(sub("^https?://", "", sources), 10, 20),
    stringsAsFactors = FALSE
  )
  
  #df_ordered <- df[,order(rownames(df)), ]
  df_ordered <- df[order(df$Package), ]
  print(df_ordered)
}

installMyPackages <- function(){
  pkgs <- c("mapview", "stars", "NLMR")
  for (pkg in pkgs) {
    message(paste("Installing and loading", pkg, "..."))
    install.packages(pkg, 
                     dependencies=TRUE, 
                     repos=c("https://predictiveecology.r-universe.dev", getOption("repos")))
    library(pkg, character.only = TRUE)
  }
  # remotes::install_github("cran/RandomFieldsUtils")
  # remotes::install_github("cran/RandomFields")
  remotes::install_github("ropensci/NLMR", dependencies = TRUE, upgrade = "never",force = TRUE)
}

uninstallAllPackages <- function (){
  pkgs <- installed.packages()
  pkgs
  user_pkgs <- pkgs[!pkgs[, "Priority"] %in% c("base", "recommended"), "Package"]
  
  # Confirm before running the next line
  remove.packages(user_pkgs)
}

clearRam <- function(){
  rm(list=ls())
}

clearRAM <- function(){
  clearRam()
}

unloadAllPackages <- function() {
  base_pkgs <- c(
    "base", "compiler", "datasets", "graphics", "grDevices",
    "grid", "methods", "parallel", "splines", "stats",
    "stats4", "tools", "utils"
  )
  
  loaded <- loadedNamespaces()
  user_pkgs <- setdiff(loaded, base_pkgs)
  
  for (pkg in user_pkgs) {
    message(pkg)
    try({
      # Detach from search path if attached
      if (paste0("package:", pkg) %in% search()) {
        message("Detaching...")
        detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
      }
      # Unload namespace
      message(paste("Unloading..."))
      
      unloadNamespace(pkg)
    }, silent = FALSE)
  }
  
  invisible(NULL)
}

###################################################
# Apply a function to every file in a folder tree
###################################################
library(purrr)
# Apply a function to every file in a folder tree
apply_to_files <- function(source_dir, regex = NULL, recursive = TRUE, fn, ...) {
  files <- list.files(
    path = source_dir,
    pattern = regex,
    recursive = recursive,
    full.names = TRUE
  )
  
  walk(files, fn, source_dir = source_dir, ...)
}

# Apply a function to every file in a folder tree in parallel
library(future.apply)
apply_to_files_p <- function(source_dir, regex = NULL, recursive = TRUE, fn, ...) {
  plan(multisession, workers = availableCores() - 1)  # cross-platform
  files <- list.files(
    path = source_dir,
    pattern = regex,
    recursive = recursive,
    full.names = TRUE
  )
  
  future_lapply(files, future.seed=NULL, FUN = fn, source_dir = source_dir, ...)
}

# get a new file name appended with "_X" if the file already exist
get_unique_filename <- function(path) {
  dir  <- dirname(path)
  name <- tools::file_path_sans_ext(basename(path))
  ext  <- tools::file_ext(path)
  
  # Build full file name
  new_path <- path
  counter <- 1
  
  while (file.exists(new_path)) {
    suffix <- paste0("_", counter)
    new_file <- paste0(name, suffix)
    
    if (nzchar(ext)) {
      new_path <- file.path(dir, paste0(new_file, ".", ext))
    } else {
      new_path <- file.path(dir, new_file)
    }
    
    counter <- counter + 1
  }
  
  new_path
}

# get a new file path in target_dir respecting the sub path of the source
new_file_path <- function(full_filepath, source_dir, target_dir = NULL, file_suffix = "", ext = NULL){
  if (is.null(target_dir)){
    target_dir <- source_dir
  }
  if (is.null(ext)){
    ext <- tools::file_ext(full_filepath)
  }
  if (!startsWith(ext, ".")){
    ext <- paste0(".", ext)
  }
  sub_dir <- dirname(sub(paste0("^", source_dir), "", full_filepath))
  filename <- tools::file_path_sans_ext(basename(full_filepath))
  filename <- paste0(filename, file_suffix, ext)
  newpath <- file.path(target_dir, sub_dir, filename)
  return(get_unique_filename(newpath))
}
# test
# undebug(new_file_path)
# new_file_path("c:/temp/sub/a.txt", "c:/temp")
# new_file_path("c:/temp/sub/a.txt", "c:/temp", "c:/temp2/")

# 
# testfnc <- function(found_file, source_dir, target_dir = NULL, file_suffix = NULL, ext = NULL) {
#   new_full_path <- new_file_path(found_file, source_dir, target_dir, file_suffix, ext)
#   # create the target dir if it does not exist
#   dir_name <- dirname(new_full_path)
#   if (!dir.exists(dir_name)) {
#     dir.create(dir_name, recursive = TRUE, showWarnings = FALSE)
#   }
#   
#   file.create(new_full_path)
# }
# test testfnc
# undebug(testfnc)
# testfnc("c:/temp/a/a1.txt", "c:/temp")
# testfnc("c:/temp/a/a1.txt", "c:/temp", "c:/temp2")

# test apply_to_files
# undebug(apply_to_files)
# apply_to_files(source_dir = "c:/temp", regex = "^[ab]1.*", fn = testfnc, target_dir = "c:/temp2", file_suffix = "_next")

# apply_to_files(source_dir = "G:/Home/FromMaria", regex = "^pixelGroupMap.*", fn = testfnc, target_dir = "G:/Home/FromMariaOutput", file_suffix = "_next", ext="tif")

parallel_chunks <- function(data, FUN, workers = NULL, nchunks = NULL, ...) {
  browser()
  # Detect workers if not specified
  if (is.null(workers)) {
    workers <- availableCores()
  }
  plan(multisession, workers = workers)
  
  n <- if (is.matrix(data) || is.data.frame(data)) {
    nrow(data)
  } else {
    length(data)
  }
  
  # Decide number of chunks (default = #workers)
  if (is.null(nchunks)) nchunks <- workers
  
  # Create index chunks (lightweight!)
  idx_chunks <- split(seq_len(n), cut(seq_len(n), nchunks, labels = FALSE))
  
  # Run in parallel
  res <- future_lapply(idx_chunks, function(idx, data=data) {
    FUN(data[idx, , drop = FALSE], ...)
  },
  future.globals = FALSE)
  
  return(res)
}
######################################
# List all the sub-object of an object with their classes
######################################

ls_with_class <- function(x, class_filter = NULL) {
  out <- sapply(x, function(obj) paste(class(obj), collapse = ", "))
  
  if (!is.null(class_filter)) {
    class_filter <- tolower(class_filter)
    
    keep <- vapply(x, function(obj) {
      any(tolower(class(obj)) %in% class_filter)
    }, logical(1))
    
    out <- out[keep]
  }
  
  cat(paste0(names(out), " (", out, ")"), sep = "\n")
}

######################################
# SpaDES management
######################################
resetSpades <- function (){
  #rm(list = ls())
  #unloadAllPackages()
  uninstallAllPackages()
  install.packages(c("SpaDES", "reproducible", "LandR", "SpaDES.project"), 
                   repos=c("predictiveecology.r-universe.dev", getOption("repos")), 
                   dependencies=TRUE)
  library(SpaDES)
  library(SpaDES.core)
  library(SpaDES.project)
  library(LandR)
  #  library(reproducible)
}

setBasePath <- function(basePath){
  setPaths(modulePath = file.path(basePath, "modules"),
           cachePath = file.path(basePath, "cache"),
           inputPath = file.path(basePath, "input"),
           outputPath = file.path(basePath, "output"),
           rasterPath = file.path(basePath, "raster"),
           scratchPath = file.path(basePath, "scratch"),
           terraPath = file.path(basePath, "terra"))
}


