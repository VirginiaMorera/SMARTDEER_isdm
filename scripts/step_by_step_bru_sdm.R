##---------------------------------------##
#### Load objects needed for modelling ####
##---------------------------------------##

in_bound <- readRDS("data/inner_boundary.RDS")
mesh0 <- readRDS("data/mesh.RDS")
covar_stack <- readRDS("data/covar_stack.RDS")
PO_data <- read.csv("data/PO_data_RD.csv", row.names = NULL)
PA_data <- read.csv("data/PA_data_RD.csv", row.names = NULL)

library(inlabru)
library(INLA)

##-----------------------------------##
#### source to load all parameters ####
##-----------------------------------##
source("scripts/params_bru_sdm.R")

##-----------------------------------------------------##
#### Initial checks that everything we need is there ####
##-----------------------------------------------------##

if (is.null(spatialcovariates)) stop("Spatial covariates not provided.")

if (is.null(poresp) | is.null(paresp)) stop("Either the precense only or the precense absence response is null.")

if (poresp == paresp) stop("Please provide different names for presence only and presence absence datasets.")

if (length(paresp) > 1 | length(poresp) > 1) stop("More than one name given for presences column.")

if (length(trialname) > 1) stop("More than one name given for number of trials column.")

if (length(coords) != 2) stop("Coordinates must have two components.")

if(ncol(spatialcovariates) > 1 & is.null(tol)) stop("Tolerance parameter not provided.")

if (as.character(proj@projargs) == "+proj=longlat +ellps=WGS84 +no_defs") warning("Default CRS is being used. Please change if incorrect.")

if (!is.null(poformula)) {
  
  if (as.character(poformula[2]) != poresp) stop(paste("Response variable for presence only datasets should be: ", poresp,'.', sep = ""))
  
}

if(!is.null(paformula)) {
  
  if (as.character(paformula[2]) != paresp) stop(paste("Response variable for presence absence datasets should be: ", paresp,'.', sep = ""))    
}

if (class(spatialcovariates) == 'data.frame' & is.null(tol)) stop('Please provide a tolerance parameter to convert the spatial covariates to a SpatialPixelsDataFrame.')

if (class(proj) != 'CRS') stop("Proj needs to be a CRS object.")

if (is.null(mesh) & is.null(meshpars)) stop("Either a mesh or mesh parameters need to be provided.")

if (is.null(mesh)) {
  
  if (sum(names(meshpars)%in%c( "cutoff", "max.edge", "offset")) < 3)
    
    stop("Meshpars requires three items in the list: cut.off, max.edge and offset.")
  
}

if (!marks) {
 # if there are no marks, this will make empty (NULL) objects named names_marks and data_marks  
  names_marks <- NULL 
  data_marks <- NULL
  multinom_vars <- NULL
  marksspatial <- FALSE
  
}

if (!marks & !is.null(inclmarks)) {
  
  warning('Marks to include is non null but include marks is set to FALSE.\nMarks are now being set to TRUE.')
  marks <- TRUE
  
}

if (!is.null(residuals)) {
  if (!residuals%in%c('response','pearson','deviance')) {
    stop("Residuals needs to be one of: 'response', 'pearson' or 'deviance'.")
  }
}

##--------------------------------------------##
#### provide data and prepare for modelling ####
##--------------------------------------------##

# datasets = list(...)

# po and pa datasets, hopefully order is irrelevant? 
datasets <- list(PA_data, PO_data)

# check that all datasets are either data.frames or SpatialPoints* object
datasets_class = sapply(datasets, class)

if (any(!datasets_class%in%c('SpatialPointsDataFrame','SpatialPoints', 'data.frame'))) {
  
  stop('Datasets need to be either a SpatialPoints* object or a data frame.')
  
}

# are the coords names (X,Y) in the objects in "datasets"? Boolean vector with TRUE TRUE TRUE TRUE
coords_in = unlist(lapply(datasets, function(dat) {
  
  if (class(dat) == 'data.frame') coords%in%names(dat)
  else
    if (inherits(dat, 'Spatial')) {
      x_coord <- colnames(dat@coords)[1]
      y_coord <- colnames(dat@coords)[2]
      coords%in%c(x_coord,y_coord)
    }
  
}))

if (!all(coords_in)) stop("At least one dataset does not have coordinates in it.\nEither check your datasets or change your coordinates argument.")


# This only works when run inside the function, but its output is a character vector with the names of the datasets in the order they were entered in the function call, so we'll reproduce it here


# data_names <- setdiff(as.character(match.call(expand.dots=TRUE)), 
#                       as.character(match.call(expand.dots=FALSE)))

data_names <- c("PA_data", "PO_data")


#Separate PO and PA data by inclusion/exclusion of 'trialname'.

data_attributes <- lapply(datasets,function(dat) {
  if (inherits(dat,"Spatial")) {
    if (class(dat) == "SpatialPoints") {

      dat <- sp::SpatialPointsDataFrame(coords = sp::coordinates(dat),
                                        data = data.frame(resp = rep(1,nrow(coordinates(dat)))),
                                        proj4string = proj)
      names(dat) <- poresp
      attr(dat,'family') <- 'poisson'
      attr(dat,'data_type') <- 'Present only'
      dat

    }
    else #if class == SpatialPointsDataFrame
      if (paresp%in%colnames(dat@data)){

        dat <- sp::SpatialPointsDataFrame(coords = dat@coords,
                                          data = as.data.frame(dat@data),
                                          proj4string = proj)
        dat@data[,paresp] <- as.numeric(dat@data[,paresp])
        if (!is.null(trialname)) {
          if (trialname%in%colnames(dat@data)) attr(dat,'Ntrials') <- dat@data[,trialname]
          else attr(dat,'Ntrials') <- 1
        }
        attr(dat,'family') <- 'binomial'
        attr(dat,'data_type') <- 'Present absent'
        dat

      }
    else {

      dat <- sp::SpatialPointsDataFrame(coords = dat@coords,
                                        data = as.data.frame(dat@data),
                                        proj4string = proj)
      if (!poresp%in%colnames(dat@data)) {
        dat@data[,poresp] <- 1

      }
      attr(dat,'family') <- 'poisson'
      attr(dat,'data_type') <- 'Present only'
      dat

    }

  }

  
  else # if data is not a spatial object
    if (class(dat) == 'data.frame') {
      if (paresp%in%colnames(dat)) { # this will happen if paresp ("PA") is a column in the dataset, so it will only happen for PA_data
        
        ## check that the dataset has enough columns (at least 3, coords and PA)
          if (ncol(dat) == 1) stop("Only one column provided in data frame.\nNeed two columns for coordinates and one for presence name\nfor presence/absence data.")
        else
          if (ncol(dat) == 2) stop("Only two columns provided for presence/absence data. Either coordinates is of length one or presence name not given.")
        else {
          
          # names of columns that are not coordinates
          names <- names(dat)[!names(dat)%in%c(coords)]
          
          # turn dataset into SpatialPointsDataFrame 
          # all the many steps seem awkward but they're probably necessary for the function to be generalisable
          
          dat <- sp::SpatialPointsDataFrame(coords = dat[,coords],
                                            data = as.data.frame(dat[,names]),
                                            # data = as.data.frame(dat[,!names(dat)%in%coords]), # see, we change this line to the much simpler line above and it still works
                                            proj4string = proj) #dat[,!names(x)%in%coords]
          colnames(dat@data) <- names
          
          # turn PA column to numeric
          dat@data[,paresp] <- as.numeric(dat@data[,paresp])
          
          # this bit has to do with repeated trials (I think), and we don't have repeated sampling, so commenting out for clarity 
          if (!is.null(trialname)) {
            if (trialname%in%colnames(dat@data)) attr(dat,'Ntrials') <- dat@data[,trialname]
            else attr(dat,'Ntrials') <- 1
          }
          
          # No idea what this is doing 
          attr(dat,'family') <- 'binomial'
          attr(dat,'data_type') <- 'Present absent'
          dat
          
        }
      }
      else { # so, the following will be done to the PO dataaset
        # in this case we'd only need 2 columns, one for each coordinate, since this is Presence Only
        if (ncol(dat) == 1) stop('Data Frame provided only has one column. Coordinates require two columns.')
        else
          if (ncol(dat) == 2){ # this will happen if the po dataset only has two columns. ours has more
            
            dat <- sp::SpatialPointsDataFrame(coords = dat[,coords],
                                              data = data.frame(resp = rep(1,nrow(dat))),
                                              proj4string = proj)
            names(dat) <- poresp
            attr(dat,'family') <- 'poisson'
            attr(dat,'data_type') <- 'Present only'
            dat
            
          }
        else  { # if the PO dataset has some extra columns apart from X and Y
          
          names <- names(dat)[!names(dat)%in%c(coords)]
          dat <- sp::SpatialPointsDataFrame(coords = dat[,coords],
                                            data = as.data.frame(dat[,!names(dat)%in%coords]),
                                            proj4string = proj)
          colnames(dat@data) <- names
          if (!poresp%in%colnames(dat@data)) {
            dat[,poresp] <- 1
          }
          ### ??? hopefully giving it attributes that we'll use later but... who knows
          attr(dat,'family') <- 'poisson'
          attr(dat,'data_type') <- 'Present only'
          dat
          
        }
        
      }
    }
})

names(data_attributes) <- data_names


if (inclcoords) { # we have selected FALSE in this so all this is not happening
  ##Should I include?
  for (i in 1:length(data_attributes)) {

    data_attributes[[i]]@data[,coords] <- data_attributes[[i]]@coords

  }

}


if (marks) {
  
  data_marks = list()
  #Make a unique SpatialPointsDataframe for each mark (to be run on spatial covariates).
  #Incorporate standardized list of marks/covariates for data.
  #Include non numeric marks as well
  for (i in 1:length(data_attributes)) {
    
    ind <- 1 + length(data_marks)
    
    if (class(data_attributes[[i]]) == 'SpatialPoints') data_marks[[ind]] <- FALSE
    else {
      
      names = names(data_attributes[[i]])[!names(data_attributes[[i]])%in%c(poresp,paresp,coords,trialname)]
      #Variable for class of variable:
      #if numeric run family as user specified
      #Else if character or factor run as multinomial
      
      class_marks <- sapply(data_attributes[[i]]@data[names], class)
      
      if (!is.null(inclmarks)) names <- names[names%in%inclmarks]
      
      if (is.null(names) | identical(names,character(0))) data_marks[[ind]] <- FALSE
      
      else
        for(j in 1:length(names)) {
          
          index <- ind + j - 1
          
          if (class_marks[j] == 'character'| class_marks[j] == 'factor') {
            
            if (attributes(data_attributes[[i]])$family == 'poisson')  mark_response <- data_attributes[[i]]@data[,poresp]
            
            else mark_response <- data_attributes[[i]]@data[,paresp]
            
            mark <- sp::SpatialPointsDataFrame(coords = coordinates(data_attributes[[i]]),
                                               data = data.frame(factor((data_attributes[[i]]@data[,names[j]]))),
                                               proj4string = proj)
            colnames(mark@data) <- names[j]
            mark@data[,paste0(names[j],'_phi')] <- rep(1,nrow(mark@coords))
            mark@data[,paste0(names[j],'_response')] <- mark_response ##How do we run the response for marks below??
            
            mark@data[,'mark_response_weights'] <- mark_response
            n_species <- sum( mark@data[,'mark_response_weights'])
            #FOR count data do this:
            weights = as(mark@data,'data.table')
            weights = weights[, .(weight = n_species/sum(mark_response_weights)), by = eval(names[j])]
            mark@data[,'weights'] <- weights[as.numeric(mark@data[,names[j]])]$weight
            #weights <- nrow(mark@coords)/(as.numeric(table(mark@data[,names[j]])))
            #Weights wont work for count data
            #mark@data[,'weights'] <- weights[as.numeric(mark@data[,names[j]])]
            attr(mark,'family') <- 'poisson'
            attr(mark,'data_type') <- 'Multinomial mark'
            ##Add phi and factor_variable names as attributes
            attr(mark,'mark_name') <- names[j]
            attr(mark, 'phi') <- paste0(names[j],'_phi')
            attr(mark,'weights') <- TRUE
            ##Then when adding them to component joint say unique(phi) etc... do avoid duplicates
            data_marks[[index]] <- mark
            names(data_marks)[[index]] <- paste0(names(data_attributes)[i],'_',names[j])
          }
          else
            if (class_marks[j] == 'numeric' | class_marks[j] == 'integer')
            {
              
              mark <- sp::SpatialPointsDataFrame(coords = coordinates(data_attributes[[i]]),
                                                 data = as.data.frame(data_attributes[[i]]@data[,names[j]]),
                                                 proj4string = proj)
              colnames(mark@data) <- names[j] #paste0(names(data_attributes)[i],'_',names[j]) #Should I do this? Would we not want group effects for the marks?#But then Names marks is not the same?
              attr(mark,'family') <- markfamily
              capital_markfamily <- gsub("^(\\w)(\\w+)", "\\U\\1\\L\\2", 
                                         markfamily, perl = TRUE)
              attr(mark,'data_type') <- paste0(capital_markfamily,' mark')
              attr(mark,'mark_name') <- names[j]
              attr(mark,'phi') <- NA
              attr(mark, 'weights') <- FALSE
              data_marks[[index]] <- mark
              names(data_marks)[[index]] <- paste0(names(data_attributes)[i],'_',names[j])
              
            }
          #else FALSE
        }
    }
  }
  
  data_marks[sapply(data_marks,is.logical)] <- NULL
  
  if (length(data_marks) == 0) stop("Either marks have been set to TRUE and no datasets contain marks, or marks to include only contains marks not present in any dataset.")
  
  names_marks <- sapply(data_marks, function(mark) attributes(mark)$mark_name)
  
  multinom_incl <- sapply(data_marks, function(mark) attributes(mark)$data_type == 'Multinomial mark')
  
  if (any(multinom_incl)) {
    
    multinom_vars <- unique(unlist(sapply(data_marks, function(mark) {
      
      if(attributes(mark)$data_type == 'Multinomial mark') attributes(mark)$mark_name
      
    })))
    
    data_attributes <- lapply(data_attributes, function(dat){
      
      if (any(multinom_vars%in%names(dat))) {
        
        dat@data[,multinom_vars] <- NULL
        dat
        
      }
      else dat
      
    })
    
  }
}


##---------------------------------##
#### mesh and integration points ####
##---------------------------------##

if (is.null(mesh)) { # if mesh is not provided, but we have (lines 492 - 532)
  
  warning("Mesh not provided. Will try to create own mesh.")
  
  #Make mesh same way as PointedSDMs
  if (is.null(bdry)) {
    
    if (inherits(spatialcovariates, "Spatial")) spatcoords <- sp::SpatialPoints(coords = spatialcovariates@coords,
                                                                                proj4string = proj)
    else spatcoords <- sp::SpatialPoints(coords = spatialcovariates[,coords],
                                         proj4string  = proj)
    
    bstart <- min(c(diff(sort(unique(spatcoords@coords[,1]))), diff(sort(unique(spatcoords@coords[,2])))))
    
    poly.tmp <- rgeos::gBuffer(spatcoords, width=bstart, byid=TRUE)
    
    bdry <- rgeos::gBuffer(rgeos::gUnaryUnion(poly.tmp), width=bstart)
    
  }
  
  else {
    if (class(bdry)!="SpatialPolygons") {
      
      bdry <- sp::SpatialPolygons(Srl=list(Polygons(srl=list(bdry), ID="eek")))
      
    } 
    else {
      if (!is.projected(bdry)) bdry <- spTransform(bdry, CRSobj = proj)
    }
    
  }
  
  region.bdry <- inla.sp2segment(bdry)
  
  #Is there a nice way to add other parameters to inla.mesh.2d?
  mesh <- inla.mesh.2d(boundary=region.bdry, 
                       cutoff=meshpars$cutoff,
                       max.edge=meshpars$max.edge, 
                       offset=meshpars$offset)
  
}

if (is.null(ips)) { # if integration points were not provided, the function will calculate them
  
  warning('Integration points not provided. Will try to create own points')
  
  ips <- ipoints(samplers = bdry,
                 domain = mesh)
  
}

#Should I do the same for raster data??
#Easiest way to fix is by spatdata <- as(rasterdata, 'SpatialPixelsDataFrame')
#Is there loss in data this way??
#Does this do the same thing as 'GetNearestCovariate'?
#When inlabru update comes, change SpatialPointsDataFrame part
#SpatialGridDataFrame?
#Remove the if ncol == 1,

##------------------------##
#### spatial covariates ####
##------------------------##

if (inherits(spatialcovariates,'Spatial')) { # if spatialcovariates is a Spatial* object
  if (ncol(spatialcovariates) == 1) { # if we only have one spatial covariate (lines 559 - 578)
    if(class(spatialcovariates) == 'SpatialPixelsDataFrame') { # if spatialcovariates is a SpatialPixelDataframe
      
      proj4string(spatialcovariates) <- proj
      spatnames <- names(spatialcovariates@data)
      assign(spatnames, spatialcovariates)}
    
    else { # if it's a SpatialPointsDataFrame (although according to the function it could be any other Spatial* object, but it would fail with anything other than SpatialPointsDataFrame)
      
      spatnames <- names(spatialcovariates)
      spatcoords <- spatialcovariates@coords
      spatdata <- spatialcovariates@data[,!colnames(spatialcovariates@data)%in%coords]
      spatpix <- sp::SpatialPixelsDataFrame(points = spatcoords,
                                            data = data.frame(spatdata), 
                                            tolerance = tol,
                                            proj4string = proj)
      names(spatpix@data) <- spatnames
      assign(names(spatpix@data),spatpix)
    }
  }
  else { # if there's more than 1 spatial covariate (not making difference between 1 spatial covariate or >1!)
    
    # get coords, data, and names of vairables separately
    spatcoords <- spatialcovariates@coords
    spatdata <-  spatialcovariates@data[,!colnames(spatialcovariates@data)%in%coords]
    spatnames <- names(spatdata)
    
    
    # turn a spatialpixeldataframe of N covariates into N spatialpixelsdataframe of 1 covariate each. 
    for (i in 1:ncol(spatdata)) {
      # i = 1
      spatpix <- sp::SpatialPixelsDataFrame(points = spatcoords,
                                            grid = spatialcovariates@grid,
                                            data = data.frame(spatdata[,i]), 
                                            tolerance = tol,
                                            proj4string = proj)
      colnames(spatpix@data) = colnames(spatdata)[i]
      assign(names(spatpix@data),spatpix)
      
    }
  }
}  else # if the spatialcovariates is NOT a Spatial* object (ours is, so we ignore this!)
  if (class(spatialcovariates) == 'data.frame') {
    
    warning("Spatialcovariates is of class 'data.frame'.\nWill convert it to a SpatialPixelsDataFrame.")
    spatcoords <- spatialcovariates[,coords]
    spatdata <- as.data.frame(spatialcovariates[,!colnames(spatialcovariates)%in%coords])
    spatnames <- names(spatialcovariates)[!colnames(spatialcovariates)%in%coords]
    
    for (i in 1:ncol(spatdata)) {
      spatpix <- sp::SpatialPixelsDataFrame(points = spatcoords,
                                            data = data.frame(spatdata[,i]), 
                                            tolerance = tol,
                                            proj4string = proj)
      colnames(spatpix@data) = spatnames[i]
      
      assign(names(spatpix@data),spatpix)
      
    }
    
  }

##-------------------------##
#### inla specifications ####
##-------------------------##

## SPDE ##
spde2 <- inla.spde2.matern(mesh)

##Construct joint components for the likelihoods.
##Will need to change with inclusion of separate covariates.


## Formulas ##
# If we haven't specified formulas, it'll construct the same for both models, apparently, and including all covariates in spatialcovariates [Lines 637 - 653]

if (is.null(poformula) | is.null(paformula)) { 
  
  components_joint <- formula(paste(c('~ 0',paste0(spatnames,'(main = ',spatnames,', model = "linear")')), collapse = '+'))
  
  if (inclcoords) { # if coords should be added as fixed effects (in our case "no" to start with)
    
    components_joint <- update(components_joint, paste0('~ . +',coords, collapse = '+'))
  }
  
  if (intercept) { # if a joint intercept should be added ("no" in our case)
    
    components_joint <- update(components_joint, ~ . + Intercept(1))
    
  }
  
  
} 

## Model parameters ##

likelihoods = list()

family <- unlist(sapply(data_attributes, function(x) attributes(x)$family))

trials <- sapply(data_attributes, function(x){ # if trials is null it will just be a named vector of 1s (one for each dataset)
  
  if (!is.null(attributes(x)$Ntrials)) attributes(x)$Ntrials
  else 1
  
}) 

# family uis a named character with the distribution families for PA_data and PO_data

E_param <- sapply(family, function(x) {
  if (x == 'poisson') 0
  else
    if (x == 'binomial') 1
  
})

#E_param is a named numeric vector, with a "0" for the PA_data and a 1 for the PO_data

## Create separated formulas for each likelihood ##

# This will take the common formula above and generate a list with the formula for PO data (containing the intercept of the PA likelihood as a variable) and viceversa
##Take out any brackets from 'components_joint'.
##I.e (for now) run coordinates only on spatial covariates (and optional others).

form_elements <- gsub(" *\\(.*?\\) *", "",components_joint)

formula <- mapply(function(fam,ind) {
  if (!is.null(poformula) & fam == 'poisson') {
    
    formula <- poformula
  }
  else
    if (is.null(poformula) & fam == 'poisson') {
      
      formula <- formula(paste0(c(poresp,'~', form_elements[2]),collapse = " ")) 
      
    }
  
  else 
    if(!is.null(paformula) & fam == 'binomial') {
      
      formula <- paformula
      
    }
  else
    if(is.null(paformula) & fam == 'binomial'){
      
      formula <- formula(paste0(c(paresp,"~", form_elements[2]),collapse = " "))
      
    }
  
  if (indivintercepts) {
    
    formula <- update(formula,paste0(' ~ . +', paste0(data_names[ind],'_intercept'), collapse = ' + '))
    
  }
  else formula
  
  if (pointsspatial) {
    
    formula <- update(formula, paste0('~ . +',data_names[[ind]],'_spde'))
    
  }
  else formula
  
}, fam = family, ind = 1:length(family))

## formula is a list of the two formulas, one for PA_data and one for PO_data

## this generates a list of lists. each sublist contains all the necessary objects for one likelihood (PA and PO)
for (i in 1:1) {
  
  lhoods <- inlabru::like(formula = formula[[i]], ##Add tag to this likelihood somehow?
                          family = family[i],
                          data = data_attributes[[i]],
                          mesh = mesh,
                          ips = ips,
                          Ntrials = trials[i],
                          E_param[i])
  likelihoods <- like_list(lhoods)
  
  
  if (length(family) > 1) { #Better way of doing this??
    for (j in 2:length(family)) {
      
      lhoods <- inlabru::like(formula = formula[[j]],
                              family = family[j],
                              data = data_attributes[[j]],
                              mesh = mesh,
                              ips = ips,
                              Ntrials = trials[j],
                              E_param[j])
      likelihoods[[j]] <- lhoods
      
    }
  }
  
  likelihoods
} ## likelihoods is a bru_like list with everything needed for tthe PA_data and the PO_data models 


if (marks) {
  
  family_marks <- sapply(data_marks, function(x) attributes(x)$family)
  formula_marks <- list()
  likelihoods_marks <- list()
  
  mark_weights <- lapply(data_marks, function(x){
    
    if (attributes(x)$weights) x@data[,'weights']
    else 1
    
    
  })
  
  for (i in 1:length(family_marks)) {
    
    formula_marks[[i]] <- formula(paste0(c(names_marks[i],'~',form_elements[2]),collapse = " "))
    
    if (marksspatial) {
      
      formula_marks[[i]] <- update(formula_marks[[i]], paste0(" . ~ . +", names_marks[i],'_spde'))
      
    }
    
    if (indivintercepts) { #probably fix something here? No indiv intercepts for multinomial response, but indiv intercepts for marks
      
      if (attributes(data_marks[[i]])$data_type != 'Multinomial mark'){
        
        formula_marks[[i]] <- update(formula_marks[[i]],paste0('. ~ . +', paste0(names_marks[i],'_intercept'), collapse = ' + '))
        
      }
    }
    
    if (attributes(data_marks[[i]])$data_type == 'Multinomial mark') {
      
      formula_marks[[i]] <- update(formula_marks[[i]], paste0(paste0(attributes(data_marks[[i]])$mark_name,'_response'), ' ~ . + ', paste(attributes(data_marks[[i]])$mark_name, attributes(data_marks[[i]])$phi, sep = ' + ')))
      
    }
    
  }
  
  for (k in 1:length(family_marks)) {
    ##Need to add exposure parameter here
    ## So probably need to add a new sapply if weights in data attributes
    ## otherwise E = 0
    lhoods <- inlabru::like(formula = formula_marks[[k]],
                            family = family_marks[k],
                            data = data_marks[[k]],
                            mesh = mesh,
                            ips = ips,
                            E = mark_weights[[k]])
    likelihoods_marks[[k]] <- lhoods
    
    
  }
  n <- length(likelihoods)
  for (l in 1:length(likelihoods_marks)) {
    
    #Better way to do this?
    likelihoods[[l + n]] <- likelihoods_marks[[l]]
    
  }
  
}


ips$int_resp <- 0
#ips <- spTransform(ips, proj) <- doesn't work if ips is not projected
proj4string(ips) <- proj # <- is this fine? # it's the one it already has

## Likelihood for Integration Points????? ## 

#Run integration points only on spatialcovariates?
like_ip = like(formula = formula(paste0(c('int_resp ~ 0', c(spatnames)) ,collapse = '+')), #formula(paste0(c("resp ~", form_elements[2]),collapse = " ")),
               family = 'poisson',
               mesh = mesh,
               E = ips$weight,
               data = ips)

likelihoods[[length(likelihoods) + 1]] = like_ip

# names(likelihoods) <- c(data_names,names_marks, species_names, 'like_ip')

if (indivintercepts) {
  
  components_joint <- update(components_joint, paste0(' ~ . +', paste0(c(data_names),'_intercept(1)'), collapse = ' + '))
  
  for (i in 1:length(data_marks)) {
    if (marks) {
      if (attributes(data_marks[[i]])$data_type != "Multinomial mark") {
        
        components_joint <- update(components_joint, paste0(' ~ . +', paste0(c(names_marks[i]),'_intercept(1)'), collapse = ' + '))
        
      }
    }
    
  }
  
}

if (pointsspatial) {
  
  components_joint <- update(components_joint, paste('. ~ . +',paste0(data_names,'_spde(main = coordinates, model = spde2)',collapse = ' + ')))
  
}

if (marksspatial) {
  # if (!is.null(names_marks)) { #I.e. if marks is null but species is non null. Should I add a seperate random for species?
  components_joint <- update(components_joint, paste('. ~ . +',paste0(names_marks,'_spde(main = coordinates, model = spde2)',collapse = ' + ')))
  #}
}

if (marks) {
  if (any(multinom_incl)) {
    
    factor_vars <- sapply(data_marks, function(name) attributes(name)$mark_name)
    factor_vars <- unique(factor_vars[multinom_incl])
    components_joint <- update(components_joint, paste(' . ~ . + ', paste0(factor_vars,'(main = ', factor_vars, ', model = "iid",constr = FALSE, fixed=TRUE)', collapse = ' + ')))
    
    phi_vars <- sapply(data_marks, function(name) attributes(name)$phi)
    phi_vars <- unique(phi_vars[multinom_incl])
    components_joint <- update(components_joint, paste(' . ~ . +', paste0(phi_vars, '(main = ',phi_vars, ', model = "iid", initial = -10, fixed = TRUE)', collapse = ' + ')))
    
  }
}


## assign cloglog link to the PA data
for (i in 1:length(likelihoods)) {
  
  if (likelihoods[[i]]$response == paresp) options[['control.family']][[i]] <- list(link = 'cloglog')
  
  else options[['control.family']][[i]] <- list(link = 'default')
  
}

model_joint <- bru(components = components_joint,
                   likelihoods, options = options)

if (!is.null(residuals)) {
  
  name_resp <- c()
  
  for (i in 1:(length(likelihoods) - 1)) {
    
    name_resp[i] <- gsub("\\(|\\)","",likelihoods[[i]]$formula[2])
    
  }
  
  fitted_residuals = list()
  for (i in 1:(length(likelihoods) - 1)) {
    res <- c()
    for (j in (model_joint$bru_iinla$inla_stack$data$index[[i]][1]):tail(model_joint$bru_iinla$inla_stack$data$index[[i]], n = 1)) {
      
      res[j] <- INLA::inla.emarginal(function(x) x, model_joint$marginals.fitted.values[[j]])
      
      
    }
    
    fitted_residuals[[i]] <- as.vector(na.omit(res))
    
  }
  
  calc_residuals = list()
  if (residuals == 'response') {
    for (k in 1:length(fitted_residuals)) {
      
      calc_residuals[[k]] <- likelihoods[[k]]$data@data[,name_resp[k]] - fitted_residuals[[k]]
      
    }
  }
  else
    if (residuals == 'pearson') {
      for (k in 1:length(fitted_residuals)) {
        stop('FIX THIS')
        calc_residuals[[k]] <- (likelihoods[[k]]$data@data[,name_resp[k]] - fitted_residuals[[k]])/sqrt(fitted_residuals[[k]])
        
      }
      
    } 
  else
    if (residuals == 'deviance') {
      for (k in 1:length(fitted_residuals)) {
        stop('FIX THIS')
        calc_residuals[[k]] <- sign(likelihoods[[k]]$data@data[,name_resp[k]] - fitted_residuals[[k]]) * sqrt(2 * likelihoods[[k]]$data@data[,name_resp[k]] * log(likelihoods[[k]]$data@data[,name_resp[k]]/fitted_residuals[[k]]) - (likelihoods[[k]]$data@data[,name_resp[k]] - fitted_residuals[[k]]))
        
      }
    }
  
  names(calc_residuals) <- c(data_names,names_marks)
  
  model_joint[['model_residuals']] = calc_residuals
  
  
}

data_type <- sapply(c(data_attributes,data_marks), function(x) attributes(x)[['data_type']])
names(data_type) <- c(data_names,names_marks)
model_joint[['data_type']] <- data_type

if (!is.null(multinom_vars)) { 
  
  model_joint[['multinom_vars']] <- multinom_vars
  
}

class(model_joint) <- c('bru_sdm',class(model_joint))
return(model_joint)