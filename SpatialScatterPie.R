# Create an isolated environment that have the helper functions to use scatterpie
#' @import ggplot2
CreateScatterPieEnv <- function() {
  require(scatterpie)
  helper.env <- new.env()
  with(helper.env, {
    # Take an array representing an image and plot it with ggplot2
    #' @import ggplot2
    .plot_image <- function(x) {
        # Check necessary packages are installed and if not STOP
        #.test_installed(c("grid", "ggplot2"))
        x <- grid::rasterGrob(x,
            interpolate = FALSE,
            width = grid::unit(1, "npc"),
            height = grid::unit(1, "npc"))
        ggplot() +
            annotation_custom(
                grob = x,
                xmin = 0,
                xmax = ncol(x$raster),
                ymin = 0,
                ymax = nrow(x$raster)) +
            coord_fixed(
                xlim = c(0, ncol(x$raster)),
                ylim = c(0, nrow(x$raster))) +
            theme_void()
    }

    #Extract image and convert it to array from allowed classes
    .extract_image <- function(x, slice = NULL) {
        # Iterate over all the accepted classes and convert the image to array
        if (is.character(x)) {
            # Check if the file exists
            stopifnot(file.exists(x))
            # Check the file is in the right format
            typ <- c("jpg", "jpeg", "png")
            pat <- paste0(".", typ, "$")
            idx <- vapply(pat, grepl, x = x, logical(1), ignore.case = TRUE)
            if (!any(idx)) {
                stop("'x' should be of file type JPG, JPEG or PNG")
            }
            # Read file
            x <- switch(typ[idx],
                png = png::readPNG(x),
                jpeg::readJPEG(x))
        } else if (is(x, "Seurat")) {
            # Stop if there are no images or the name selected doesn't exist
            stopifnot(
                !is.null(SeuratObject::Images(x)),
                slice %in% SeuratObject::Images(x))
            # If image is null use the first slice
            if (is.null(slice)) 
                slice <- SeuratObject::Images(x)[1]
            # Extract Image in raster format
            x <- SeuratObject::GetImage(x, image = slice, mode = "raster")
            # Conver to matrix
            x <- as.matrix(x)
        } else {
            stop("Couldn't extract image, See ?plotImage for valid image inputs.")
        }
        return(x)
    }

    plotImage <- function(x, slice = NULL) {
        # check validity of input arguments
        stopifnot(
            # Check for valid x classes
            is.matrix(x) | is.character(x) | is.array(x) | is(x, "rastergrob") | 
                is(x, "Seurat") | is(x, "SpatialExperiment"),
            # Check for valid slice classes
            is.null(slice) | is.character(slice))
        
        if (!is.array(x))
            x <- .extract_image(x)
            
        # Plot image
        plt <- .plot_image(x)
    }

    plotSpatialScatterpie <- function(
        x,
        y,
        cell_types = colnames(y),
        img = FALSE,
        slice = NULL,
        scatterpie_alpha = 1,
        pie_scale = 0.4,
        degrees = NULL,
        axis = NULL,
        ...) {
        # Check necessary packages are installed and if not STOP
        
        # Class checks
        stopifnot(
            # Check x inputs
            is.matrix(x) | is.data.frame(x) |
                is(x, "Seurat") | is(x, "SpatialExperiment"),
            # Check y inputs
            is.matrix(y) | is.data.frame(y),
            # cell_types needs to be a character with max length = ncol(y)
            is.character(cell_types) & length(cell_types) <= ncol(y),
            # Check img
            # img not checked since its checked in plotImage()
            # Check slice name
            is.character(slice) | is.null(slice),
            # Check plotting parameters are numeric
            is.numeric(scatterpie_alpha),
            is.numeric(pie_scale),
            is.numeric(degrees) | is.null(degrees),
            axis %in% c("h", "v") | is.null(axis)
        )
        
        # If image is passed add it as the base layer, if not, no image
        # Need to use isFALSE bc img can have many different inputs
        # Set ymax to overlap image and piecharts
        if (isFALSE(img)) {
            p <- ggplot() +
                coord_fixed()
            ymax <- 0
        } else {
            # Extract image from Seurat or SE objects when img is TRUE
            # If image is not TRUE and not FALSE an acceptable class for plotImage
            # has been passed
            if (is(x, "Seurat") | is(x, "SpatialExperiment") & isTRUE(img)) {
                img <- .extract_image(x, slice)
            }
            
            p <- plotImage(x = img)
            ymax <- max(p$coordinates$limits$y)
        }
        
        # Extract coordinate matrix from x
        if (!is.matrix(x))
            x <- .extract_coord(x = x, slice = slice, img = img)
        
        # Check colnames
        x <- .x_cnames(x)
        
        # Convert y to matrix format
        if (!is.matrix(x)) {
            y <- as.matrix(x)
        }
        
        # Stop if x and y don't have the same number of columns or if the
        # rownames are not common between them
        stopifnot(
            nrow(x) == nrow(y),
            all(rownames(x) %in% rownames(y)))
        
        # merge by row names (by=0 or by="row.names")
        df <- merge(x, y, by = 0, all = TRUE)
        # make y negative
        df$coord_y_i <- abs(df$coord_y - ymax)
        
        # Plot
        p + scatterpie::geom_scatterpie(
            data = df,
            aes_string(
                x = "coord_x",
                y = "coord_y_i"
            ),
            cols = cell_types,
            color = NA,
            alpha = scatterpie_alpha,
            pie_scale = pie_scale,
            ...) +
            # Below not needed bc comes from plotImage
            # coord_fixed() +
            theme_void() + 
            theme(legend.key.size = unit(0.5, "lines"))
        }

    .x_cnames <- function(x) {
        # If the column names of x aren't right fix them
        cnames <- c("coord_y", "coord_x")
        if (!all(colnames(x) %in% cnames)) {
            colnames(x) <- cnames
        }
        x
    }

    #Coordinates and return a matrix object where each row is a spot and the
    #columns are the x and y coordinates
    .extract_coord <- function(x, slice, img) {
        # Iterate over all the accepted classes and return spot coordinates
        if (is.data.frame(x)) {
            # Convert to matrix
            x <- as.matrix(x)
        } else if (is(x, "Seurat")) {
            # Stop if there are no images or the name selected doesn't exist
            stopifnot(
                # Stop if there are no images
                !is.null(SeuratObject::Images(x)),
                # Stop if the image doesn't exist
                is.null(slice) | slice %in% SeuratObject::Images(x))
            
            # If image is null use the first slice
            if (is.null(slice))
                slice <- SeuratObject::Images(x)[1]
            
            # Extract spatial coordinates
            x <- as.matrix(SeuratObject::GetTissueCoordinates(x, image = slice))
        } else {
            stop("Couldn't extract image coordinates.
                Please check class(x) is SpatialExperiment, Seurat,
                dataframe or matrix")
        }
        return(x)
    }
  })
  return(helper.env)
}

# Take an array representing an image and plot it with ggplot2
#' @import ggplot2
.plot_image <- function(x) {
    # Check necessary packages are installed and if not STOP
    #.test_installed(c("grid", "ggplot2"))
    
    x <- grid::rasterGrob(x,
        interpolate = FALSE,
        width = grid::unit(1, "npc"),
        height = grid::unit(1, "npc"))
    
    ggplot() +
        annotation_custom(
            grob = x,
            xmin = 0,
            xmax = ncol(x$raster),
            ymin = 0,
            ymax = nrow(x$raster)) + 
        coord_fixed(
            xlim = c(0, ncol(x$raster)),
            ylim = c(0, nrow(x$raster))) + 
        theme_void()
        # theme_classic()
}

# Extract image and convert it to array from allowed classes
.extract_image <- function(x, slice = NULL) {
    # Iterate over all the accepted classes and convert the image to array
    if (is.character(x)) {
        #.test_installed(c("jpeg", "png"))
        
        # Check if the file exists
        stopifnot(file.exists(x))
        
        # Check the file is in the right format
        typ <- c("jpg", "jpeg", "png")
        pat <- paste0(".", typ, "$")
        idx <- vapply(pat, grepl, x = x, logical(1), ignore.case = TRUE)
        if (!any(idx)) {
            stop("'x' should be of file type JPG, JPEG or PNG")
        }
        
        # Read file
        x <- switch(typ[idx],
            png = png::readPNG(x),
            jpeg::readJPEG(x))
        
    } else if (is(x, "Seurat")) {
        #.test_installed(c("SeuratObject"))
        # Stop if there are no images or the name selected doesn't exist
        stopifnot(
            !is.null(SeuratObject::Images(x)),
            slice %in% SeuratObject::Images(x))
        
        # If image is null use the first slice
        if (is.null(slice)) 
            slice <- SeuratObject::Images(x)[1]
        
        # Extract Image in raster format
        x <- SeuratObject::GetImage(x, image = slice, mode = "raster")
        # Conver to matrix
        x <- as.matrix(x)
        
    } else {
        stop("Couldn't extract image, See ?plotImage for valid image inputs.")
    }
    return(x)
}

plotImage <- function(x, slice = NULL) {
    # check validity of input arguments
    stopifnot(
        # Check for valid x classes
        is.matrix(x) | is.character(x) | is.array(x) | is(x, "rastergrob") | 
            is(x, "Seurat") | is(x, "SpatialExperiment"),
        # Check for valid slice classes
        is.null(slice) | is.character(slice))
    
    if (!is.array(x))
        x <- .extract_image(x)
        
    # Plot image
    plt <- .plot_image(x)
}

plotSpatialScatterpie_bt <- function(
    x,
    y,
    cell_types = colnames(y),
    img = FALSE,
    slice = NULL,
    scatterpie_alpha = 1,
    pie_scale = 0.4,
    degrees = NULL,
    axis = NULL,
    ...) {
    # Check necessary packages are installed and if not STOP
    
    # Class checks
    stopifnot(
        # Check x inputs
        is.matrix(x) | is.data.frame(x) |
            is(x, "Seurat") | is(x, "SpatialExperiment"),
        # Check y inputs
        is.matrix(y) | is.data.frame(y),
        # cell_types needs to be a character with max length = ncol(y)
        is.character(cell_types) & length(cell_types) <= ncol(y),
        # Check img
        # img not checked since its checked in plotImage()
        # Check slice name
        is.character(slice) | is.null(slice),
        # Check plotting parameters are numeric
        is.numeric(scatterpie_alpha),
        is.numeric(pie_scale),
        is.numeric(degrees) | is.null(degrees),
        axis %in% c("h", "v") | is.null(axis)
    )
    
    # If image is passed add it as the base layer, if not, no image
    # Need to use isFALSE bc img can have many different inputs
    # Set ymax to overlap image and piecharts
    if (isFALSE(img)) {
        p <- ggplot() +
            coord_fixed()
        ymax <- 0
    } else {
        # Extract image from Seurat or SE objects when img is TRUE
        # If image is not TRUE and not FALSE an acceptable class for plotImage
        # has been passed
        if (is(x, "Seurat") | is(x, "SpatialExperiment") & isTRUE(img)) {
            img <- .extract_image(x, slice)
        }
        
        p <- plotImage(x = img)
        ymax <- max(p$coordinates$limits$y)
    }
    
    # Extract coordinate matrix from x
    if (!is.matrix(x))
        x <- .extract_coord(x = x, slice = slice, img = img)
    
    # Check colnames
    x <- .x_cnames(x)
    
    # Convert y to matrix format
    if (!is.matrix(x)) {
        y <- as.matrix(x)
    }
    
    # Stop if x and y don't have the same number of columns or if the
    # rownames are not common between them
    stopifnot(
        nrow(x) == nrow(y),
        all(rownames(x) %in% rownames(y)))
    
    # merge by row names (by=0 or by="row.names")
    df <- merge(x, y, by = 0, all = TRUE)
    # make y negative
    df$coord_y_i <- abs(df$coord_y - ymax)
    
    # Plot
    p + scatterpie::geom_scatterpie(
        data = df,
        aes_string(
            x = "coord_x",
            y = "coord_y_i"
        ),
        cols = cell_types,
        color = NA,
        alpha = scatterpie_alpha,
        pie_scale = pie_scale,
        ...) +
        # Below not needed bc comes from plotImage
        # coord_fixed() +
        theme_void() + 
        theme(legend.key.size = unit(0.5, "lines"))
    }

.x_cnames <- function(x) {
    # If the column names of x aren't right fix them
    cnames <- c("coord_y", "coord_x")
    if (!all(colnames(x) %in% cnames)) {
        colnames(x) <- cnames
    }
    x
}

# Coordinates and return a matrix object where each row is a spot and the
# columns are the x and y coordinates
.extract_coord <- function(x, slice, img) {
    # Iterate over all the accepted classes and return spot coordinates
    if (is.data.frame(x)) {
        # Convert to matrix
        x <- as.matrix(x)
    } else if (is(x, "Seurat")) {
        # Stop if there are no images or the name selected doesn't exist
        stopifnot(
            # Stop if there are no images
            !is.null(SeuratObject::Images(x)),
            # Stop if the image doesn't exist
            is.null(slice) | slice %in% SeuratObject::Images(x))
        
        # If image is null use the first slice
        if (is.null(slice))
            slice <- SeuratObject::Images(x)[1]
        
        # Extract spatial coordinates
        x <- as.matrix(SeuratObject::GetTissueCoordinates(x, image = slice))
    } else {
        stop("Couldn't extract image coordinates.
            Please check class(x) is SpatialExperiment, Seurat,
            dataframe or matrix")
    }
    return(x)
    
}
