hex2rgb <- function(img.data) {
  height <- dim(img.data)[1]
  width <- dim(img.data)[2]
  red.result <- list()
  green.result <- list()
  blue.result <- list()
  
  apply(img.data, 2, function(row) {
    rgb <- col2rgb(row)
    red.result <<- append(red.result, list(rgb[1, ]))
    green.result <<- append(green.result, list(rgb[2, ]))
    blue.result <<- append(blue.result, list(rgb[3, ]))
    return(TRUE)
  })
  
  red.result <- matrix(unlist(red.result), ncol = width, nrow = height)
  green.result <- matrix(unlist(green.result), ncol = width, nrow = height)
  blue.result <- matrix(unlist(blue.result), ncol = width, nrow = height)
  
  return(abind:::abind(red.result/255, green.result/255, blue.result/255, along=3))
}

rgb2hex <- function(img.data) {
  hex.result <- list()
  height <- dim(img.data)[1]
  width <- dim(img.data)[2]
  sapply(1:nrow(img.data), function(row.index) {
    red.channel <- img.data[, , 1]
    blue.channel <- img.data[, , 2]
    green.channel <- img.data[, , 3]
    
    merged.channels <- abind:::abind(red.channel[row.index, ],
                                     blue.channel[row.index, ],
                                     green.channel[row.index, ],
                                     along=2)
    merged.hex <- apply(merged.channels, 1, function(x) rgb(x[1], x[2], x[3], 1))
    hex.result <<- append(hex.result, list(merged.hex))
  })
  hex.result <- matrix(unlist(hex.result), ncol = width, nrow = height)
  return(hex.result)
}
