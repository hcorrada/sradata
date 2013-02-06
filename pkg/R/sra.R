setClassUnion("listOrNULL", c("list", "NULL"))

setClass("SRATable", representation(filePath = "character",
                                    ePtr = "externalptr",
                                    info = "listOrNULL"))

SRATable <- function(filePath) {
  new("SRATable", filePath)
}

setMethod("initialize", "SRATable", function(.Object, filePath, ...) {
  if (missing(filePath))
    return(.Object)

  .Object@ePtr <- .Call("sra_open", filePath, package = "SRAdata");
  .Object@filePath <- filePath
  .Object@info <- NULL
  return(.Object)
})

setGeneric("SRAGetTableInfo", function(sraTable, ...) {
  standardGeneric("SRAGetTableInfo")
})

setMethod("SRAGetTableInfo", c("SRATable"), function(sraTable, ...) {
  callres <- .Call("sra_get_info", sraTable@ePtr, PACKAGE="SRAdata")
  names(callres) <- c("maxSpotId", "nlanes", "lanes", "nTiles", "starts", "rlens")

  res <- list()
  res$maxSpotId <- callres$maxSpotId
  res$nlanes <- callres$nlanes
  res$lanes <- callres$lanes[1:res$nlanes]
  res$ntiles <- callres$nTiles[1:res$nlanes]
  res$rlens <- callres$rlens[1:res$nlanes]

  res$spotids <- vector("list", res$nlanes)
  cs <- c(1,cumsum(res$ntiles))
  for (i in 1:res$nlanes) {
    ii <- cs[i]+(0:(cs[i+1]-1))
    res$spotids[[i]] <- callres$starts[ii]
  }
  return(res)
})

setMethod("show", "SRATable", function(object) {
  cat("SRATable Object:\n")
  cat("  filePath: ", object@filePath, "\n")
  
  if (is.null(object@info)) {
    cat("  info: none\n");
  } else {
    cat("  info: maxSpotId = ", object@info$maxSpotId, "\n")

    nreads <- c(sapply(object@info$spotids, "[", 1), object@info$maxSpotId)
    nreads <- nreads[-1] - nreads[-length(nreads)] + 1

    for (i in seq_along(object@info$lanes)) {
      cat(sprintf("        lane %d: nTiles = %d, nReads = %d, readLength=%d\n",
                  object@info$lanes[i], object@info$ntiles[i], as.integer(nreads)[i], object@info$rlens[i]))
    }
  }
})

setGeneric("SRATableReadFastq", function(sraTable, ...) {
  standardGeneric("SRATableReadFastq")
})

setMethod("SRATableReadFastq", "SRATable", function(sraTable, lane=NULL, tile=NULL, max.nreads=NULL, ...) {
  require(ShortRead)
  
  if (is.null(sraTable@info)) {
    stop("Error in SRATableReadFastq: info slot on SRATable object is NULL. Call SRATableGetInfo() first\n");
  }

  if (is.null(lane)) {
    lane <- sraTable@info$lanes[1]
  }
  
  if (length(lane) > 1) {
    stop("Will only do one lane at a time")
  }
    
  m <- match(lane, sraTable@info$lanes)
  if (is.null(m)) {
    stop("Requested lane ", lane, " not found in SRA table")
  }

  if (!is.null(tile)) {
    # tile is given, make sure only one requested
    if (length(tile) > 1) {
      stop("Will only do one tile at a time")
    }
    if (tile < 1 || tile > sraTable@info$ntiles[m]) {
      stop("Requested tile ", tile, " is out of range")
    }

    startSpotId = sraTable@info$spotids[[m]][tile]
    # check if this is the last tile in the lane
    if (tile == sraTable@info$ntiles[m]) {
      # this is the last tile in the lane
      # check if this is the last lane
      if (m == sraTable@info$nlanes) {
        # this is the last lane: stop at the last spot
        endSpotId <- sraTable@info$maxSpotId
      } else {
        # this is not the last lane: stop before start of next lane
        endSpotId <- sraTable@info$spotids[[m+1]][1] - 1
      }
    } else {
      # this is not the last tile in the lane: stop before start of next tile
      endSpotId <- sraTable@info$spotids[[m]][tile+1] - 1
    }
  } else {
    # no tile given, return full lane
    startSpotId <- sraTable@info$spotids[[m]][1]

    # check if this is the last lane
    if (m == sraTable@info$nlanes) {
      # this is the last lane: stop at the last spot
      endSpotId <- sraTable@info$maxSpotId
    } else {
      # this is not the last lane: stop before start of the next lane
      endSpotId <- sraTable@info$spotids[[m+1]][1]-1
    }
  }

  cur.nreads = endSpotId - startSpotId + 1
  if (!is.null(max.nreads) && (cur.nreads > max.nreads)) {
    endSpotId = startSpotId + max.nreads - 1
  }

  rlen <- sraTable@info$rlens[m]

  sread <- DNAStringSet(.Call("sra_getBases", sraTable@ePtr, as.integer(startSpotId), as.integer(endSpotId), as.integer(rlen), PACKAGE="SRAdata"))
  id <- BStringSet(.Call("sra_getIds", sraTable@ePtr, as.integer(startSpotId), as.integer(endSpotId), PACKAGE="SRAdata"))
  quals <- .Call("sra_getQuals", sraTable@ePtr, as.integer(startSpotId), as.integer(endSpotId), as.integer(rlen), PACKAGE="SRAdata")
  quals <- FastqQuality(BStringSet(quals, start=1, width=rlen))

  res <- new("ShortReadQ", sread=sread, quality=quals, id=id)
  return(res)
})

setGeneric("SRATableReadIntensities", function(sraTable, ...) {
  standardGeneric("SRATableReadIntensities")
})

setMethod("SRATableReadIntensities", "SRATable", function(sraTable, lane=NULL, tile=NULL, max.nreads=NULL, ...) {
  require(ShortRead)
  
  if (is.null(sraTable@info)) {
    stop("Error in SRATableReadFastq: info slot on SRATable object is NULL. Call SRATableGetInfo() first\n");
  }

  if (is.null(lane)) {
    lane <- sraTable@info$lanes[1]
  }
  
  if (length(lane) > 1) {
    stop("Will only do one lane at a time")
  }
    
  m <- match(lane, sraTable@info$lanes)
  if (is.na(m)) {
    stop("Requested lane ", lane, " not found in SRA table")
  }

  if (!is.null(tile)) {
    # tile is given, make sure only one requested
    if (length(tile) > 1) {
      stop("Will only do one tile at a time")
    }
    if (tile < 1 || tile > sraTable@info$ntiles[m]) {
      stop("Requested tile ", tile, " is out of range")
    }

    startSpotId = sraTable@info$spotids[[m]][tile]
    # check if this is the last tile in the lane
    if (tile == sraTable@info$ntiles[m]) {
      # this is the last tile in the lane
      # check if this is the last lane
      if (m == sraTable@info$nlanes) {
        # this is the last lane: stop at the last spot
        endSpotId <- sraTable@info$maxSpotId
      } else {
        # this is not the last lane: stop before start of next lane
        endSpotId <- sraTable@info$spotids[[m+1]][1] - 1
      }
    } else {
      # this is not the last tile in the lane: stop before start of next tile
      endSpotId <- sraTable@info$spotids[[m]][tile+1] - 1
    }
  } else {
    # no tile given, return full lane
    startSpotId <- sraTable@info$spotids[[m]][1]

    # check if this is the last lane
    if (m == sraTable@info$nlanes) {
      # this is the last lane: stop at the last spot
      endSpotId <- sraTable@info$maxSpotId
    } else {
      # this is not the last lane: stop before start of the next lane
      endSpotId <- sraTable@info$spotids[[m+1]][1]-1
    }
  }

  cur.nreads = endSpotId - startSpotId + 1
  if (!is.null(max.nreads) && (cur.nreads > max.nreads)) {
    endSpotId = startSpotId + max.nreads - 1
  }

  rlen <- sraTable@info$rlens[m]
  intInfo <- .Call("sra_getReadInfo", sraTable@ePtr, as.integer(startSpotId), as.integer(endSpotId), PACKAGE="SRAdata");
  intInfo <- SolexaIntensityInfo(intInfo[[1]], intInfo[[2]], intInfo[[3]], intInfo[[4]]);

  int <- .Call("sra_getIntensity", sraTable@ePtr, as.integer(startSpotId), as.integer(endSpotId), as.integer(rlen), PACKAGE="SRAdata");
  nreads <- endSpotId - startSpotId + 1

  int <- array(int, dim=c(4L, rlen, nreads), dimnames=list(c("A","C","G","T"), NULL, NULL));
  int <- aperm(int, c(3,1,2));

  res <- SolexaIntensity(int, readInfo=intInfo)
  return(res)
})
          
