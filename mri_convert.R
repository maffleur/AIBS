suppressMessages(require(oro.nifti))

### converts mni coordinates to voxels ###

getTransform <- function(niifile){
  nim <- readNIfTI(niifile)
  T <- rbind(nim@srow_x, nim@srow_y,nim@srow_z, c(0, 0, 0, 1))
  return(T)

}

mni2vox <- function(mni, T){

  invT <- solve(T)
  if (is.null(dim(mni))){
    mni <- matrix(mni, nrow=1)
  }
  mni <- cbind(mni, 1)
  res <- round(mni %*% t(invT))
  res[res < 1] <- 1
  res <- res[,1:3]
  return(res)


}

mni2vox_file <- function(mni, niifile){
  T <- getTransform(niifile)
  res <- mni2vox(mni, T)
  return(res)
}

vox2mni <- function(vox, T){
  vox <- matrix(c(vox,1), nrow=1)
  mni <- T %*% t(vox)
  return(mni[1:3])
}


#### kind of unrelated ###
voxel2content <- function(vox, nim){

  if (sum(vox <= 0))
    return(0)
  if (sum(vox > dim(nim)))
    return(0)
  return(nim[vox[1], vox[2], vox[3]])
}


checkNeighborhood <- function(vox, nsize, nim, reject=0){

  #as long as the return value equals reject we search the
  #simplistic neighborhood of the voxel.
  result=reject

  #create
  ttest <- vox
  for(i in 1:length(vox)){
    tmp <- vox[i] + nsize
    ttest <- rbind(ttest, tmp)
    tmp <- vox[i] - nsize
    ttest <- rbind(ttest, tmp)
  }
  mapping <- apply(ttest, 1, voxel2content, nim)
  idx <- min(which(mapping != reject), na.rm=T)

  if (length(idx) > 0)
    return(mapping[idx])

  return(0)

}

voxel2content_cube <- function(vox, nim){

  #generate cube (-1,0,+1) of vox
  voxels <- c()
  offset <- c(-1,0,1)
  for(i in offset){
    for(j in offset){
      for(k in offset){
        voxels <- rbind(voxels, vox + c(i,j,k))
      }
    }
  }

  res <- apply(voxels, 1, voxel2content, nim)
  res[res == 0]    <- NA
  res[is.nan(res)] <- NA
  mmm <- mean(res, na.rm=T)
  if (is.nan(mmm))
    mmm <- 0
  #sss <- sd(res, na.rm=T)
  return(mmm)
}


