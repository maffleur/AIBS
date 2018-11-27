rm(list=ls())
setwd("~/Desktop/TIG_research_project/AIBS_studies/")
source("mri_convert.R")

aibs_raw <- "."
samples <- c("H0351.1009","H0351.1012","H0351.1015","H0351.1016","H0351.2001","H0351.2002")

atlas.path <- "./maps"
tmaps <- list.files(atlas.path)
atlas <- paste(atlas.path, tmaps[1], sep="/")

ofname <- "NODDIandPCAmaps_sample_info.csv"

### gene expression info
mainfile <- "MicroarrayExpression.csv"
pafile   <- "PACall.csv"

merged.dat <- c();
pa.dat     <- c();
sampleID   <- c();

for (s in samples){
  message(s)

  fname <- paste(aibs_raw, s, mainfile, sep="/")
  dat   <- read.csv(fname, head=F, row.names=1)
  merged.dat <- rbind(merged.dat, t(dat))

  fname <- paste(aibs_raw, s, pafile, sep="/")
  dat   <- read.csv(fname, head=F, row.names=1)
  pa.dat <- rbind(pa.dat, t(dat==1))

  sampleID <- c(sampleID, rep(s, ncol(dat)))
}


##sample info
sample.info <- c()
for (s in samples){
  message(s)
  fname <- paste(aibs_raw, s, "SampleAnnot.csv", sep="/")
  sinfo <- read.csv(fname, as.is=T)
  #Mirror right into left for the last two donors, and change the structure name field accordingly
  if (s %in% c("H0351.2001", "H0351.2002")){
    R_idx <- grep('right', sinfo$structure_name)
    sinfo[R_idx, "mni_x"] <- -1*sinfo[R_idx, "mni_x"]
    sinfo[R_idx, "structure_name"] <- gsub(pattern = "right", replacement = 'left', x = sinfo[R_idx, "structure_name"])
  }
  sample.info <- rbind(sample.info, sinfo)
}
##add brain iD
sample.info <- cbind(sample.info, sampleID)

##rename columns of merged.dat with probe names
probe.info <- read.csv(paste(aibs_raw,s,"Probes.csv",sep="/"), as.is=T)

#be on the safe side
probeid2name <- paste(probe.info[,"probe_name"])
names(probeid2name) <- probe.info[,"probe_id"]

colnames(merged.dat) <- probeid2name[colnames(merged.dat)]
colnames(pa.dat) <- colnames(merged.dat)

merged.dat <- data.frame(merged.dat)
pa.dat <- data.frame(pa.dat)

##get transform
Tra <- getTransform(atlas)

##mni_coords
mni_coords <- sample.info[,paste("mni",c("x","y","z"),sep="_")]

vox_coords <- t(apply(mni_coords, 1, function(x){
  mni2vox(x, Tra)
}))

info.new <- c()
for (mmap in tmaps){
  message(mmap)
  atlas.tmp <- readNIfTI(paste(atlas.path, mmap, sep="/"))
  probe2atr <- apply(vox_coords, 1, voxel2content_cube, atlas.tmp)
  info.new <- cbind(info.new, probe2atr)
}
colnames(info.new) <- c("Bingham_DAb", "Bingham_ODIp", "Bingham_ODIs", "Bingham_ODItot","Bingham_Vin","Bingham_Viso","FA","MD", "NODDI_ODI","NODDI_Vin","NODDI_Viso","PCA")
sample.info.new <- cbind(sample.info, info.new)

write.csv(sample.info.new, ofname, row.names=F)

#save RData file
save(merged.dat, pa.dat, sample.info.new, file="TestExpressionData.RData")

#######################################################################################################################
#######################################################################################################################