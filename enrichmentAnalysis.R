# =============================================================================.
# This is a workflow script to be executed secondly (after dataPreparation.sh)
# when running the demo analysis. The script chains normalization, probes
# filtering and multi-resolution analysis of the 4C enrichments.
# =============================================================================.

# =============================================================================.
# Load R package
# -----------------------------------------------------------------------------.
library("MRA.TA")

# =============================================================================.
# Import microarray data
# -----------------------------------------------------------------------------.
# Read predefined experimental setup
experiment.design <- read.delim(
  "Experiment_Design/GSE23887_dm6.txt", stringsAsFactors=F
)

# Read raw data (two replicates)
array.data.format <- file.formats$nimblegen.pair
r1.ct <- readData(experiment.design$SAMPLE_PATH[1], array.data.format) # control
r1.4C <- readData(experiment.design$SAMPLE_PATH[2], array.data.format) # 4C
r2.ct <- readData(experiment.design$SAMPLE_PATH[3], array.data.format) # control
r2.4C <- readData(experiment.design$SAMPLE_PATH[4], array.data.format) # 4C

# Read updated and cleaned array design annotations
array.design.path <- unique(experiment.design$DESIGN_PATH)
array.design.format <- file.formats$nimblegen.updated
array.design <- readData(array.design.path, array.design.format)

# Match array data to updated design using probe identifiers
r1.4C <- r1.4C[match(array.design$PROBE_ID, r1.4C$PROBE_ID),]
r1.ct <- r1.ct[match(array.design$PROBE_ID, r1.ct$PROBE_ID),]
r2.4C <- r2.4C[match(array.design$PROBE_ID, r2.4C$PROBE_ID),]
r2.ct <- r2.ct[match(array.design$PROBE_ID, r2.ct$PROBE_ID),]

# =============================================================================.
# Import coordinates of restriction fragments and match with probe coordinates
# -----------------------------------------------------------------------------.
# Read precomputed restriction fragments
fragments <- readData(
  "Processed_Data/RestrictionMap/dm6_DpnII_Fragments.txt.gz",
  file.formats$restriction.fragment
)

# Define restriction fragments’ start and end (which depends on the enzyme)
fr.len <- with(fragments, RS3.START - RS5.END + 1)
fr.start <- fragments$RS5.END
fr.end <- fragments$RS3.START

# Make genomic intervals with probe and restriction fragment coordinates
probes.grg <- GRanges(
  seqnames=as.character(array.design$CHR),
  ranges=IRanges(
    start=array.design$START, end=array.design$END,
    names=as.character(array.design$PROBE_ID)
  ),
  strand=array.design$STRAND
)
fragments.grg <- GRanges(
  seqnames=as.character(fragments$CHR),
  ranges=IRanges(
    start=fr.start, end=fr.end,
    names=as.character(fragments$RFID)
  )
)

# Match probes to fragments using respective genomic intervals
probes.grg <- matchProbesToFragments(probes.grg, fragments.grg)

# Reorder data consistently with matched probes 
idx <- match(names(probes.grg), array.design$PROBE_ID)
array.design <- array.design[idx,] # For consistency with data and probes.grg
r1.ct <- r1.ct[idx,]; r1.4C <- r1.4C[idx,]
r2.ct <- r2.ct[idx,]; r2.4C <- r2.4C[idx,]

# =============================================================================.
# Normalization procedure
# -----------------------------------------------------------------------------.
# Compute raw A and M values in each replicate
r1.A <- (log2(r1.4C$PM) + log2(r1.ct$PM))/2 # average
r1.M <- (log2(r1.4C$PM) - log2(r1.ct$PM))   # log2 ratio
r2.A <- (log2(r2.4C$PM) + log2(r2.ct$PM))/2 # average
r2.M <- (log2(r2.4C$PM) - log2(r2.ct$PM))   # log2 ratio

# Apply background bias correction and lowess normalization to each replicate
# Replicate 1
res <- normalizeArrayData(
  r1.A, r1.M, name="r1_norm", plots=TRUE
)
r1.A <- res$A; r1.M <- res$M # Update A and M with normalized values
# Replicate 2
res <- normalizeArrayData(
  r2.A, r2.M, name="r2_norm", plots=TRUE
)
r2.A <- res$A; r2.M <- res$M # Update A and M with normalized values

# =============================================================================.
# Normalization procedure (Step by step alternative)
# -----------------------------------------------------------------------------.
if(T) {
  # Compute raw A and M values in each replicate
  r1.A <- (log2(r1.4C$PM) + log2(r1.ct$PM))/2
  r1.M <- (log2(r1.4C$PM) - log2(r1.ct$PM))
  r2.A <- (log2(r2.4C$PM) + log2(r2.ct$PM))/2
  r2.M <- (log2(r2.4C$PM) - log2(r2.ct$PM))
  
  # Estimate background bias
  pdf("backgroundBiasEstimation.pdf", width=5, height=5)
  bb.r1 <- backgroundBiasEstimation(r1.A, r1.M, plots = T)
  bb.r2 <- backgroundBiasEstimation(r2.A, r2.M, plots = T)
  dev.off()

  # Correct background bias
  res <- backgroundBiasCorrection(r1.A, r1.M, theta=bb.r1)
  r1.A <- res$x; r1.M <- res$y
  res <- backgroundBiasCorrection(r2.A, r2.M, theta=bb.r2)
  r2.A <- res$x; r2.M <- res$y

  # Apply lowess correction
  pdf("lowessCorrection.pdf", width=5, height=5)
  res <- lowessCorrection(r1.A, r1.M, lowess.f=0.2, plots = T)
  r1.A <- res$x; r1.M <- res$y
  res <- lowessCorrection(r2.A, r2.M, lowess.f=0.2, plots = T)
  r2.A <- res$x; r2.M <- res$y
  dev.off()
}

# =============================================================================.
# Probe selection
# -----------------------------------------------------------------------------.
# Identify probes likely associated with best/problematic signals with raw data
pdf("probesFitering.r1.pdf", width=5, height=5)
r1.EQ <- enrichmentQuality(r1.ct$PM, r1.4C$PM, plots=T) # Replicate 1
dev.off()
pdf("probesFitering.r2.pdf", width=5, height=5)
r2.EQ <- enrichmentQuality(r2.ct$PM, r2.4C$PM, plots=T) # Replicate 2
dev.off()

# Make control plot of normalized 4C signal versus probe distance to site
pdf("probeDistanceControls_r1.M.pdf", width=5, height=5)
plotProbeDistanceControls(r1.M, rnk=probes.grg$RF_RANK, dis=probes.grg$RF_DIST, dlim=c(-1500, 1500), QF=r1.EQ)
dev.off()
pdf("probeDistanceControls_r2.M.pdf", width=5, height=5)
plotProbeDistanceControls(r2.M, rnk=probes.grg$RF_RANK, dis=probes.grg$RF_DIST, dlim=c(-1500, 1500), QF=r2.EQ)
dev.off()

# -----------------------------------------------------------------------------.
# Filtering option 1:
# fragment too short OR no fragment assigned OR distance to site too large
reject <- with(probes.grg, RF_LEN < 50 | is.na(RF_ID) | RF_DIST>250)

# -----------------------------------------------------------------------------.
# Filtering option 2:
# Reject when fragment too short OR no fragment assigned OR rank to site too large
reject <- with(probes.grg, RF_LEN < 50 | is.na(RF_ID) | RF_RANK>2)

# -----------------------------------------------------------------------------.
# Always filter out probes associated to a problematic signal in both replicates
reject <- reject | (r1.EQ$is.worst & r2.EQ$is.worst)

# Make control plot of selected probes (4C signal versus probe distance to site)
pdf("probeSelection_r1.pdf", width=5, height=5)
plotSelectedProbes(r1.A, r1.M, dis=probes.grg$RF_DIST, sel = ! reject)
dev.off()
pdf("probeSelection_r2.pdf", width=5, height=5)
plotSelectedProbes(r2.A, r2.M, dis=probes.grg$RF_DIST, sel = ! reject)
dev.off()

# Update probes and normalized data to retain accepted probes only
probes.grg <- probes.grg[! reject]
r1.A <- r1.A[! reject] # For consistency with r1.M and probes.grg
r1.M <- r1.M[! reject]
r2.A <- r2.A[! reject] # For consistency with r2.M and probes.grg
r2.M <- r2.M[! reject]

# =============================================================================.
# 4C enrichement scores
# -----------------------------------------------------------------------------.
# Calculate the number of accepted probes per half-fragment (5’ and 3’ ends)
probes.grg <- countProbesPerFragment(probes.grg) 

# Average 4C enrichments in each half-fragment
y1 <- combineByFragment(r1.M, probes.grg, FUN=mean) # Replicate 1
y2 <- combineByFragment(r2.M, probes.grg, FUN=mean) # Replicate 2

# Compute statistical scores
y1$SCORE <- enrichmentScore(y1$VALUE)
y2$SCORE <- enrichmentScore(y2$VALUE)

# Pool replicates into combined 4C enrichment scores
Yi <- y1$SCORE + y2$SCORE
Yi <- sapply(-2*Yi, pchisq, df=2*2, lower.tail=FALSE, log.p=TRUE)

# =============================================================================.
# Multi-resolution visualization of 4C enrichement profiles (domainogram)
# -----------------------------------------------------------------------------.
# Compute genomic coordinates for visualization
vis.coords <- visualizationCoordinates(y1$SITE, y1$SITE)

# 4C-anchoring (bait) location
anchors <- read.table(
  "4C_Bait_Sequences/4C_Primers_Aligned_dm6.txt", sep='\t', header=FALSE
)
x.bait <- anchors$V4

# Define the range of genomic locations to be represented
xlim <- c(min(vis.coords$start), max(vis.coords$end))

# Define the range of resolutions (window sizes) to be represented
wlim <- c(1, 2^16)
w2y <- wSize2yAxis(wmax=wlim[2], logscale=T)

# Plot domainogram
png("4C_Fragment-based_3R.png", width=4096, height=1024)
domainogram(Yi, vis.coords$start, vis.coords$end, w2y, xlim, wlim)
abline(v=x.bait, col=rgb(1,0,1), lwd=2) # highlight 4C-anchoring location
dev.off()

# =============================================================================.
# 4C-enriched chromatin domains
# -----------------------------------------------------------------------------.
# Multi-resolution segmentation
sgm <- segmentation(
  Yi, name="4C_segmentation", wmax=2^12, wmin=5, gamma=1E-4, Tw=-10
)

# Load segmentation results
opts <- read.delim(sgm$file.segments, stringsAsFactors=F, skip=1)
doms <- read.delim(sgm$file.domains, stringsAsFactors=F)
doms.mr <- read.delim(sgm$file.maxresolution, stringsAsFactors=F)
doms.ms <- read.delim(sgm$file.maxscale, stringsAsFactors=F)

# Visualize segmentation results
png("4C-enriched_Segments_3R.png", width=4096, height=1024)
plotOptimalSegments(opts, vis.coords$start, vis.coords$end, w2y, xlim, wlim)
dev.off()
png("4C-enriched_ChromatinDomains_3R.png", width=4096, height=1024)
plotDomains(doms, vis.coords$start, vis.coords$end, w2y, xlim, wlim, col="black")
dev.off()

# Visualize max. resolution and max. scale domains
png("4C-enriched_FilteredDomains_3R.png", width=4096, height=1024)
plotDomains(doms.mr, vis.coords$start, vis.coords$end, w2y, xlim, wlim, col=rgb(0,1,0,0.5), border=rgb(0,0,0,0), lwd=2, lty=1)
plotDomains(doms.ms, vis.coords$start, vis.coords$end, w2y, xlim, wlim, col=rgb(1,0,0,0.5), border=rgb(0,0,0,0), lwd=2, lty=1, add=T)
dev.off()

# Retreive genomic coordinates of 4C-enriched chromatin domains
doms.mr$RS_START <- y1$SITE[doms.mr$start]; doms.mr$RS_END <- y1$SITE[doms.mr$end]
doms.ms$RS_START <- y1$SITE[doms.ms$start]; doms.ms$RS_END <- y1$SITE[doms.ms$end]
