##################
#  ABIDE data reader (paper Sec. VIII)
#  Reads in the NYU ROI time series of ABIDE I and generates n*p*T arrays h[i,j,t],
#  one for the ASD group (h.1) and one for the control group (h.2).
#
#  Data: ABIDE I, Preprocessed Connectomes Project (PCP), pipeline cpac, strategy
#  nofilt_noglobal, derivative rois_aal (AAL atlas, p = 116 ROIs), site NYU
#  (T = 176 time points). Each subject is one file <FILE_ID>_rois_aal.1D, e.g.
#  NYU_0050952_rois_aal.1D, where FILE_ID is the column of the phenotypic file
#    https://s3.amazonaws.com/fcp-indi/data/Projects/ABIDE_Initiative/Phenotypic_V1_0b_preprocessed1.csv
#  (rows with SITE_ID == "NYU"). The files can be downloaded from
#    https://s3.amazonaws.com/fcp-indi/data/Projects/ABIDE_Initiative/Outputs/cpac/nofilt_noglobal/rois_aal/<FILE_ID>_rois_aal.1D
#  or with PCP's download_abide_preproc.py (-d rois_aal -p cpac -s nofilt_noglobal -t NYU -o <dir>).
#  PCP provides 175 NYU subjects; the paper uses 171 of them (73 ASD, 98 controls).
#  The other four (SUB_ID 50952, 51030, 51127, 51159) were not part of the data set
#  used in the paper and are skipped below if present.
#
#  Usage:   Rscript real-data/01_read_abide.R   (from the repository root)
#  Inputs:  data/abide/Phenotypic_V1_0b_preprocessed1.csv
#           data/abide/rois_aal/<FILE_ID>_rois_aal.1D (FILE_ID = NYU_<7-digit SUB_ID>)
#  Outputs: data/time.series.ABIDE.RData (h.1: 73 x 116 x 176, h.2: 98 x 116 x 176)
##################

if (!file.exists("functions/bootstrap_sample.R")) stop("Run this script from the root of the GA-Experiments repository.")

abide.path <- "data/abide/"
data.path <- "data/abide/rois_aal/"
file.list <- list.files(data.path, "^NYU_[0-9]+_rois_aal\\.1D$")
if(length(file.list) == 0) stop("No NYU .1D files found in ", data.path)
phenotypic <- read.csv(paste(abide.path,"Phenotypic_V1_0b_preprocessed1.csv",sep=""))

# extract SUB_ID
sub_id <- as.numeric(sapply(strsplit(file.list,"_"),"[[",2))

# skip the four NYU subjects that are not in the data set of the paper
not.used <- sub_id %in% c(50952, 51030, 51127, 51159)
file.list <- file.list[!not.used]
sub_id <- sub_id[!not.used]

f1 <- read.delim(paste(data.path, file.list[1], sep=""))
f1 <- as.matrix(f1)

n <- length(file.list)
p <- dim(f1)[2]
tau <- dim(f1)[1]

h <- array(0, c(n,p,tau))
dx_group <- rep(0,n)

for(i in 1:n){
  f1 <- read.delim(paste(data.path, file.list[i], sep=""))
  f1 <- as.matrix(f1)
  # every NYU file is 176*116 (time points * ROIs)
  h[i,,] <- t(f1)

  # finding diagnosis label (DX_GROUP: 1 = ASD, 2 = control)
  dx_group[i] <- phenotypic$DX_GROUP[which(phenotypic$SUB_ID == sub_id[i])]
}


# 73 subjects have DX_GROUP 1 (ASD) and 98 have DX_GROUP 2 (control), 171 overall.
if(sum(dx_group==1) != 73 || sum(dx_group==2) != 98) warning("Expected 73 ASD and 98 control subjects; the subject removals in functions/load_abide.R assume this data set.")
h.1 <- h[which(dx_group==1),,]
h.2 <- h[which(dx_group==2),,]

dir.create("data", recursive = TRUE, showWarnings = FALSE)
save(h.1, h.2, file="data/time.series.ABIDE.RData")
