# Convert DNAStringSet to dataframe
dss2df = function(dss) data.frame(PeakID=names(dss), seq=as.character(dss), width=width(dss))