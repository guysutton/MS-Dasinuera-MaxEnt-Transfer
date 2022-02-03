# Custom function to convert DMS to decimal degrees
# I can't remember where I found this function online.
# - Please contact me if you are the author or know who the author is 
#   so that I can give you the appropriate credit. 

# Converting coordinates from DMS or DdM formats to decimal
# DMS = "degrees minutes seconds"; DdM = "degrees decimal minutes" 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dg2dec <- function(varb, Dg=NA, Min=NA, Sec=NA, SW.Hemisphere="S|W") {
  # Dg=decimal, Min=minutes and Sec=seconds; 
  # NOTE 1 - if the format is "degrees decimal minutes - DdM" (e.g. 40?? 26.767?? N) and not 
  # "degrees minutes seconds - DMS" (e.g. 40?? 26?? 46?? N), then call the function only with 
  # Dg and Min arguments, like dg2dec(varb, Dg="??", Min="??N"). 
  # Same applies when there is no seconds symbol (e.g. 45??12'7.38).
  # Note that one should not use blank spaces in Dg, Min or Sec arguments (will return NA).
  # For more details on formats see: 
  # https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#Coordinate_format_conversion
  
  # Use paste0("[", Dg, Min, Sec, "]") to build regex [] pattern
  # therefore, strsplit() will split string "varb" by what symbols you give to Dg, Min, Sec
  DMS <- sapply(strsplit(varb, paste0('[', Dg, Min, Sec, ']')), as.numeric)
  
  # DMS is a matrix; first row contains degrees; second - minutes; third - seconds.
  # If the format is "degrees decimal minutes" (e.g. 40?? 26.767?? N) and not 
  # "degrees minutes seconds" (e.g. 40?? 26?? 46?? N), then the matrix has only two valid rows: 
  # first row contains degrees; the second - minutes;
  # therefore, compute conversion for seconds only if there are more than 2 rows in DMS 
  # and Sec is different from NA (if there are seconds in the DMS format)
  decdg <- abs(DMS[1, ]) + DMS[2, ]/60 + ifelse(dim(DMS)[1] > 2  & !is.na(Sec), DMS[3, ]/3600, 0)
  
  # all cordinates from Southern or Western Hemispheres become negative in their decimal format
  SW <- grepl(pattern = SW.Hemisphere, x = varb, ignore.case = TRUE)
  return(ifelse(SW, -1, 1) * decdg)
}

#x <- read.table(text="DdM     DMS_lat      DMS_lon        DMS
#                45??12.123'  45??12'7.38''S  45??12'7.38''W  45??12'7.38
#                31??29.17'   31??29'10.2''N  31??29'10.2''E  31??29'10.2 ",
#                header=TRUE, stringsAsFactors=FALSE)
#x

# case of "degrees decimal minutes" (DdM)
#dg2dec(varb=x$DdM, Dg="??", Min="'")


