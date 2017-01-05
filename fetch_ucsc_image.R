library(RCurl)
library(png)

#retrieve a formatted ucsc png
#chrm - UCSC formatted chromosome, like chr1, chr2, etc.
#start, end - integer range of coordinates to view
#base_url - tells UCSC what to display and how it should appear.  default is for hg38.
fetch_ucsc_image = function(chrm, start, end, base_url = NULL){
  if(is.null(base_url)){
    # warning("Using default UCSC URL.  Only appropriate for hg38 and might not work at all.")
    # warning("To make your own URL, select the appropriate genome and desired display settings in UCSC genome browser. Append everything after 'hgsid=' to the following URL:\nhttp://genome.ucsc.edu/cgi-bin/hgRenderTracks?position=CHR:START-END&pubs=pack%27&refGene=pack%27&hgsid=")
    base_url = "http://genome.ucsc.edu/cgi-bin/hgRenderTracks?position=CHR:START-END&pubs=pack%27&refGene=pack%27&hgsid="
  }  
  sym_url = base_url
  sym_url = sub("CHR", chrm, sym_url)
  sym_url = sub("START", start, sym_url)
  sym_url = sub("END", end, sym_url)
  my_image <-  readPNG(getURLContent(sym_url))
  return(my_image)
}

#draw the UCSC image returned by fetch_ucsc_image() in existiong plot
#my_image - matrix returned by fetch_ucsc_image()
#lab_width - how many pixels to trim off the left side of image.  Used to remove default track name area.
#xmin, xmax, ymin, ymax - area in existing plot where UCSC image should appear.
plot_ucsc_image = function(my_image, lab_width = 117, xmin = 0, xmax = 1, ymin = 0, ymax = 1){
  img_height = nrow(my_image)
  img_width = ncol(my_image)
  #   par(mai = rep(0, 4))
  #   plot(x = c(xmin, xmax), y = c(0,1), axes = F, type = "n", xlab = "", ylab = "")
  #   axis(1, at = c(xmin, xmax))
  w = xmax - xmin
  fudge = w * lab_width / (800 - lab_width)
  xmin = xmin -  fudge
  #   rect(xleft = xmin, ybottom = ymin, xright = xmax, ytop = ymax)
  rasterImage(my_image, xleft = xmin, ybottom = ymin, xright = xmax, ytop = ymax)
  return(img_height)
}
