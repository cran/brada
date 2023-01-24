## ---- include = FALSE---------------------------------------------------------
old <- options() 
options(rmarkdown.html_vignette.check_title = FALSE)
# chunk option dev="png" is the default raster graphics format for HTML output
knitr::opts_chunk$set(dev="png")
options(mc.cores=2)

## ---- warning=FALSE, message=FALSE, class.source="chunks"---------------------
library(brada)

## ---- class.source="chunks"---------------------------------------------------
design = brada(Nmax = 30, batchsize = 5, nInit = 10, 
                   p_true = 0.4, p0 = 0.4, p1 = 0.4, 
                   nsim = 100, 
                   theta_T = 0.90, theta_L = 0.1, theta_U = 1, 
                   method = "PP",
                   cores = 2)

## ---- class.source="chunks", class.output="output"----------------------------
monitor(design, obs = c(0,1,0,0,0,0,0,1,0,0))

## ---- include = FALSE---------------------------------------------------------
options(old)

