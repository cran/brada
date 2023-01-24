## ---- include = FALSE---------------------------------------------------------
old <- options()
options(rmarkdown.html_vignette.check_title = FALSE)
# chunk option dev="png" is the default raster graphics format for HTML output
knitr::opts_chunk$set(dev="png")
options(mc.cores=2)

## ----fig.align='center', echo=FALSE, fig.cap="The predictive probability design for binary endpoints", out.width = "700px", dpi = 100----
knitr::include_graphics("pp-min.png")

## ---- warning=FALSE, message=FALSE, class.source="chunks", class.output="bg-warning"----
library(brada)

## ---- warning=FALSE, message=FALSE, fig.align="center", class.source="chunks"----
design = brada(Nmax = 40, batchsize = 5, nInit = 10, 
               p_true = 0.2 , p0 = 0.2, p1 = 0.2, 
               nsim = 100, 
               cores = 2)

## ---- warning=FALSE, message=FALSE, fig.align="center", fig.height=6, class.source="chunks"----
design = brada(Nmax = 40, batchsize = 5, nInit = 10, 
               p_true = 0.2 , p0 = 0.2, p1 = 0.2, 
               nsim = 100,
               a0 = 1, b0 = 1, 
               theta_T = 0.90, theta_L = 0.1, theta_U = 1, 
               method = "PP",
               cores = 2)

## ---- class.source="chunks", class.output="output"----------------------------
new_design = brada(Nmax = 40, batchsize = 5, nInit = 10, 
                   p_true = 0.2 , p0 = 0.2, p1 = 0.2, 
                   nsim = 100,
                   a0 = 0.5, b0 = 0.5, 
                   theta_T = 0.90, theta_L = 0.05, theta_U = 0.95, 
                   method = "PP",
                   cores = 2)

## ---- class.source="chunks", class.output="output"----------------------------
summary(design)

## ---- class.source="chunks", class.output="output"----------------------------
summary(new_design)

## ----fig.align='center', fig.width = 7, fig.height = 5, out.width = "600", out.height = "450", class.source="chunks", dpi=100----
plot(design)

## ---- class.source="chunks", class.output="output"----------------------------
trial_data = generateData(p = 0.2, Nmax = 40, nsim = 100, seed = 420)
library(DT)
datatable(trial_data, rownames = FALSE, filter="top", 
          options = list(pageLength = 5, scrollX=T),
          caption = 'Simulated trial data for a phase IIA trial with a binary endpoint')

## ---- class.source="chunks", class.output="output"----------------------------
head(new_design$trials)

## ---- class.source="chunks", class.output="output"----------------------------
sum(new_design$trials[,6] < 0.05)/100

## ---- include = FALSE---------------------------------------------------------
options(old)

