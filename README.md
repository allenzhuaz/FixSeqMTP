# FixSeqMTP

## Fixed Sequence Multiple Testing Procedures

[![](https://www.r-pkg.org/badges/version/FixSeqMTP?color=orange)](https://cran.r-project.org/package=FixSeqMTP) [![](http://cranlogs.r-pkg.org/badges/grand-total/FixSeqMTP?color=blue)](https://cran.r-project.org/package=FixSeqMTP) [![](https://img.shields.io/badge/lifecycle-stable-freshgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)

Overview
--------

This package implements several FWER, FDR and Directional Error (mdFWER) controlling procedures for fixed sequence hypotheses.

The methodology is described in the following papers:

[Lynch, G., Guo, W., Sarkar, S. K., & Finner, H. (2017). The control of the false discovery rate in fixed sequence multiple testing. Electronic Journal of Statistics, 11(2), 4649-4673.](https://projecteuclid.org/euclid.ejs/1510974129)

[Qiu, Z., Guo, W., & Lynch, G. (2015). On generalized fixed sequence procedures for controlling the FWER. Statistics in medicine, 34(30), 3968-3983.](https://onlinelibrary.wiley.com/doi/full/10.1002/sim.6603)

[Grandhi, A., Guo, W., & Romano, J. P. (2019). CONTROL OF DIRECTIONAL ERRORS IN FIXED SEQUENCE MULTIPLE TESTING. Statistica Sinica, 29, 1047-1064.](https://arxiv.org/abs/1602.02345)

Installation
------------

Open R console, install the pacakge directly from [CRAN](https://cran.r-project.org/web/packages/FixSeqMTP/index.html):

```r
install.packages("FixSeqMTP")
library(FixSeqMTP)
```

Or install the development version from GitHub, first make sure to install the `devtools` package:

```r
# install.packages("devtools")
devtools::install_github("allenzhuaz/FixSeqMTP")
```
