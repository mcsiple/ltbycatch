# Package ltbycatch
master: <a href="https://badge.fury.io/gh/mcsiple%2Fltbycatch"><img src="https://badge.fury.io/gh/mcsiple%2Fltbycatch.svg" alt="GitHub version" height="18"></a>

ltbycatch is an R package that generates marine mammal population projections based on starting abundance, life history, and bycatch rates, based on the BALEEN II population dynamics model.

## Authors
Margaret Siple  
André Punt


Contents
--------

-   [details](#details)
-   [Install](#install)
-   [Contributing](#contributing)
-   [References](#references)
<!-- end toc -->

Details
-------------------------------
This R package contains the functions used in the Marine Mammal Bycatch Impacts Exploration Tool (mmBIET), a Shiny app built by Margaret Siple and André Punt for the Ocean Modeling Forum's [Marine Mammal Bycatch Working Group](https://oceanmodelingforum.org/working-groups/marine-mammal-bycatch-working-group/). The app is available [here](https://msiple.shinyapps.io/mammaltool/). The functions in this package, and the app, are both intended to be used in cases where data on bycatch and/or population status are sparse or unavailable. 

Our target audience is stakeholders interested in projecting marine mammal populations to examine the impacts of bycatch. Those code could also be used as a teaching tool, or for anyone who is more familiar with R than FORTRAN and wants to use some components of the BALEEN II model (Punt 1999). 

## Install
This package can be downloaded directly from GitHub:
```
library(devtools)
install_github("mcsiple/ltbycatch")
```



## Contributing [![contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/dwyl/esta/issues)


## References
Punt, A. E. 1999. Annex R: A full description of the standard Baleen II model and some variants thereof. Division of Marine Research, CSIRO Marine Laboratories, Hobart, Australia. Available from https://duwamish.lib.washington.edu/uwnetid/illiad.dll?Action=10&Form=75&Value=1651729 (accessed August 7, 2018).

