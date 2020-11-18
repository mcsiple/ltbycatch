# ltbycatch: Projecting long-term marine mammal abundance with bycatch
master: <a href="https://badge.fury.io/gh/mcsiple%2Fltbycatch"><img src="https://badge.fury.io/gh/mcsiple%2Fltbycatch.svg" alt="GitHub version" height="18"></a>

ltbycatch is an R package that generates marine mammal population projections based on starting abundance, life history, and bycatch rates, based on the BALEEN II population dynamics model.

## Authors
Margaret Siple  
André Punt  
Tessa Francis  
Phil S. Hammond  
Dennis Heinemann  
Kristy J. Long  
Jeff Moore  
Maritza Sepulveda  
Randall R. Reeves  
Guðjón Már Sigurðsson  
Gísli Víkingsson  
Paul R. Wade  
Rob Williams  
Alexandre N. Zerbini  


## Contents
-   [Details](#details)
-   [Install](#install)
-   [Contributing](#contributing)
-   [References](#references)
<!-- end toc -->

## Details
This R package contains the functions used in the Marine Mammal Bycatch Impacts Exploration Tool (mmBIET), a Shiny app built by Margaret Siple and André Punt for the Ocean Modeling Forum's [Marine Mammal Bycatch Working Group](https://oceanmodelingforum.org/working-groups/marine-mammal-bycatch-working-group/). The app is available [here](https://msiple.shinyapps.io/mammaltool/). The functions in this package, and the app, are both intended to be used in cases where data on bycatch and/or population status are sparse or unavailable. 

Our target audience is stakeholders interested in projecting marine mammal populations to examine the impacts of bycatch. Those code could also be used as a teaching tool, or for anyone who is more familiar with R than FORTRAN and wants to use some components of the BALEEN II model (Punt 1999). 

## Install
This package can be downloaded directly from GitHub:
```
library(devtools)
install_github("mcsiple/ltbycatch")
```


## Contributing [![contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/dwyl/esta/issues)
We would like this package to be sustainable in the long term and welcome contributions. If you encounter a bug, please leave a note on the Issues page. You can also leave comments there about additional functionality. If you are interested in contributing, we direct you to the R package [contribution advice](http://r-pkgs.had.co.nz/git.html) from Hadley Wickham.

## Accessing the mmBIET Shiny app
The functions in this package can also be accessed through the Shiny app for this project, which is located online [here](https://msiple.shinyapps.io/mammaltool/). The app provides an easy way to explore outcomes. 

<img src="https://github.com/mcsiple/ltbycatch/blob/master/docs/screenshot1.png" alt="screenshot1" width="400">

## References
Punt, A. E. 1999. Annex R: A full description of the standard Baleen II model and some variants thereof. Division of Marine Research, CSIRO Marine Laboratories, Hobart, Australia. Available from https://duwamish.lib.washington.edu/uwnetid/illiad.dll?Action=10&Form=75&Value=1651729 (accessed August 7, 2018).

## Citation (t.b.d.)
Margaret C. Siple, André E. Punt, Tessa B. Francis, Phil S. Hammond, Dennis Heinemann, Kristy J. Long, Jeff Moore, Maritza Sepulveda, Randall R. Reeves, Guðjón Már Sigurðsson, Gísli Víkingsson, Paul R. Wade, Rob Williams, and Alexandre N. Zerbini (t.b.d.). ltbycatch: Projecting long-term marine mammal abundance with bycatch. R package version 1.0.0. url: https://github.com/mcsiple/ltbycatch
