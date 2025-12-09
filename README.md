## <font color="#6A5ACD">PBSwicked: Wicked Wonders from Another Universe</font> ##

<font color="red">&copy; Fisheries and Oceans Canada (2025)</font>

**PBSwicked** provides pseudo-R code (largely based on the Hadley Wickham R-verse that parallels our universe) to bedazzle fisheries stock assessment    scientists. It brings to mind the Thomas Dolby 1982 song *She Blinded Me with Science*. The functions herein are kept separate from other PBS packages, largely because there are so many R package dependencies.

**PBSwicked** depends on **PBStools**, **sdmTMB**, **sp**, **ggplot2**, **dplyr**, and **akima**. Lord knows what packages these depend on; however, the mildly curious can navigate to CRAN or GitHub and find out for themselves. 

<font color="red"><h3>Installation</h3></font>

Although **PBSwicked** is not available on <a href="https://cran.r-project.org/">CRAN</a> (Comprehensive R Archive Network), the source code appears on <a href="https://github.com/pbs-software/pbs-tools">GitHub</a> and can be built in R using:

`devtools::install_github("pbs-software/pbs-wicked/PBSwicked")`

However, at this point the package is not constructible. At some point it might be, but don't hold your breath. Once package formation has been achieved, the source code will be checked using CRAN's `R CMD check --as-cran` routine using a recent R-devel installation on a **Windows 11** 64-bit system.

<font color="red"><h3>Disclaimer</h3></font>

"Fisheries and Oceans Canada (DFO) GitHub project code is provided on an 'as is' basis and the user assumes responsibility for its use. DFO relinquishes control of the information and assumes no responsibility to protect the integrity, confidentiality, or availability of the information. Any claims against DFO stemming from the use of its GitHub project will be governed by all applicable Canadian Federal laws. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favouring by DFO. The Fisheries and Oceans Canada seal and logo, or the seal and logo of a DFO bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DFO or the Canadian Government.”

As with any freely available product, there is no warranty or promise that **PBSwicked** will perform adequately for all circumstances. Additionally, coding errors are inevitable, and users should contact the package maintainer if bugs are detected. Merci au revoir.

Maintainer: <a href="mailto:rowan.haigh@dfo-mpo.gc.ca">Rowan Haigh</a>


<!---<p align="right"><img src="DFOlogo_small.jpg" alt="DFO logo" style="height:30px;"></p>-->

<!---<img src="Uranus.png" alt="Description" width="300" height="300" style="opacity: 0.75;">-->

<img src="Uranus.jpg" alt="Uranus" style="width:5%" align="right"  hspace="5" />
<img src="DFOlogo_small.jpg" alt="dfo" style="width:50%" align="right" />
