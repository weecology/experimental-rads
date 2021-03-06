experimental-rads
=================

Research on the form of the rank-abundance distribution (RAD) using data from terrestrial animal taxa. 
Contains code for reproducing the results of the analyses using a dataset compiled by Sarah R. Supp from the published literature. 
Collaborators on this project include Sarah R. Supp and S. K. Morgan Ernest. Code was written by Sarah R. Supp.

The code and data in this repository allow for the analyses and figures to be fully replicated using a global-span 
local-scale database of maniuplated community data. 

The project and code in this repository are still under development. 

Requirements: R 3.x, 
R packages `vegan`, `BiodiversityR`, `plotrix`, `graphics`, `CCA`, `VGAM`, `nlme`, `lme4`, `languageR`, `poilog`, `scatterplot3d`, `hydroGOF`, and 
`VennDiagram`, `knitr`, `ggplot2`,  and the file containing functions specific to this code, `expRADsFunctions.R`.

The analyses can be replicated by changing the working directory at the top of the file `expRAD_md_script.R` to the location on 
your computer where you have stored the `.R` and `.csv` files.
If you wish to make your own version of the manuscript, change the pathname also at the top of the `expRAD_knitr.Rmd` and run.

Code should take less than 5 minutes to run start to finish. 
Figures should output as pdfs to your working directory.

Data use: Data is provided in this supplement for the purposes of replication. 
If you wish to use the data for additional research, they should be obtained from Sarah R. Supp (sarah@weecology.org).

Included Files: 
* `expRAD_md_script.R` script -- cleans up the data, runs the statistical analyses, and outputs figures.
* `ExpRADsFunctions.R` script -- holds the relevant functions for executing the script.
* `ExpRADs_knitr.Rmd` markdown doc -- used to recreate the manuscript and figures. Sources the previous two functions. Note that this is not the most updated version of the manuscript.
* `expRAD_dev.R` script -- holds some earlier code that I wanted to keep, not necessary for replicating the published project.
* `data` folder -- holds the 5 data tables described in the appendix to the paper (references, sites, experiments, communities, and comparisons)
    * `ref_analsyis_data.csv` - table with references for the data
    * `sites_analysis_data.csv` - table with site-level information
    * `experiments_analysis_data.csv` - table with experimental details for each site
    * `community_analysis_data.csv` - table with the species abundance community-level data for each site
    * `orderedcomparisons.csv` - table that shows which sites can be compared (as paired control-manipulation) based on the original references for each study


License This code is available under a BSD 2-Clause License.

Copyright (c) 2013 Weecology. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. 
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the 
documentation and/or other materials provided with the distribution. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS 
AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Contact information Sarah Supp's email: sarah@weecology.org and sarah.supp@stonybrook.edu

Sarah's website: http://weecology.org/people/sarahsupp/Sarah_Supp/
