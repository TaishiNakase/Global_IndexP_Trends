# A retrospective analysis of climate-based dengue virus transmission suitability and population growth over the last four decades

This study uses estimates of dengue transmission suitability based on historical temperature and humidity data to examined how shifts in climate and human population growth have contributed to the change in the geographical distribution and size of the global population living in areas with high climate suitability from 1979 to 2022.

- Details on the methodology are described in: \
  [TO BE UPDATED]

## Version History 

**June 29, 2024**: Uploaded first version of code. 

## Data
- **Climate-based dengue virus transmission suitability data**: Gridded spatiotemporal data of climate-based dengue virus transmission suitability from 1979 to 2022 at a time resolution of 1 month and a spatial resolution of 360 arcseconds are publicly available on figshare (https://doi.org/10.6084/m9.figshare.21502614.v5).

## Contents

`indexP_trend_estimation.R`: This R script estimates long-term climate suitability trends using the seasonal Mann-Kendall test on the monthly Index P time series from 1979 to 2022. Sen's slope is used to estimate the linear rate of change for each pixel. Autocorrelation artifacts are corrected for using the prewhitening algorithm proposed by Yue et al (https://doi.org/10.1002/hyp.1095). 

`indexP_summarization.R`: This R script summarizes climate suitability for dengue virus transmission by calculating the average Index P over a specified time period. 

Details of specific R packages required for each step are described in the respective R scripts. 

## Support 

Please direct question or bug reports to tnakase@stanford.edu
