# RockCoast10Be
Latest version: 2.0 (2025-07-01)

<em>RockCoast10Be.m</em> is a backward geometry-based model developped in MATLAB to explore cosmogenic <sup>10</sup>Be concentrations across an active shore platform as a function of cliff retreat and shore platform down-wearing scenarios. It may be used to either explore the expected concentrations or to find a best-fit scenario if 10Be data are available.

Please refer to the original paper when using this work: 
Swirad, Z. M. et al. 2020. Cosmogenic exposure dating reveals limited long-term variability in erosion of a rocky coastline. Nature Communcations  11: 3804. https://doi.org/10.1038/s41467-020-17611-9

The original code is provided as a supplement to Swirad et al. (2020). The GitHub repository contains a simplified and updated model. On the top of clearer nomenclature, the differences include:
1. The empirical model of shore platform erosion at Staithes, UK was removed because of its limited utility elsewhere;
2. The require tidal information was limited to a single value of tidal range and the tidal-dependent highest and lowest shore platform elevation (~highest and lowest astronomical tide);
3. Topographic shielding is based on cliff geometry (inclination and subtended azimuth zngles);
4. Best-fit scenario can be identified based on mean squared difference. If no <sup>10</sup>Be data are availabe, please comment lines #24 and #508:end.
