# Community_Energy_Storage
Code used for analysis in Community Energy Storage paper

This repo contains ipynb files used in the analysis from the paper "Community Energy Storage: a smart choice for the smart grid?"

PREP: To run the analysis, first unzip the input data, creating a folder called "InputData".
Folders called "IntermediateData" and "Results" should also be created and will be used.

1. First run the file "CreateRoadNetworkMicrogrids.ipynb". This will create the communities joined by the road network.
NOTE: The seeding for the microgrids has already been done via k-means clustering, so is included in the input data.

2. Run the file "RadPlotter.ipynb" to generate the translation between the solar irradiance data in Austin and in Boston

3. Run "PecanStreetInspect-15min.ipynb". 
This file takes the input data from Pecan and creates demand and generation profiles in Cambridge.
Each household is a node in step two and is part of a community. Each demand and potential generation simulated.

4. Run "CommunityEnergyStorage96_Month.ipynb". This performs the analysis using the intermediate data generated in step 3.
The solar represents the average month (though this can easily be changed).
Batteries are simulated for all households with PV and all communities.

5. Run "CommunityEnergyStorage96_SolarSensitivity.ipynb".
This performs the same analysis with different solar months.
