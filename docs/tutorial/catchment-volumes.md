## Runoff Depths vs Catchment Volumes

To perform river routing on a vector watershed definition, you need to catchment level volume time series. Runoff data is most commonly generated on
regular grids from a land surface process model. Runoff grid files may have a single file containing multiple time steps or multiple files which each
contain a single time step. You might also have multiple realizations of those data for each member in an ensemble simulation.

## Weight Table

You can use GIS methods to aggregate the distributed runoff depths to catchment level volumes. The steps are (1) intersect the grid cell boundaries
with the catchment boundaries to create polygons of each grid cell and catchment combination, then (2) multiply the runoff depth value of each polygon
times the area of that polygon (i.e. convert depth to volume). Repeat this for each time step. Repeat again for each ensemble member, if applicable.

Intersections of grid cells and catchment boundaries are computationally expensive and do not change. One method to avoid repeating the expensive
intersection step (and therefore increase speed) is to cache the intersection results in a "weight table" that describes which grid cells contribute
to which catchments. The cached results are unique to the grid cell geometry (and catchment boundaries). You need to recreate the cached intersections
if you use runoff grids with different resolutions or extents or if you change your watershed definition (catchment boundaries).
