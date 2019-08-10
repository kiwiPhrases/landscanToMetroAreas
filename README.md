# Landscan to metropolitan areas

For this project, I use 2015 [Landscan](https://landscan.ornl.gov/) data to find metropolitan areas in Peru and compute De La Roca and Puga Exposure rates for each metro area. Metropolitan areas are determined in R while post-processing is performed in Python (thank you GeoPandas and Fiona). I am glad to have been able to find a way to render both arcGIS and qGIS dispensable for computational purposes (of course, they're still great for mapping). 

There are a few additional layers to the project. I used [GADM](https://gadm.org/) data for Peruvian administrative borders so I can geocode otherwise nameless metropolitan areas. I also used world borders from [Thematic Maps](http://thematicmapping.org/downloads/world_borders.php) to clip the world LandScan data down to Peru. The starting points for density and population thresholds and metro-finding algorythms for sattelite-image came from [Definition Matters](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3298761) though there are a few other references that you can find my notes on here. 

The project is under the behest of Professor Jorge De La Roca.
