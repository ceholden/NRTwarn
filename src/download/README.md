MODIS Tile Storage Location
==========================

# Data
## MOD09GQ / MYD09GQ
Daily 250m SR bands 1 - 2 (Red / NIR) from Terra and Aqua

https://lpdaac.usgs.gov/products/modis_products_table/mod09gq
https://lpdaac.usgs.gov/products/modis_products_table/myd09gq

## MOD09GA / MYD09GA
Daily 500m/1km SR from Terra and Aqua

https://lpdaac.usgs.gov/products/modis_products_table/mod09ga
https://lpdaac.usgs.gov/products/modis_products_table/myd09ga

# Download Utilities
The download script "fetch_modis.sh" relies on "get_modis.py" from:

    https://github.com/jgomezdans/get_modis

to be in your PATH.

# Example

An example of how to download all data required for 2000 - 2014:

    > for y in $(/usr/bin/seq 2000 2014); do 
        qsub -N dl_${y} fetch_modis.sh h10v08 $y
    done

This submits a QSUB job for every year to download data from the "h10v08" tile.
