# RdNapTransCpp
C++ 14 implementation of RDNAPTRANS™  
Meant as a first stab at modernizing the code.  
Modeled according to Java RdNapTrans (https://github.com/PDOK/rdnaptrans-java  )  
Converts Spatial Coordinates from RD_New (EPSG:28992 https://spatialreference.org/ref/epsg/amersfoort-rd-new/) to ETRS89 (EPSG:4258 https://spatialreference.org/ref/epsg/4258/) and back.  
Original (C) source has been retrieved via:
https://zakelijk.kadaster.nl/transformatie-van-coordinaten  
The makefile should compile on Windows and Linux (only tested on Windows to be honest).  
Make sure you have got the following grid files in the same directory as the executable: **nlgeo04.grd, x2c.grd, y2c.grd**  
(If you want the code to run as part of another application, add everything except from Runner.cpp to your code base and call the methods in the Transformer class directly from your code. * *Just make sure that the three mentioned files are in the same folder as your executable)* *.
