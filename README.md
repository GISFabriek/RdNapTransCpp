# RdNapTransCpp
C++ 14 implementation of RDNAPTRANS™  
Meant as a first stab at modernizing the code.  
Modeled according to Java RdNapTrans (https://github.com/PDOK/rdnaptrans-java)  
Converts Spatial Coordinates from RD_New (EPSG:28992 https://spatialreference.org/ref/epsg/amersfoort-rd-new/) to ETRS89 (EPSG:4258 https://spatialreference.org/ref/epsg/4258/) and back.  
Original (C) source has been retrieved via:
https://zakelijk.kadaster.nl/transformatie-van-coordinaten  
The makefile project should compile on Windows, Linux and macOS
(Tested on Windows using both MSVC and Clang, on Linux and on macOS using GCC).
The three correction grid files **(nlgeo04.b64, x2c.b64, y2c.b64)** needed for the grid interpolation (see https://nl.wikipedia.org/wiki/Rijksdriehoekscoördinaten) are now included as base64 encoded string resources (making use of this resource compiler: https://github.com/vector-of-bool/cmrc).
Thus, we no longer need the three grid files to exist in the same directory as the executable (as was the case in the original C code).   
*If you want the code to run as part of another application, add everything except from Runner.cpp to your code base and call the methods in the Transformer class directly from your code.*
