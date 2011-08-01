# Moment Tensor Code
*Contributors: Frank Vernon (FLV), Rob Newman (RLN), Juan Reyes (JCR), Gert-Jan van den Hazel (GJV)*

### Moment tensor integration into BRTT Antelope Datascope using the Python interface to Antelope

### Dependencies:
You must have the following CSS3.0 schema extensions installed:
* moment_tensor_greensfuncs
* moment_tensor_images
These are located in $ANTELOPE_CONTRIB/data/schemas/css3.0.ext
where $ANTELOPE_CONTRIB is your path to the contributed code repository.

### History:

#### Fall 2010
* Initial code review and breakup of tasks. 
* GJV writes dbmoment
* JCR writes Green's functions in Python from Fortran originals
* RLN creates schema extensions to Datascope for Green's functions

#### Spring 2011
* GJV and JCR continue work in Amsterdam on their component parts

#### Summer 2011
* RLN rewrite of dbmoment
* JCR finalize Fortran translation of Green's functions into Python
* Integration of dbmoment with native Python Green's functions

### Comments:
#### 2011-07-29
* dbmoment.py is not currently tethered to the moment tensor Datascope schema extension __moment_tensor_greensfuncs__
* dbmoment.py uses 4 local files for testing: data1, data2, data3 (the three components) and green (the Greens function, aka synthetic)
* These 4 data files need to be copied to your realtime system directly and placed in the db/data directory

#### 2011-08-01
* dbmoment.py __is now__ tethered to the moment tensor Datascope schema extensions
