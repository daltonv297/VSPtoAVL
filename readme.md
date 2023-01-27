# OpenVSP to AVL Converter
### v1.3


This script creates an AVL file from the lifting geometry of an OpenVSP model.

## Instructions:
All desired geometry in OpenVSP must be of type Wing. To specify the surface from which to extract
reference dimensions, rename the geometry to "Wing_ref". Additionally, all wing geometry must have unique names. To
exclude wing geometry from the conversion, add an asterisk to its name.

To add a control surface, add a sub-surface to the wing geometry and name it one of the following: "aileron",
"elevator", "rudder", or "flap".

When specifying an airfoil from a file in OpenVSP, make sure that the first line in the file is the filename without the
extension. Example: the first line of SD7062.dat should be SD7062