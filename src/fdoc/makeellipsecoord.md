# MakeEllipseCoord

Compute the coordinates of an ellipse placed at a given latitude and longitude.

# Usage

call MakeEllipseCoord (`coord`, `lat`, `lon`, `dec`, `a_theta`, `b_theta`, `cinterval`, `cnum`, `exitstatus`)

# Parameters

`coord` : output, real*8, dimension(360/`cinterval`, 2)
:   The latitude (:,1) and longitude (:,2) coordinates of the ellipse in degrees.

`lat` : input, real\*8
:   The latitude of the center of the ellipse in degrees.

`lon` : input, real\*8
:   The longitude of the center of the ellipse in degrees.

`dec` : input, real\*8
:   Rotation angle of the semi-major axis of the ellipse in degrees with respect to local north.

`a_theta` : input, real\*8
:   The angular radius of the semi-major axis of the ellipse in degrees.

`b_theta` : input, real\*8
:   The angular radius of the semi-minor axis of the ellipse in degrees.

`cinterval` : optional, input, real\*8, default = 1
:   Angular spacing in degrees of the output latitude and longitude points. If not present, the default is 1.

`cnum` : optional, output, integer
:   Number of elements in the output arrays.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`MakeEllipseCoord` will calculate the (lat, long) coordinates of an ellipse placed on a sphere at position (`lat`, `lon`). The semi-major and semi-minor axes, expressed in angular radii in degrees, are given by `a_theta` and `b_theta`, respectively. The semimajor axis of the ellipse is initially directed due north, and it is then rotated clockwise by the angle `dec`. This is useful for plotting ellipses on geographic maps. The first index in the output vectors corresponds to the northern rotated semimajor axis, and subsequent points are arranged in a clockwise manner. Input and output units are in degrees.

# See also

[makecirclecoord](makecirclecoord.html)
