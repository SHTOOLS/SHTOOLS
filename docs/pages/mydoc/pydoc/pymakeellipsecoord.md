---
title: MakeEllipseCoord (Python)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: pymakeellipsecoord.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Compute the coordinates of an ellipse placed at a given latitude and longitude.

## Usage

`coord` = MakeEllipseCoord (`lat`, `lon`, `dec`, `a_theta`, `b_theta`, [`cinterval`])

## Returns

`coord` :float, dimension(360/`cinterval`, 2)
:   The latitude [:,0] and longitude [:,1] coordinates of the ellipse in degrees.

## Parameters

`lat` : float
:   The latitude of the center of the ellipse in degrees.

`lon` : float
:   The longitude of the center of the ellipse in degrees.

`dec` : float
:   Rotation angle of the semi-major axis of the ellipse in degrees with respect to local north.

`a_theta` : float
:   The angular radius of the semi-major axis of the ellipse in degrees.

`b_theta` : float
:   The angular radius of the semi-minor axis of the ellipse in degrees.

`cinterval` : optional, float, default = 1
:   Angular spacing in degrees of the output latitude and longitude points. If not present, the default is 1.

## Description

`MakeEllipseCoord` will calculate the (lat, long) coordinates of an ellipse placed on a sphere at position (`lat`, `lon`). The semi-major and semi-minor axes, expressed in angular radii in degrees, are given by `a_theta` and `b_theta`, respectively. The semimajor axis of the ellipse is initially directed due north, and it is then rotated clockwise by the angle `dec`. This is useful for plotting ellipses on geographic maps. The first index in the output vectors corresponds to the northern rotated semimajor axis, and subsequent points are arranged in a clockwise manner. Input and output units are in degrees.

## See also

[makecirclecoord](pymakecirclecoord.html)
