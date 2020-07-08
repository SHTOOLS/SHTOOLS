---
title: MakeCircleCoord()
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: pymakecirclecoord.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Compute the coordinates of a circle placed at a given latitude and longitude.

## Usage

```python
`coord` = MakeCircleCoord (`lat`, `lon`, `theta0`, [`cinterval`])
```

## Returns

`coord` : float, dimension(360/`cinterval`, 2)
:   The latitude [:,0] and longitude [:,1] coordinates of the circle in degrees.

## Parameters

`lat` : float
:   The latitude of the center of the circle in degrees.

`lon` : float
:   The longitude of the center of the circle in degrees.

`theta0` : float
:   The angular radius of the circle in degrees.

`cinterval` : optional, float, default = 1
:   Angular spacing in degrees of the output latitude and longitude points. If not present, the default is 1.

## Description

`MakeCircleCoord` will calculate the latitude and longitude coordinates of a circle of angular radius `theta0` placed on a sphere at position (`lat`, `lon`). This is useful for plotting circles on geographic maps. The first index in the output vectors corresponds to the point directly north of the cirlce origin, and subsequent points are arranged in a clockwise manner. Input and output units are in degrees.
