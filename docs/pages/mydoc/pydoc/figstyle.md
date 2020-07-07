---
title: figstyle()
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: figstyle.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Set matplotlib parameters for creating publication quality graphics.

## Usage

figstyle([rel_width, screen_dpi, aspect_ratio, max_width, figsize, units])

## Parameters

**rel_width : float, optional, default = 0.75**
:   The relative width of the plot (from 0 to 1) wih respect to width.

**screen_dpi : int, optional, default = 114**
:   The screen resolution of the display in dpi, which determines the size of the plot on the display.

**aspect_ratio : float, optional, default = 4/3**
:   The aspect ratio (horizontal/vertical) of the plot.

**max_width : float, optional, default = 7.48031**
:   The maximum width of the usable area of a journal page in inches.

**figsize : list, optional, default = None**
:   A list containing the width and height of the figure.

**units : str, optional, default = 'i'**
:   The measurement units of max_width and figsize, either 'i' for inches or 'c' or 'cm' for centimeters.

## Notes

This function sets a variety of matplotlib parameters for creating
publication quality graphics. The default parameters are tailored to
AGU/Wiley-Blackwell journals that accept relative widths of 0.5, 0.75,
and 1. To reset the maplotlib parameters to their default values, use

matplotlib.pyplot.style.use('default')
    