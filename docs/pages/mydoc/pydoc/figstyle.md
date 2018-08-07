---
title: figstyle
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: figstyle.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

<style>
table:nth-of-type(n) {
    display:table;
    width:100%;
}
table:nth-of-type(n) th:nth-of-type(2) {
    width:75%;
}
</style>

A collection of matplotlib style sheets for creating publication quality graphics.

## Usage

```python
import matplotlib.pyplot as plt
import pyshtools.utils.figstyle as figstyle
plt.style.use(figstyle.shtools) # use the shtools style file
plt.style.use([figstyle.shtools, figstyle.half]) # combine multiple style files
```

## Styles

| Style | Description |
| ----- | ----------- |
| `shtools` | Core style file for pyshtools. |
| `full` | Set the figure width and height for figures that span the entire width of a page. |
| `threequarters` | Set the figure width and height for figures that span three-quarters of a page (default for `shtools` style). |
| `half` | Set the figure width and height for figures that span half a page. |
| `map` | Set parameters for plotting global maps. |
