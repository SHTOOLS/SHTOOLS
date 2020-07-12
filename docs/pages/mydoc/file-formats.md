---
title: "File formats"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: file-formats.html
summary: pyshtools can read spherical harmonic coefficients using several standard file formats.
toc: true
folder: mydoc
---

<style>
table:nth-of-type(n) {
    display:table;
    width:100%;
}
</style>

## Supported file formats

It often arises that you are given a file containing spherical harmonic coefficients, and that you need to import these into your program. pyshtools supports several standard file formats for spherical harmonic coefficients, and also provides support for skipping or reading header lines. In most cases, spherical harmonic coefficients can be read from a file simply by specifying the optional argument `format` as in these examples:

```python
hlm = pysh.SHCoeffs.from_file('file.sh', format='shtools')
glm = pysh.SHMagCoeffs.from_file('file.shc', format='dov')
clm = pysh.SHGravCoeffs.from_file('file.gfc', format='icgem')
```

### shtools

The 'shtools' file format is the default. This text-based format contains comment lines that start with '#', up to two header lines, and lines with the degree, order, and cosine and sine coefficients. Optionally, the file can also contain the error coefficients. The beginning of the file should look something like this:

```bash
# Comment lines start with '#' and are ignored.
# The next line is a header that will be read and stored. Up to two headers can be included.
0.17380000E+04, 0.49028002E+04, 0.73047000E-05,  900,  900
1,    0,  0.00000000E+00,  0.00000000E+00, 0.00000000E+00, 0.00000000E+00
1,    1,  0.00000000E+00,  0.00000000E+00, 0.00000000E+00, 0.00000000E+00
2,    0, -0.90881248E-04,  0.00000000E+00, 0.14809068E-11, 0.00000000E+00
2,    1,  0.23517464E-09,  0.10444740E-08, 0.16629736E-11, 0.17906454E-11
2,    2,  0.34674086E-04, -0.13397766E-09, 0.27246809E-11, 0.27256909E-11
```

The spherical harmonic coefficients in the file should be formatted as

```bash
l, m, coeffs[0, l, m], coeffs[1, l, m]
```

where *l* and *m* are the spherical harmonic degree and order, respectively. If the errors are included, the line should be formatted as

```bash
l, m, coeffs[0, l, m], coeffs[1, l, m], errors[0, l, m], errors[1, l, m]
```

For each value of increasing degree, all the angular orders are listed in increasing order. The list of spherical harmonic coefficients can start at any degree (such as degree 1 in the above example). For gravity and magnetic field data, the header lines often contain the GM and reference radius of the model.

### dov

The 'dov' format (for degree, order, value) is nearly identical to the 'shtools' format. Instead of providing both the cosine and sine coefficients on a single line, each line of these files contains only a single coefficient. The cosine coefficient has a positive order, whereas the sine coefficient has a negative order. The beginning of the file should look something like this:

```bash
Reference radius:  0.339350E+04
  1          0        -0.1776756860E+01
  1          1        -0.3727294260E+00
  1         -1        -0.3445717196E+00
  2          0        -0.1093718332E+00
  2          1         0.1045706028E+01
  2         -1         0.4747853109E-01
  2          2         0.9616234618E+00
  2         -2         0.8172537989E+00
```

The spherical harmonic coefficients in the file should be formatted as

```bash
l, m, coeffs[0, l, m]
l, -m, coeffs[1, l, m]
```

If the errors are to be read, the line should be formatted as

```bash
l, m, coeffs[0, l, m], errors[0, l, m]
l, -m, coeffs[1, l, m], errors[1, l, m]
```

### bshc

The 'bshc' format (for binary spherical harmonic coefficients) is a binary format that does not include any metadata or errors. The file is composed solely of (little endian) 8-byte floats, starting with the minimum and maximum degree, and followed by the cosine coefficients and then sine coefficients (with all orders being listed, one degree at a time). For a 100 degree file, the contents would be the following:

```bshc
0 100
C(0,0), C(1,0), C(1,1), C(2,0), C(2,1), ... C(100,99), C(100,100)
S(0,0), S(1,0), S(1,1), S(2,0), S(2,1), ... S(100,99), S(100,100)
```

where *C(l,m)* and *S(l,m)* are the cosine and sine coefficients, respectively.

### icgem

The 'icgem' format is a file format for gravitational potential coefficients used by the [International Centre of Global Earth Models](http://icgem.gfz-potsdam.de/home). These files contain a multiline header that provides metadata such as GM, the reference radius, and the type of errors. Following the header, each line contains spherical harmonic coefficients and errors as described [here](http://icgem.gfz-potsdam.de/ICGEM-Format-2011.pdf).