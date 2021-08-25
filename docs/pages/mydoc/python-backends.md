---
title: pysh.backends
keywords: pyshtools, backends, shtools, ducc
sidebar: mydoc_sidebar
permalink: python-backends.html
summary: This module defines functions for accessing and setting optional parameters of the backends used for the spherical harmonic transforms in pyshtools.
toc: false
folder: mydoc
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

| Function name | Description |
| ------------- | ----------- |
| [select_preferred_backend](select_preferred_backend.html) | Set the preferred backend module and options. |
| [preferred_backend](preferred_backend.html) | Return the name of the current preferred backend. |
| [preferred_backend_module](preferred_backend_module.html) | Return a reference to the preferred backend module. |
| [backend_module](backend_module.html) | Return a reference to the specified backend module. |
