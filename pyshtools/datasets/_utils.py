"""Utility functions for pyshtools datasets."""

from dataclasses import dataclass
import re as _re

import pooch as _pooch
import numpy as _np


def _choose_sh_model(
    archive: _pooch.Pooch,
    user_lmax: int,
    downloader=_pooch.DOIDownloader(progressbar=True),
) -> tuple[str, int]:
    """FOR PYSHTOOLS DEVELOPER USE ONLY. When a model has multiple files with
    various resolutions (e.g. `Mars.MOLA_shape` has 5759, 2879, 1439, and 719),
    this function selects the most appropriate file to use based on the
    requested lmax. Specifically, if a high-res model with sufficient lmax
    already exists in the cache, it will be used. Otherwise, download the
    lowest-res model with sufficient lmax.

    Parameters
    ----------
    archive : _pooch.Pooch
        A Pooch object with a registry of files. The filenames must be
        identical except for a single integer (one or more digits, any
        position), which is taken to be the degree of the model.

    user_lmax : int
        Spherical harmonic degree requested by the user. If this value is
        negative or exceeds the highest available degree, the highest available
        degree will be used.

    downloader : _pooch.DOIDownloader | _pooch.HTTPDownloader, optional
        The downloader to use for fetching files. Defaults to
        `_pooch.DOIDownloader` with a progress bar.


    Returns
    -------
    fpath : str
        Path to an existing or downloaded file.
    target_lmax : int
        The max degree of the selected model (not necessarily the same as
        `user_lmax`).

    Raises
    ------
    ValueError
        If the filenames in the archive do not contain exactly one varying
        integer.

    Notes
    -----
    If `_choose_sh_model` is defined in `pyshtools/datasets/__init__.py`, we
    get a circular import error. Therefore we define it in a new file,
    `pyshtools/datasets/_utils.py`, and import directly from dataset files
    (e.g. `pyshtools/datasets/Ceres.py`).
    """

    # Find all integers in filenames
    def _extract_ints_from_fnames(fnames: tuple[str]) -> _np.ndarray:
        """Output is always 2D numpy array of integers.

        EXAMPLES:
            In:
                ["Mars_MOLA_shape_5759.bshc.gz",
                 "Mars_MOLA_shape_2879.bshc.gz"]
            Out:
                [ [5759],
                  [2879] ]

            In:
                ["Moon_LDEM128_shape_pa_1439.sh.gz",
                "Moon_LDEM128_shape_pa_719.sh.gz" ]
            Out:
                [ [128, 1439],
                  [128, 719] ]
        """
        regex: str = r"(\d+)"
        ints_in_fnames: list[list[str]] = []
        for fname in fnames:
            matches: list[str] = _re.findall(regex, fname)
            ints_in_fnames.append(matches)
        return _np.array(ints_in_fnames, dtype=int)

    fnames = tuple(archive.registry.keys())
    ints_in_fnames: _np.ndarray = _extract_ints_from_fnames(
        fnames
    )  # numpy array (2d, ints)

    # Ignore any integers that are constant across all filenames
    def _extract_sh_degrees(ints2d: _np.ndarray) -> _np.ndarray:
        """
        EXAMPLES:
            In:
                [ [128, 1439],
                  [128, 719] ]
            Out:
                [1439, 719]
        """
        # boolean mask of "varies relative to the first row"
        varying: _np.ndarray = (ints2d != ints2d[0]).any(axis=0)
        if varying.sum() != 1:
            raise ValueError(
                # TODO: ask maintainers if this is an appropriate error
                # message. It should really only be seen by developers.
                f"[DEVELOPER ERROR] Expected exactly one varying column, "
                f"but found {varying.sum()} in\n{ints2d}"
            )
        return ints2d[:, varying].ravel()

    sh_degrees: _np.ndarray = _extract_sh_degrees(
        ints_in_fnames
    )  # numpy array (1d, ints)

    # Sort filnames, hashes, and degrees by degree
    @dataclass(frozen=True)
    class FileInfo:
        fname: str
        sha256: str
        lmax: int

    arr_fileinfo: list[FileInfo] = [
        FileInfo(fname, archive.registry[fname], int(deg))
        for fname, deg in zip(fnames, sh_degrees)
    ]

    arr_fileinfo = sorted(
        arr_fileinfo,
        key=lambda x: x.lmax,
    )

    # Toss any models with lmax less than user-requested
    target_lmax: int
    if user_lmax < 0:
        target_lmax = arr_fileinfo[-1].lmax
    else:
        target_lmax = min(user_lmax, arr_fileinfo[-1].lmax)

    for i, fileinfo in enumerate(arr_fileinfo):
        if fileinfo.lmax >= target_lmax:
            arr_fileinfo = arr_fileinfo[i:]
            break

    file_to_fetch: FileInfo
    # First, check if a model with the desired lmax or greater already exists
    # in the cache.
    for fileinfo in arr_fileinfo:
        if (archive.path / fileinfo.fname).exists():
            file_to_fetch = fileinfo
            break
    # Otherwise, choose the smallest model which satisfies the desired lmax.
    else:
        file_to_fetch = arr_fileinfo[0]

    # Download the file
    fpath: str = archive.fetch(
        file_to_fetch.fname,
        downloader=downloader,
    )

    return fpath, target_lmax
