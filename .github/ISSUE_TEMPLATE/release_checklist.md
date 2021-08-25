---
name: SHTOOLS release checklist
about: Checklist for a new SHTOOLS release.
title: 'Release SHTOOLS x.x.x'
labels: 'maintenance'
assignees: ''

---

**Scheduled Date:** YYYY/MM/DD

### Before release ###
Make all changes on the branch `develop`. Verify that the version numbers and other metadata are up to date in the following files:
- [ ] `Makefile` : update shtools version number (for use with fortran man page documentation only; not required for maintenance releases x.x.>0)
- [ ] `docs/_data/sidebars/mydoc_sidebar.yml` : update pyshtools version number for web documentation
- [ ] `docs/_data/sidebars/fortran_sidebar.yml` : update shtools version number for web documentation
- [ ] `docs/pages/mydoc/release-notes-v4.md` : update release notes
- [ ] `docs/pages/fortran/fortran-release-notes-v4.md` : update release notes
- [ ] `requirements.txt` : update version numbers of the python dependencies, if necessary
- [ ] `requirements-dev.txt` : update version numbers of the python developer dependencies, if necessary
- [ ] `environment.yml` : update version numbers of the conda environment, if necesseary
- [ ] `binder/environment.yml` : update version number of pyshtools and other dependencies
- [ ] `fpm.toml` : update shtools version number
- [ ] `AUTHORS.md`, `docs/pages/mydoc/contributors.md`, `docs/pages/fortran/fortran-contributors.md` : Add new contributors, if necessary

Update the documentation files and man pages
- [ ] `cd docs; bundle update; cd ..` : update the Gemfile for the jekyll web documentation
- [ ] `make remove-doc` : this ensures that the correct version number will be written to the fortran man pages
- [ ] `make doc` : make the fortran man pages, create markdown files from the python docstrings, and create web documentation

### Release ###
- [ ] Commit all changes to the `develop` branch and then merge all changes to the `master` branch.
- [ ] Go to [GitHub Release](https://github.com/SHTOOLS/SHTOOLS/releases), create a tag of the form `vX.X`, and draft a new release.
- [ ] Update the master branch on your personal repo, along with the newly created tag, using `git pull shtools master --tags`, where shtools is the name of the remote repo on github.

### Verify workflow execution ###
- [ ] Creation and upload of a zipped archive to [Zenodo](https://doi.org/10.5281/zenodo.592762)
- [ ] Upload of a source tarball as a release asset.
- [ ] Upload of the source distribution to [pypi](https://pypi.org/project/pyshtools/).
- [ ] Creation and upload of macOS and Linux binary wheels to [pypi](https://pypi.org/project/pyshtools/).
- [ ] Manually trigger Appveyor on the [build-shtools](https://github.com/SHTOOLS/build-shtools) repo to create and upload a windows wheel.

### Update Homebrew ###
Update the homebrew installation by editing the file `shtools.rb` in the [homebrew-core](https://github.com/Homebrew/homebrew-core) repo.
- [ ] Change "url" to point to the new version (the link to the tar.gz archive can be found on the release page).
- [ ] Update the sha256 hash of the tar.gz pypi upload by using `shasum -a 256 filename`.
- [ ] Commit and push changes.

### Update MacPorts ###
Update the MacPorts installation by editing the file `Portfile` in the [macports-ports](https://github.com/macports/macports-ports) repo.
- [ ] Change the version number in `github.setup`.
- [ ] Update the sha256, rmd160 and file size of the release asset.
- [ ] Commit and push changes.

### Update conda-forge ###
- [ ] Go to [pyshtools-feedstock](https://github.com/conda-forge/pyshtools-feedstock) and check that the `recipe/meta.yaml` file has been updated. This is usually done automatically by conda-forge's bot. If necessary, merge changes into the feedstock.

### Update web documentation ###
Update the web documentation at shtools.oca.eu.
- [ ] Push the branch `master` to the private repo at gitlab.oca.eu.
- [ ] Execute `gitlab-runner run` locally to build and upload the web documentation.

### Post release ###
- [ ] Build a new Binder image by running one of the tutorials.
- [ ] Advertise on Twitter and Element/Matrix.
