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
- [ ] `docs/pages/mydoc/fortran-release-notes-v4.md` : update release notes
- [ ] `requirements.txt` : update version numbers of the python dependencies, if necessary
- [ ] `requirements-dev.txt` : update version numbers of the python developer dependencies, if necessary
- [ ] `environment.yml` : update version numbers of the conda environment, if necesseary
- [ ] `binder/environment.yml` : update version number of pyshtools and other dependencies
- [ ] `fpm.toml` : update shtools version number

Update the documentation files and man pages
- [ ] `cd docs; bundle update; cd ..` : update the Gemfile for the jekyll web documentation
- [ ] `make remove-doc` : this ensures that the correct version number will be written to the fortran man pages
- [ ] `make doc` : make the fortran man pages, create markdown files from the python docstrings, and create web documentation

### Release ###
- [ ] Commit all changes to the `develop` branch and then merge all changes to the `master` branch.
- [ ] Go to [GitHub Release](https://github.com/SHTOOLS/SHTOOLS/releases), create a tag of the form `vX.X`, and draft a new release. After this is done, a zipped archive will be sent to [Zenodo](https://doi.org/10.5281/zenodo.592762), which will create a doi for citation.
- [ ] Update the master branch on your personal repo, along with the newly created tag, using `git pull shtools master --tags`, where shtools is the name of the remote repo on github.

### Update pypi ###
For the next steps to work, the file ```.pypirc``` with the username and password needs to be set (see [this link](https://packaging.python.org/guides/migrating-to-pypi-org/#uploading)). ```pandoc``` needs to be installed with either ```conda install -c conda-forge pypandoc``` or ```pip install pypandoc```, and ```twine``` neeeds to be installed by ```pip install twine```.
- [ ] A pypi upload can only be done once for a given version. It is therefore recommended to test it first on pypitest.
    ```
    python3 setup.py sdist
    gpg --detach-sign -a dist/pyshtools-x.x.tar.gz
    twine upload dist/* -r pypitest
    ```
- [ ] Inspect the pypi page at https://test.pypi.org/project/pyshtools/x.x/ for errors in the project description or metadata. Then install pyshtools in a directory that is different from your local pyshtools repo, and run a couple of quick tests.
    ```
    pip3 uninstall pyshtools
    pip3 install -i https://test.pypi.org/simple pyshtools --no-binary pyshtools
    pip3 uninstall pyshtools
    ```
- [ ] Upload to pypi:
    ```
    python3 setup.py sdist
    gpg --detach-sign -a dist/pyshtools-x.x.tar.gz
    twine upload dist/* -r pypi
    pip3 install pyshtools  # check that it works!
    ```
    As a sanity check, inspect the pypi project page at https://pypi.org/project/pyshtools/.

### Build wheels ###
- [ ] Build wheels for linux, macOS and windows and upload to pypi. This is done using multibuild in combination with Travis and Appveyor.
    ```
    git clone https://github.com/shtools/build-shtools.git # only necessary the first time.
    cd build-shtools
    git submodule update --remote  # add the option --init the first time you use this command.
    git commit -a -m "Update shtools and multibuild to master"
    git push
    ```
- [ ] Inspect the pypi project page at https://pypi.org/project/pyshtools/ to ensure that the wheels were uploaded.

### Update Homebrew ###
Update the homebrew installation by editing the file `shtools.rb` in the homebrew-shtools repo.
- [ ] Change "url" to point to the new version (the link to the tar.gz archive can be found on the release page).
- [ ] Update the sha256 hash of the tar.gz pypi upload (either from the pypi files, or by using `shasum -a 256 filename`).
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
