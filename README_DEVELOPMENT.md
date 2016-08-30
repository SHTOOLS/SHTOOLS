# Development Instructions

### Participation
You are welcome to participate and report issues on
[GITHUB](www.github.com/SHTOOLS/SHTOOLS).

### Versioning

SHTOOLS uses the following convention for version numbers:
* the version number is a *float* (e.g. 3.4) and specified in the file
  VERSION in the root directory of SHTOOLS
* a release has to be tagged with the VERSION number in git
* ISRELEASED in the file setup.py has to be *True* if the version indicated in
  VERSION has already been released. In this case (if the version is not tagged in git)
  .post0+commit is appended to the version number (e.g. pyshtools-3.3.post0+f9fc65e).
  If the version that is indicated in the VERSION file is not yet released .dev0+commit is appended (e.g. pyshtools-3.4.dev0+f9fc65e)
