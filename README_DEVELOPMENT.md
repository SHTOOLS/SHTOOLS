# Development Instructions

### Participation
You are welcome to participate and report issues on
[GITHUB](www.github.com/SHTOOLS/SHTOOLS).

### Versioning

SHTOOLS uses the folling convention for versioning:
* the version number is a *float* (e.g. 3.4) and specified in the file
  VERSION in the root directory of SHTOOLS
* a release has to be tagged with the VERSION number in git
* if the version that is indicated in the VERSION file has already been
  released, the variable ISRELEASED in the file setup.py has to be *True*.
  In this case, if the version is not tagged in git, .post0+commit is appended.
  If the version that is indicated in the VERSION file is not yet released,
  the variable ISRELEASED in the file setup.py has to be *False*.
  In this case, if the version is not tagged in git, .dev0+commit is appended.
  For exemple: 3.4.dev0+392123 indicates that it is a development version that
  will be released in the future as version 3.4 . 3.4.post0+392123 indicates
  that this version contains changes to the already released version 3.4 .
