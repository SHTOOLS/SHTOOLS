# Development Instructions

### Participation
You are welcome to participate and report issues on
[GitHub](www.github.com/SHTOOLS/SHTOOLS).

If you would like to contribute to the development process, please follow these guidelines:

* The SHTOOLS repository consists primarily of the branches `master` and `develop`. `master` is used only for tagged releases, whereas `develop` is the main branch for development.

* All changes to the code base are made by *pull requests*. Start by forking the SHTOOLS repo into your GitHub account, and then clone that repo onto your local computer. Commit changes to your personal repo, and when you think that your changes are ready to be incorporated into the main code, make a pull request from your account on GitHub.

* All pull requests will be subjected to discussion. During this phase, you may be asked to make additional changes. When the lead developers are happy with the modifications, they will merge them into the `develop` branch.

### Versioning

SHTOOLS uses the following convention for version numbers:

* The version number is a *float* (e.g. 3.4) and specified in the file
  VERSION in the root directory of SHTOOLS.
* A release has to be tagged with the VERSION number in git
* ISRELEASED in the file setup.py has to be *True* if the version indicated in
  VERSION has already been released. In this case (if the version is not tagged in git)
  .post0+commit is appended to the version number (e.g. pyshtools-3.3.post0+f9fc65e).
  If the version that is indicated in the VERSION file is not yet released .dev0+commit is appended (e.g. pyshtools-3.4.dev0+f9fc65e).
