"""
pyshtools datasets.

To load a dataset, call the relevant method as in this example:

    hlm = pysh.datasets.Venus.VenusTopo719()

When accessing a pyshtools dataset, the file will first be downloaded from the
orginal source and stored in the user's cache directory (if it had not been
done previously). The file hash will be verified to ensure that it has not been
modified, and the file will then be used to instantiate and return an SHCoeffs,
SHGravCoeffs or SHMagCoeffs class instance.

For datasets of spherical harmonic coefficients, the coefficients can be read
up to a maximum specified degree by providing the optional variable lmax. For
magnetic field data, the coefficients can be returned in Tesla or nT by
use of the variable nt.
"""
from . import Mercury
from . import Venus
from . import Earth
from . import Moon
from . import Mars
from . import Eros
from . import Vesta
from . import Ceres
from . import Io
from . import Europa
from . import Ganymede
from . import Callisto
from . import Saturn
from . import Titan
from . import Enceladus
from . import Uranus
from . import Neptune


__all__ = ['Mercury', 'Venus', 'Earth', 'Moon', 'Mars', 'Eros', 'Vesta',
           'Ceres', 'Io', 'Europa', 'Ganymede', 'Callisto', 'Saturn', 'Titan',
           'Enceladus', 'Uranus', 'Neptune']
