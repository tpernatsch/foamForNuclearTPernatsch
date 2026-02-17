# -*- coding: utf-8 -*-
"""
foamForNuclear:

    /*--------------------------------------------------------------------------*\
    |       ______ ______ _   __    |                                            |
    |      / ____// ____// | / /    | foamForNuclear - Python API                |
    |     / /_   / /_   /  |/ /     | Website: https://gitlab.com/foamForNuclear |
    |    / __/  / __/  / /|  /      |                                            |
    |   /_/    /_/    /_/ |_/       |                                            |
    |                                                                            |
    |  Built on OpenFOAM v2512                                                   |
    |  Copyright 2011-2016 OpenFOAM Foundation, 2017-2025 OpenCFD Ltd.           |
    \*--------------------------------------------------------------------------*/

License
    This file is part of foamForNuclear.

    foamForNuclear is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    foamForNuclear is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    This offering is not approved or endorsed by the OpenFOAM Foundation nor
    OpenCFD Limited, producer and distributor of the OpenFOAM(R)software via
    www.openfoam.com, and owner of the OPENFOAM(R) and OpenCFD(R) trademarks.

    This particular snippet of code is developed according to the developer's
    knowledge and experience in OpenFOAM. The users should be aware that
    there is a chance of bugs in the code, though we've thoroughly test it.
    The source code may not be in the OpenFOAM coding style, and it might not
    be making use of inheritance of classes to full extent.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
"""

from . import _attrs_tools
from . import boundaryConditions

# To be checked how to group
from foamForNuclear.checkvalue import (check_value, check_positive, check_type, CheckedList)
from . import common

from . import control

# To be checked how to group
from . import coupling
from foamForNuclear.executor import run, allclean, run_preprocessing, run_reconstruction, copyFolder
from . import executor
from foamForNuclear.externalSource import ExternalSource

from . import fields

# # To be checked what this import does
# from foamForNuclear.fmi import *

from . import functions

# To be checked how to group
from . import nuclearData

from . import case
from . import mesh
from . import numerics
from . import offbeat_lib
from . import openmcTools

# To be checked how to group
from foamForNuclear.openfoamFile import OpenFOAMFile

from . import porous_medium
from . import preprocessing
from . import profiles

# To be checked how to group
from . import quadratureSet

from . import solvers

# To be checked how to group
from . import timeFolder
from . import timeProfile

from . import thermo
from . import transport
from . import turbulence