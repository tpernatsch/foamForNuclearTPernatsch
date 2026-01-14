=================
How to contribute
=================

-------------
Merge request
-------------

Contributions are accepted using standard GitLab procedures.

---------
Tutorials
---------

Tutorials should always include:

- ``Allrun`` or ``Allrun_parallel``
- ``Allclean``
- ``README.md``
- ``0.orig`` or ``0`` time folder without additional time steps unless
  necessary (e.g for restart)
- Optional: ``Alltest``

-----------------
C++ documentation
-----------------

All contributed code should follow the guidelines below.

For code styling refer to `OpenFOAM styling <https://openfoam.org/dev/coding-style-guide/>`_.


.H header
=========

The header of the ``.H`` file requires after the license:

- the name of the class
- a minimum description
- the name of the authors of the class
- the source files


.. code :: cpp

  <code banner with copyright>

  License
      <license content>

  Class
      Foam::<className>

  Description
      <minimal description of the class>

  \mainauthor
      <author1 name> <author1 email>, <institution1 name> (<country1>);
      <author2 name> <author2 email>, <institution2 name> (<country2>)

  SourceFiles
      <className>.C


**Important notes:**

- Be careful with the use of ``//-`` in the code. This only refer to doxygen
  inline documentation. Do not use for regular comment in the code.
- Include comments for all members and member functions.


.C header
=========

The header of the ``.C`` only requires the banner and the license. It doesn't
require description or duplicated information from the ``.H`` file.

.. code :: cpp

    /*--------------------------------------------------------------------------*\
    |       ______ ______ _   __    |                                            |
    |      / ____// ____// | / /    | foamForNuclear                             |
    |     / /_   / /_   /  |/ /     | Website: https://gitlab.com/foamForNuclear |
    |    / __/  / __/  / /|  /      |                                            |
    |   /_/    /_/    /_/ |_/       |                                            |
    |                                                                            |
    |  Built on OpenFOAM v2506                                                   |
    |  Copyright 2011-2016 OpenFOAM Foundation, 2017-2025 OpenCFD Ltd.           |
    ------------------------------------------------------------------------------
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

    \*---------------------------------------------------------------------------*/



YAML Doc file
=============

Provide a ``.yaml`` file for each model. A new merge request with a new class
without a ``.yaml`` will be denied.

We provide below a template for this ``.yaml`` file.

**Important notes:**

- The indentation matters


.. code :: yaml

    description: |
      <add your description>
      <over multiple lines>

      <add a inline formula $R = \Sigma \phi$>


      <add a formula as a block>

      $$
      R = \Sigma \phi
      $$


      <add code block use the verbatim keyword>

      \verbatim
      someDict
      {
          keyA        varA;
          keyB        varB;
      }
      \endverbatim


      <add a table, can be used for internal variables of interest in a formula>

      \table
          GeN-Foam Parameter  | Param    | Req'd | Default | Notation                 | Description
          fuelFraction        | Constant | No    | 1       | $\alpha_{fuel}$          | Volume of fuel per lattice volume
          IV                  | Constant | Yes   | -       | $1/v_g$                  | Inverse velocity
          D                   | Variable | Yes   | -       | $D_g$                    | Diffusion coefficient
          nuSigmaEff          | Variable | Yes   | -       | $\nu \Sigma_{f,eff, g}$  | Fission neutrons cross-section
      \endtable

    admonitions:
      - kind: warning
        body: |
          <example of warning>
          <over multiple lines>
      - kind: note
        body: |
          <example of note>
          <over multiple lines>

    options:
      - key: key1
        type: scalar
        required: false
        default: 1.0
        description: <description of key1>
      - key: key2
        type: int
        required: false
        default: 1
        description: <description of key1>

    usage:
      - comment: |
          <comment before usage1>
        snippet: |
          feature1
          {
              key1  2.0;
              key2  1;
          }

      - comment: |
          <comment before usage2>
        snippet: |
          feature1
          {
              key1  2.0;
          }