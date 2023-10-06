/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2306                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2022 OpenCFD Ltd.         |
-------------------------------------------------------------------------------
License
    This file is part of GeN-Foam.

    GeN-Foam is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    GeN-Foam is distributed in the hope that it will be useful, but WITHOUT
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
Class
    Foam::CHFModel::lookUpTableCHF
    
Description
    Compute the Critical Heat Flux with a Look Up Table.

Author
    - Yuhang Niu (niu.yuhang@epfl.ch)
    - A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE 
    DE LAUSANNE, Switzerland, Laboratory for Reactor Physics and Systems 
    Behaviour)
    - C. Fiorina - Texas A&M University (TAMU)


SourceFiles
    lookUpTableCHF.C

\*---------------------------------------------------------------------------*/

#ifndef lookUpTableCHF_H
#define lookUpTableCHF_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CHFModel.H"
#include "InterpolateTables.H"
#include "InterpolateTable.H"
#include "scalarFieldFieldINew.H"
#include "PtrListScalarFieldFieldINew.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class FFPair;

namespace criticalHeatFluxModels
{

/*---------------------------------------------------------------------------*\
            Class lookUpTableCHF Declaration
\*---------------------------------------------------------------------------*/

class lookUpTableCHF
:
    public CHFModel
{
protected:

    //- Ref to flow quality
    const volScalarField& mass_flow_quality_;

    //- Ref to pressure 
    const volScalarField& p_;

    //- Ref to 
    const volScalarField& rhoL_;

    //- Ref to rho vapour
    const volScalarField& rhoV_;

    //- Ref to liquid velocity 
    const volScalarField& uL_;

    //- Ref to vapour velocity 
    const volScalarField& uV_;

    //- Ref to alpha normalized liquid 
    const volScalarField& normalizedL_;

    //- Ptr to FFPair (breaking encapsulation as usual)
    mutable const FFPair* FFPairPtr_;

    //- Ref to Hydraulic diameter
    const volScalarField& Dh_;

    // -pitch-to-diameter ratio
    scalar PitchToDiameter_;

    //- Tabulated pressure values (direction z)
    scalarField pressureValues_;
    
    //- Tabulated mass flow rate values (direction x)
    scalarField massFlowRateValues_;
    
    //- Tabulated quality values (direction y)
    scalarField qualityValues_;
    
    //- Tabulated ptrList of 2D tables (i.e. raw data for the 3D interpolation
    // table)
    PtrList<FieldField<Field, scalar>> data_;

    //- Pressure interpolation method
    interpolateTableBase::interpolationMethod pMethod;
    
    //- Mass flow rate interpolation method
    interpolateTableBase::interpolationMethod gMethod_;
    
    //- Quality interpolation method
    interpolateTableBase::interpolationMethod xeMethod_;
    
    //- 3D Interpolation table (pressure as the main direction)
    scalarFieldFieldInterpolateTable pTable_;


    
public:

    //- Runtime type information
    TypeName("lookUpTableCHF");

    // Constructors

        lookUpTableCHF
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~lookUpTableCHF(){}

    // Member Functions

        //- Return CHF
        virtual scalar value
        (
            const label& celli
        ) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace criticalHeatFluxModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
