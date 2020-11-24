/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::twoPhaseDragMultiplierModels::LottesFlinn

Description
    Model in the spirit of Lottes-Flinn, i.e.: 

    phi2 = 1/(alpha^exp) 
    
    with alpha being the normalized phase fraction (by the void fraction) of
    the multiplierFluid and exp defaulting to 2 (which is the actual
    Lottes-Flinn model). The exp can be user-provided (it is know from e.g.
    TRACE development that this exponent can be slightly smaller for bubbly
    flows, e.g. in the 1.7-1.8 range). The optional mode keyword is 
    used to specify whether the actual phase-fraction of the vapour should
    be used of if alpha should be computed via the correlation by Nguyen,
    H. Nguyen "A finite element method for one-dimensional sodium boiling
    simulation" in Third International Topical Meeting on Reactor 
    Thermal-hydraulics, p.17C

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    LottesFlinn.C

DictEntryExample

    twoPhaseDragMultiplierModel
    {
        type            LottesFlinn;
        multiplierFluid liquid;
        exp             2;         //- Not required, defaults to 2  
        mode            alpha;     //- Not required, defaults to alpha
    }


\*---------------------------------------------------------------------------*/

#ifndef LottesFlinn_H
#define LottesFlinn_H

#include "twoPhaseDragMultiplierModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace twoPhaseDragMultiplierModels
{

/*---------------------------------------------------------------------------*\
                           Class LottesFlinn Declaration
\*---------------------------------------------------------------------------*/

class LottesFlinn
:
    public twoPhaseDragMultiplierModel
{
public:

    enum mode
    {
        alpha,
        XLM
    };

    static const Enum<mode> modeNames_;

private:

    //- exp in phi2 = 1/(alpha^exp)
    scalar exp_;

    //- If set to XLM, the Nguyen correlation for alpha will be used instead of
    //  the phase fraction field available from the solver
    const mode mode_;

public:

    //- Runtime type information
    TypeName("LottesFlinn");

    //- Constructors

        LottesFlinn
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const fvMesh& mesh
        );

    //- Destructor
    virtual ~LottesFlinn();

    //- Member Functions

        //- Calculate phi2_;
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace twoPhaseDragMultiplierModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //