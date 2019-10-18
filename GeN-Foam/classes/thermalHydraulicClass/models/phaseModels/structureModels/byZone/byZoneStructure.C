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

\*---------------------------------------------------------------------------*/

#include "byZoneStructure.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcDiv.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace structureModels
{
    defineTypeNameAndDebug(byZone, 0);
    addToRunTimeSelectionTable
    (
        structureModel,
        byZone,
        structureModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::structureModels::byZone::byZone
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    structureModel
    ( 
        dict,
        mesh
    )
{
    forAll(dict.toc(), i)
    {
        //- Read cellZone and subDicts
        word zoneName(dict.toc()[i]);
        if (zoneName == "type") continue;
        const dictionary& zoneDict(dict.subDict(zoneName));
        const dictionary& actDict
        (
            zoneDict.subDict("powerModel")
        );
        const dictionary& pasDict
        (
            zoneDict.subDict("passiveProperties")
        );
        const labelList& zoneCellList(mesh.cellZones()[zoneName]);
        
        //- Adjust fields
        regions_.append(zoneName);
        cellLists_.insert
        (
            zoneName,
            zoneCellList
        );
        scalarField zoneCellField(mesh.cells().size(), 0.0);
        forAll(zoneCellList, i)
        {
            zoneCellField[zoneCellList[i]] = 1.0;
        }
        cellFields_.insert
        (
            zoneName,
            zoneCellField
        );
        scalar alpha(zoneDict.lookupType<scalar>("volumeFraction"));
        scalar Dh(zoneDict.lookupType<scalar>("Dh"));
        scalar Tpas(pasDict.lookupOrDefault<scalar>("T", 298.15));
        scalar iApas(pasDict.lookupOrDefault<scalar>("iA", 0.0));
        scalar rhopas(pasDict.lookupOrDefault<scalar>("rho", 1e-9));
        scalar Cppas(pasDict.lookupOrDefault<scalar>("Cp", 1e-9));
        rhopas = std::max(rhopas, 1e-9);
        Cppas = std::max(Cppas, 1e-9);
        forAll(zoneCellList, i)
        {   
            label celli(zoneCellList[i]);
            (*this)[celli] = alpha;
            Dh_[celli] = Dh;
            Tpas_[celli] = Tpas;
            iApas_[celli] = iApas;
            rhopas_[celli] = rhopas;
            Cppas_[celli] = Cppas;
        }
        this->correctBoundaryConditions();
        Dh_.correctBoundaryConditions();
        Tpas_.correctBoundaryConditions();
        iApas_.correctBoundaryConditions();
        rhopas_.correctBoundaryConditions();
        Cppas_.correctBoundaryConditions();
        volScalarField zoneAlphaField(*this);
        zoneAlphaField.primitiveFieldRef() *= zoneCellField;
        alphaFields_.insert
        (
            zoneName,
            zoneAlphaField
        );

        //- Now for the rotation matrices to move from the global to the local
        //  reference frame
        vector localX
        (
            zoneDict.lookupOrDefault<vector>("localX", vector(1,0,0))
        );
        localX /= mag(localX);
        vector localZ
        (
            zoneDict.lookupOrDefault<vector>("localZ", vector(0,0,1))
        );
        localZ /= mag(localZ);
        //- Absorb non-orthogonalities in X
        localX -= (localX&localZ)*localZ/mag(localZ);
        localX /= mag(localX);
        //- Compute third axis
        vector localY(localZ ^ localX);
        localY /= mag(localY);

        //- Construct matrices
        forAll(zoneCellList, i)
        {
            label celli(zoneCellList[i]);
            Rg2l_[celli][0] = localX[0];
            Rg2l_[celli][1] = localY[0];
            Rg2l_[celli][2] = localZ[0];
            Rg2l_[celli][3] = localX[1];
            Rg2l_[celli][4] = localY[1];
            Rg2l_[celli][5] = localZ[1];
            Rg2l_[celli][6] = localX[2];
            Rg2l_[celli][7] = localY[2];
            Rg2l_[celli][8] = localZ[2];
        }
        Rg2l_.correctBoundaryConditions();

        //- Construct lDh_ (for isotropic structures each component of
        //  lDh_ is equal to Dh cell by cell)
        vector lDhAnisotropy
        (
            zoneDict.lookupOrDefault<vector>("localDhAnisotropy", vector::one)
        );
        forAll(zoneCellList, i)
        {
            label celli(zoneCellList[i]);
            lDh_[celli][0] = lDhAnisotropy[0]*Dh_[celli];
            lDh_[celli][1] = lDhAnisotropy[1]*Dh_[celli];
            lDh_[celli][2] = lDhAnisotropy[2]*Dh_[celli];
        }
        lDh_.correctBoundaryConditions();

        //- Construct global tortuosity tensor by first constructing the local
        //  one and rotating in back to the global reference frame
        vector tortuosityV
        (
            zoneDict.lookupOrDefault<vector>("localTortuosity", vector::one)
        );
        tensor R
        (
            localX[0], localY[0], localZ[0],
            localX[1], localY[1], localZ[1],
            localX[2], localY[2], localZ[2]
        );
        tensor tortuosity(tensor::zero);
        tortuosity[0] = tortuosityV[0];
        tortuosity[4] = tortuosityV[1];
        tortuosity[8] = tortuosityV[2];
        forAll(zoneCellList, i)
        {
            label celli(zoneCellList[i]);
            tortuosity_[celli] = 
                R.T() & tortuosity & R;
        }
        tortuosity_.correctBoundaryConditions();

        //- Construct powerModel
        powerModels_.insert
        (
            zoneName,
            powerModel::New
            (
                *this,
                actDict,
                zoneName
            )
        );
    }

    //- Adjust iAact that was left out as it is a member of
    //  powerModel
    forAllIter
    (
        powerModelTable,
        powerModels_,
        iter
    )
    {
        iAact_ += iter()->iA();
    }
    iAact_.correctBoundaryConditions();

    //- Inverse rotation matrix to move from local to global is the transpose
    //  of the one to move from global to local
    Rl2g_ = Rg2l_.T();
    Rl2g_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// ************************************************************************* //
