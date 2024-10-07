/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "userParameters.H"
#include "IFstream.H"
#include "demandDrivenData.H"
#include "messageStream.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

HashTable<dimensionedScalar, word>* userParametersRegistryPtr_(nullptr);
HashTable<dimensionedScalar, word>* userParametersPtr_(nullptr);


HashTable<dimensionedScalar, word>& userParametersRegistry()
{
    if (!userParametersRegistryPtr_)
    {
        userParametersRegistryPtr_ = new HashTable<dimensionedScalar, word>();
    }

    return *userParametersRegistryPtr_;
}


void readFromDict(const dictionary& groupDict, const word& group)
{
    HashTable<dimensionedScalar, word>& table = userParameters();
                
    wordList names(groupDict.sortedToc());
    
    forAll(names, j)
    {
        word name = names[j];
        word key = group + "::" + names[j];
        
        if(!groupDict.isDict(name))
        {
            dimensionedScalar value(name, groupDict.lookup(name));
            table.insert(key, value);
        }
        else
        {
            const dictionary& groupSubDict = groupDict.subDict(name);

            readFromDict(groupSubDict, key);
        }
    }
}


HashTable<dimensionedScalar, word>& userParameters()
{
    if (!userParametersPtr_)
    {
        userParametersPtr_ = new HashTable<dimensionedScalar, word>();
        
        IFstream ifs("system/userParameters");
        
        if (ifs.good())
        {
            dictionary dict(ifs);
            
            wordList groups(dict.sortedToc());
            
            forAll(groups, i)
            {
                word group = groups[i];
                
                const dictionary& groupDict = dict.subDict(group);

                readFromDict(groupDict, group);
            }
            
            Info << endl;
        }
    }

    return *userParametersPtr_;
}


void registerUserParameter
(
    const char* group,
    const char* name,
    const dimensionedScalar& defaultValue
)
{
    HashTable<dimensionedScalar, word>& table = userParametersRegistry();
    word key = word(group) + "::" + word(name);

    if (table.found(key))
    {
        FatalErrorInFunction << "User parameter \"" << key 
            << "\" already registered." << abort(FatalError);
        
        return;
    }
    else
    {
        table.insert(key, defaultValue);
    }
}


dimensionedScalar userParameter
(
    const char* group,
    const char* name,
    const dimensionSet& dimensions,
    const scalar defaultValue
)
{
    word key = word(group) + "::" + word(name);
    dimensionedScalar value(key, dimensions, defaultValue);
    
    registerUserParameter(group, name, value);
    
    HashTable<dimensionedScalar, word>& table = userParameters();

    if (table.found(key))
    {
        return table[key];
    }
    else
    {
        return value;
    }
}


dimensionedScalar userParameter
(
    const char* group,
    const char* name
)
{
    HashTable<dimensionedScalar, word>& table = userParameters();
    word key = word(group) + "::" + word(name);

    if (table.found(key))
    {
        return table[key];
    }
    else
    {
        return registeredUserParameter(group, name);
    }
}


bool isInUserParameters
(
    const char* group,
    const char* name
)
{
    HashTable<dimensionedScalar, word>& table = userParameters();
    word key = word(group) + "::" + word(name);

    if (table.found(key))
    {
        return true;
    }
    else
    {
        return false;
    }
}


dimensionedScalar registeredUserParameter
(
    const char* group,
    const char* name
)
{
    const HashTable<dimensionedScalar, word>& table = userParametersRegistry();
    word key = word(group) + "::" + word(name);

    if (!table.found(key))
    {
        FatalErrorInFunction << "User parameter \"" << key
            << "\" has not been registered. This is likely a developer error."
            << abort(FatalError);
    }
    
    return table[key];
}


void listUserParameters()
{
    Info<< "List of user-supplied parameters:" << endl;
    
    const HashTable<dimensionedScalar, word>& table = userParameters();
    
    wordList keys(table.sortedToc());
    
    forAll(keys, i)
    {
        word key = keys[i];
        Info<< tab << key << tab << table[key] << endl;
    }
    
    Info << endl;
}


void listRegisteredUserParameters()
{
    Info<< "List of all registered (default) user parameters:" << endl;
    
    const HashTable<dimensionedScalar, word>& table = userParametersRegistry();
    
    wordList keys(table.sortedToc());
    
    forAll(keys, i)
    {
        word key = keys[i];
        Info<< tab << key << tab << table[key] << endl;
    }
    
    Info << endl;
}


void checkUserParameters()
{
    const HashTable<dimensionedScalar, word>& registry = userParametersRegistry();
    const HashTable<dimensionedScalar, word>& table = userParameters();
    
    wordList keys(table.sortedToc());
    
    forAll(keys, i)
    {
        const word key = keys[i];
        
        if (!registry.found(key))
        {
            Ostream& os = FatalErrorInFunction;
            os  << "User parameter \"" << key
                << "\" has not been registered. " << nl
                << "The following user parameters are currently registered:"
                << nl;
            
            wordList keys(registry.sortedToc());
            
            forAll(keys, i)
            {
                word key = keys[i];
                os<< tab << key << tab << registry[key] << nl;
            }
            
            os << abort(FatalError);
        }
    }
}


// To ensure the pointer are deleted at the end of the run
class deleteUserParametersPtr
{
public:

    deleteUserParametersPtr()
    {}

    ~deleteUserParametersPtr()
    {
        deleteDemandDrivenData(userParametersPtr_);
        deleteDemandDrivenData(userParametersRegistryPtr_);
    }
};

deleteUserParametersPtr deleteUserParametersPtr_;

}

// ************************************************************************* //
