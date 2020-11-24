#include "myOps.H"

namespace Foam
{
namespace myOps
{
    bool relax(const fvMesh& mesh, word name)
    {
        bool finalIter
        (
            mesh.data::template getOrDefault<bool>
            (
                "finalIteration",
                false
            ) 
        );
        if (finalIter)
            name += "Final";

        return mesh.relaxField(name);
    }

    scalar relaxationFactor(const fvMesh& mesh, word name)
    {
        bool finalIter
        (
            mesh.data::template getOrDefault<bool>
            (
                "finalIteration",
                false
            ) 
        );
        if (finalIter)
            name += "Final";

        return 
            (mesh.relaxField(name)) ?
            mesh.fieldRelaxationFactor(name) :
            1.0;
    }
}
}