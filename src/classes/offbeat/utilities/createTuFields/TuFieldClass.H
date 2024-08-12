// I created this class because I need to build an hash table of TuField objects
// That have a correspoding OFFBEAT field name, the dimension of the field and
// The multiplier to convert the Tu units to the OFFBEAT units for a specific 
// field.
// Using an hash table is convenient because it is very easy to query to look 
// for the specific info relative to a specific field. 

class TuField
{
public:
    
    // Public data

        //- OFFBEAT name
        word nameOFFBEAT_;
        
        //- Dimensions
        dimensionSet dimension_;
        
        //- Conversion multiplier for unit of measure
        scalar unitMultiplier_;

        // - Disallow default bitwise copy construct
        // TuField(const TuField&);

        //- Disallow default bitwise assignment
        // void operator=(const TuField&);

    // Constructors

        //- Construct from tabulated values
        TuField
        (
            word nameOFFBEAT,
            dimensionSet dimension,
            scalar unitMultiplier
        )
        :
            nameOFFBEAT_(nameOFFBEAT),
            dimension_(dimension),
            unitMultiplier_(unitMultiplier)
        {}


    // Destructor
    ~TuField() {}


    
    // Member Operators 
};
