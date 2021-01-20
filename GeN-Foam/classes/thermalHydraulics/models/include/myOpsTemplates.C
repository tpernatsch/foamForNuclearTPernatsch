namespace Foam
{
namespace myOps
{
    template<class Type>
    List<Type> split(const Type& compoundItem, char d)
    {
        List<Type> items(1, "");
        int i = 0;
        for (char const &c : compoundItem)
        {
            if (c == d)
            {
                i += 1;
                items.append("");
            }
            else items[i] += c;
        }
        return items;
    }

    template<class Type>
    List<Type> split(const Type& compoundItem, List<char> ds)
    {
        List<Type> items(1, "");
        int i = 0;
        int j = 0;
        int compoundItemSize(0);
        for (char const &c : compoundItem)
        {
            (void) c;
            compoundItemSize += 1;
        }
        for (char const &c : compoundItem)
        {
            bool addChar = true;
            for (char const &d : ds)
            {
                if (c == d)
                {
                    if 
                    (
                        items[items.size()-1] != ""
                    and j != compoundItemSize-1)
                    {
                        i += 1;
                        items.append("");
                    }
                    addChar = false;
                    break;
                }
            }
            j += 1;
            if (addChar)
                items[i] += c;
        }
        return items;
    }
}
}
