#pragma once

#include <unordered_map>
#include <functional>

/**
 * @todo write docs
 */
template<class PARAMETER, class RESULT>
class FunctionTableCache
{
private:
    std::unordered_map<PARAMETER, RESULT> _cache;
    std::function<RESULT (PARAMETER)> _function;

public:
    FunctionTableCache(const auto & function)
    : _cache(), _function(function)
    {}

    ~FunctionTableCache()
    {}

    RESULT operator() (const PARAMETER parameter)
    {
        if (GetCache().contains(parameter))
        {
            return GetCache()[parameter];
        }
        else
        {
             auto result = _function(parameter);
             GetCache()[parameter] = result;
             return result;
        }
    }

    std::unordered_map<PARAMETER, RESULT> &GetCache()
    {
        return _cache;
    }
};
