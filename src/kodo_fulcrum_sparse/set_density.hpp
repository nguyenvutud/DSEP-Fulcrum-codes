// Copyright Steinwurf ApS 2011.
// Distributed under the "STEINWURF RESEARCH LICENSE 1.0".
// See accompanying file LICENSE.rst or
// http://www.steinwurf.com/licensing

#pragma once

#include <cassert>
#include <memory>

#include "has_set_density.hpp"

namespace kodo_fulcrum
{
/// @ingroup generic_api
/// @copydoc set_systematic_off(T&)
template
<
    class T,
    typename std::enable_if<has_set_density<T>::value, uint8_t>::type = 0
>
inline void set_density(T& t, double density)
{
    t.set_density(density);
}

/// @ingroup generic_api
/// @copydoc set_systematic_off(T&)
template
<
    class T,
    typename std::enable_if<
        !has_set_density<T>::value, uint8_t>::type = 0
>
inline void set_density(T& t, double density)
{
    (void) t;
    (void) density;

    // We do the assert here - to make sure that this call is not
    // silently ignored in cases where the stack does not have the
    // set_systematic_off() function. However, this assert can
    // be avoided by using has_set_systematic_off.
    assert(0 && "Missing function: set_density");
}
}
