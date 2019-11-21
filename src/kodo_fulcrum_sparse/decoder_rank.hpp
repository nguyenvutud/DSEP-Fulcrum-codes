// Copyright Steinwurf ApS 2011.
// Distributed under the "STEINWURF RESEARCH LICENSE 1.0".
// See accompanying file LICENSE.rst or
// http://www.steinwurf.com/licensing

#pragma once

#include <cassert>
#include <cstdint>
#include <type_traits>

#include "has_decoder_rank.hpp"

namespace kodo_fulcrum
{
/// @ingroup generic_api
/// @copydoc write_feedback(const T&)
template
<
    class T,
    typename std::enable_if<has_decoder_rank<T>::value, uint8_t>::type=0
>
inline uint32_t decoder_rank(const T& t)
{
    return t.decoder_rank();
}

/// @ingroup generic_api
/// @copydoc write_feedback(const T&)
template
<
    class T,
    typename std::enable_if<!has_decoder_rank<T>::value, uint8_t>::type=0
>
inline uint32_t decoder_rank(const T& t)
{
    (void) t;


    // We do the assert here - to make sure that this call is not
    // silently ignored in cases where the stack does not have the
    // write_feedback() function. However, this assert can
    // be avoided by using the has_write_feedback.hpp functionality
    assert(0  && "Missing function: decoder_rank");

    return 0;
}
}
