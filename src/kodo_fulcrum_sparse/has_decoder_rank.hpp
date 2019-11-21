// Copyright Steinwurf ApS 2011.
// Distributed under the "STEINWURF RESEARCH LICENSE 1.0".
// See accompanying file LICENSE.rst or
// http://www.steinwurf.com/licensing

#pragma once

#include <type_traits>
#include <utility>

namespace kodo_fulcrum
{
/// @ingroup type_traits
///
/// Type trait helper allows compile time detection of whether a
/// codec contains a layer with the member function
/// set_symbol_size(uint32_t)
///
/// Example:
///
/// using encoder_t kodo_core::full_rlnc8_encoder;
///
/// if (kodo_core::has_set_symbol_size<encoder_t>::value)
/// {
///     // Do something here
/// }
///
template<typename T>
struct has_decoder_rank
{
private:
    using yes = std::true_type;
    using no = std::false_type;

    template<typename U>
    static auto test(int) ->
        decltype(std::declval<U>().decoder_rank(0), yes());

    template<typename> static no test(...);

public:

    static const bool value = std::is_same<decltype(test<T>(0)),yes>::value;
};
}
