// Copyright Steinwurf ApS 2014.
// Distributed under the "STEINWURF RESEARCH LICENSE 1.0".
// See accompanying file LICENSE.rst or
// http://www.steinwurf.com/licensing

#pragma once

#include <kodo_core/coefficient_info.hpp>
#include <kodo_core/final_layer.hpp>
#include <kodo_core/finite_field_layers.hpp>
#include <kodo_core/finite_field_counter.hpp>
#include <kodo_core/payload_info.hpp>
#include <kodo_core/plain_symbol_id_reader.hpp>
#include <kodo_core/plain_symbol_id_size.hpp>
#include <kodo_core/pool_factory.hpp>
#include <kodo_core/proxy_layer.hpp>
#include <kodo_core/storage_block_info.hpp>
#include <kodo_core/symbol_id_decoder.hpp>
#include <kodo_core/systematic_decoder_layers.hpp>

namespace kodo_fulcrum
{
/// @ingroup utility fulcrum
///
/// @brief Decodes the payload portion of the inner code in a
///        fulcrum decoder.
///
/// When used this layer will consume layer::read_payload(uint8_t*)
/// calls and forward them to layer::read_symbol(uint8_t*,
/// uint8_t*) or layer::read_uncoded_symbol(uint8_t*, uint32_t) in the
/// main stack.
///
/// This is implemented as a secondary stack since the inner and outer
/// codes do typically not use the same finite field. Having
/// fulcrum_payload_decoder nested means that the user will only see
/// the outer code. Which is typically what we want.
///
/// @tparam MainStack The type of the "main stack" where calls not
///         implemented in this stack will be forwarded curtecy of
///         the proxy_layer.
///
template<class MainStack>
class fulcrum_payload_decoder_counter_high_field : public
    // Payload API
    kodo_core::payload_info<
    // Codec Header API
    kodo_core::systematic_decoder_layers<
    kodo_core::symbol_id_decoder<
    // Symbol ID API
    kodo_core::plain_symbol_id_reader<
    kodo_core::plain_symbol_id_size<
    // Coefficient Storage API
    kodo_core::coefficient_info<
    // Storage API
    kodo_core::storage_block_info<
    // Finite Field API
    kodo_core::finite_field_counter<
    kodo_core::finite_field_layers<fifi::binary4,
    // Proxy
    kodo_core::proxy_layer<MainStack,
    // Final layer
    kodo_core::final_layer
    >>>>>>>>>>
{
public:
    using factory = kodo_core::pool_factory<fulcrum_payload_decoder_counter_high_field>;
};
}
