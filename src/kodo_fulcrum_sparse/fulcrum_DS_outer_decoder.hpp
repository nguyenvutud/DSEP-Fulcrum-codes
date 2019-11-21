// Copyright Steinwurf ApS 2014.
// Distributed under the "STEINWURF RESEARCH LICENSE 1.0".
// See accompanying file LICENSE.rst or
// http://www.steinwurf.com/licensing

#pragma once

#include <cstdint>
#include <type_traits>

#include <fifi/prime2325.hpp>

#include <kodo_core/coefficient_storage_layers.hpp>
#include <kodo_core/common_decoder_layers.hpp>
#include <kodo_core/final_layer.hpp>
#include <kodo_core/find_enable_trace.hpp>
#include <kodo_core/finite_field_layers.hpp>
#include <kodo_core/finite_field_counter.hpp>
#include <kodo_core/nested_payload_size.hpp>
#include <kodo_core/nested_read_payload.hpp>
#include <kodo_core/pool_factory.hpp>
#include <kodo_core/proxy_args.hpp>
#include <kodo_core/select_storage_layers.hpp>
#include <kodo_core/deep_storage_layers.hpp>
#include <kodo_core/trace_layer.hpp>
#include <kodo_core/uniform_generator.hpp>

#include <kodo_fulcrum/systematic_coefficient_mapper.hpp>
#include <kodo_fulcrum/fulcrum_info.hpp>
#include <kodo_fulcrum/fulcrum_outer_symbol_mapper.hpp>
#include <kodo_fulcrum/fulcrum_payload_decoder.hpp>
#include <kodo_fulcrum/fulcrum_proxy_stack.hpp>

namespace kodo_fulcrum
{
/// @ingroup fec_stacks
///
/// @brief The fulcrum outer decoder maps all incoming symbols
///        before starting the decoding process.
///
/// For a detailed description of the fulcrum codec see the
/// following paper on arxiv: http://arxiv.org/abs/1404.6620 by
/// Lucani et. al.
///
/// @tparam Field @copydoc layer_types::Field
///
template
<
    class Field,
    class Features = meta::typelist<>
>
class fulcrum_DS_outer_decoder : public
    // Payload API
    kodo_core::nested_read_payload<
    kodo_core::nested_payload_size<
    fulcrum_proxy_stack<kodo_core::proxy_args<>, fulcrum_payload_decoder_counter,
    // Decoder API
    fulcrum_outer_symbol_mapper<
    systematic_coefficient_mapper<

	///Outer Generation
    // Coefficient Generator API
	kodo_fulcrum::specific_sparse_uniform_index_generator<
//	kodo_core::sparse_uniform_index_generator<
//	kodo_core::uniform_generator<

    kodo_core::common_decoder_layers<Features,
    // Coefficient Storage API
    kodo_core::coefficient_storage_layers<
    // Storage API
    kodo_core::select_storage_layers<
    kodo_core::deep_storage_layers, Features,
    // Finite Field API
	kodo_core::finite_field_counter<
    kodo_core::finite_field_layers<Field,
    // Fulcrum API
    fulcrum_info<
    std::integral_constant<uint32_t,20>, // MaxExpansion
    std::integral_constant<uint32_t,4>,  // DefaultExpansion
    // Trace Layer
    kodo_core::trace_layer<kodo_core::find_enable_trace<Features>,
    // Final Layer
    kodo_core::final_layer
    >>>>>>>>>>>>>
{
public:

    static_assert(!std::is_same<Field, fifi::prime2325>::value,
                  "The mapping between inner and outer code requires "
                  "that both are binary extension fields");

    // Define the nested factory type
    using factory = kodo_core::pool_factory<fulcrum_DS_outer_decoder>;
};
}
