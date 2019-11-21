// Copyright Steinwurf ApS 2011.
// Distributed under the "STEINWURF RESEARCH LICENSE 1.0".
// See accompanying file LICENSE.rst or
// http://www.steinwurf.com/licensing

#pragma once

#include <kodo_core/basic_proxy_stack.hpp>
#include <kodo_core/block_decoder.hpp>
#include <kodo_core/coefficient_storage_layers.hpp>
#include <kodo_core/common_decoder_layers.hpp>
#include <kodo_core/final_layer.hpp>
#include <kodo_core/find_enable_trace.hpp>
#include <kodo_core/finite_field_layers.hpp>
#include <kodo_core/finite_field_counter.hpp>
#include <kodo_core/nested_write_payload.hpp>
#include <kodo_core/payload_info.hpp>
#include <kodo_core/plain_symbol_id_reader.hpp>
#include <kodo_core/plain_symbol_id_size.hpp>
#include <kodo_core/select_storage_layers.hpp>
#include <kodo_core/deep_storage_layers.hpp>
#include <kodo_core/symbol_id_decoder.hpp>
#include <kodo_core/systematic_decoder_layers.hpp>
#include <kodo_core/trace_layer.hpp>
#include <kodo_core/select_symbol_id_reader_layers.hpp>
#include <kodo_core/plain_symbol_id_reader_layers.hpp>

#include <kodo_rlnc/full_vector_recoding_stack.hpp>

namespace kodo_rlnc
{
    /// @ingroup fec_stacks
    ///
    /// @brief Implementation of a complete RLNC decoder
    ///
    /// This configuration adds the following features (including those
    /// described for the encoder):
    /// - Recoding using the recoding_stack
    /// - Linear block decoder using Gauss-Jordan elimination.
    template
    <
        class Field,
        class Features = meta::typelist<>
    >
    class full_vector_decoder_counter : public
        // Payload API
        kodo_core::nested_write_payload<
        kodo_core::basic_proxy_stack<
            kodo_core::proxy_args<Features>, full_vector_recoding_stack,
        kodo_core::payload_info<
        // Block Coder API
        kodo_core::block_decoder<
        // Codec Header API
        kodo_core::systematic_decoder_layers<
        kodo_core::symbol_id_decoder<
        // Symbol ID API
        kodo_core::select_symbol_id_reader_layers<
            kodo_core::plain_symbol_id_reader_layers, Features,
        // Decoder API
        kodo_core::common_decoder_layers<Features,
        // Coefficient Storage API
        kodo_core::coefficient_storage_layers<
        // Storage API
        kodo_core::select_storage_layers<
            kodo_core::deep_storage_layers, Features,
        // Finite Field API
		kodo_core::finite_field_counter<
        kodo_core::finite_field_layers<Field,
        // Trace Layer
        kodo_core::trace_layer<kodo_core::find_enable_trace<Features>,
        // Final Layer
        kodo_core::final_layer
        > > > > > > > > > > > > >
    {
    public:
        using factory = kodo_core::pool_factory<full_vector_decoder_counter>;
    };
}
