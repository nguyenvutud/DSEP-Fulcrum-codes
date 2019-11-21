// Copyright Steinwurf ApS 2011.
// Distributed under the "STEINWURF RESEARCH LICENSE 1.0".
// See accompanying file LICENSE.rst or
// http://www.steinwurf.com/licensing

#pragma once

#include <kodo_core/block_encoder.hpp>
#include <kodo_core/coefficient_info.hpp>
#include <kodo_core/coefficient_storage.hpp>
#include <kodo_core/coefficient_value_access.hpp>
#include <kodo_core/common_encoder_layers.hpp>
#include <kodo_core/default_on_systematic_encoder.hpp>
#include <kodo_core/final_layer.hpp>
#include <kodo_core/find_enable_trace.hpp>
#include <kodo_core/finite_field_layers.hpp>
#include <kodo_core/finite_field_counter.hpp>
#include <kodo_core/linear_block_encoder.hpp>
#include <kodo_core/payload_info.hpp>
#include <kodo_core/plain_symbol_id_writer.hpp>
#include <kodo_core/select_storage_layers.hpp>
#include <kodo_core/deep_storage_layers.hpp>
#include <kodo_core/storage_aware_encoder.hpp>
#include <kodo_core/symbol_id_encoder.hpp>
#include <kodo_core/trace_layer.hpp>
#include <kodo_core/uniform_generator_layers.hpp>
#include <kodo_core/write_symbol_tracker.hpp>
#include <kodo_core/zero_symbol_encoder.hpp>
#include <kodo_core/select_generator_layers.hpp>
#include <kodo_core/select_symbol_id_writer_layers.hpp>
#include <kodo_core/plain_symbol_id_writer_layers.hpp>

namespace kodo_rlnc
{
    /// @ingroup fec_stacks
    ///
    /// @brief Complete stack implementing a RLNC encoder.
    ///
    /// The key features of this configuration is the following:
    /// - Systematic encoding (uncoded symbols produced before switching
    ///   to coding)
    /// - Full encoding vectors, this stack uses the plain_symbol_id_writer
    ///   which sends the full encoding vector with every encoded symbol.
    ///   Encoding vectors are generated using a random uniform generator.
    /// - Deep symbol storage which makes the encoder allocate its own
    ///   internal memory.
    template
    <
        class Field,
        class Features = meta::typelist<>
    >
    class full_vector_encoder_counter : public
        // Payload Codec API
        kodo_core::payload_info<
        // Block Coder API
        kodo_core::block_encoder<
        // Codec Header API
        kodo_core::default_on_systematic_encoder<
        kodo_core::symbol_id_encoder<
        // Symbol ID API
        kodo_core::select_symbol_id_writer_layers<
            kodo_core::plain_symbol_id_writer_layers, Features,
        // Coefficient Generator API
        kodo_core::select_generator_layers<
            kodo_core::uniform_generator_layers, Features,
        // Codec API
        kodo_core::common_encoder_layers<Features,
        // Coefficient Storage API
        kodo_core::coefficient_value_access<
        kodo_core::coefficient_info<
        // Symbol Storage API
        kodo_core::select_storage_layers<
            kodo_core::deep_storage_layers, Features,
        // Finite Field API
		kodo_core::finite_field_counter<
        kodo_core::finite_field_layers<Field,
        // Trace layer
        kodo_core::trace_layer<kodo_core::find_enable_trace<Features>,
        // Final Layer
        kodo_core::final_layer
        > > > > > > > > > > > > >
    {
    public:
        using factory = kodo_core::pool_factory<full_vector_encoder_counter>;
    };
}
