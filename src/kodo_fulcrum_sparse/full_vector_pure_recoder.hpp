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

#include <kodo_rlnc/on_the_fly_recoder_type.hpp>
#include <kodo_rlnc/recoder_read_uncoded_symbol.hpp>
#include <kodo_rlnc/recoder_mix_in.hpp>
#include <kodo_rlnc/recoder_mix_out.hpp>
#include <kodo_rlnc/recoder_factory.hpp>
#include <kodo_rlnc/recoder_coefficient_storage.hpp>
#include <kodo_rlnc/recoder_coefficients_generator.hpp>
#include <kodo_rlnc/recoder_symbol_storage.hpp>
#include <kodo_rlnc/recoder_symbols_info.hpp>
#include <kodo_rlnc/recoder_symbol_id_writer.hpp>
#include <kodo_rlnc/trace_recoder.hpp>

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
    class full_vector_pure_recoder : public
        // Payload API
        kodo_core::nested_write_payload<
        kodo_core::basic_proxy_stack<
            kodo_core::proxy_args<Features>, full_vector_recoding_stack,
        kodo_core::payload_info<
		// Codec API
		kodo_core::non_systematic_encoder<
		kodo_core::systematic_decoder<
		kodo_core::systematic_base_coder<
		kodo_core::symbol_id_encoder<
		kodo_core::symbol_id_decoder<
		// Select Recoding Layers
		trace_recoder<kodo_core::find_enable_trace<Features>,
		recoder_symbol_id_writer<
		recoder_mix_out<
		recoder_read_uncoded_symbol<
		recoder_mix_in<
		recoder_coefficients_generator<
		recoder_symbol_storage<
		recoder_coefficient_storage<
		recoder_factory<
		on_the_fly_recoder_type<Features,
		recoder_symbols_info<
		kodo_core::coefficient_value_access<
		//select_recoding_layers<recode_in_recode_out_layers, Features,
		// Symbol ID API
		kodo_core::plain_symbol_id_reader<
		kodo_core::plain_symbol_id_size<
		// Coefficient Storage API
//		kodo_core::coefficient_info<
//		// Storage API
//		kodo_core::storage_block_size<

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
		> > > > > > > > > > > > > > > > > > > > > > > > > > >
    {
    public:
        using factory = kodo_core::pool_factory<full_vector_pure_recoder>;
    };
}
