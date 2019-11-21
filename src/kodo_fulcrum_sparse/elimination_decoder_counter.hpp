// Copyright Steinwurf ApS 2011.
// Distributed under the "STEINWURF RESEARCH LICENSE 1.0".
// See accompanying file LICENSE.rst or
// http://www.steinwurf.com/licensing

#pragma once

#include <fifi/default_field.hpp>

#include <kodo_core/symbol_id_decoder.hpp>
#include <kodo_core/plain_symbol_id_reader.hpp>
#include <kodo_core/forward_linear_block_decoder.hpp>
#include <kodo_core/symbol_decoding_status_tracker.hpp>
#include <kodo_core/symbol_decoding_status_counter.hpp>
#include <kodo_core/coefficient_storage.hpp>
#include <kodo_core/deep_symbol_storage.hpp>
#include <kodo_core/storage_bytes_used.hpp>
#include <kodo_core/storage_block_size.hpp>
#include <kodo_core/storage_block_length.hpp>
#include <kodo_core/final_layer.hpp>
#include <kodo_core/finite_field_math.hpp>
#include <kodo_core/finite_field.hpp>
#include <kodo_core/finite_field_counter.hpp>
#include <kodo_core/coefficient_info.hpp>
#include <kodo_core/elimination_coefficient_info.hpp>
#include <kodo_core/elimination_coefficient_value_access.hpp>
#include <kodo_core/elimination_coefficient_offset.hpp>
#include <kodo_core/trace_linear_block_decoder.hpp>
#include <kodo_core/trace_read_symbol.hpp>
#include <kodo_core/trace_symbol.hpp>
#include <kodo_core/trace_decoder_rank.hpp>
#include <kodo_core/trace_layer.hpp>
#include <kodo_core/select_storage_layers.hpp>
#include <kodo_core/deep_storage_layers.hpp>

namespace kodo_core
{
    /// @ingroup fec_stacks
    ///
    /// @brief The elimination decoder uses a non-square decoding
    ///        matrix to ignore certain coding coding
    ///        coefficients. One typical use for this is create an
    ///        elimination decoer for a specific part of encoded
    ///        symbols when the elimination decoder reaches full rank
    ///        it will effectively be able to eliminate that part for
    ///        the encoding, while still producing valid encoded
    ///        symbols.
    template
    <
        class Field,
        class Features = meta::typelist<>
    >
    class elimination_decoder_counter : public
        // Decoder API
        trace_decoder_rank<find_enable_trace<Features>,
        trace_read_symbol<find_enable_trace<Features>,
        trace_symbol<find_enable_trace<Features>,
        trace_linear_block_decoder<find_enable_trace<Features>,
        forward_linear_block_decoder<
        symbol_decoding_status_counter<
        symbol_decoding_status_tracker<
        // Coefficient Storage API
        coefficient_storage<
        elimination_coefficient_value_access<
        elimination_coefficient_info<
        elimination_coefficient_offset<
        coefficient_value_access<
        coefficient_info<
        // Storage API
        select_storage_layers<deep_storage_layers, Features,
        // Finite Field API
		kodo_core::finite_field_counter<
        finite_field_math<typename fifi::default_field<Field>::type,
        finite_field<Field,
        // Trace Layer
        trace_layer<find_enable_trace<Features>,
        // Final Layer
        final_layer
        >>>>>>>>>>>>>>>>>>
    {
    public:
        using factory = pool_factory<elimination_decoder_counter>;
    };
}
