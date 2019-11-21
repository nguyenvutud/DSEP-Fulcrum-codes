// Copyright Steinwurf ApS 2011.
// Distributed under the "STEINWURF RESEARCH LICENSE 1.0".
// See accompanying file LICENSE.rst or
// http://www.steinwurf.com/licensing

#pragma once

#include <kodo_core/coefficient_info.hpp>
#include <kodo_core/coefficient_value_access.hpp>
#include <kodo_core/common_encoder_layers.hpp>
#include <kodo_core/default_on_systematic_encoder.hpp>
#include <kodo_core/final_layer.hpp>
#include <kodo_core/finite_field_layers.hpp>
#include <kodo_core/finite_field_counter.hpp>
#include <kodo_core/linear_block_encoder.hpp>
#include <kodo_core/payload_info.hpp>
#include <kodo_core/plain_symbol_id_writer.hpp>
#include <kodo_core/plain_symbol_id_size.hpp>
#include <kodo_core/pool_factory.hpp>
#include <kodo_core/seed_generator.hpp>
#include <kodo_core/select_storage_layers.hpp>
#include <kodo_core/deep_storage_layers.hpp>
#include <kodo_core/sparse_uniform_index_generator.hpp>
#include <kodo_core/sparse_uniform_generator.hpp>
#include <kodo_core/storage_aware_encoder.hpp>
#include <kodo_core/symbol_id_encoder.hpp>
#include <kodo_core/trace_layer.hpp>
#include <kodo_core/write_symbol_tracker.hpp>
#include <kodo_core/zero_symbol_encoder.hpp>

namespace kodo_rlnc
{
/// @ingroup fec_stacks
///
/// @brief Complete stack implementing a sparse RLNC encoder.
///
/// The encoder is identical to the full_vector_encoder except for
/// the fact that it uses a density based random generator, which can be
/// used to control the density i.e. the number of non-zero elements in
/// the encoding vector.
template
<
    class Field,
    class Features = meta::typelist<>
>
class fulcrum_sparse_full_vector_encoder : public
    // Payload Codec API
    kodo_core::payload_info<
    // Codec Header API
    kodo_core::default_on_systematic_encoder<
    kodo_core::symbol_id_encoder<
    // Symbol ID API
    kodo_core::plain_symbol_id_writer<
    kodo_core::plain_symbol_id_size<
    // Coefficient Generator API
    kodo_core::seed_generator<
    kodo_core::sparse_uniform_index_generator<
    // Encoder API
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
    // Trace Layer
    kodo_core::trace_layer<kodo_core::find_enable_trace<Features>,
    // Final Layer
    kodo_core::final_layer
    > > > > > > > > > > > > > >
{
public:
    using factory = kodo_core::pool_factory<fulcrum_sparse_full_vector_encoder>;
};
}
