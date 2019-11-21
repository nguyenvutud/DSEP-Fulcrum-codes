// Copyright Steinwurf ApS 2011.
// Distributed under the "STEINWURF RESEARCH LICENSE 1.0".
// See accompanying file LICENSE.rst or
// http://www.steinwurf.com/licensing

#pragma once

#include <kodo_fulcrum_sparse/fulcrum_sparse_codes.hpp>

namespace kodo_fulcrum
{
template<class Field>
class full_vector_encoder_count : public
    // Payload Codec API
    kodo_core::payload_info<
    // Codec Header API
    kodo_core::default_on_systematic_encoder<
    kodo_core::symbol_id_encoder<
    // Symbol ID API
    kodo_core::plain_symbol_id_writer<
    kodo_core::plain_symbol_id_size<
    // Coefficient Generator API
    kodo_core::uniform_generator<
    // Encoder API
    kodo_core::write_symbol_tracker<
    kodo_core::zero_symbol_encoder<
    kodo_core::linear_block_encoder<
    kodo_core::storage_aware_encoder<
    // Coefficient Storage API
    kodo_core::coefficient_value_access<
    kodo_core::coefficient_info<
    // Symbol Storage API
    kodo_core::const_partial_shallow_symbol_storage<
    kodo_core::const_shallow_symbol_storage<
    kodo_core::shallow_symbol_storage<storage::const_storage,
    kodo_core::storage_bytes_used<
    kodo_core::storage_block_length<
    kodo_core::storage_block_size<
    // Finite Field API
    kodo_core::finite_field_counter<
    kodo_core::finite_field_math<typename fifi::default_field<Field>::type,
    kodo_core::finite_field<Field,
    // Final Layer
    kodo_core::final_layer
    > > > > > > > > > > > > > > > > > > > > >
{
public:
    using factory = kodo_core::pool_factory<full_vector_encoder_count>;
};

template<class Field>
class full_vector_decoder_count : public
    // Payload API
    kodo_core::payload_info<
    // Codec Header API
    kodo_core::systematic_decoder_layers<
    kodo_core::symbol_id_decoder<
    // Symbol ID API
    kodo_core::plain_symbol_id_reader<
    kodo_core::plain_symbol_id_size<
    // Decoder API
    kodo_core::forward_linear_block_decoder<
    kodo_core::symbol_decoding_status_counter<
    kodo_core::symbol_decoding_status_tracker<
    // Coefficient Storage API
    kodo_core::coefficient_value_access<
    kodo_core::coefficient_storage<
    kodo_core::coefficient_info<
    // Storage API
    kodo_core::deep_symbol_storage<
    kodo_core::storage_bytes_used<
    kodo_core::storage_block_length<
    kodo_core::storage_block_size<
    // Finite Field API
    kodo_core::finite_field_counter<
    kodo_core::finite_field_math<typename fifi::default_field<Field>::type,
    kodo_core::finite_field<Field,
    // Final Layer
    kodo_core::final_layer
    > > > > > > > > > > > > > > > > > >
{
public:
    using factory = kodo_core::pool_factory<full_vector_decoder_count>;
};

template<class Field>
class full_delayed_rlnc_decoder_count : public
    // Payload API
    kodo_core::payload_info<
    // Codec Header API
    kodo_core::systematic_decoder_layers<
    kodo_core::symbol_id_decoder<
    // Symbol ID API
    kodo_core::plain_symbol_id_reader<
    kodo_core::plain_symbol_id_size<
    // Decoder API
    kodo_core::linear_block_decoder_delayed<
    kodo_core::forward_linear_block_decoder<
    kodo_core::symbol_decoding_status_counter<
    kodo_core::symbol_decoding_status_tracker<
    // Coefficient Storage API
    kodo_core::coefficient_value_access<
    kodo_core::coefficient_storage<
    kodo_core::coefficient_info<
    // Storage API
    kodo_core::deep_symbol_storage<
    kodo_core::storage_bytes_used<
    kodo_core::storage_block_length<
    kodo_core::storage_block_size<
    // Finite Field API
    kodo_core::finite_field_counter<
    kodo_core::finite_field_math<typename fifi::default_field<Field>::type,
    kodo_core::finite_field<Field,
    // Final Layer
    kodo_core::final_layer
    > > > > > > > > > > > > > > > > > > >
{
public:
    using factory =
        kodo_core::pool_factory<full_delayed_rlnc_decoder_count>;
};
}
