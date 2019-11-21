// Copyright Steinwurf ApS 2011.
// Distributed under the "STEINWURF RESEARCH LICENSE 1.0".
// See accompanying file LICENSE.rst or
// http://www.steinwurf.com/licensing

#pragma once

#include <kodo_core/mutable_shallow_storage_layers.hpp>

#include <kodo_fulcrum_sparse/full_vector_decoder_counter.hpp>

namespace kodo_rlnc
{
    /// @ingroup fec_stacks
    ///
    /// @brief Complete stack implementing a shallow storage RLNC decoder.
    ///
    /// The decoder is identical to the full_vector_decoder except for
    /// the fact that is uses a shallow storage layer.
    template
    <
        class Field,
        class Features = meta::typelist<>
    >
    using shallow_full_vector_decoder_counter = full_vector_decoder_counter<
        Field, meta::typelist<kodo_core::mutable_shallow_storage_layers>::
            extend<Features>>;
}
