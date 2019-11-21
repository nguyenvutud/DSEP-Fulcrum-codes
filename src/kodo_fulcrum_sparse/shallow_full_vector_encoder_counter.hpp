// Copyright Steinwurf ApS 2011.
// Distributed under the "STEINWURF RESEARCH LICENSE 1.0".
// See accompanying file LICENSE.rst or
// http://www.steinwurf.com/licensing

#pragma once

#include <kodo_core/const_shallow_storage_layers.hpp>
#include <kodo_fulcrum_sparse/full_vector_encoder_counter.hpp>
//#include <kodo_rlnc/sparse_full_vector_encoder.hpp>

namespace kodo_rlnc
{
    /// @ingroup fec_stacks
    ///
    /// @brief Complete stack implementing a shallow storage RLNC encoder.
    ///
    /// The encoder is identical to the full_vector_encoder except for
    /// the fact that is uses a shallow storage layer.
    template
    <
        class Field,
        class Features = meta::typelist<>
    >
    using shallow_full_vector_encoder_counter = full_vector_encoder_counter<
        Field,
        meta::typelist<
            kodo_core::const_shallow_storage_layers>::
            extend<Features>>;
}
