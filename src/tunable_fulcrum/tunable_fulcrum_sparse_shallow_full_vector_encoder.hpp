// Copyright Steinwurf ApS 2011.
// Distributed under the "STEINWURF RESEARCH LICENSE 1.0".
// See accompanying file LICENSE.rst or
// http://www.steinwurf.com/licensing

#pragma once

#include <kodo_core/const_shallow_storage_layers.hpp>
//#include "full_vector_encoder.hpp"
//#include <kodo_rlnc/sparse_full_vector_encoder.hpp>
#include <tunable_fulcrum/tunable_sparse_full_vector_encoder.hpp>

namespace kodo_fulcrum
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
    using tunable_fulcrum_sparse_shallow_full_vector_encoder = kodo_fulcrum::tunable_sparse_full_vector_encoder<
        Field,
        meta::typelist<
            kodo_core::const_shallow_storage_layers>::
            extend<Features>>;
}
