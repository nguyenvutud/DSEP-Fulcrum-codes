// Copyright Steinwurf ApS 2013.
// Distributed under the "STEINWURF RESEARCH LICENSE 1.0".
// See accompanying file LICENSE.rst or
// http://www.steinwurf.com/licensing

#pragma once

#include <cstdint>

#include <kodo_core/finite_field_info.hpp>

#include <fifi/binary.hpp>

namespace kodo_fulcrum
{
/// @ingroup fulcrum
///
/// @brief Keeps track of the expansion factor used in the anycode
///        implementation.
///
/// @tparam MaxExpansion The maximum number of expansion symbols
///         that will be supported. Specified e.g. using
///         std::integral_constant.
///
/// @tparam DefaultExpansion The default number of expansion
///         symbols used (can be changed at run time with the
///         fulcrum::set_expansion(uint32_t) function.
///
/// @tparam SuperCoder @copydoc layer_types::SuperCoder
///
template
<
	class Field,
    class MaxExpansion,
    class DefaultExpansion,
    class SuperCoder
>
class fulcrum_sparse_info : public SuperCoder
{
public:

    /// The finite field used by the inner code. Currently we only support
    /// binary, in the future this might change.

	using inner_field = kodo_core::finite_field_info<Field>;  //Dinh nghia inner field
//	using inner_field = kodo_core::finite_field_info<fifi::binary4>;  //Dinh nghia inner field

public:

    /// @ingroup config_layers
    /// The factory layer associated with this coder.
    class config : public SuperCoder::config
    {
    public:

        /// @copydoc layer::factory::factory(uint32_t,uint32_t)
        config(uint32_t max_symbols, uint32_t max_symbol_size) :
            SuperCoder::config(max_symbols, max_symbol_size),
            m_expansion(DefaultExpansion::value),
            m_max_inner_symbols(max_symbols + MaxExpansion::value)
        {
            assert(m_max_inner_symbols);
        }

        /// @return The maximum expansion supported
        uint32_t max_expansion() const
        {
            return MaxExpansion::value;
        }

        /// @return The expansion factor used. The expansion factor
        ///         denotes the number of additional symbols
        ///         created by the outer code.
        uint32_t expansion() const
        {
            return m_expansion;
        }

        /// Sets the number of expansion symbols
        /// @param expansion The number of expansion symbols to use
        void set_expansion(uint32_t expansion)
        {
            assert(expansion <= MaxExpansion::value);
            m_expansion = expansion;
        }

        /// @return the maximum number of symbols in the inner code
        uint32_t max_inner_symbols() const
        {
            return m_max_inner_symbols;
        }

    private:

        /// The number of symbols in the inner code
        uint32_t m_expansion;

        /// The total number of symbols in the inner code
        uint32_t m_max_inner_symbols;
    };

public:

    /// @copydoc layer::initialize(Factory&)
    template<class Factory>
    void initialize(Factory& the_factory)
    {
        SuperCoder::initialize(the_factory);
        m_expansion = the_factory.expansion();
        m_inner_symbols = m_expansion + the_factory.symbols();
    }

    /// @return The maximum expansion supported
    uint32_t max_expansion() const
    {
        return MaxExpansion::value;
    }

    /// @return The expansion factor used. The expansion factor
    ///         denotes the number of additional symbols
    ///         created by the outer code.
    uint32_t expansion() const
    {
        return m_expansion;
    }

    /// @return the number of symbols in the inner code
    uint32_t inner_symbols() const
    {
        return m_inner_symbols;
    }

private:

    /// The number of symbols added by the inner code
    uint32_t m_expansion;

    /// The total number of symbols in the inner code
    uint32_t m_inner_symbols;
};
}
