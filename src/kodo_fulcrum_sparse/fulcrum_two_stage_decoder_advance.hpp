// Copyright Steinwurf ApS 2013.
// Distributed under the "STEINWURF RESEARCH LICENSE 1.0".
// See accompanying file LICENSE.rst or
// http://www.steinwurf.com/licensing

#pragma once

#include <algorithm>

#include <kodo_core/has_deep_symbol_storage.hpp>

namespace kodo_fulcrum
{
/// @ingroup fulcrum
///
/// @brief In this layer we implement the two stage decoder used by the
///        fulcrum combined decoder.
///
/// The basic concept in fulcrum is to create an expansion of the
/// original data (this is the purpose of the outer code). The
/// expansion represents some linear combination of the original
/// data. For example:
///
///                        +------+------ The three coefficients represent
///                        |      |       the encoding vector in the
///                        |      |       outer code.
///                        v      v
///
///             symbol 1:  01 00 00
///             symbol 2:  00 01 00
///             symbol 3:  00 00 01
///   expansion symbol 1:  c1 c2 c3
///   expansion symbol 2:  c4 c5 c6
///
/// In the above example we have three original data symbols and we
/// create two expansion symbols (which is some linear combination of
/// the original symbols).
///
/// As an example expansion symbol 1 is the following linear combination:
///
///   symbol 1 * c1 + symbol 2 * c2 + symbol 3 * c3
///
/// The inner code is then linear combinations of the five symbols
/// (i.e. the orignal symbols and the expansion symbols). As an example
/// the following encoding vector:
///
///   1 0 0 0 1 0
///
/// Would be a linear combination of "symbol 1" and "expansion symbol 1".
///
/// What the two stage decoder does is to create two decoders with the
/// following purpose:
///
/// The stage one decoder: is responsible for decoding the expansion
///   created by the outer encoder. Since all symbols belonging to the
///   expansion are consumed by the stage one decoder we can eliminate
///   them before passing the encoded symbol to the stage two decoder.
///
/// The stage two decoder: is responsible for decoding the original
///   data part. This also means that when the stage two decoder
///   reaches full rank decoding would be complete.
///
/// We try to decode as much as possible in the stage two decoder since
/// this maps to "trivial" encoding vectors in the outer code. The
/// stage one decoder represents the dense high field encoded expansion
/// symbols.
///
/// Decoding is possible when the combined rank of the stage one and
/// stage two decoders matches the number of original symbols, not
/// including the expansion.
///
/// @todo Currently both the stage one and two decoders work on the
///       entire inner coding vector. This might not be needed,
///       something to investigate.
///
/// @tparam StageOneDecoder The type of the stage one decoder. It needs
///         to support decoding from an offset in the encoding
///         vector. Such that we can decode the expansion only.
///
/// @tparam StageTwoDecoder The type of the stage two decoder. It works
///         on the original data part of the encoding vector.
///
/// @tparam SuperCoder @copydoc layer_types::SuperCoder
///
template<class StageOneDecoder, class StageTwoDecoder, class SuperCoder>
class fulcrum_two_stage_decoder_advance : public SuperCoder
{
public:

    /// The stage one decoder
    using stage_one_decoder_type = StageOneDecoder;

    /// The stage one decoder factory
    using stage_one_factory = typename stage_one_decoder_type::factory;

    /// The stage two decoder
    using stage_two_decoder_type = StageTwoDecoder;

    /// The stage two decoder factory
    using stage_two_factory = typename stage_two_decoder_type::factory;

    /// Check storage type
    static_assert(
        kodo_core::has_deep_symbol_storage<StageOneDecoder>::value,
        "Should be a deep storage decoder");

    /// Check storage type
    static_assert(
        kodo_core::has_deep_symbol_storage<StageTwoDecoder>::value,
        "Should be a deep storage decoder");

public:

    /// @ingroup config_layers
    ///
    /// The factory layer associated with this coder.
    class config : public SuperCoder::config
    {
    public:

        /// @copydoc layer::config::config(uint32_t,uint32_t)
        config(uint32_t max_symbols, uint32_t max_symbol_size) :
            SuperCoder::config(max_symbols, max_symbol_size),
            m_stage_one_factory(SuperCoder::config::max_expansion() +
                                max_symbols,
                                max_symbol_size),
            m_stage_two_factory(SuperCoder::config::max_expansion() +
                                max_symbols,
                                max_symbol_size)
        { }

        /// Build and return the stage one decoder
        /// @return The stage one decoder
        typename stage_one_factory::pointer build_stage_one()
        {
            uint32_t symbol_size = SuperCoder::config::symbol_size();
            uint32_t symbols = SuperCoder::config::expansion();
            uint32_t offset = SuperCoder::config::symbols();

            m_stage_one_factory.set_symbol_size(symbol_size);
            m_stage_one_factory.set_symbols(symbols);
            m_stage_one_factory.set_elimination_offset(offset);

            return m_stage_one_factory.build();
        }

        /// Build and return the stage two decoder
        /// @return The stage two decoder
        typename stage_two_factory::pointer build_stage_two()
        {
            uint32_t symbol_size = SuperCoder::config::symbol_size();
            uint32_t symbols = SuperCoder::config::symbols() +
                               SuperCoder::config::expansion();

            m_stage_two_factory.set_symbol_size(symbol_size);
            m_stage_two_factory.set_symbols(symbols);

            return m_stage_two_factory.build();
        }

    protected:

        /// The stage one factory
        stage_one_factory m_stage_one_factory;

        /// The stage two factory
        stage_two_factory m_stage_two_factory;
    };

public:

    /// @copydoc layer::construct(Factory&)
    template<class Factory>
    void construct(Factory& the_factory)
    {
        SuperCoder::construct(the_factory);

        m_stage_one_decoder_copied.resize(
            the_factory.max_expansion());

        m_stage_two_decoder_copied.resize(
            the_factory.max_expansion() + the_factory.max_symbols());

        m_outer_coefficients.resize(
            the_factory.max_coefficient_vector_size());

        m_outer_symbol.resize(
            the_factory.max_symbol_size());
    }

    /// @copydoc layer::initialize(Factory&)
    template<class Factory>
    void initialize(Factory& the_factory)
    {
        SuperCoder::initialize(the_factory);

        m_stage_one_decoder = the_factory.build_stage_one();
        m_stage_two_decoder = the_factory.build_stage_two();

        assert(m_stage_one_decoder->symbols() ==
               SuperCoder::expansion());

        assert(m_stage_two_decoder->symbols() ==
               SuperCoder::expansion() + SuperCoder::symbols());

        assert(m_stage_one_decoder->symbol_size() ==
               m_stage_two_decoder->symbol_size());

        assert(m_stage_one_decoder->symbol_size() ==
               SuperCoder::symbol_size());

        assert(m_stage_one_decoder->coefficient_vector_size() ==
               m_stage_two_decoder->coefficient_vector_size());

        assert(m_stage_one_decoder->symbols() <=
               m_stage_one_decoder_copied.size());

        assert(m_stage_two_decoder->symbols() <=
               m_stage_two_decoder_copied.size());

        std::fill_n(m_stage_two_decoder_copied.begin(),
                    m_stage_two_decoder->symbols(), false);

        std::fill_n(m_stage_one_decoder_copied.begin(),
                    m_stage_one_decoder->symbols(), false);
    }

    /// @copydoc layer::read_uncoded_symbol(uint8_t*,uint32_t)
    void read_uncoded_symbol(uint8_t* symbol_data, uint32_t symbol_index)
    {
        if (SuperCoder::is_outer_systematic())
        {
            if (symbol_index < SuperCoder::symbols())
            {
                // In this case the symbol represents original
                // data we will pass it first the "stage two"
                // decoder to make sure the symbol is available
                // for later elimination of coded symbols. After
                // that we immediately copy it ot the main decoder
                // (SuperCoder) and mark it as copied.
                m_stage_two_decoder->read_uncoded_symbol(
                    symbol_data, symbol_index);

                assert(symbol_index < m_stage_two_decoder_copied.size());
                m_stage_two_decoder_copied[symbol_index] = true;

                SuperCoder::read_uncoded_symbol(symbol_data, symbol_index);
            }
            else
            {
                assert(SuperCoder::symbols() <= symbol_index);

                assert(symbol_index <
                       SuperCoder::symbols() + SuperCoder::expansion());

                // In this case the systematic symbol represents
                // the extension of the outer code. So we pass it
                // to the "stage one" decoder and check whether
                // decoding is possible.
                uint32_t stage_one_index =
                    symbol_index - SuperCoder::symbols();

                m_stage_one_decoder->read_uncoded_symbol(
                    symbol_data, stage_one_index);
            }

            check_combined_rank();
        }
        else
        {
            // We have not yet implemented the case where the
            // outer code is not systematic
            assert(0);
        }
    }

    /// @copydoc layer::read_symbol(uint8_t*,uint8_t*)
    void read_symbol(uint8_t* symbol_data, uint8_t* coefficients)
    {
        uint32_t rank_stage_one = m_stage_one_decoder->rank();
        if(rank_stage_one < m_stage_one_decoder->symbols()){

        // We always pass the packet to the stage one decoder it
        // is responsible for eliminating the extension.
        m_stage_one_decoder->read_symbol(symbol_data, coefficients);

        if (m_stage_one_decoder->rank() > rank_stage_one)
        {
            // If the rank increased in the stage one coder we
            // check whether we have achieved full rank. We do not
            // need to go to the stage two decoder since we have
            // accumulated the additional degree of freedom
            // contained in the encoded symbol.
            check_combined_rank();
            return;
        }
        }

        // Getting here means that the stage one decoder did not
        // increase its rank, next step is to pass the reduced
        // symbol to the stage two decoder
        uint32_t rank_stage_two = m_stage_two_decoder->rank();

        m_stage_two_decoder->read_symbol(symbol_data, coefficients);

        if (m_stage_two_decoder->rank() > rank_stage_two)
        {
            // If the stage two decoder increased the rank we
            // check whether we can map anything to the outer
            // decoder
            check_combined_rank();
        }
    }

    /// This function is checking the combined rank of the stage one and
    /// stage two decoders. If the rank is sufficient we may be able to
    /// decode the outer code we map the symbols, that have not previously,
    /// been mapped to the outer code.
    void check_combined_rank()
    {
        uint32_t combined_rank = m_stage_one_decoder->rank() +
                                 m_stage_two_decoder->rank();

        assert(combined_rank <= SuperCoder::symbols() +
               SuperCoder::expansion());

        if (combined_rank < SuperCoder::symbols())
            return;

        // We map from stage two first. The is advantageous if the outer
        // encoder has produced systematic symbols. Since in that case we
        // are able to map symbols directly between the outer and inner
        // decoders.
        for (uint32_t i = 0; i < m_stage_two_decoder->symbols(); ++i)
        {
            bool is_pivot = m_stage_two_decoder->is_symbol_pivot(i);
            bool is_copied = m_stage_two_decoder_copied[i];

            if (is_pivot && !is_copied)
            {
                m_stage_two_decoder_copied[i] = true;
                map_stage_two(i);
            }
        }


        for (uint32_t i = 0; i < m_stage_one_decoder->symbols(); ++i)
        {
            bool is_pivot = m_stage_one_decoder->is_symbol_pivot(i);
            bool is_copied = m_stage_one_decoder_copied[i];

            if (is_pivot && !is_copied)
            {
                m_stage_one_decoder_copied[i] = true;
                map_stage_one(i);
            }
        }
    }

    /// Takes a specific symbol from the stage one decoder and
    /// maps it to the outer decoder.
    ///
    /// @param symbol_index The index of the symbol to map
    void map_stage_one(uint32_t symbol_index)
    {
        const uint8_t* symbol_data =
            m_stage_one_decoder->const_symbol(symbol_index);

        const uint8_t* coefficients =
            m_stage_one_decoder->coefficient_vector_data(symbol_index);

        map_symbol(symbol_data, coefficients);
    }

    /// Takes a specific symbol from the stage two decoder and
    /// maps it to the outer decoder.
    ///
    /// @param symbol_index The index of the symbol to map
    void map_stage_two(uint32_t symbol_index)
    {
        const uint8_t* symbol_data =
            m_stage_two_decoder->const_symbol(symbol_index);

        const uint8_t* coefficients =
            m_stage_two_decoder->coefficient_vector_data(symbol_index);

        map_symbol(symbol_data, coefficients);
    }

    /// Maps the coding coefficients from the inner code to the
    /// outer and copy the symbol for decoding in the outer code.
    ///
    /// @param symbol_index The index of the symbol to map
    void map_symbol(const uint8_t* symbol_data,
                    const uint8_t* coefficients)
    {
        // Copy the pivot symbol to a modifiable buffer, this copy
        // might not be needed - if we could guarentee
        // decoding. But we cannot right now.
        std::copy_n(symbol_data, SuperCoder::symbol_size(),
                    m_outer_symbol.data());

        SuperCoder::map_to_outer(
            coefficients, m_outer_coefficients.data());

        SuperCoder::read_symbol(
            m_outer_symbol.data(), m_outer_coefficients.data());
    }

protected:

    /// @return The stage one decoder
    stage_one_decoder_type& stage_one_decoder()
    {
        return *m_stage_one_decoder;
    }

    /// @return The stage two decoder
    stage_two_decoder_type& stage_two_decoder()
    {
        return *m_stage_two_decoder;
    }

protected:

    /// The stage one decoder
    typename stage_one_factory::pointer m_stage_one_decoder;

    /// The stage two decoder
    typename stage_two_factory::pointer m_stage_two_decoder;

    /// Keeps track of which symbols have been copied from the
    /// stage one decoder to the inner decoder
    std::vector<bool> m_stage_one_decoder_copied;

    /// Keeps track of which symbols have been copied from the
    /// stage two decoder to the inner decoder
    std::vector<bool> m_stage_two_decoder_copied;

    /// Buffer for the outer coding coefficients
    std::vector<uint8_t> m_outer_coefficients;

    /// Buffer for the outer coding symbol data
    std::vector<uint8_t> m_outer_symbol;
};
}
