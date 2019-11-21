// Copyright Steinwurf ApS 2013.
// Distributed under the "STEINWURF RESEARCH LICENSE 1.0".
// See accompanying file LICENSE.rst or
// http://www.steinwurf.com/licensing

#pragma once

#include <cstdint>
#include <algorithm>
#include <vector>
#include <cassert>

#include <sak/aligned_allocator.hpp>

namespace kodo_fulcrum
{
/// @ingroup utility fulcrum
///
/// @brief The role of the coefficient mapper is to produce the
///        mapping between a symbol in the inner code and a symbol
///        in the outer code.
///
/// This layer maps coding coefficients from the inner code to the
/// outer. The coding matrix of the outer code is known i.e. lets
/// say we have 3 original symbols and 2 expansion symbols then
/// the outer encoding vectors could look something like:
///
///             symbol 1:  01 00 00
///             symbol 2:  00 01 00
///             symbol 3:  00 00 01
///   expansion symbol 1:  c1 c2 c3
///   expansion symbol 2:  c4 c5 c6
///
/// Here the coefficients of the outer code is determined by the
/// outer field. In this case {c1,...,c6} represents random values
/// in the outer field.
///
/// An encoding vector of the inner code would then be some linear
/// combination of these 5 encoding vectors.
///
/// Notice that the outer code is systematic so the first 3
/// encoding vectors simply refers to the original symbol (this
/// layer only works with systematic outer codes).
///
/// An example of a inner encoding vector could be:
///
///    1 1 0 0 0
///
/// Lets see how this would map to the outer code:
///
///     1 * (01 00 00)
///          +
///     1 * (00 01 00)
///          +
///     0 * (00 00 01)
///          +
///     0 * (c1 c2 c3)
///          +
///     0 * (c4 c5 c6)
///    ---------------
///          01 01 00
///
/// So the corresponding outer encoding vector is 01 01 00
template<class SuperCoder>
class systematic_coefficient_mapper_high_field : public SuperCoder
{
public:

    /// @copydoc layer::field_type
    using field_type = typename SuperCoder::field_type;

    /// @copydoc layer::value_type
    using value_type = typename SuperCoder::value_type;

public:

    /// @copydoc layer::construct(Factory&)
    template<class Factory>
    void construct(Factory& the_factory)
    {
        SuperCoder::construct(the_factory);

        // Prepare buffer for the coefficients of the outer code
        m_coefficients_lookup.resize(the_factory.max_expansion());

        for (uint32_t i = 0; i < m_coefficients_lookup.size(); ++i)
        {
            auto& v = m_coefficients_lookup[i];
            v.resize(the_factory.max_coefficient_vector_size());
        }
    }

    /// @copydoc layer::initialize(Factory&)
    template<class Factory>
    void initialize(Factory& the_factory)
    {
        SuperCoder::initialize(the_factory);

        for (uint32_t i = 0; i < m_coefficients_lookup.size(); ++i)
        {
            auto& v = m_coefficients_lookup[i];
            SuperCoder::set_seed(i);
            SuperCoder::generate(v.data());
        }
    }

    /// Given an encoded symbol in the inner code map this to the
    /// corresponding coding coefficients in the outer code.
    ///
    /// @param inner_coefficients The coding coefficients of the inner code
    ///
    /// @param outer_coefficients The coding coefficients of the outer code
    void map_to_outer(const uint8_t* inner_coefficients,
                      uint8_t* outer_coefficients)
    {
        assert(inner_coefficients);
        assert(outer_coefficients);

        // Zero the destination coefficients buffer
        std::fill_n(outer_coefficients,
                    SuperCoder::coefficient_vector_size(), 0U);

        // For the systematic part
        for (uint32_t i = 0; i < SuperCoder::symbols(); ++i)
        {
            auto coefficient = SuperCoder::inner_field::get_value(
                inner_coefficients, i);

            if (!coefficient)
                continue;

//            SuperCoder::field::set_value(outer_coefficients, i, 1U);
            SuperCoder::field::set_value(outer_coefficients, i, coefficient);
        }

        // Loop over the inner coding coefficients and create the outer
        // coding coefficients
        for (uint32_t i = 0; i < SuperCoder::expansion(); ++i)
        {
            auto coefficient_index = i + SuperCoder::symbols();
            assert(i < SuperCoder::inner_symbols());

            auto coefficient = SuperCoder::inner_field::get_value(
                inner_coefficients, coefficient_index);

            if (!coefficient)
                continue;

            assert(i < m_coefficients_lookup.size());

            const auto& v = m_coefficients_lookup[i];
            SuperCoder::add(outer_coefficients, v.data(),
                            SuperCoder::coefficient_vector_size());

        }
    }

    /// Given an encoded symbol in the inner code map this to the
       /// corresponding coding coefficients in the outer code.
       ///
       /// @param inner_coefficients The coding coefficients of the inner code
       ///
       /// @param outer_coefficients The coding coefficients of the outer code
       void map_to_outer_high_field(uint8_t* symbol_data, const uint8_t* inner_coefficients,
                         uint8_t* outer_coefficients)
       {
           assert(inner_coefficients);
           assert(outer_coefficients);

           // Zero the destination coefficients buffer
           std::fill_n(outer_coefficients,
                       SuperCoder::coefficient_vector_size(), 0U);

           // For the systematic part
           for (uint32_t i = 0; i < SuperCoder::symbols(); ++i)
           {
               auto coefficient = SuperCoder::inner_field::get_value(
                   inner_coefficients, i);

               if (!coefficient)
                   continue;

   //            SuperCoder::field::set_value(outer_coefficients, i, 1U);
               SuperCoder::field::set_value(outer_coefficients, i, coefficient);
//               SuperCoder::field::get_value();
           }

           // Loop over the inner coding coefficients and create the outer
           // coding coefficients
           for (uint32_t i = 0; i < SuperCoder::expansion(); ++i)
           {
               auto coefficient_index = i + SuperCoder::symbols();
               assert(i < SuperCoder::inner_symbols());

               auto coefficient = SuperCoder::inner_field::get_value(
                   inner_coefficients, coefficient_index);

               if (!coefficient)
                   continue;

               assert(i < m_coefficients_lookup.size());

               const auto& v = m_coefficients_lookup[i];

//               SuperCoder::multiply_add(outer_coefficients, v.data(), coefficient,
//            		   SuperCoder::coefficient_vector_size());
//
               SuperCoder::add(outer_coefficients, v.data(),
                                         SuperCoder::coefficient_vector_size());

               std::cout<<"Coefficient:" <<unsigned(coefficient)<<std::endl;

//               auto temp = SuperCoder::invert(coefficient);
//               std::cout<<"Invert:" <<unsigned(temp)<<std::endl;
//
//               SuperCoder::multiply(outer_coefficients, temp, SuperCoder::coefficient_vector_size());
//               SuperCoder::add(outer_coefficients, v.data(),
//                                           SuperCoder::coefficient_vector_size());
//               auto data_index = SuperCoder::inner_field::get_value(
//                                 symbol_data, coefficient_index);
//               uint8_t* symbol_i = SuperCoder::mutable_symbol(i);
//               std::cout<<"Data symbol i from mutable symbol:" <<unsigned(*symbol_i)<<std::endl;
//               std::cout<<"Symbol data:" <<unsigned(*symbol_data)<<std::endl;
//               std::cout<<"Data symbol i:" <<unsigned(data_index)<<std::endl;
//               SuperCoder::multiply(symbol_data, temp, SuperCoder::coefficient_vector_size());
           }
       }


    /// Given a systematic symbol from the inner code what will be
    /// the coefficients from the outer code
    ///
    /// @param inner_symbol The index of the inner symbol
    ///
    /// @param outer_coefficients The coding coefficients use in
    ///        the outer code.
    void map_uncoded_to_outer(uint32_t inner_symbol,
                              uint8_t* outer_coefficients)
    {
        assert(inner_symbol < SuperCoder::inner_symbols());
        assert(outer_coefficients);

        if (inner_symbol < SuperCoder::symbols())
        {
            // Zero the destination coefficients buffer
            std::fill_n(outer_coefficients,
                        SuperCoder::coefficient_vector_size(), 0U);

            SuperCoder::field::set_value(
                outer_coefficients, inner_symbol, 1U);
        }
        else
        {
            // Get the offset into the extension
            uint32_t offset = inner_symbol - SuperCoder::symbols();

            assert(offset < SuperCoder::expansion());
            assert(offset < m_coefficients_lookup.size());

            const auto& src = m_coefficients_lookup[offset];

            std::copy_n(src.data(), SuperCoder::coefficient_vector_size(),
                        outer_coefficients);
        }
    }

    /// @return True because this coefficient mapper only works
    ///         for a systematic outer code.
    bool is_outer_systematic() const
    {
        return true;
    }

private:

    /// The buffer type
    using aligned_vector =
        std::vector<uint8_t, sak::aligned_allocator<uint8_t>>;

    /// Coefficients use for the inner code
    std::vector<aligned_vector> m_coefficients_lookup;
};
}
