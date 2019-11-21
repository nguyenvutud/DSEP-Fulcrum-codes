// Copyright Steinwurf ApS 2015.
// Distributed under the "STEINWURF RESEARCH LICENSE 1.0".
// See accompanying file LICENSE.rst or
// http://www.steinwurf.com/licensing

#pragma once

#include <boost/random/linear_congruential.hpp>
#include <boost/random/binomial_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include <cstdint>
#include <cassert>

namespace kodo_fulcrum
{
/// @ingroup coefficient_generator_layers
///
/// @brief Generate uniformly distributed coefficients with a specific
///        density.
///
/// If we are using a sparse code with a very low density, we increase
/// the probability of generating all zero encoding vectors and in turn
/// produce encoded symbols containing no useful information at all.
///
/// This layer implements an alternative approach to generate sparse coding
/// vectors that cannot be all zero. It generates a binomial random number
/// between 1 and SuperCoder::symbols() that decides the number nonzero
/// coefficients to set in each coding vector. The coefficients are set one
/// by one and inserted in the coding vector on a random index that is zero.
template<class SuperCoder>
class tunable_expansion_uniform_index_generator : public SuperCoder
{
public:

    /// @copydoc layer::value_type
    using value_type = typename SuperCoder::value_type;

    /// The random generator used
    using generator_type = boost::random::rand48;

    /// @copydoc layer::seed_type
    using seed_type = generator_type::result_type;

public:

    /// Constructor
    tunable_expansion_uniform_index_generator() :
        m_value_distribution(1, SuperCoder::field::max_value())
    {
    	m_non_zero_expansion = 0;
	m_number_nonzero_coefficient = 0;
    }

    /// @copydoc layer::initialize(Factory& factory)
    template<class Factory>
    void initialize(Factory& factory)
    {
        SuperCoder::initialize(factory);

        assert(SuperCoder::symbols() > 0);
        m_symbol_selection =
            symbol_selection_distribution(0, SuperCoder::symbols() - 1);

        // We can reset the binomial distribution by resetting the density.
        // This should be done since the generation size (symbols) may have
        // changed
        set_density(density());
    }
    void set_coefficient_value_at_index(uint8_t* coefficients, uint32_t index)
    {
    	if (SuperCoder::field::is_binary())
		{
			SuperCoder::set_coefficient_value(coefficients, index, 1U);
		}
		else
		{
			auto coefficient = m_value_distribution(m_random_generator);
			SuperCoder::set_coefficient_value(coefficients, index,
											  coefficient);
		}
    }
    /// @copydoc layer::generate(uint8_t*)
    void generate(uint8_t* coefficients)
    {
        assert(coefficients != nullptr);

        // Since we will not set all coefficients we should ensure
        // that the non-specified ones are zero
        std::fill_n(coefficients, SuperCoder::coefficient_vector_size(), 0);
        //NVu: changed
//        m_non_zero_count =
//                   std::max<uint32_t>(1U,  m_number_nonzero_coefficient);
        //uint32_t expansion = SuperCoder::expansion();
        uint32_t symbols = SuperCoder::symbols();
//        std::cout<<"Symbols:"<<symbols<<std::endl;

        // The number of non-zero coefficients we should set
	m_non_zero_count = m_number_nonzero_coefficient;
	if (m_number_nonzero_coefficient == 0)
	{
		m_non_zero_count =
		    std::max<uint32_t>(1U,  m_binomial(m_random_generator));
		// To avoid generating m_non_zero_count > real symbols, while num_non_zero coef = 0
		while (m_non_zero_count > symbols - m_expansion)
			m_non_zero_count =
			            std::max<uint32_t>(1U,  m_binomial(m_random_generator));
	}


//        std::cout<<"Number of nonzero in generate function:"<<m_non_zero_count<<std::endl;
        // The number of non zero coefficients we have set
        uint32_t non_zero_count = 0;
        uint32_t start_exp_position = symbols;
        bool all_extra_coef_zero = false;
        //uint32_t upper_bound = symbols - m_expansion + m_non_zero_expansion;
        uint32_t real_symbol = symbols - m_expansion;
        uint32_t upper_bound = real_symbol - 1;

        do
        {
        	uint32_t index = 0;
        	//check whether all_non_zero expansion coefficient is set or not
        	if (m_all_non_zero)
        	{
        		upper_bound = real_symbol - 1;
        		m_normal = false;

        		for (index = real_symbol; index < real_symbol + m_non_zero_expansion; index = index + 1)
        			if (non_zero_count < m_non_zero_count)
        			{
        				set_coefficient_value_at_index(coefficients, index);
        				non_zero_count +=1;
        			}
        			else
					{
        				break;
					}
        		m_all_non_zero = false;

        	}
        	else if (m_first_increment)
        	{
        		upper_bound = real_symbol - 1;
        		m_normal = false;

        		set_coefficient_value_at_index(coefficients, real_symbol + m_non_zero_expansion - 1);
        		non_zero_count +=1;
        		m_first_increment = false;
        	}
        	else{
        		if (m_normal)
        			upper_bound = real_symbol + m_non_zero_expansion - 1;
        		else
        			upper_bound = real_symbol - 1;

        		index = m_symbol_selection(m_random_generator);

             	if (index > upper_bound)
					continue;

             	// Don't overwrite a coefficient that is already set
             	// This condition was missing in TEP paper
				if (SuperCoder::coefficient_value(coefficients, index) != 0)
				{
					continue;
				}


             	set_coefficient_value_at_index(coefficients, index);

             	++non_zero_count;

        	}
        }
        while (non_zero_count < m_non_zero_count);

        m_normal = true;
    }

    /// @copydoc layer::generate_partial(uint8_t*)
    void generate_partial(uint8_t* coefficients)
    {
        assert(coefficients != nullptr);

        // There should be at least 1 symbol available
        assert(SuperCoder::rank() > 0);

        // Since we will not set all coefficients we should ensure
        // that the non-specified ones are zero
        std::fill_n(coefficients, SuperCoder::coefficient_vector_size(), 0);

        // changed
        m_non_zero_count = std::max(1U, std::min<uint32_t>(
                         m_number_nonzero_coefficient,
                         SuperCoder::rank()));
       //        std::cout<<"Number of nonzero in generate_partial() function:"<<m_non_zero_count<<std::endl;

        // The number of non-zero coefficients we should set
//        m_non_zero_count = std::max(1U, std::min<uint32_t>(
//            m_binomial(m_random_generator),
//            SuperCoder::rank()));


        // The number of non zero coefficients we have set
        uint32_t non_zero_count = 0;
        do
        {
        	std::cout<<"PERFORMING GENERATE PARTIAL\n";
            // Generate coefficient index that we should write nonzero
            // values to
            uint32_t index = m_symbol_selection(m_random_generator);

            // Retry if we have not received the symbol
            if (!SuperCoder::is_symbol_pivot(index))
            {
                continue;
            }

            // Don't overwrite a coefficient that is already set
            if (SuperCoder::coefficient_value(coefficients, index) != 0)
            {
                continue;
            }

            // Set coefficient
            if (SuperCoder::field::is_binary())
            {
                SuperCoder::set_coefficient_value(coefficients, index, 1U);
            }
            else
            {
                auto coefficient = m_value_distribution(m_random_generator);
                SuperCoder::set_coefficient_value(coefficients, index,
                                                  coefficient);
            }

            ++non_zero_count;

        }
        while (non_zero_count < m_non_zero_count);
    }

    /// @copydoc layer::set_seed(seed_type)
    void set_seed(seed_type seed_value)
    {
        m_random_generator.seed(seed_value); //Nvu: se tao bo sinh cho coefficient_lookup dua vao Index (row) cua no
    }

    void set_num_nonzero_expansion(uint32_t num_nonzero_expansion, uint32_t exp,
    		const bool first_increment = false, const bool all_non_zero = false,
			const bool normal = true)
    {
    	m_non_zero_expansion = num_nonzero_expansion;
    	m_expansion = exp;
    	m_first_increment = first_increment;
    	m_all_non_zero = all_non_zero;
    	m_normal = normal;
    }

    /// Set the density of the coefficients generated
    ///
    /// @param density coefficients density
    void set_density(double density)
    {
        assert(density > 0.0);

        if (SuperCoder::field::is_binary())
        {
            // If binary, the density should be less than 1
            assert(density < 1.0 && "The density must be less than 1.0 "
                   "for the binary field to avoid generating all-one "
                   "coefficient vectors");
        }
        else
        {
            // If not binary, the density should be less than or equal to 1
            assert(density <= 1.0 && "The density must be less than or "
                   "equal to 1.0 for all fields except the binary field");
        }

        m_binomial = binomial_distribution(SuperCoder::symbols(), density);
    }

    /// Set the average number of nonzero symbols
    ///
    /// @param symbols the average number of nonzero symbols
    void set_average_nonzero_symbols(double symbols)
    {
        set_density(symbols / SuperCoder::symbols());
    }

    /// Get the density of the coefficients generated
    ///
    /// @return the density of the generator
    double density() const
    {
        return m_binomial.p();
    }

    /// @return The number of nonzero coefficient generated
    uint32_t nonzeros_generated() const
    {
        return m_non_zero_count;
    }
    void set_number_nonzero_coefficient(uint32_t symbols)
    {
    	m_number_nonzero_coefficient = symbols;
    }
    uint32_t get_number_nonzero_coefficient() const
    {
       		return m_number_nonzero_coefficient;
    }
private:

    /// The random generator
    generator_type m_random_generator;

    /// Non-zero counter - i.e. the number of non zero coefficients
    /// generated
    uint32_t m_non_zero_count, m_expansion = 0;

    uint32_t m_non_zero_expansion = 0;

    uint32_t m_number_nonzero_coefficient = 0;

    bool m_first_increment = false;

    bool m_all_non_zero = false;

    bool m_normal = true;

    /// The type of the binomial distribution
    using binomial_distribution =
        boost::random::binomial_distribution<>;

    /// The distribution controlling the density of the coefficients
    binomial_distribution m_binomial;

    /// The type of the symbol_selection distribution
    using symbol_selection_distribution =
        boost::random::uniform_int_distribution<uint32_t>;

    /// Distribution that generates random indexes
    symbol_selection_distribution m_symbol_selection;

    /// The type of the value_type distribution
    using value_type_distribution =
        boost::random::uniform_int_distribution<value_type>;

    /// Distribution that generates random values from a finite field
    value_type_distribution m_value_distribution;
};
}
