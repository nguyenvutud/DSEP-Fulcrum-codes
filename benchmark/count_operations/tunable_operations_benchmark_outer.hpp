// Copyright Steinwurf ApS 2011.
// Distributed under the "STEINWURF RESEARCH LICENSE 1.0".
// See accompanying file LICENSE.rst or
// http://www.steinwurf.com/licensing

#include <cassert>
#include <ctime>
#include <memory>
#include <cmath>

#include <gauge/console_printer.hpp>
#include <gauge/csv_printer.hpp>
#include <gauge/gauge.hpp>
#include <gauge/python_printer.hpp>


#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <tables/table.hpp>

#include <kodo_core/set_systematic_off.hpp>

#include <kodo_fulcrum/fulcrum_codes.hpp>
#include <kodo_core/operations_counter.hpp>

#include "time_measurement.hpp"

template<class Encoder, class Recoder, class Decoder>
class tunable_operations_benchmark_outer : public gauge::benchmark
{
public:

    using encoder_factory = typename Encoder::factory;
    using encoder_ptr = typename Encoder::factory::pointer;

    using decoder_factory = typename Decoder::factory;
    using decoder_ptr = typename Decoder::factory::pointer;

    using recoder_factory = typename Recoder::factory;
    using recoder_ptr = typename Recoder::factory::pointer;

    tunable_operations_benchmark_outer()
    {
        // Seed the random generator controlling the erasures
        m_random_generator.seed((uint32_t)time(0));
    }
    /// Starts a measurement and saves the counter
    void start()
    {
    	tm.start();
        m_encoder->reset_operations_counter();
        m_decoder->reset_operations_counter();
        m_recoder->reset_operations_counter();
        m_encoder->nested()->reset_operations_counter();
        m_recoder->nested()->reset_operations_counter();
        m_decoder->nested()->reset_operations_counter();

    }

    /// Stops a measurement and saves the counter
    void stop()
    {
        gauge::config_set cs = get_current_configuration();

        std::string type = cs.get_value<std::string>("type");
        if (type == "encoder")
        {
            m_counter_nested = m_encoder->nested()->get_operations_counter();
            m_counter = m_encoder->get_operations_counter();
        }
        else if (type == "decoder")
        {
	    m_counter_nested = kodo_core::operations_counter();

            m_counter = m_decoder->get_operations_counter();

        }
        else if(type =="recoder")
        {
        	m_counter_nested = m_recoder->nested()->get_operations_counter();
        	m_counter = m_recoder->get_operations_counter();
        }
        else
        {
            assert(0);
        }

        tm.stop();
    }

    void store_run(tables::table& results)
    {
        if (!results.has_column("add"))
            results.add_column("add");

        results.set_value("add", m_counter.m_add+m_counter_nested.m_add);

        if (!results.has_column("sub"))
            results.add_column("sub");

        results.set_value("sub", m_counter.m_subtract+m_counter_nested.m_subtract);

        if (!results.has_column("mul"))
            results.add_column("mul");

        results.set_value("mul", m_counter.m_multiply+m_counter_nested.m_multiply);

        if (!results.has_column("add_mul"))
            results.add_column("add_mul");

        results.set_value("add_mul",
                          m_counter.m_multiply_add + m_counter_nested.m_multiply_add);

        if (!results.has_column("sub_mul"))
            results.add_column("sub_mul");

        results.set_value("sub_mul",
                          m_counter.m_multiply_subtract + m_counter_nested.m_multiply_subtract);

        if (!results.has_column("invert"))
            results.add_column("invert");
        results.set_value("invert", m_counter.m_invert + m_counter_nested.m_invert);

        if (!results.has_column("symbols_used"))
            results.add_column("symbols_used");
        results.set_value("symbols_used", symbols_used);

        if (!results.has_column("sent_symbols"))
		   results.add_column("sent_symbols");
	   results.set_value("sent_symbols", sent_symbols);

		if (!results.has_column("time"))
			results.add_column("time");
		results.set_value("time", tm.measurement());
        if (!results.has_column("num_regions"))
            results.add_column("num_regions");
        results.set_value("num_regions", num_regions);

    	if (!results.has_column("ranks"))
    			results.add_column("ranks");
		results.set_value("ranks", m_ranks);

		if (!results.has_column("add_operations"))
			results.add_column("add_operations");
		results.set_value("add_operations", m_add_operations);

		if (!results.has_column("mul_operations"))
				results.add_column("mul_operations");
		results.set_value("mul_operations", m_mul_operations);

		if (!results.has_column("operations"))
			results.add_column("operations");
		results.set_value("operations", m_operations);

		if (!results.has_column("expansions"))
			results.add_column("expansions");
		results.set_value("expansions", m_expansions);

    }

    //NVu: add more
    void setup_factories()
    {
        //Super::setup_factories();
        gauge::config_set cs = get_current_configuration();

        uint32_t expansion = cs.get_value<uint32_t>("expansion");

        m_decoder_factory->set_expansion(expansion);
        m_encoder_factory->set_expansion(expansion);
        m_recoder_factory->set_expansion(expansion);
       

    }

    /// Prepares the measurement for every run
    void setup()
    {
	
        gauge::config_set cs = get_current_configuration();

        uint32_t symbols = cs.get_value<uint32_t>("symbols");
        uint32_t symbol_size = cs.get_value<uint32_t>("symbol_size");

        double density = cs.get_value<double>("density");
        double erasure = cs.get_value<double>("erasure");

        m_distribution = boost::random::bernoulli_distribution<>(erasure);

        m_decoder_factory = std::make_shared<decoder_factory>(
            symbols, symbol_size);

        m_encoder_factory = std::make_shared<encoder_factory>(
            symbols, symbol_size);

        m_recoder_factory = std::make_shared<recoder_factory>(
                  symbols, symbol_size);


        m_decoder_factory->set_symbols(symbols);
        m_decoder_factory->set_symbol_size(symbol_size);

        m_encoder_factory->set_symbols(symbols);
        m_encoder_factory->set_symbol_size(symbol_size);

        m_recoder_factory->set_symbols(symbols);
        m_recoder_factory->set_symbol_size(symbol_size);

        setup_factories(); //Nvu: add

        m_encoder = m_encoder_factory->build();
        //m_encoder->nested()->set_density(density);
        m_decoder = m_decoder_factory->build();
        m_recoder = m_recoder_factory->build();

        m_payload_buffer.resize(m_encoder->payload_size(), 0);
        m_encoded_data.resize(m_encoder->block_size(), 'x');
        m_encoder->set_const_symbols(storage::storage(m_encoded_data));
        symbols_used = 0;
        sent_symbols = 0;

        m_ranks.resize(symbols);
        std::fill(m_ranks.begin(), m_ranks.end(), 0);

        m_operations.resize(symbols);
		std::fill(m_operations.begin(), m_operations.end(), 0);

		 m_add_operations.resize(symbols);
				std::fill(m_add_operations.begin(), m_add_operations.end(), 0);
		 m_mul_operations.resize(symbols);
						std::fill(m_mul_operations.begin(), m_mul_operations.end(), 0);

		m_expansions.resize(symbols);
		std::fill(m_expansions.begin(), m_expansions.end(), 0);


    }
    double calculate_region(uint32_t symbols, uint32_t num_regions, const uint32_t exp = 0)
    {
    	return (symbols) * (pow(2, num_regions) - 1)/ pow(2, num_regions);
    }

    uint32_t get_add_operations()
            {
          	m_counter_nested = kodo_core::operations_counter();

          	m_counter = m_decoder->get_operations_counter();

                uint32_t sum_operations = 0;

                sum_operations = m_counter.m_add+m_counter_nested.m_add +
                				 m_counter.m_subtract+m_counter_nested.m_subtract;

                return sum_operations;
            }
       uint32_t get_mul_operations()
            {
   			m_counter_nested = kodo_core::operations_counter();

   			m_counter = m_decoder->get_operations_counter();

                uint32_t sum_operations = 0;

                sum_operations = m_counter.m_multiply+m_counter_nested.m_multiply +
        						 m_counter.m_multiply_add + m_counter_nested.m_multiply_add +
        						 m_counter.m_multiply_subtract + m_counter_nested.m_multiply_subtract+
        						 m_counter.m_invert + m_counter_nested.m_invert;

                return sum_operations;
            }
       /// Where the actual measurement is performed

    uint32_t get_operations()
      {
    	return get_add_operations() + get_mul_operations();
      }

   void run_fulcrum()
    {
        // Check with the result memory whether we already have
        // results for this configuration
        gauge::config_set cs = get_current_configuration();
        uint32_t expansion = cs.get_value<uint32_t>("expansion");
        uint32_t symbols = cs.get_value<uint32_t>("symbols");
        // We switch any systematic operations off so we code
        // symbols from the beginning
        if (kodo_core::has_set_systematic_off<Encoder>::value)
            kodo_core::set_systematic_off(*m_encoder);

        m_encoder->nested()->set_num_nonzero_expansion(expansion, expansion, false, false, true);

        RUN
        {
        	uint32_t pre_rank = 0;

            while (!m_decoder->is_complete())
            {

                // Encode a packet into the payload buffer
                m_encoder->write_payload(m_payload_buffer.data());
                sent_symbols +=1;

                if(m_distribution(m_random_generator))
				{
						continue;
				}

                m_decoder->read_payload(m_payload_buffer.data());
//                std::cout<<"Symbol_used:"<<symbols_used<<std::endl;
                uint32_t r = m_decoder->rank();
				if(r > symbols) r = symbols;
				++m_ranks[r-1];
				m_operations[r-1] = get_operations();
				m_add_operations[r-1] = get_add_operations();
				m_mul_operations[r-1] = get_mul_operations();
            	symbols_used++;

            }

        }
    }

    /// Where the actual measurement is performed
    void run_dynamic_tep()
         {
             // Check with the result memory whether we already have
             // results for this configuration
             gauge::config_set cs = get_current_configuration();
             uint32_t symbols = cs.get_value<uint32_t>("symbols");
             uint32_t expansion = cs.get_value<uint32_t>("expansion");
             double density = cs.get_value<double>("density");

             // We switch any systematic operations off so we code
             // symbols from the beginning
             if (kodo_core::has_set_systematic_off<Encoder>::value)
                 kodo_core::set_systematic_off(*m_encoder);


             RUN
             {
             		uint32_t cur_r = 0;
     			num_regions = 1;
     			double size_cur_region = 0.0;
     			m_encoder->nested()->set_num_nonzero_expansion(cur_r, expansion);

     			m_encoder->nested()->set_density(density);

     			size_cur_region = std::min(double(symbols-expansion + cur_r), calculate_region(symbols, num_regions));
    // 			std::cout<<"new region:"<<num_regions<<" with size:"<<size_cur_region<<std::endl;
     			//std::cout<<"new region:"<<num_regions<<" with size:"<<size_cur_region<<std::endl;
     			while (!m_decoder->is_complete())
     			{


     				if ((sent_symbols >= size_cur_region) && (cur_r < expansion))
    				{


    					if(cur_r < expansion)
    					{
    						cur_r += 1;
    						m_encoder->nested()->set_num_nonzero_expansion(cur_r, expansion, true);
    //						std::cout<<"introduce new coefficient:"<<cur_r<<std::endl;
    //						std::cout<<"current coded packets:"<<sent_symbols<<std::endl;

    					}
    					num_regions += 1;
    					size_cur_region = std::min(double(symbols-expansion + cur_r),
    							calculate_region(symbols, num_regions));
    					//std::cout<<"new region:"<<num_regions<<" with size:"<<size_cur_region<<std::endl;


    				}
//     				if (sent_symbols >= symbols)
//     					m_encoder->nested()->set_num_nonzero_expansion(cur_r, expansion, false, true);

     				// Encode a packet into the payload buffer
     				m_encoder->write_payload(m_payload_buffer.data());
     				sent_symbols += 1;


     				if(m_distribution(m_random_generator))
					{
							continue;
					}


     				m_decoder->read_payload(m_payload_buffer.data());


     				uint32_t r = m_decoder->rank();
     				if(r > symbols) r = symbols;
     				++m_ranks[r-1];
     				m_operations[r-1] = get_operations();
     				m_add_operations[r-1] = get_add_operations();
					m_mul_operations[r-1] = get_mul_operations();
     				m_expansions[r-1] = cur_r;
     				symbols_used++;

                 }

             }
         }
        void run_static_tep()
        {
            // Check with the result memory whether we already have
            // results for this configuration
            gauge::config_set cs = get_current_configuration();
            uint32_t symbols = cs.get_value<uint32_t>("symbols");
            uint32_t expansion = cs.get_value<uint32_t>("expansion");
            uint32_t threshold = cs.get_value<uint32_t>("threshold");
            double density = cs.get_value<double>("density");
            // We switch any systematic operations off so we code
            // symbols from the beginning
            if (kodo_core::has_set_systematic_off<Encoder>::value)
                kodo_core::set_systematic_off(*m_encoder);


			//if (symbols > 64) threshold = 6;

            RUN
            {

				uint32_t cur_r = 0;

				m_encoder->nested()->set_density(density);

    			m_encoder->nested()->set_num_nonzero_expansion(cur_r, expansion);


    			while (!m_decoder->is_complete())
    			{

    				//using like-systematic the extra coefficient, random the one when send extra packets
    				if (cur_r < expansion)
					{
						if (symbols - threshold <= sent_symbols)
						{
							cur_r = expansion;
							m_encoder->nested()->set_num_nonzero_expansion(cur_r, expansion, true);
							//std::cout<<"r = "<<cur_r<<"--> symbols sent:"<<sent_symbols<<std::endl;

						}
						else{
							//to avoid cur_r being never increased when symbols - threshold - expansion < 0
							int sub = symbols - threshold - expansion;
							if ((sub < 0) || (symbols - threshold - expansion <= sent_symbols))
							{
								cur_r += 1;
								m_encoder->nested()->set_num_nonzero_expansion(cur_r, expansion, true);
//								std::cout<<"r = "<<cur_r<<"--> symbols sent:"<<sent_symbols<<std::endl;

							}
						}


					}
//    				if (sent_symbols >= symbols)
//    					m_encoder->nested()->set_num_nonzero_expansion(cur_r, expansion, false, true);

    				// Encode a packet into the payload buffer
    				m_encoder->write_payload(m_payload_buffer.data());
    				sent_symbols += 1;

                    if(m_distribution(m_random_generator))
    				{
    						continue;
    				}

    				m_decoder->read_payload(m_payload_buffer.data());

     			 	
     			 	uint32_t r = m_decoder->rank();
    				if(r > symbols) r = symbols;
    				++m_ranks[r-1];
    				m_operations[r-1] = get_operations();
    				m_add_operations[r-1] = get_add_operations();
					m_mul_operations[r-1] = get_mul_operations();
    				m_expansions[r-1] = cur_r;

    				symbols_used++;

                }

            }
        }


	void get_options(gauge::po::variables_map& options)
	{
		auto symbols = options["symbols"].as<std::vector<uint32_t> >();
		auto symbol_size = options["symbol_size"].as<std::vector<uint32_t> >();
		auto types = options["type"].as<std::vector<std::string> >();

		auto expansion = options["expansion"].as<std::vector<uint32_t> >();

		auto erasure = options["erasure"].as<std::vector<double> >();
		auto density = options["density"].as<std::vector<double> >();

		auto threshold = options["threshold"].as<std::vector<uint32_t> >();

		assert(symbols.size() > 0);
		assert(symbol_size.size() > 0);
		assert(types.size() > 0);

		assert(expansion.size() > 0);

		assert(erasure.size() > 0);
		assert(density.size() > 0);


		for (uint32_t i = 0; i < symbols.size(); ++i)
		{
			for (uint32_t j = 0; j < symbol_size.size(); ++j)
			{
				for (uint32_t u = 0; u < types.size(); ++u)
				{
						for (const auto& e: expansion)
						{
							for (const auto& er: erasure)
							{
								for (const auto& den: density)
								{
									for (const auto& th: threshold)
									{
										gauge::config_set cs;
										cs.set_value<uint32_t>("symbols", symbols[i]);
										cs.set_value<uint32_t>("symbol_size", symbol_size[j]);
										cs.set_value<std::string>("type", types[u]);
										cs.set_value<uint32_t>("expansion", e);
										cs.set_value<uint32_t>("threshold", th);

										cs.set_value<double>("erasure", er);
										cs.set_value<double>("density", den);

										add_configuration(cs);
									}
								}

							}
						}

				}
			}
		}
	}

    /// The unit we measure in
    std::string unit_text() const
    { return "operations per symbol"; }

protected:

    /// The decoder factory
    std::shared_ptr<decoder_factory> m_decoder_factory;


    /// The recoder factory
    std::shared_ptr<recoder_factory> m_recoder_factory;

    /// The encoder factory
    std::shared_ptr<encoder_factory> m_encoder_factory;

    /// The decoder
    decoder_ptr m_decoder;

    /// The Encoder
    encoder_ptr m_encoder;

    /// The Recoder
    recoder_ptr m_recoder;

    /// The payload buffer
    std::vector<uint8_t> m_payload_buffer;

    /// The payload buffer
    std::vector<uint8_t> m_encoded_data;

    /// The counter containing the measurement results
    kodo_core::operations_counter m_counter;
    kodo_core::operations_counter m_counter_nested;

    // The distribution wrapping the random generator
//    boost::random::bernoulli_distribution<> m_distribution1;
//
//    // The distribution wrapping the random generator
//    boost::random::bernoulli_distribution<> m_distribution2;
    uint32_t symbols_used = 0;
    uint32_t sent_symbols = 0;

    uint32_t num_regions = 0;
    time_measurement tm;

    std::vector<uint32_t> m_operations;
    std::vector<uint32_t> m_add_operations;
    std::vector<uint32_t> m_mul_operations;
	std::vector<uint32_t> m_ranks;
	std::vector<uint32_t> m_expansions;

	boost::random::mt19937 m_random_generator;

	boost::random::bernoulli_distribution<> m_distribution;


};
