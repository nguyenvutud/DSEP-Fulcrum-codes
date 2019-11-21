// Copyright Steinwurf ApS 2011.
// Distributed under the "STEINWURF RESEARCH LICENSE 1.0".
// See accompanying file LICENSE.rst or
// http://www.steinwurf.com/licensing

#pragma once

#include <cstdint>
#include <cassert>
#include <vector>
#include <string>
#include <type_traits>

#include <gauge/gauge.hpp>
#include <tables/table.hpp>

#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <fifi/is_prime2325.hpp>
#include <fifi/prime2325_binary_search.hpp>
#include <fifi/prime2325_apply_prefix.hpp>

#include <kodo_core/has_deep_symbol_storage.hpp>
#include <kodo_core/has_set_mutable_symbols.hpp>
#include <kodo_core/set_systematic_off.hpp>
#include <kodo_core/set_mutable_symbols.hpp>
#include <kodo_core/read_payloads.hpp>
#include <kodo_core/write_payloads.hpp>

/// A test block represents an encoder and decoder pair
template<class Encoder, class Decoder>
struct dynamic_sparse_tunable_throughput_benchmark : public gauge::time_benchmark
{
    using encoder_factory = typename Encoder::factory;
    using encoder_ptr = typename Encoder::factory::pointer;

    using decoder_factory = typename Decoder::factory;
    using decoder_ptr = typename Decoder::factory::pointer;
    dynamic_sparse_tunable_throughput_benchmark(){

        // Seed the random generator controlling the erasures
        m_random_generator.seed((uint32_t)time(0));
    }
    void init()
    {
        m_factor = 1;
        gauge::time_benchmark::init();
    }

    void start()
    {
        m_encoded_symbols = 0;
 
        m_recovered_symbols = 0;
        m_processed_symbols = 0;
        gauge::time_benchmark::start();
    }

    void stop()
    {
        gauge::time_benchmark::stop();
    }

    double calculate_throughput(bool goodput = false)
    {
        // Get the time spent per iteration
        m_time = gauge::time_benchmark::measurement();

        gauge::config_set cs = get_current_configuration();
        std::string type = cs.get_value<std::string>("type");
        //uint32_t symbols = cs.get_value<uint32_t>("symbols");
        uint32_t symbol_size = cs.get_value<uint32_t>("symbol_size");

        // The number of bytes {en|de}coded
        uint64_t total_bytes = 0;

        if (type == "decoder")
        {
            if (goodput)
                total_bytes = m_recovered_symbols * symbol_size;
            else
                total_bytes = m_processed_symbols * symbol_size;
        }
        else if (type == "encoder")
        {
            total_bytes = m_encoded_symbols * symbol_size;
        }
        else
        {
            assert(0 && "Unknown benchmark type");
        }

        // The bytes per iteration
        uint64_t bytes =
            total_bytes / gauge::time_benchmark::iteration_count();

        return bytes / m_time; // MB/s for each iteration
    }

    void store_run(tables::table& results)
    {
        if (!results.has_column("throughput"))
            results.add_column("throughput");

        if (!results.has_column("goodput"))
            results.add_column("goodput");

        results.set_value("throughput", calculate_throughput(false));
        results.set_value("goodput", calculate_throughput(true));
	
	 if (!results.has_column("used"))
            results.add_column("used");
        results.set_value("used", m_used);
//        if (!results.has_column("time"))
//                   results.add_column("time");
//               results.set_value("time", m_time);

    }

    bool accept_measurement()
    {
        gauge::config_set cs = get_current_configuration();

        std::string type = cs.get_value<std::string>("type");

        if (type == "decoder")
        {
            // If we are benchmarking a decoder we only accept
            // the measurement if the decoding was successful
            if (!m_decoder->is_complete())
            {
                // We did not generate enough payloads to decode successfully,
                // so we will generate more payloads for next run
                ++m_factor;

                return false;
            }

            // At this point, the output data should be equal to the input
            // data
            assert(m_data_out == m_data_in);
        }

        return gauge::time_benchmark::accept_measurement();
    }

    std::string unit_text() const
    {
        return "MB/s";
    }

    void get_options(gauge::po::variables_map& options)
    {
        auto symbols = options["symbols"].as<std::vector<uint32_t>>();
        auto symbol_size = options["symbol_size"].as<std::vector<uint32_t>>();
        auto types = options["type"].as<std::vector<std::string>>();
        auto expansion = options["expansion"].as<std::vector<uint32_t> >();
        auto erasure = options["erasure"].as<std::vector<double> >();

        auto redundant = options["redundant"].as<std::vector<uint32_t> >();

        assert(symbols.size() > 0);
        assert(symbol_size.size() > 0);
        assert(types.size() > 0);
        assert(expansion.size() > 0);

        assert(erasure.size() > 0);

        for (const auto& s : symbols)
        {
            for (const auto& p : symbol_size)
            {
                for (const auto& t : types)
                {

                    for (const auto& e: expansion)
                    {
					   for (const auto& er: erasure)
						{
						    for (const auto& red: redundant)
						    {
								gauge::config_set cs;
								cs.set_value<uint32_t>("symbols", s);
								cs.set_value<uint32_t>("symbol_size", p);
								cs.set_value<std::string>("type", t);
								cs.set_value<uint32_t>("expansion", e);
								cs.set_value<uint32_t>("redundant", red);

								cs.set_value<double>("erasure", er);

								add_configuration(cs);
						    }
						}
                    }

                }
            }
        }
    }

    void setup()
    {
    	gauge::config_set cs = get_current_configuration();

		uint32_t symbols = cs.get_value<uint32_t>("symbols");
		uint32_t symbol_size = cs.get_value<uint32_t>("symbol_size");
		//uint32_t sparse = cs.get_value<uint32_t>("sparse");

//		double density = cs.get_value<double>("density");

		m_decoder_factory = std::make_shared<decoder_factory>(
			symbols, symbol_size);

		m_encoder_factory = std::make_shared<encoder_factory>(
			symbols, symbol_size);


		m_decoder_factory->set_symbols(symbols);
		m_decoder_factory->set_symbol_size(symbol_size);

		m_encoder_factory->set_symbols(symbols);
		m_encoder_factory->set_symbol_size(symbol_size);


		setup_factories(); //Nvu: add

		m_encoder = m_encoder_factory->build();


		m_decoder = m_decoder_factory->build();

		m_payload_buffer.resize(m_encoder->payload_size(), 0);

		m_data_in.resize(m_encoder->block_size());

	    std::generate_n(m_data_in.begin(), m_data_in.size(), rand);

		m_encoder->set_const_symbols(storage::storage(m_data_in));

		m_data_out.resize(m_encoder->block_size());

		// Make sure that the encoder and decoder payload sizes are equal
		assert(m_encoder->payload_size() == m_decoder->payload_size());

		// Create the payload buffer
		uint32_t payload_size = m_encoder->payload_size();

		m_temp_payload.resize(payload_size);

		// Prepare storage to the encoded payloads
		uint32_t payload_count = symbols * m_factor;

		assert(payload_count > 0);

		// Allocate contiguous payload buffer and store payload pointers
		m_payload_buffer.resize(payload_count * payload_size);
		m_payloads.resize(payload_count);

		for (uint32_t i = 0; i < payload_count; ++i)
		{
			m_payloads[i] = &m_payload_buffer[i * payload_size];
		}

		m_used = 0;
    }

    /// Setup the factories - we re-factored this into its own function
    /// since some of the codecs like sparse, perpetual and fulcrum
    /// customize the setup of the factory.
    void setup_factories()
      {
          //Super::setup_factories();
          gauge::config_set cs = get_current_configuration();

          uint32_t expansion = cs.get_value<uint32_t>("expansion");

          m_decoder_factory->set_expansion(expansion);
          m_encoder_factory->set_expansion(expansion);


      }

    int calculate_region(uint32_t symbols, uint32_t num_regions, uint32_t red)
     {
     	return (int) ((symbols) * (pow(2, num_regions) - 1)/ pow(2, num_regions));
     }


    void encode_payloads_dynamic(uint32_t symbols, uint32_t expansion, uint32_t redundant)
            {
                // If we are working in the prime2325 finite field we have to
                // "map" the data first. This cost is included in the encoder
                // throughput (we do it with the clock running), just as it
                // would be in a real application.
                if (fifi::is_prime2325<typename Encoder::field_type>::value)
                {
                    uint32_t block_length =
                        fifi::size_to_length<fifi::prime2325>(m_encoder->block_size());

                    fifi::prime2325_binary_search search(block_length);
                    m_prefix = search.find_prefix(storage::storage(m_data_in));

                    // Apply the negated prefix
                    fifi::apply_prefix(storage::storage(m_data_in), ~m_prefix);
                }

                m_encoder->set_const_symbols(storage::storage(m_data_in));

                // We switch any systematic operations off so we code
                // symbols from the beginning
                if (kodo_core::has_set_systematic_off<Encoder>::value)
                    kodo_core::set_systematic_off(*m_encoder);


                uint32_t cur_sparse = 5, received_symbols = 0;
				num_regions = 1;
				int size_cur_region = 0;
//				m_encoder->nested()->set_num_nonzero_expansion(cur_r, expansion);

				size_cur_region =  calculate_region(symbols + expansion, num_regions, redundant);

//				std::cout<<"size of region:"<<size_cur_region<<std::endl;

				m_encoder->nested()-> set_average_nonzero_symbols(K[redundant/5-1][0]);
//				std::cout<<"size of region:"<<size_cur_region<<std::endl;
//				std::cout<<"K = "<<K[redundant/5-1][0]<<std::endl;

				uint32_t sum = 0;

				for (uint8_t* payload : m_payloads)
				{
					if ((m_encoded_symbols >= size_cur_region))
					{
						num_regions += 1;

						if(num_regions < 9){
							size_cur_region = calculate_region(symbols+expansion, num_regions, redundant);
//							std::cout<<"size of region:"<<size_cur_region<<std::endl;

							m_encoder->nested()-> set_average_nonzero_symbols(K[redundant/5-1][num_regions-1]);
//							std::cout<<"K = "<<K[redundant%5][num_regions-1]<<std::endl;
						}

					}

					++m_encoded_symbols;

					m_encoder->write_payload(payload);
//					std::cout<<"Number of nonzero generated:"<<m_encoder->nested()->nonzeros_generated()<<std::endl;
//					std::cout<<"Encoded packets:"<<m_encoded_symbols<<std::endl;

					sum = sum + m_encoder->nested()->nonzeros_generated();
				}
//				std::cout<<"Avg of nonzeros:"<<sum/(m_encoded_symbols*1.0)<<std::endl;



			if (fifi::is_prime2325<typename Encoder::field_type>::value)
			{
				// Apply the negated prefix
				fifi::apply_prefix(storage::storage(m_data_in), ~m_prefix);
			}
            }



    void decode_payloads()
    {
        // If the decoder uses shallow storage, we have to initialize
        // its decoding buffer
        if (kodo_core::has_set_mutable_symbols<Decoder>::value)
        {
            kodo_core::set_mutable_symbols(*m_decoder, storage::storage(m_data_out));
        }

            for (const uint8_t* payload : m_payloads)
            {
         		std::copy_n(payload, m_decoder->payload_size(),
                            m_temp_payload.data());

                m_decoder->read_payload(m_temp_payload.data());

                ++m_processed_symbols;
                ++m_used;

                //std::cout<<"Decoder data:" << m_processed_symbols <<std::endl;

                if (m_decoder->is_complete())
                {
                	break;
                }
            }


        if (m_decoder->is_complete())
        {
            // If the decoder uses deep storage, we have to copy
            // its decoding buffers for verification.
            // This would also be necessary in real applications, so it
            // should be part of the timed benchmark.
            if (kodo_core::has_deep_symbol_storage<Decoder>::value)
            {
                m_decoder->copy_from_symbols(storage::storage(m_data_out));
            }

            if (fifi::is_prime2325<typename Decoder::field_type>::value)
            {
                // Now we have to apply the negated prefix to the decoded data
                fifi::apply_prefix(storage::storage(m_data_out), ~m_prefix);
            }

            m_recovered_symbols += m_decoder->symbols();

            //std::cout<<"Decode completely:" <<m_recovered_symbols <<std::endl;
        }
    }

    ///FULCRUM DYNAMIC SPARSE
    /// Run the encoder
    void run_fulcrum_encode_dynamic()
    {
  	  gauge::config_set cs = get_current_configuration();
  	  uint32_t symbols = cs.get_value<uint32_t>("symbols");

  	  uint32_t expansion = cs.get_value<uint32_t>("expansion");

  	  uint32_t redundant = cs.get_value<uint32_t>("redundant");

  	  if (kodo_core::has_set_systematic_off<Encoder>::value)
  		  kodo_core::set_systematic_off(*m_encoder);

        // The clock is running
        RUN
        {
            // We have to make sure the encoder is in a "clean" state
            m_encoder->initialize(*m_encoder_factory);

//		    m_encoder->nested()->set_number_nonzero_coefficient(sparse);
            //m_encoder->nested()->set_average_nonzero_symbols(sparse);
            m_encoder->nested()->set_num_nonzero_expansion(expansion, expansion, false, false, true);

            encode_payloads_dynamic(symbols, expansion, redundant);
        }
    }

    /// Run the decoder
    void run_fulcrum_decode_dynamic()
    {
    	gauge::config_set cs = get_current_configuration();
		uint32_t symbols = cs.get_value<uint32_t>("symbols");
		uint32_t expansion = cs.get_value<uint32_t>("expansion");
		uint32_t redundant = cs.get_value<uint32_t>("redundant");

		m_encoder->initialize(*m_encoder_factory);

		//m_encoder->nested()->set_average_nonzero_symbols(sparse);
		m_encoder->nested()->set_num_nonzero_expansion(expansion, expansion, false, false, true);

		if (kodo_core::has_set_systematic_off<Encoder>::value)
		  kodo_core::set_systematic_off(*m_encoder);
        // Encode some data
        encode_payloads_dynamic(symbols, expansion, redundant);

        // Zero the data buffer for the decoder
        std::fill_n(m_data_out.begin(), m_data_out.size(), 0);

        // The clock is running
        RUN
        {
            // We have to make sure the decoder is in a "clean" state
            // i.e. no symbols already decoded.
            m_decoder->initialize(*m_decoder_factory);

            // Decode the payloads
            decode_payloads();
        }
    }

    void run_body_fulcrum_dynamic()
    {
        gauge::config_set cs = get_current_configuration();

        std::string type = cs.get_value<std::string>("type");

        if (type == "encoder")
        {
        	run_fulcrum_encode_dynamic();
        }
        else if (type == "decoder")
        {
            run_fulcrum_decode_dynamic();
        }
        else
        {
            assert(0 && "Unknown benchmark type");
        }
    }



protected:

    /// The decoder factory
    std::shared_ptr<decoder_factory> m_decoder_factory;

    /// The encoder factory
    std::shared_ptr<encoder_factory> m_encoder_factory;


    /// The encoder to use
    encoder_ptr m_encoder;

    /// The number of encoded symbols
    uint32_t m_encoded_symbols;


    /// The decoder to use
    decoder_ptr m_decoder;

    /// The number of symbols recovered by the decoder
    uint32_t m_recovered_symbols;

    /// The number of symbols processed by the decoder
    uint32_t m_processed_symbols;

    /// The input data
    std::vector<uint8_t> m_data_in;

    /// The output data
    std::vector<uint8_t> m_data_out;

    /// Buffer storing the encoded payloads before adding them to the
    /// decoder. This is necessary since the decoder will "work on"
    /// the encoded payloads directly. Therefore if we want to be able
    /// to run multiple iterations with the same encoded paylaods we
    /// have to copy them before injecting them into the decoder. This
    /// of course has a negative impact on the decoding throughput.
    std::vector<uint8_t> m_temp_payload;


    /// Contiguous buffer for coded payloads
    std::vector<uint8_t> m_payload_buffer;


    /// Pointers to each payload in the payload buffer
    std::vector<uint8_t*> m_payloads;

    /// Multiplication factor for payload_count
    uint32_t m_factor;

    /// Prefix used when testing with the prime2325 finite field
    uint32_t m_prefix;

    boost::random::mt19937 m_random_generator;

    boost::random::bernoulli_distribution<> m_distribution;

    uint32_t m_used = 0;

    uint32_t num_regions = 0;

    double m_time;
//    uint32_t K[4][8] ={{6, 10, 19, 33, 48, 67, 54, 54},
//    				   {4, 8, 15, 25, 40, 58, 57, 57},
//    				   {3, 7, 13, 21, 35, 54, 59, 59},
//    				   {3, 6, 11, 19, 32, 47, 62, 62}};
//New data: row 0: min, row 1, 3: estimate, row 2: max
    uint32_t K[4][8] ={	{2, 5, 9, 18, 31, 47, 52, 52},
    					{3, 6, 12, 22, 35, 52, 52, 52},
    					{5, 9, 17, 27, 40, 52, 52, 52},
			   	   	   {6, 10, 19, 33, 48, 67, 54, 54}
      				   };
};
