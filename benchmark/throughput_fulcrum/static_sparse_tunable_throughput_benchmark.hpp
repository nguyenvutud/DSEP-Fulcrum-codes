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
struct static_sparse_tunable_throughput_benchmark : public gauge::time_benchmark
{
    using encoder_factory = typename Encoder::factory;
    using encoder_ptr = typename Encoder::factory::pointer;

    using decoder_factory = typename Decoder::factory;
    using decoder_ptr = typename Decoder::factory::pointer;
    static_sparse_tunable_throughput_benchmark(){

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
        //std::cout<<"START COUTING\n";
    }

    void stop()
    {
        gauge::time_benchmark::stop();
        //std::cout<<"STOP COUTING\n";
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

	if (!results.has_column("num_trans"))
			   results.add_column("num_trans");
		   results.set_value("num_trans", m_transmissions);

     if (!results.has_column("dependency"))
			results.add_column("dependency");
	 results.set_value("dependency", m_dependency);

	 if (!results.has_column("nonzero_coefficient"))
			results.add_column("nonzero_coefficient");
	 results.set_value("nonzero_coefficient", m_nonzero_coefficient);
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

        auto sparse = options["sparse"].as<std::vector<uint32_t> >();
        auto density = options["density"].as<std::vector<double> >();

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
                	for (const auto& d: density){
                    for (const auto& e: expansion)
                    {
					   for (const auto& er: erasure)
						{
						    for (const auto& sp: sparse)
						    {
								gauge::config_set cs;
								cs.set_value<uint32_t>("symbols", s);
								cs.set_value<uint32_t>("symbol_size", p);
								cs.set_value<std::string>("type", t);
								cs.set_value<uint32_t>("expansion", e);
								cs.set_value<uint32_t>("sparse", sp);
								cs.set_value<double>("density", d);
								cs.set_value<double>("erasure", er);

								add_configuration(cs);
						    }
						}
                    }

                	}
                }
            }
        }
    }

    void setup()
    {
    	//std::cout<<"START SETUP\n";
    	gauge::config_set cs = get_current_configuration();

		uint32_t symbols = cs.get_value<uint32_t>("symbols");
		uint32_t symbol_size = cs.get_value<uint32_t>("symbol_size");
		uint32_t expansion = cs.get_value<uint32_t>("expansion");
		//uint32_t sparse = cs.get_value<uint32_t>("sparse");

		double density = cs.get_value<double>("density");
		double erasure = cs.get_value<double>("erasure");

		m_distribution = boost::random::bernoulli_distribution<>(erasure);

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
		//set density
//		m_encoder->nested()->set_density(density);

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

		m_dependency.resize(symbols + 10 + expansion+1);
		std::fill(m_dependency.begin(), m_dependency.end(), 0);



		// Prepare storage to the encoded payloads
		uint32_t payload_count = symbols * m_factor;

		assert(payload_count > 0);

		// Allocate contiguous payload buffer and store payload pointers
		m_payload_buffer.resize(payload_count * payload_size);
		m_payloads.resize(payload_count);

		m_nonzero_coefficient.resize(symbols+expansion);
		std::fill(m_nonzero_coefficient.begin(), m_nonzero_coefficient.end(), 0);


		for (uint32_t i = 0; i < payload_count; ++i)
		{
			m_payloads[i] = &m_payload_buffer[i * payload_size];
		}

		m_used = 0;
		m_transmissions = 0;
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
    void encode_payloads_with_step(uint32_t symbols, uint32_t expansion)
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

            //use for normal fulcrum decoder
           m_encoder->nested()->set_num_nonzero_expansion(0, expansion, false, false, true);

            // We switch any systematic operations off so we code
            // symbols from the beginning
            if (kodo_core::has_set_systematic_off<Encoder>::value)
                kodo_core::set_systematic_off(*m_encoder);

            uint32_t g = 0;

    		for (uint8_t* payload : m_payloads)
    		{
    			//aim to set the same extra packets to other schemes
    			if (m_encoded_symbols >= symbols - 10)
    				m_encoder->nested()->set_num_nonzero_expansion(expansion, expansion, false, false, true);

    			m_encoder->write_payload(payload);

    			if (g < symbols + expansion){
    				m_nonzero_coefficient[g] = m_encoder->nested()->nonzeros_generated();
//    				std::cout<<"nonzero coefficient at "<<g<<" is: "<<m_nonzero_coefficient[g]<<std::endl;
    				g++;
    			}
    			++m_encoded_symbols;

    		}

            if (fifi::is_prime2325<typename Encoder::field_type>::value)
            {
                // Apply the negated prefix
                fifi::apply_prefix(storage::storage(m_data_in), ~m_prefix);
            }
        }
    void encode_payloads(uint32_t symbols, uint32_t expansion)
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

        //use for normal fulcrum decoder
       m_encoder->nested()->set_num_nonzero_expansion(expansion, expansion, false, false, true);

        // We switch any systematic operations off so we code
        // symbols from the beginning
        if (kodo_core::has_set_systematic_off<Encoder>::value)
            kodo_core::set_systematic_off(*m_encoder);

        uint32_t g = 0;

		for (uint8_t* payload : m_payloads)
		{
			m_encoder->write_payload(payload);

			if (g < symbols + expansion){
				m_nonzero_coefficient[g] = m_encoder->nested()->nonzeros_generated();
//				std::cout<<"nonzero coefficient at "<<g<<" is: "<<m_nonzero_coefficient[g]<<std::endl;
				g++;
			}
			++m_encoded_symbols;

		}

        if (fifi::is_prime2325<typename Encoder::field_type>::value)
        {
            // Apply the negated prefix
            fifi::apply_prefix(storage::storage(m_data_in), ~m_prefix);
        }
        //std::cout<<"FINISH ENCODING:"<<m_payloads.size()<<std::endl;
    }
    double calculate_region(uint32_t symbols, uint32_t num_regions)
    {
    	return symbols * (pow(2, num_regions) - 1)/ pow(2, num_regions);
    }



    void decode_payloads()
    {
    	gauge::config_set cs = get_current_configuration();
    	uint32_t expansion = cs.get_value<uint32_t>("expansion");

        // If the decoder uses shallow storage, we have to initialize
        // its decoding buffer
        if (kodo_core::has_set_mutable_symbols<Decoder>::value)
        {
            kodo_core::set_mutable_symbols(*m_decoder, storage::storage(m_data_out));
        }
        uint32_t dependent_counter = 0, pre_rank = 0;
		for (const uint8_t* payload : m_payloads)
		{
			std::copy_n(payload, m_decoder->payload_size(),
						m_temp_payload.data());
			++m_transmissions;
			if (m_distribution(m_random_generator))
			{
//				std::cout<<"LOST HAPPPED\n";
				continue;
			}

			m_decoder->read_payload(m_temp_payload.data());

			++m_processed_symbols;
			++m_used;

			uint32_t symbols = m_decoder->symbols();

			///count the dependency following rank
//			uint32_t r = m_decoder->rank();
//			if(m_decoder->rank() == pre_rank)
//				++m_dependency[r];
//			else pre_rank = m_decoder->rank();

			///count the dependency following transmissions this use to calculate prob of innovative coded packets
			dependent_counter +=1;
			if(dependent_counter < symbols + expansion + 10)
			{
				if(m_decoder->rank() == pre_rank)
					++m_dependency[dependent_counter];
				else pre_rank = m_decoder->rank();
//				std::cout<<"symbol used: " << m_used <<" is: "<<m_dependency[dependent_counter]<<std::endl;
			}



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
    ///FULCRUM STATIC SPARSE
     /// Run the encoder
     void run_fulcrum_encode_static()
     {
   	  gauge::config_set cs = get_current_configuration();
   	  uint32_t symbols = cs.get_value<uint32_t>("symbols");

   	  uint32_t expansion = cs.get_value<uint32_t>("expansion");

   	  uint32_t sparse = cs.get_value<uint32_t>("sparse");

   	  if (kodo_core::has_set_systematic_off<Encoder>::value)
   		  kodo_core::set_systematic_off(*m_encoder);

         // The clock is running
         RUN
         {
             // We have to make sure the encoder is in a "clean" state
             m_encoder->initialize(*m_encoder_factory);

 		    m_encoder->nested()->set_number_nonzero_coefficient(sparse);
             encode_payloads(symbols, expansion);
         }
     }

     /// Run the decoder
     void run_fulcrum_decode_static()
     {
     	gauge::config_set cs = get_current_configuration();
 		uint32_t symbols = cs.get_value<uint32_t>("symbols");
 		uint32_t expansion = cs.get_value<uint32_t>("expansion");
 		uint32_t sparse = cs.get_value<uint32_t>("sparse");
 		m_encoder->initialize(*m_encoder_factory);

 		 m_encoder->nested()->set_number_nonzero_coefficient(sparse);

 		if (kodo_core::has_set_systematic_off<Encoder>::value)
 		  kodo_core::set_systematic_off(*m_encoder);
         // Encode some data
         encode_payloads(symbols, expansion);

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

     void run_body_fulcrum_static()
     {
         gauge::config_set cs = get_current_configuration();

         std::string type = cs.get_value<std::string>("type");

         if (type == "encoder")
         {
         	run_fulcrum_encode_static();
         }
         else if (type == "decoder")
         {
             run_fulcrum_decode_static();
         }
         else
         {
             assert(0 && "Unknown benchmark type");
         }
     }

    ///FULCRUM AVERAGE SPARSE
    /// Run the encoder
    void run_fulcrum_encode_average()
    {
  	  gauge::config_set cs = get_current_configuration();
  	  uint32_t symbols = cs.get_value<uint32_t>("symbols");

  	  uint32_t expansion = cs.get_value<uint32_t>("expansion");

  	  uint32_t sparse = cs.get_value<uint32_t>("sparse");

  	  if (kodo_core::has_set_systematic_off<Encoder>::value)
  		  kodo_core::set_systematic_off(*m_encoder);

        // The clock is running
        RUN
        {
            // We have to make sure the encoder is in a "clean" state
            m_encoder->initialize(*m_encoder_factory);

//		    m_encoder->nested()->set_number_nonzero_coefficient(sparse);
            m_encoder->nested()->set_average_nonzero_symbols(sparse);
            m_encoder->nested()->set_num_nonzero_expansion(expansion, expansion, false, false, true);

            encode_payloads(symbols, expansion);
        }
    }

    /// Run the decoder
    void run_fulcrum_decode_average()
    {
    	gauge::config_set cs = get_current_configuration();
		uint32_t symbols = cs.get_value<uint32_t>("symbols");
		uint32_t expansion = cs.get_value<uint32_t>("expansion");
		uint32_t sparse = cs.get_value<uint32_t>("sparse");

		m_encoder->initialize(*m_encoder_factory);

		m_encoder->nested()->set_average_nonzero_symbols(sparse);
		m_encoder->nested()->set_num_nonzero_expansion(expansion, expansion, false, false, true);

		if (kodo_core::has_set_systematic_off<Encoder>::value)
		  kodo_core::set_systematic_off(*m_encoder);
        // Encode some data
        encode_payloads(symbols, expansion);

//        linear_dependency.resize(symbols+expansion, 0);

        // Zero the data buffer for the decoder
        std::fill_n(m_data_out.begin(), m_data_out.size(), 0);

        // The clock is running
        uint32_t num_runs = 0;
        RUN
        {
            // We have to make sure the decoder is in a "clean" state
            // i.e. no symbols already decoded.
            m_decoder->initialize(*m_decoder_factory);

            // Decode the payloads
            decode_payloads();
//            num_runs += 1;

        }
//        for(int i = 0; i< linear_dependency.size(); i++)
//			std::cout<<linear_dependency[i]<<" ";
//		std::cout<<std::endl;
//		std::cout<<"Runs:"<<num_runs<<std::endl;
    }

    void run_body_fulcrum_average()
    {
        gauge::config_set cs = get_current_configuration();

        std::string type = cs.get_value<std::string>("type");

        if (type == "encoder")
        {
        	run_fulcrum_encode_average();
        }
        else if (type == "decoder")
        {
            run_fulcrum_decode_average();
        }
        else
        {
            assert(0 && "Unknown benchmark type");
        }
    }

    ///ORIGINAL FULCRUM
    /// Run the encoder
     void run_fulcrum_encode()
     {
   	  gauge::config_set cs = get_current_configuration();
   	  uint32_t symbols = cs.get_value<uint32_t>("symbols");

   	  uint32_t expansion = cs.get_value<uint32_t>("expansion");

   	  double density = cs.get_value<double>("density");


   	  if (kodo_core::has_set_systematic_off<Encoder>::value)
   		  kodo_core::set_systematic_off(*m_encoder);

         // The clock is running
         RUN
         {
   		  	  //std::cout<<"START RUNNING ENCODER\n";
             // We have to make sure the encoder is in a "clean" state
             m_encoder->initialize(*m_encoder_factory);

             m_encoder->nested()->set_density(density);

             encode_payloads(symbols, expansion);
         }
     }

     /// Run the decoder
     void run_fulcrum_decode()
     {
     	gauge::config_set cs = get_current_configuration();
 		uint32_t symbols = cs.get_value<uint32_t>("symbols");
 		uint32_t expansion = cs.get_value<uint32_t>("expansion");

 		m_encoder->initialize(*m_encoder_factory);

 		double density = cs.get_value<double>("density");
        m_encoder->nested()->set_density(density);

 		if (kodo_core::has_set_systematic_off<Encoder>::value)
 		  kodo_core::set_systematic_off(*m_encoder);
         // Encode some data
         encode_payloads(symbols, expansion);

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

     void run_body_fulcrum()
     {
         gauge::config_set cs = get_current_configuration();

         std::string type = cs.get_value<std::string>("type");

         if (type == "encoder")
         {
         	run_fulcrum_encode();
         }
         else if (type == "decoder")
         {
             run_fulcrum_decode();
         }
         else
         {
             assert(0 && "Unknown benchmark type");
         }
     }
     void run_fulcrum_encode_with_step()
     {
   	  gauge::config_set cs = get_current_configuration();
   	  uint32_t symbols = cs.get_value<uint32_t>("symbols");

   	  uint32_t expansion = cs.get_value<uint32_t>("expansion");

   	  double density = cs.get_value<double>("density");


   	  if (kodo_core::has_set_systematic_off<Encoder>::value)
   		  kodo_core::set_systematic_off(*m_encoder);

         // The clock is running
         RUN
         {
             // We have to make sure the encoder is in a "clean" state
             m_encoder->initialize(*m_encoder_factory);

             m_encoder->nested()->set_density(density);

             encode_payloads(symbols, expansion);
         }
     }

     void run_fulcrum_decode_with_step()
     {
     	gauge::config_set cs = get_current_configuration();
 		uint32_t symbols = cs.get_value<uint32_t>("symbols");
 		uint32_t expansion = cs.get_value<uint32_t>("expansion");

 		m_encoder->initialize(*m_encoder_factory);

 		double density = cs.get_value<double>("density");
        m_encoder->nested()->set_density(density);

 		if (kodo_core::has_set_systematic_off<Encoder>::value)
 		  kodo_core::set_systematic_off(*m_encoder);

 		// Encode some data
         encode_payloads_with_step(symbols, expansion);

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

     void run_body_fulcrum_with_step()
     {
         gauge::config_set cs = get_current_configuration();

         std::string type = cs.get_value<std::string>("type");

         if (type == "encoder")
         {
         	run_fulcrum_encode_with_step();
         }
         else if (type == "decoder")
         {
             run_fulcrum_decode_with_step();
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

    uint32_t m_transmissions = 0;

    double m_time;
    std::vector<uint32_t> m_dependency;
    std::vector<uint32_t> m_nonzero_coefficient;
};
