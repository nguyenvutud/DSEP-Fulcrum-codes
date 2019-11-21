// Copyright Steinwurf ApS 2011.
// Distributed under the "STEINWURF RESEARCH LICENSE 1.0".
// See accompanying file LICENSE.rst or
// http://www.steinwurf.com/licensing

#include <cassert>
#include <ctime>
#include <memory>

#include <gauge/console_printer.hpp>
#include <gauge/csv_printer.hpp>
#include <gauge/gauge.hpp>
#include <gauge/python_printer.hpp>


#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <tables/table.hpp>

#include <kodo_core/set_systematic_off.hpp>

#include <kodo_fulcrum/fulcrum_codes.hpp>
//#include <kodo_fulcrum_sparse/set_density.hpp>
//#include <kodo_fulcrum_sparse/set_number_nonzero_coefficient.hpp>
//#include <kodo_fulcrum_sparse/has_set_number_nonzero_coefficient.hpp>
//#include <kodo_fulcrum_sparse/nested_set_number_nonzero_coefficient.hpp>

#include <kodo_core/operations_counter.hpp>

template<class Encoder, class Recoder, class Decoder>
class operations_benchmark_inner : public gauge::benchmark
{
public:

    using encoder_factory = typename Encoder::factory;
    using encoder_ptr = typename Encoder::factory::pointer;

    using decoder_factory = typename Decoder::factory;
    using decoder_ptr = typename Decoder::factory::pointer;

    using recoder_factory = typename Recoder::factory;
    using recoder_ptr = typename Recoder::factory::pointer;

    /// Starts a measurement and saves the counter
    void start()
    {
        m_encoder->reset_operations_counter();
        m_decoder->reset_operations_counter();
        m_recoder->reset_operations_counter();
        m_encoder->nested()->reset_operations_counter();
        m_decoder->nested()->reset_operations_counter();
        m_recoder->nested()->reset_operations_counter();
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
            m_counter_nested = m_decoder->nested()->get_operations_counter();
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
    }

    void store_run(tables::table& results)
    {
        if (!results.has_column("addition"))
            results.add_column("addition");

        results.set_value("addition", m_counter.m_add+m_counter_nested.m_add);

        if (!results.has_column("subtraction"))
            results.add_column("subtraction");

        results.set_value("subtraction", m_counter.m_subtract+m_counter_nested.m_subtract);

        if (!results.has_column("multiplication"))
            results.add_column("multiplication");

        results.set_value("multiplication", m_counter.m_multiply+m_counter_nested.m_multiply);

        if (!results.has_column("addition_multiplication"))
            results.add_column("addition_multiplication");

        results.set_value("addition_multiplication",
                          m_counter.m_multiply_add + m_counter_nested.m_multiply_add);

        if (!results.has_column("subtraction_multiplication"))
            results.add_column("subtraction_multiplication");

        results.set_value("subtraction_multiplication",
                          m_counter.m_multiply_subtract + m_counter_nested.m_multiply_subtract);

        if (!results.has_column("invert"))
            results.add_column("invert");
        results.set_value("invert", m_counter.m_invert + m_counter_nested.m_invert);
        if (!results.has_column("symbols_used"))
                  results.add_column("symbols_used");
        results.set_value("symbols_used", symbols_used);
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
        m_decoder = m_decoder_factory->build();
        m_recoder = m_recoder_factory->build();

        m_payload_buffer.resize(m_encoder->payload_size(), 0);
        m_encoded_data.resize(m_encoder->block_size(), 'x');
        m_encoder->set_const_symbols(storage::storage(m_encoded_data));
        symbols_used = 0;
    }

    /// Where the actual measurement is performed
    void run()
    {
        // Check with the result memory whether we already have
        // results for this configuration
        gauge::config_set cs = get_current_configuration();

        // We switch any systematic operations off so we code
        // symbols from the beginning
        if (kodo_core::has_set_systematic_off<Encoder>::value)
            kodo_core::set_systematic_off(*m_encoder);

	
        RUN
        {


            while (!m_decoder->is_complete())
            {
            	symbols_used++;
                // Encode a packet into the payload buffer
                m_encoder->write_payload(m_payload_buffer.data());

                m_recoder->read_payload(m_payload_buffer.data());
                m_recoder->nested()->write_payload(m_payload_buffer.data());
                // Pass that packet to the decoder
                m_decoder->read_payload(m_payload_buffer.data());

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

		assert(symbols.size() > 0);
		assert(symbol_size.size() > 0);
		assert(types.size() > 0);

		assert(expansion.size() > 0);

		assert(erasure.size() > 0);

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
								gauge::config_set cs;
								cs.set_value<uint32_t>("symbols", symbols[i]);
								cs.set_value<uint32_t>("symbol_size", symbol_size[j]);
								cs.set_value<std::string>("type", types[u]);
								cs.set_value<uint32_t>("expansion", e);

								cs.set_value<double>("erasure", er);

								add_configuration(cs);

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
};
