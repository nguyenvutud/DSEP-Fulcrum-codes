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
#include <cmath>

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
struct dynamic_both_throughput_benchmark: public gauge::time_benchmark {
	using encoder_factory = typename Encoder::factory;
	using encoder_ptr = typename Encoder::factory::pointer;

	using decoder_factory = typename Decoder::factory;
	using decoder_ptr = typename Decoder::factory::pointer;
	dynamic_both_throughput_benchmark() {

		// Seed the random generator controlling the erasures
		m_random_generator.seed((uint32_t) time(0));
	}
	void init() {
		m_factor = 1;
		gauge::time_benchmark::init();
	}

	void start() {
		m_encoded_symbols = 0;

		m_recovered_symbols = 0;
		m_processed_symbols = 0;
		gauge::time_benchmark::start();
	}

	void stop() {
		gauge::time_benchmark::stop();
	}

	double calculate_throughput(bool goodput = false) {
		// Get the time spent per iteration
		m_time = gauge::time_benchmark::measurement();

		gauge::config_set cs = get_current_configuration();
		std::string type = cs.get_value < std::string > ("type");
		//uint32_t symbols = cs.get_value<uint32_t>("symbols");
		uint32_t symbol_size = cs.get_value < uint32_t > ("symbol_size");

		// The number of bytes {en|de}coded
		uint64_t total_bytes = 0;

		if (type == "decoder") {
			if (goodput)
				total_bytes = m_recovered_symbols * symbol_size;
			else
				total_bytes = m_processed_symbols * symbol_size;
		} else if (type == "encoder") {
			total_bytes = m_encoded_symbols * symbol_size;
		} else {
			assert(0 && "Unknown benchmark type");
		}

		// The bytes per iteration
		uint64_t bytes = total_bytes / gauge::time_benchmark::iteration_count();

		return bytes / m_time; // MB/s for each iteration
	}

	void store_run(tables::table& results) {
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
//        if (!results.has_column("time"))
//                   results.add_column("time");
//               results.set_value("time", m_time);
		if (!results.has_column("dependency"))
			results.add_column("dependency");
		results.set_value("dependency", m_dependency);

		if (!results.has_column("nonzero_coefficient"))
			results.add_column("nonzero_coefficient");
		results.set_value("nonzero_coefficient", m_nonzero_coefficient);

	}

	bool accept_measurement() {
		gauge::config_set cs = get_current_configuration();

		std::string type = cs.get_value < std::string > ("type");

		if (type == "decoder") {
			// If we are benchmarking a decoder we only accept
			// the measurement if the decoding was successful
			if (!m_decoder->is_complete()) {
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

	std::string unit_text() const {
		return "MB/s";
	}

	void get_options(gauge::po::variables_map& options) {
		auto symbols = options["symbols"].as<std::vector<uint32_t>>();
		auto symbol_size = options["symbol_size"].as<std::vector<uint32_t>>();
		auto types = options["type"].as<std::vector<std::string>>();
		auto expansion = options["expansion"].as<std::vector<uint32_t> >();
		auto erasure = options["erasure"].as<std::vector<double> >();
		auto threshold = options["threshold"].as<std::vector<uint32_t> >();
		auto redundant = options["redundant"].as<std::vector<uint32_t> >();

		assert(symbols.size() > 0);
		assert(symbol_size.size() > 0);
		assert(types.size() > 0);
		assert(expansion.size() > 0);

		assert(erasure.size() > 0);

		for (const auto& s : symbols) {
			for (const auto& p : symbol_size) {
				for (const auto& t : types) {
					for (const auto& e : expansion) {
						for (const auto& er : erasure) {
							for (const auto& red : redundant) {
								for (const auto& thr : threshold) {
									gauge::config_set cs;
									cs.set_value < uint32_t > ("symbols", s);
									cs.set_value < uint32_t
											> ("symbol_size", p);
									cs.set_value < std::string > ("type", t);
									cs.set_value < uint32_t > ("expansion", e);
									cs.set_value < uint32_t
											> ("redundant", red);
									cs.set_value < uint32_t
											> ("threshold", thr);
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

	uint32_t calculate_region(uint32_t symbols, uint32_t num_regions,
			uint32_t red) {
		return (uint32_t)(
				(symbols) * (pow(2, num_regions) - 1) / pow(2, num_regions));
	}

	double count_sparsity(uint32_t mu, uint32_t r, uint32_t n, uint32_t delta,
			uint32_t rx) {
//    	double root = pow(1 - n/ (1.0*(n + delta)), 1/(1.0*(n - rx))); //Eq. (13) main 11. tex
		double root = pow(1 - (n + r) / (1.0 * (n + r + delta)),
				1 / (1.0 * (n + mu - rx))); //Eq. (12) main 11. tex
		return (1 - root) * (n + mu);
	}

	void setup() {
		gauge::config_set cs = get_current_configuration();

		uint32_t symbols = cs.get_value < uint32_t > ("symbols");
		uint32_t symbol_size = cs.get_value < uint32_t > ("symbol_size");
		uint32_t expansion = cs.get_value < uint32_t > ("expansion");
		uint32_t redundant = cs.get_value < uint32_t > ("redundant");

		double erasure = cs.get_value<double>("erasure");

		m_distribution = boost::random::bernoulli_distribution<>(erasure);

		m_decoder_factory = std::make_shared < decoder_factory
				> (symbols, symbol_size);

		m_encoder_factory = std::make_shared < encoder_factory
				> (symbols, symbol_size);

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

		m_temp_payload1.resize(payload_size);

		m_dependency.resize(symbols + expansion + 1);
		std::fill(m_dependency.begin(), m_dependency.end(), 0);

		// Prepare storage to the encoded payloads
		uint32_t payload_count = symbols * m_factor;

		assert(payload_count > 0);

		// Allocate contiguous payload buffer and store payload pointers
		m_payload_buffer.resize(payload_count * payload_size);
		m_payloads.resize(payload_count);

		m_nonzero_coefficient.resize(symbols + expansion);
		std::fill(m_nonzero_coefficient.begin(), m_nonzero_coefficient.end(),
				0);

		for (uint32_t i = 0; i < payload_count; ++i) {
			m_payloads[i] = &m_payload_buffer[i * payload_size];
		}

		m_used = 0;

		m_transmissions = 0;

		m_sparsity.resize(symbols);
		uint32_t cur_r = 0;
		num_regions = 1;
		uint32_t size_cur_region = std::min(symbols - expansion + cur_r,
				calculate_region(symbols, num_regions, redundant));

		for (uint32_t i = 0; i < symbols; i++) {
			if ((i >= size_cur_region)) {

				num_regions += 1;
				size_cur_region = std::min(symbols - expansion + cur_r,
						calculate_region(symbols, num_regions, redundant));
				if (cur_r < expansion)
					cur_r += 1;
			}
			uint32_t s = (uint32_t) round(
					count_sparsity(cur_r, expansion, symbols, redundant, i));

			if ((s == 0) & (i <= symbols / 2))
				s = 1;

			m_sparsity[i] = std::min(symbols / 2, s);
//			std::cout<<m_sparsity[i]<<" ";
		}
	}

	/// Setup the factories - we re-factored this into its own function
	/// since some of the codecs like sparse, perpetual and fulcrum
	/// customize the setup of the factory.
	void setup_factories() {
		//Super::setup_factories();
		gauge::config_set cs = get_current_configuration();

		uint32_t expansion = cs.get_value < uint32_t > ("expansion");

		m_decoder_factory->set_expansion(expansion);
		m_encoder_factory->set_expansion(expansion);

	}

	void encode_payloads_DTEP(uint32_t symbols, uint32_t expansion,
			uint32_t redundant) {
		// If we are working in the prime2325 finite field we have to
		// "map" the data first. This cost is included in the encoder
		// throughput (we do it with the clock running), just as it
		// would be in a real application.
		if (fifi::is_prime2325<typename Encoder::field_type>::value) {
			uint32_t block_length = fifi::size_to_length < fifi::prime2325
					> (m_encoder->block_size());

			fifi::prime2325_binary_search search(block_length);
			m_prefix = search.find_prefix(storage::storage(m_data_in));

			// Apply the negated prefix
			fifi::apply_prefix(storage::storage(m_data_in), ~m_prefix);
		}

		m_encoder->set_const_symbols(storage::storage(m_data_in));

		// We switch any systematic operations off so we code
		// symbols from the beginning
		if (kodo_core::has_set_systematic_off < Encoder > ::value)
			kodo_core::set_systematic_off(*m_encoder);

		uint32_t cur_sparse = 5, cur_r = 0, received_symbols = 0;
		num_regions = 1;
		int size_cur_region = 0;

		//size_cur_region =  calculate_region(symbols + expansion, num_regions, redundant);
		size_cur_region = std::min(symbols - expansion + cur_r,
				calculate_region(symbols, num_regions, redundant));
		//std::cout<<"size of region:"<<size_cur_region<<std::endl;

		//initialize the expansion
		m_encoder->nested()->set_num_nonzero_expansion(cur_r, expansion);

		//initialize the first sparsity
//				m_encoder->nested()-> set_average_nonzero_symbols(K[redundant / 5 - 1][0]);

//				std::cout<<"size of region:"<<size_cur_region<<std::endl;
//				std::cout<<"K = "<<K[redundant / 5 - 1][0]<<std::endl;

		uint32_t sum = 0, pre_sparsity = 0;

		bool normal_mode = false;
		uint32_t g = 0, m_sent = 0;

		for (uint8_t* payload : m_payloads) {

			if ((m_sent >= size_cur_region)) {

				num_regions += 1;
				size_cur_region = std::min(symbols - expansion + cur_r,
						calculate_region(symbols, num_regions, redundant));

				//check to change the expansion packets
				if (cur_r < expansion) {
					cur_r += 1;
					m_encoder->nested()->set_num_nonzero_expansion(cur_r,
							expansion, false, false, true);

//							std::cout<<"Changing Current expansion to:"<<cur_r<<" rx: "<<m_encoded_symbols<<std::endl;
				} else if (!normal_mode) {
					//set normal mode: randomly generate nonzero coeff from 1->n+r
					m_encoder->nested()->set_num_nonzero_expansion(expansion,
							expansion, false, false, true);
					normal_mode = true;
				}

//						std::cout<<"Number encoded:"<<m_encoded_symbols<<" Size of region: "<<size_cur_region<<std::endl;

			}
			//count the current sparsity
//					std::cout<<"n:"<<symbols<<" r:"<<expansion<<" mu:"<<cur_r<<" rx:"<<m_encoded_symbols<<std::endl;

			uint32_t w = 1;
			if (m_sent > symbols - 1)
				w = symbols / 2;
			else
				w = m_sparsity[m_sent];
//					std::cout<<"Sparsity:"<< w <<" rx: "<<m_encoded_symbols<<std::endl;
//					std::cout<< w <<std::endl;

			m_encoder->nested()->set_number_nonzero_coefficient(w);
			++m_encoded_symbols;
			++m_sent;

			m_encoder->write_payload(payload);

			if (g < symbols + expansion) {
				m_nonzero_coefficient[g] =
						m_encoder->nested()->nonzeros_generated();
				//				std::cout<<"nonzero coefficient at "<<g<<" is: "<<m_nonzero_coefficient[g]<<std::endl;
				g++;
			}

//					sum = sum + m_encoder->nested()->nonzeros_generated();

//					std::cout<<"Number of nonzero generated:"<<m_encoder->nested()->nonzeros_generated()<<std::endl;
		}
//				std::cout<<"Avg of nonzeros:"<<sum/(m_encoded_symbols*1.0)<<std::endl;

		if (fifi::is_prime2325<typename Encoder::field_type>::value) {
			// Apply the negated prefix
			fifi::apply_prefix(storage::storage(m_data_in), ~m_prefix);
		}
	}
	void encode_payloads_DTEP_max_sparsity(uint32_t symbols, uint32_t expansion,
			uint32_t redundant) {
		// If we are working in the prime2325 finite field we have to
		// "map" the data first. This cost is included in the encoder
		// throughput (we do it with the clock running), just as it
		// would be in a real application.
		if (fifi::is_prime2325<typename Encoder::field_type>::value) {
			uint32_t block_length = fifi::size_to_length < fifi::prime2325
					> (m_encoder->block_size());

			fifi::prime2325_binary_search search(block_length);
			m_prefix = search.find_prefix(storage::storage(m_data_in));

			// Apply the negated prefix
			fifi::apply_prefix(storage::storage(m_data_in), ~m_prefix);
		}

		m_encoder->set_const_symbols(storage::storage(m_data_in));

		// We switch any systematic operations off so we code
		// symbols from the beginning
		if (kodo_core::has_set_systematic_off < Encoder > ::value)
			kodo_core::set_systematic_off(*m_encoder);

		uint32_t cur_sparse = 1, cur_r = 0, received_symbols = 0;
		num_regions = 1;
		int size_cur_region = 0;

		//size_cur_region =  calculate_region(symbols + expansion, num_regions, redundant);
		size_cur_region = std::min(symbols - expansion + cur_r,
				calculate_region(symbols, num_regions, redundant));
		//std::cout<<"size of region:"<<size_cur_region<<std::endl;

		cur_sparse = std::min(symbols / 2,
				(uint32_t) count_sparsity(cur_r, expansion, symbols, redundant,
						size_cur_region));

//				std::cout<<"current sparse:"<<cur_sparse<<" in region: "<<num_regions<< " with size:"<<size_cur_region<<std::endl;

		m_encoder->nested()->set_number_nonzero_coefficient(cur_sparse);

		//initialize the expansion
		m_encoder->nested()->set_num_nonzero_expansion(cur_r, expansion);

		//initialize the first sparsity
//				m_encoder->nested()-> set_average_nonzero_symbols(K[redundant / 5 - 1][0]);

//				std::cout<<"size of region:"<<size_cur_region<<std::endl;
//				std::cout<<"K = "<<K[redundant / 5 - 1][0]<<std::endl;

		uint32_t sum = 0, pre_sparsity = 0;

		bool normal_mode = false;
		uint32_t g = 0;

		for (uint8_t* payload : m_payloads) {

			if ((m_encoded_symbols >= size_cur_region)) {

				num_regions += 1;
				size_cur_region = std::min(symbols - expansion + cur_r,
						calculate_region(symbols, num_regions, redundant));

				cur_sparse = std::min(symbols / 2,
						(uint32_t) count_sparsity(cur_r, expansion, symbols,
								redundant, size_cur_region));

//						std::cout<<"current sparse:"<<cur_sparse<<" in region: "<<num_regions<< " with size:"<<size_cur_region<<std::endl;

				m_encoder->nested()->set_number_nonzero_coefficient(cur_sparse);

				//check to change the expansion packets
				if (cur_r < expansion) {
					cur_r += 1;
					m_encoder->nested()->set_num_nonzero_expansion(cur_r,
							expansion, true);

//							std::cout<<"Current expansion:"<<cur_r<<std::endl;
				} else if (!normal_mode) {

					//set normal mode: randomly generate nonzero coeff from 1->n+r
					m_encoder->nested()->set_num_nonzero_expansion(expansion,
							expansion, false, false, true);
					normal_mode = true;
				}

//						std::cout<<"Number encoded:"<<m_encoded_symbols<<" Size of region: "<<size_cur_region<<std::endl;

			}

			++m_encoded_symbols;

			m_encoder->write_payload(payload);

			if (g < symbols + expansion) {
				m_nonzero_coefficient[g] =
						m_encoder->nested()->nonzeros_generated();
				//				std::cout<<"nonzero coefficient at "<<g<<" is: "<<m_nonzero_coefficient[g]<<std::endl;
				g++;
			}

//					sum = sum + m_encoder->nested()->nonzeros_generated();

//					std::cout<<"Number of nonzero generated:"<<m_encoder->nested()->nonzeros_generated()<<std::endl;
		}
//				std::cout<<"Avg of nonzeros:"<<sum/(m_encoded_symbols*1.0)<<std::endl;

		if (fifi::is_prime2325<typename Encoder::field_type>::value) {
			// Apply the negated prefix
			fifi::apply_prefix(storage::storage(m_data_in), ~m_prefix);
		}
	}

	void encode_payloads_STEP(uint32_t symbols, uint32_t expansion,
			uint32_t redundant, uint32_t threshold) {
		// If we are working in the prime2325 finite field we have to
		// "map" the data first. This cost is included in the encoder
		// throughput (we do it with the clock running), just as it
		// would be in a real application.
		if (fifi::is_prime2325<typename Encoder::field_type>::value) {
			uint32_t block_length = fifi::size_to_length < fifi::prime2325
					> (m_encoder->block_size());

			fifi::prime2325_binary_search search(block_length);
			m_prefix = search.find_prefix(storage::storage(m_data_in));

			// Apply the negated prefix
			fifi::apply_prefix(storage::storage(m_data_in), ~m_prefix);
		}

		m_encoder->set_const_symbols(storage::storage(m_data_in));

		// We switch any systematic operations off so we code
		// symbols from the beginning
		if (kodo_core::has_set_systematic_off < Encoder > ::value)
			kodo_core::set_systematic_off(*m_encoder);

		uint32_t cur_sparse = 0, cur_r = 0, received_symbols = 0;

		num_regions = 1;
		int size_cur_region = 0;

		//size_cur_region =  calculate_region(symbols + expansion, num_regions, redundant);
		size_cur_region = std::min(symbols - expansion + cur_r,
				calculate_region(symbols, num_regions, redundant));
		//std::cout<<"size of region:"<<size_cur_region<<std::endl;

		//initialize the expansion
		m_encoder->nested()->set_num_nonzero_expansion(0, expansion);

		uint32_t sum = 0, pre_sparsity = 1;

		bool normal_mode = false;
		uint32_t g = 0, m_sent = 0; //use m_sent to avoid the number of encoded packets is over symbols

		for (uint8_t* payload : m_payloads) {

			//set the expansion packets
			//using like-systematic the extra coefficient, random the one when send extra packets
			if (cur_r < expansion) {
				if (symbols - threshold <= m_sent) {
					cur_r = expansion;
					m_encoder->nested()->set_num_nonzero_expansion(cur_r,
							expansion, true);

					//std::cout<<"set to r_max r = "<<cur_r<<"--> symbols sent:"<<m_encoded_symbols<<std::endl;

				} else {
					//to avoid cur_r being never increased when symbols - threshold - expansion < 0
					int sub = symbols - threshold - expansion;
					if ((sub < 0)
							|| (symbols - threshold - expansion <= m_sent)) {
						cur_r += 1;
						m_encoder->nested()->set_num_nonzero_expansion(cur_r,
								expansion, true);
//								std::cout<<"r = "<<cur_r<<"--> symbols sent:"<<m_sent<<std::endl;
					}
					//else{
					//	m_encoder->nested()->set_num_nonzero_expansion(expansion, expansion, false, false, true);
					//}
				}
			} else if (!normal_mode) {
				//set normal mode: randomly generate nonzero coeff from 1->n+r
				m_encoder->nested()->set_num_nonzero_expansion(expansion,
						expansion, false, false, true);
				normal_mode = true;
			}

			///change sparsity following region
			uint32_t w = 1;
			if (m_sent > symbols - 1)
				w = symbols / 2;
			else
				w = m_sparsity[m_sent];
//					std::cout<< w <<std::endl;
			//std::cout<<"Prob innovative pcks:"<<1-pow(1 - s/(1.0*(symbols+cur_r)), symbols + cur_r - m_encoded_symbols)<<std::endl;

			m_encoder->nested()->set_number_nonzero_coefficient(w);

			++m_encoded_symbols;
			++m_sent;

			m_encoder->write_payload(payload);

			if (g < symbols + expansion) {
				m_nonzero_coefficient[g] =
						m_encoder->nested()->nonzeros_generated();
				g++;
			}
//					std::cout<<"Encoded symbols:"<<m_encoded_symbols<<
//							" --Number of nonzero generated:"
//							<<m_encoder->nested()->nonzeros_generated()
//							<<"-- sent:"<<m_sent<<std::endl;

		}
//				std::cout<<"Payload size:"<<m_payloads.size()<<std::endl;

		if (fifi::is_prime2325<typename Encoder::field_type>::value) {
			// Apply the negated prefix
			fifi::apply_prefix(storage::storage(m_data_in), ~m_prefix);
		}
	}

	//This function will create the payloads without any linearly dependent packets

	void encode_payloads_STEP_with_feedback(uint32_t symbols,
			uint32_t expansion, uint32_t redundant, uint32_t threshold) {
		// If we are working in the prime2325 finite field we have to
		// "map" the data first. This cost is included in the encoder
		// throughput (we do it with the clock running), just as it
		// would be in a real application.
		if (fifi::is_prime2325<typename Encoder::field_type>::value) {
			uint32_t block_length = fifi::size_to_length < fifi::prime2325
					> (m_encoder->block_size());

			fifi::prime2325_binary_search search(block_length);
			m_prefix = search.find_prefix(storage::storage(m_data_in));

			// Apply the negated prefix
			fifi::apply_prefix(storage::storage(m_data_in), ~m_prefix);
		}

		m_decoder->initialize(*m_decoder_factory);

		m_encoder->set_const_symbols(storage::storage(m_data_in));

		// We switch any systematic operations off so we code
		// symbols from the beginning
		if (kodo_core::has_set_systematic_off < Encoder > ::value)
			kodo_core::set_systematic_off(*m_encoder);

		uint32_t cur_sparse = 0, cur_r = 0, received_symbols = 0;

		num_regions = 1;
		int size_cur_region = 0;

		//size_cur_region =  calculate_region(symbols + expansion, num_regions, redundant);
		size_cur_region = std::min(symbols - expansion + cur_r,
				calculate_region(symbols, num_regions, redundant));
		//std::cout<<"size of region:"<<size_cur_region<<std::endl;

		//initialize the expansion
		m_encoder->nested()->set_num_nonzero_expansion(0, expansion);

		uint32_t sum = 0, pre_sparsity = 1;

		bool normal_mode = false, ldp_check = false;
		uint32_t g = 0, m_sent = 0; //use m_sent to avoid the number of encoded packets is over symbols

		uint32_t ind_payload = -1, curr_rank = 0;

		while (true) {
			//do sth before this

			//set the expansion packets
			//using like-systematic the extra coefficient, random the one when send extra packets
			if (!ldp_check) {
				if (cur_r < expansion) {
					if (symbols - threshold <= m_sent) {
						cur_r = expansion;
						m_encoder->nested()->set_num_nonzero_expansion(cur_r,
								expansion, true);

//							std::cout<<"set to r_max r = "<<cur_r<<"--> symbols sent:"<<m_encoded_symbols<<std::endl;

					} else {
						//to avoid cur_r being never increased when symbols - threshold - expansion < 0
						int sub = symbols - threshold - expansion;
						if ((sub < 0)
								|| (symbols - threshold - expansion <= m_sent)) {
							cur_r += 1;
							m_encoder->nested()->set_num_nonzero_expansion(
									cur_r, expansion, true);
//								std::cout<<"r = "<<cur_r<<"--> symbols sent:"<<m_sent<<std::endl;
						}
						//else{
						//	m_encoder->nested()->set_num_nonzero_expansion(expansion, expansion, false, false, true);
						//}
					}
				} else if (!normal_mode) {
					//set normal mode: randomly generate nonzero coeff from 1->n+r
					m_encoder->nested()->set_num_nonzero_expansion(expansion,
							expansion, false, false, true);
					normal_mode = true;
				}
				///change sparsity following region
				uint32_t w = 1;
				if (m_sent > symbols - 1)
					w = symbols / 2;
				else
					w = m_sparsity[m_sent];

				m_encoder->nested()->set_number_nonzero_coefficient(w);
			}

//					std::cout<< w <<std::endl;
			//std::cout<<"Prob innovative pcks:"<<1-pow(1 - s/(1.0*(symbols+cur_r)), symbols + cur_r - m_encoded_symbols)<<std::endl;

			++m_encoded_symbols;

			m_encoder->write_payload(m_temp_payload.data());

			//store the current coded packet

			std::copy_n(m_temp_payload.data(), m_encoder->payload_size(),
					m_temp_payload1.data());

			m_decoder->read_payload(m_temp_payload.data());

			ldp_check = false;

			//this is an independent packet
			if (m_decoder->rank() > curr_rank) {
				++m_sent;

				++ind_payload;
				uint8_t* payload;

				if (ind_payload < m_payloads.size())
					payload = m_payloads[ind_payload];

				//save an independent packet
				std::copy_n(m_temp_payload1.data(), m_decoder->payload_size(),
						payload);

				curr_rank = m_decoder->rank();

			} else {
				ldp_check = true;

			}

			if (g < symbols + expansion) {
				m_nonzero_coefficient[g] =
						m_encoder->nested()->nonzeros_generated();
				g++;
			}

			if (m_decoder->is_complete())
				break;

		}
//		std::cout<<"Payload size:"<<m_payloads.size()<<std::endl;
//		std::cout<<"Real payload size:"<<ind_payload<<std::endl;
//		std::cout<<"Total independent packets:"<<m_sent<<std::endl;
//		std::cout<<"Total encoded packets:"<<m_encoded_symbols<<std::endl;

		if (fifi::is_prime2325<typename Encoder::field_type>::value) {
			// Apply the negated prefix
			fifi::apply_prefix(storage::storage(m_data_in), ~m_prefix);
		}
	}
	void decode_payloads() {
		gauge::config_set cs = get_current_configuration();
		//  uint32_t symbols = cs.get_value<uint32_t>("symbols");

		uint32_t expansion = cs.get_value < uint32_t > ("expansion");

		// If the decoder uses shallow storage, we have to initialize
		// its decoding buffer
		if (kodo_core::has_set_mutable_symbols < Decoder > ::value) {
			kodo_core::set_mutable_symbols(*m_decoder,
					storage::storage(m_data_out));
		}
		uint32_t dependent_counter = 0, pre_rank = 0;
		for (const uint8_t* payload : m_payloads) {
			std::copy_n(payload, m_decoder->payload_size(),
					m_temp_payload.data());

			++m_transmissions;

			if (m_distribution(m_random_generator)) {
//				std::cout<<"LOST HAPPPED\n";
				continue;
			}

			m_decoder->read_payload(m_temp_payload.data());

			++m_processed_symbols;
			++m_used;
			uint32_t r = m_decoder->rank();
			uint32_t symbols = m_decoder->symbols();

//			if(r > symbols + expansion) r = symbols + expansion;

			if (r == pre_rank)
				++m_dependency[r];
			else
				pre_rank = m_decoder->rank();

			//std::cout<<"Symbols:"<<symbols<<" Decodeing rank:" << r <<" dependency:"<<m_dependency[r-1]<<std::endl;

			if (m_decoder->is_complete()) {
				break;
			}
		}

		if (m_decoder->is_complete()) {
			// If the decoder uses deep storage, we have to copy
			// its decoding buffers for verification.
			// This would also be necessary in real applications, so it
			// should be part of the timed benchmark.
			if (kodo_core::has_deep_symbol_storage < Decoder > ::value) {
				m_decoder->copy_from_symbols(storage::storage(m_data_out));
			}

			if (fifi::is_prime2325<typename Decoder::field_type>::value) {
				// Now we have to apply the negated prefix to the decoded data
				fifi::apply_prefix(storage::storage(m_data_out), ~m_prefix);
			}

			m_recovered_symbols += m_decoder->symbols();

			//std::cout<<"Decode completely:" <<m_recovered_symbols <<std::endl;
		}
	}

	void decode_payloads_with_feedback() {
		gauge::config_set cs = get_current_configuration();
		//  uint32_t symbols = cs.get_value<uint32_t>("symbols");

		uint32_t expansion = cs.get_value < uint32_t > ("expansion");

		// If the decoder uses shallow storage, we have to initialize
		// its decoding buffer
		if (kodo_core::has_set_mutable_symbols < Decoder > ::value) {
			kodo_core::set_mutable_symbols(*m_decoder,
					storage::storage(m_data_out));
		}
		uint32_t dependent_counter = 0, pre_rank = 0;
		uint32_t ind_pck = 0;

		while (ind_pck < m_payloads.size()) {

			//get coded data from payloads
			const uint8_t* payload = m_payloads[ind_pck];

			std::copy_n(payload, m_decoder->payload_size(),
					m_temp_payload.data());

			++m_transmissions;

			if (m_distribution(m_random_generator)) {
				//				std::cout<<"LOST HAPPPED\n";
				continue;
			}

			m_decoder->read_payload(m_temp_payload.data());

			++m_processed_symbols;
			++m_used;
			++ind_pck;

			uint32_t r = m_decoder->rank();
			uint32_t symbols = m_decoder->symbols();

			//			if(r > symbols + expansion) r = symbols + expansion;

			if (r == pre_rank)
				++m_dependency[r];
			else
				pre_rank = m_decoder->rank();

			//std::cout<<"Symbols:"<<symbols<<" Decodeing rank:" << r <<" dependency:"<<m_dependency[r-1]<<std::endl;

			if (m_decoder->is_complete()) {
				break;
			}
		}

		if (m_decoder->is_complete()) {
			// If the decoder uses deep storage, we have to copy
			// its decoding buffers for verification.
			// This would also be necessary in real applications, so it
			// should be part of the timed benchmark.
			if (kodo_core::has_deep_symbol_storage < Decoder > ::value) {
				m_decoder->copy_from_symbols(storage::storage(m_data_out));
			}

			if (fifi::is_prime2325<typename Decoder::field_type>::value) {
				// Now we have to apply the negated prefix to the decoded data
				fifi::apply_prefix(storage::storage(m_data_out), ~m_prefix);
			}

			m_recovered_symbols += m_decoder->symbols();

			//std::cout<<"Decode completely:" <<m_recovered_symbols <<std::endl;
		}
	}

	///FULCRUM DYNAMIC SPARSE
	/// Run the encoder
	void run_fulcrum_encode_DTEP() {
		gauge::config_set cs = get_current_configuration();
		uint32_t symbols = cs.get_value < uint32_t > ("symbols");

		uint32_t expansion = cs.get_value < uint32_t > ("expansion");

		uint32_t redundant = cs.get_value < uint32_t > ("redundant");

		if (kodo_core::has_set_systematic_off < Encoder > ::value)
			kodo_core::set_systematic_off(*m_encoder);

		// The clock is running
	RUN
	{
		// We have to make sure the encoder is in a "clean" state
		m_encoder->initialize(*m_encoder_factory);

//		    m_encoder->nested()->set_number_nonzero_coefficient(sparse);
		//m_encoder->nested()->set_average_nonzero_symbols(sparse);
		//m_encoder->nested()->set_num_nonzero_expansion(expansion, expansion, false, false, true);

		encode_payloads_DTEP(symbols, expansion, redundant);
	}
}
void run_fulcrum_encode_DTEP_max_sparsity() {
	gauge::config_set cs = get_current_configuration();
	uint32_t symbols = cs.get_value < uint32_t > ("symbols");

	uint32_t expansion = cs.get_value < uint32_t > ("expansion");

	uint32_t redundant = cs.get_value < uint32_t > ("redundant");

	if (kodo_core::has_set_systematic_off < Encoder > ::value)
		kodo_core::set_systematic_off(*m_encoder);

	// The clock is running
RUN
{
	// We have to make sure the encoder is in a "clean" state
	m_encoder->initialize(*m_encoder_factory);

	//		    m_encoder->nested()->set_number_nonzero_coefficient(sparse);
	//m_encoder->nested()->set_average_nonzero_symbols(sparse);
	//m_encoder->nested()->set_num_nonzero_expansion(expansion, expansion, false, false, true);

	encode_payloads_DTEP_max_sparsity(symbols, expansion, redundant);
}
}

	/// Run the decoder
void run_fulcrum_decode_DTEP_max_sparsity() {
gauge::config_set cs = get_current_configuration();
uint32_t symbols = cs.get_value < uint32_t > ("symbols");
uint32_t expansion = cs.get_value < uint32_t > ("expansion");
uint32_t redundant = cs.get_value < uint32_t > ("redundant");

m_encoder->initialize(*m_encoder_factory);

if (kodo_core::has_set_systematic_off < Encoder > ::value)
	kodo_core::set_systematic_off(*m_encoder);
// Encode some data
encode_payloads_DTEP_max_sparsity(symbols, expansion, redundant);

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

	/// Run the decoder
void run_fulcrum_decode_DTEP() {
gauge::config_set cs = get_current_configuration();
uint32_t symbols = cs.get_value < uint32_t > ("symbols");
uint32_t expansion = cs.get_value < uint32_t > ("expansion");
uint32_t redundant = cs.get_value < uint32_t > ("redundant");

m_encoder->initialize(*m_encoder_factory);

if (kodo_core::has_set_systematic_off < Encoder > ::value)
kodo_core::set_systematic_off(*m_encoder);
	// Encode some data
encode_payloads_DTEP(symbols, expansion, redundant);

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

void run_fulcrum_decode_DTEP_with_feedback() {
gauge::config_set cs = get_current_configuration();
uint32_t symbols = cs.get_value < uint32_t > ("symbols");
uint32_t expansion = cs.get_value < uint32_t > ("expansion");
uint32_t redundant = cs.get_value < uint32_t > ("redundant");

m_encoder->initialize(*m_encoder_factory);

if (kodo_core::has_set_systematic_off < Encoder > ::value)
kodo_core::set_systematic_off(*m_encoder);
	// Encode some data
encode_payloads_DTEP(symbols, expansion, redundant);

	// Zero the data buffer for the decoder
std::fill_n(m_data_out.begin(), m_data_out.size(), 0);

	// The clock is running
RUN
{
	// We have to make sure the decoder is in a "clean" state
	// i.e. no symbols already decoded.
m_decoder->initialize(*m_decoder_factory);

	// Decode the payloads
decode_payloads_with_feedback();
}
}


void run_body_fulcrum_DTEP() {
gauge::config_set cs = get_current_configuration();

std::string type = cs.get_value < std::string > ("type");

if (type == "encoder") {
run_fulcrum_encode_DTEP();
} else if (type == "decoder") {
run_fulcrum_decode_DTEP();
} else {
assert(0 && "Unknown benchmark type");
}
}

void run_body_fulcrum_DTEP_with_feedback() {
gauge::config_set cs = get_current_configuration();

std::string type = cs.get_value < std::string > ("type");

if (type == "encoder") {
run_fulcrum_encode_DTEP();
} else if (type == "decoder") {
run_fulcrum_decode_DTEP_with_feedback();
} else {
assert(0 && "Unknown benchmark type");
}
}

void run_body_fulcrum_DTEP_max_sparsity() {
gauge::config_set cs = get_current_configuration();

std::string type = cs.get_value < std::string > ("type");

if (type == "encoder") {
run_fulcrum_encode_DTEP_max_sparsity();
} else if (type == "decoder") {
run_fulcrum_decode_DTEP_max_sparsity();
} else {
assert(0 && "Unknown benchmark type");
}
}

	///FULCRUM STEP: DYNAMIC SPARSITY and EXPANSION
	/// Run the encoder
void run_fulcrum_encode_STEP() {
gauge::config_set cs = get_current_configuration();
uint32_t symbols = cs.get_value < uint32_t > ("symbols");

uint32_t expansion = cs.get_value < uint32_t > ("expansion");

uint32_t redundant = cs.get_value < uint32_t > ("redundant");
uint32_t threshold = cs.get_value < uint32_t > ("threshold");

if (kodo_core::has_set_systematic_off < Encoder > ::value)
kodo_core::set_systematic_off(*m_encoder);

			 // The clock is running
RUN
{
			 // We have to make sure the encoder is in a "clean" state
m_encoder->initialize(*m_encoder_factory);

//		    m_encoder->nested()->set_number_nonzero_coefficient(sparse);
//m_encoder->nested()->set_average_nonzero_symbols(sparse);
//m_encoder->nested()->set_num_nonzero_expansion(expansion, expansion, false, false, true);

encode_payloads_STEP(symbols, expansion, redundant, threshold);
}
}

			 /// Run the decoder
void run_fulcrum_decode_STEP() {
gauge::config_set cs = get_current_configuration();
uint32_t symbols = cs.get_value < uint32_t > ("symbols");
uint32_t expansion = cs.get_value < uint32_t > ("expansion");
uint32_t redundant = cs.get_value < uint32_t > ("redundant");
uint32_t threshold = cs.get_value < uint32_t > ("threshold");

m_encoder->initialize(*m_encoder_factory);

if (kodo_core::has_set_systematic_off < Encoder > ::value)
kodo_core::set_systematic_off(*m_encoder);
			 // Encode some data
encode_payloads_STEP(symbols, expansion, redundant, threshold);

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

			/// Run the decoder
void run_fulcrum_decode_STEP_feedback() {
gauge::config_set cs = get_current_configuration();
uint32_t symbols = cs.get_value < uint32_t > ("symbols");
uint32_t expansion = cs.get_value < uint32_t > ("expansion");
uint32_t redundant = cs.get_value < uint32_t > ("redundant");
uint32_t threshold = cs.get_value < uint32_t > ("threshold");

m_encoder->initialize(*m_encoder_factory);

if (kodo_core::has_set_systematic_off < Encoder > ::value)
kodo_core::set_systematic_off(*m_encoder);
			// Encode some data

//  		 encode_payloads_STEP_with_feed_back(symbols, expansion, redundant, threshold);

encode_payloads_STEP(symbols, expansion, redundant, threshold);

			// Zero the data buffer for the decoder
std::fill_n(m_data_out.begin(), m_data_out.size(), 0);

			// The clock is running
RUN
{
			// We have to make sure the decoder is in a "clean" state
			// i.e. no symbols already decoded.
m_decoder->initialize(*m_decoder_factory);

			// Decode the payloads
decode_payloads_with_feedback();
}
}

void run_body_fulcrum_STEP() {
gauge::config_set cs = get_current_configuration();

std::string type = cs.get_value < std::string > ("type");

if (type == "encoder") {
run_fulcrum_encode_STEP();
} else if (type == "decoder") {
run_fulcrum_decode_STEP();
} else {
assert(0 && "Unknown benchmark type");
}
}

void run_body_fulcrum_STEP_with_feedback() {
gauge::config_set cs = get_current_configuration();

std::string type = cs.get_value < std::string > ("type");

if (type == "encoder") {
//          	run_fulcrum_encode_STEP_with_feedback();
} else if (type == "decoder") {
run_fulcrum_decode_STEP_feedback();
} else {
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
std::vector<uint8_t> m_temp_payload, m_temp_payload1;

			 /// Contiguous buffer for coded payloads
std::vector<uint8_t> m_payload_buffer;

			 /// Pointers to each payload in the payload buffer
std::vector<uint8_t*> m_payloads;

			 /// Multiplication factor for payload_count
uint32_t m_factor = 2;

			 /// Prefix used when testing with the prime2325 finite field
uint32_t m_prefix;

boost::random::mt19937 m_random_generator;

boost::random::bernoulli_distribution<> m_distribution;

uint32_t m_used = 0;

uint32_t m_transmissions = 0;

uint32_t num_regions = 0;

double m_time;
//    uint32_t K[4][8] ={{6, 10, 19, 33, 48, 67, 54, 54},
//    				   {4, 8, 15, 25, 40, 58, 57, 57},
//    				   {3, 7, 13, 21, 35, 54, 59, 59},
//    				   {3, 6, 11, 19, 32, 47, 62, 62}};
//uint32_t K[4][8] = { { 5, 7, 8, 9, 11, 13, 15, 17 }, { 5, 9, 17, 27, 40, 52, 52,
//52 }, { 3, 7, 13, 21, 35, 54, 59, 59 }, { 3, 6, 11, 19, 32, 47, 55, 59 } };
std::vector<uint32_t> m_dependency;
std::vector<uint32_t> m_nonzero_coefficient;
std::vector<double> m_sparsity;
};
