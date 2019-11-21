// Copyright Steinwurf ApS 2011.
// Distributed under the "STEINWURF RESEARCH LICENSE 1.0".
// See accompanying file LICENSE.rst or
// http://www.steinwurf.com/licensing

#include <gauge/gauge.hpp>

#include <kodo_fulcrum/fulcrum_codes.hpp>
#include <kodo_fulcrum_sparse/fulcrum_sparse_codes.hpp>
#include <tunable_fulcrum/tunable_fulcrum_encoder.hpp>

#include "static_sparse_tunable_throughput_benchmark.hpp"
#include "dynamic_sparse_tunable_throughput_benchmark.hpp"
#include "dynamic_both_throughput_benchmark.hpp"
#include "tunable_throughput_benchmark.hpp"

template<class SuperCoder>
class fulcrum_combined_decoder: public SuperCoder {
public:
	using stage_one_decoder_type = typename SuperCoder::stage_one_decoder_type;
	using stage_two_decoder_type = typename SuperCoder::stage_two_decoder_type;

	/// @return The stage one decoder
	uint32_t stage_one_decoder_rank() {
		return SuperCoder::stage_one_decoder().rank();
	}

	/// @return The stage two decoder
	uint32_t stage_two_decoder_rank() {
		return SuperCoder::stage_two_decoder().rank();
	}

	uint32_t rank()
	{
		return stage_one_decoder_rank() + stage_two_decoder_rank();
	}


	using factory = kodo_core::pool_factory<fulcrum_combined_decoder>;
};

template<class SuperCoder>
class fulcrum_inner_decoder: public SuperCoder {
public:

	uint32_t rank()
	{
		const auto& nested = SuperCoder::nested();
		assert(nested);
		return nested->rank();
	}

	using factory = kodo_core::pool_factory<fulcrum_inner_decoder>;
};

/// Using this macro we may specify options. For specifying options
/// we use the boost program options library. So you may additional
/// details on how to do it in the manual for that library.
BENCHMARK_OPTION(throughput_options)
{
    gauge::po::options_description options;

    std::vector<uint32_t> symbols;
//    symbols.push_back(16);
//    symbols.push_back(32);
    symbols.push_back(64);
//      symbols.push_back(100);
//    symbols.push_back(128);
//    symbols.push_back(256);
//    symbols.push_back(512);
//    symbols.push_back(1024);

    auto default_symbols =
        gauge::po::value<std::vector<uint32_t> >()->default_value(
            symbols, "")->multitoken();

    std::vector<uint32_t> symbol_size;
    symbol_size.push_back(1500);

    auto default_symbol_size =
        gauge::po::value<std::vector<uint32_t> >()->default_value(
            symbol_size, "")->multitoken();

    std::vector<uint32_t> expansion;
//    expansion.push_back(1);
//    expansion.push_back(2);
//    expansion.push_back(3);
    expansion.push_back(4);
//    expansion.push_back(5);
//    expansion.push_back(6);
//    expansion.push_back(7);
//    expansion.push_back(8);
//

    auto default_expansion =
        gauge::po::value<std::vector<uint32_t> >()->default_value(
            expansion, "")->multitoken();

    std::vector<uint32_t> sparse;
//    sparse.push_back(3);
    sparse.push_back(5);
//    sparse.push_back(9);
//    sparse.push_back(15);
//    sparse.push_back(17);

//    sparse.push_back(50);

    auto default_sparse =
            gauge::po::value<std::vector<uint32_t> >()->default_value(
                sparse, "")->multitoken();

    std::vector<double> density;
//    density.push_back(0.1);
    density.push_back(0.5);
    auto default_density =
            gauge::po::value<std::vector<double> >()->default_value(
                density, "")->multitoken();

    std::vector<double> erasures;
    erasures.push_back(0);
//    erasures.push_back(0.1);
//	erasures.push_back(0.2);
//    erasures.push_back(0.3);
//    erasures.push_back(0.4);
//    erasures.push_back(0.5);
//   erasures.push_back(0.6);
//   erasures.push_back(0.7);
//   erasures.push_back(0.8);
//   erasures.push_back(0.9);

    auto default_erasure =
        gauge::po::value<std::vector<double> >()->default_value(
            erasures, "")->multitoken();

    std::vector<std::string> types;
//    types.push_back("encoder");
    types.push_back("decoder");

    auto default_types =
        gauge::po::value<std::vector<std::string> >()->default_value(
            types, "")->multitoken();

    options.add_options()
    ("symbols", default_symbols, "Set the number of symbols");

    options.add_options()
    ("symbol_size", default_symbol_size, "Set the symbol size in bytes");

    options.add_options()
    ("expansion", default_expansion,
     "Set the expansion of the fulcrum codes");

    options.add_options()
      ("sparse", default_sparse,
       "Set the sparsity of the fulcrum codes");



 options.add_options()
    ("erasure", default_erasure,
     "Set the erasure of the fulcrum codes");
 options.add_options()
     ("density", default_density,
      "Set the density of the fulcrum codes");


    options.add_options()
    ("type", default_types, "Set type [encoder|decoder]");

    gauge::runner::instance().register_options(options);
}

BENCHMARK_OPTION(redundant_options)
{
    gauge::po::options_description options;

    std::vector<uint32_t> redundancies;
//    redundancies.push_back(5);
    redundancies.push_back(10);
//    redundancies.push_back(15);
//    redundancies.push_back(20);
//    redundancies.push_back(30);


    auto default_redundancies =
        gauge::po::value<std::vector<uint32_t> >()->default_value(
            redundancies, "")->multitoken();
    options.add_options()
     ("redundant", default_redundancies, "Set redundancies");

     gauge::runner::instance().register_options(options);
}
BENCHMARK_OPTION(threshod_options)
{
    gauge::po::options_description options;

    std::vector<uint32_t> threshold;
//    threshold.push_back(0);
//    threshold.push_back(2);
    threshold.push_back(4);
//    threshold.push_back(6);
//    threshold.push_back(8);
//    threshold.push_back(10);


    auto default_threshold =
        gauge::po::value<std::vector<uint32_t> >()->default_value(
            threshold, "")->multitoken();
    options.add_options()
     ("threshold", default_threshold, "Set threshold");

     gauge::runner::instance().register_options(options);
}

//------------------------------------------------------------------
// Original Fulcrum: this measurement can run with different density: 0.1 0.5
//------------------------------------------------------------------
//using setup_fulcrum_throughput_inner =
//	static_sparse_tunable_throughput_benchmark<
//    kodo_fulcrum::tunable_fulcrum_encoder<fifi::binary8>,
//	fulcrum_inner_decoder<kodo_fulcrum::fulcrum_inner_decoder<fifi::binary>>>;
//
//BENCHMARK_F_INLINE(setup_fulcrum_throughput_inner, FulcrumInner, Binary, 3)
//{
//	run_body_fulcrum();
//}
//
//using setup_fulcrum_throughput_outer =
//	static_sparse_tunable_throughput_benchmark<
//    kodo_fulcrum::tunable_fulcrum_encoder<fifi::binary8>,
//	kodo_fulcrum::fulcrum_outer_decoder<fifi::binary8>>;
//
//BENCHMARK_F_INLINE(setup_fulcrum_throughput_outer, FulcrumOuter, Binary, 3)
//{
//	run_body_fulcrum();
//}
//
//
//using setup_fulcrum_throughput_combined =
//	static_sparse_tunable_throughput_benchmark<
//    kodo_fulcrum::tunable_fulcrum_encoder<fifi::binary8>,
//	fulcrum_combined_decoder<kodo_fulcrum::fulcrum_combined_decoder<fifi::binary8>>>;
//
//BENCHMARK_F_INLINE(setup_fulcrum_throughput_combined, FulcrumCombined, Binary, 3)
//{
//	run_body_fulcrum();
//}

//////////////////////////////////////////////////

//using setup_fulcrum_throughput_inner =
//	static_sparse_tunable_throughput_benchmark<
//    kodo_fulcrum::tunable_fulcrum_encoder<fifi::binary8>,
//	fulcrum_inner_decoder<kodo_fulcrum::fulcrum_inner_decoder<fifi::binary>>>;
//
//BENCHMARK_F_INLINE(setup_fulcrum_throughput_inner, FulcrumInner, Binary, 3)
//{
//	run_body_fulcrum_with_step();
//}

//using setup_fulcrum_throughput_outer =
//	static_sparse_tunable_throughput_benchmark<
//    kodo_fulcrum::tunable_fulcrum_encoder<fifi::binary8>,
//	kodo_fulcrum::fulcrum_outer_decoder<fifi::binary8>>;
//
//BENCHMARK_F_INLINE(setup_fulcrum_throughput_outer, FulcrumOuter, Binary, 3)
//{
//	run_body_fulcrum_with_step();
//}


//using setup_fulcrum_throughput_combined =
//	static_sparse_tunable_throughput_benchmark<
//    kodo_fulcrum::tunable_fulcrum_encoder<fifi::binary8>,
//	fulcrum_combined_decoder<kodo_fulcrum::fulcrum_combined_decoder<fifi::binary8>>>;
//
//BENCHMARK_F_INLINE(setup_fulcrum_throughput_combined, FulcrumCombined, Binary, 3)
//{
//	run_body_fulcrum_with_step();
//}
//------------------------------------------------------------------
// Average sparse Fulcrum: generate the nonzero coefficient around w (=5)
//------------------------------------------------------------------
//using setup_fulcrum_throughput_inner =
//	static_sparse_tunable_throughput_benchmark<
//    kodo_fulcrum::tunable_fulcrum_encoder<fifi::binary8>,
//	fulcrum_inner_decoder<kodo_fulcrum::fulcrum_inner_decoder<fifi::binary>>>;
//
//BENCHMARK_F_INLINE(setup_fulcrum_throughput_inner, FulcrumInner, Binary, 1)
//{
//	run_body_fulcrum_average();
//}
//
//using setup_fulcrum_throughput_outer =
//	static_sparse_tunable_throughput_benchmark<
//    kodo_fulcrum::tunable_fulcrum_encoder<fifi::binary8>,
//	kodo_fulcrum::fulcrum_outer_decoder<fifi::binary8>>;
//
//BENCHMARK_F_INLINE(setup_fulcrum_throughput_outer, FulcrumOuter, Binary, 1)
//{
//	run_body_fulcrum_average();
//}
//
//using setup_fulcrum_throughput_combined =
//	static_sparse_tunable_throughput_benchmark<
//    kodo_fulcrum::tunable_fulcrum_encoder<fifi::binary8>,
//	fulcrum_combined_decoder<kodo_fulcrum::fulcrum_combined_decoder<fifi::binary8>>>;
//
//BENCHMARK_F_INLINE(setup_fulcrum_throughput_combined, FulcrumCombined, Binary, 3)
//{
//	run_body_fulcrum_average();
//}

//------------------------------------------------------------------
// Static sparse Fulcrum (fix sparsity each coding vector)
//------------------------------------------------------------------
//using setup_fulcrum_throughput_inner =
//	static_sparse_tunable_throughput_benchmark<
//    kodo_fulcrum::fulcrum_SD_encoder<fifi::binary8>,
//	kodo_fulcrum::fulcrum_sparse_inner_decoder<fifi::binary>>;
//
//BENCHMARK_F_INLINE(setup_fulcrum_throughput_inner, FulcrumInner, Binary, 3)
//{
//	run_body_fulcrum_static();
//}
//
//using setup_fulcrum_throughput_outer =
//	static_sparse_tunable_throughput_benchmark<
//    kodo_fulcrum::fulcrum_SD_encoder<fifi::binary8>,
//	kodo_fulcrum::fulcrum_outer_decoder<fifi::binary8>>;
//
//BENCHMARK_F_INLINE(setup_fulcrum_throughput_outer, FulcrumOuter, Binary, 3)
//{
//	run_body_fulcrum_static();
//}
//
//using setup_fulcrum_throughput_combined =
//	static_sparse_tunable_throughput_benchmark<
//    kodo_fulcrum::fulcrum_SD_encoder<fifi::binary8>,
//	kodo_fulcrum::fulcrum_combined_decoder<fifi::binary8>>;
//
//BENCHMARK_F_INLINE(setup_fulcrum_throughput_combined, FulcrumCombined, Binary, 3)
//{
//	run_body_fulcrum_static();
//}

//------------------------------------------------------------------
// Dynamic sparse Fulcrum
//------------------------------------------------------------------
//using setup_fulcrum_throughput_inner =
//	dynamic_sparse_tunable_throughput_benchmark<
//    kodo_fulcrum::tunable_fulcrum_encoder<fifi::binary8>,
//	kodo_fulcrum::fulcrum_inner_decoder<fifi::binary>>;
//
//BENCHMARK_F_INLINE(setup_fulcrum_throughput_inner, FulcrumInner, Binary, 1)
//{
//	run_body_fulcrum_dynamic();
//}
//
//using setup_fulcrum_throughput_outer =
//	dynamic_sparse_tunable_throughput_benchmark<
//    kodo_fulcrum::tunable_fulcrum_encoder<fifi::binary8>,
//	kodo_fulcrum::fulcrum_outer_decoder<fifi::binary8>>;
//
//BENCHMARK_F_INLINE(setup_fulcrum_throughput_outer, FulcrumOuter, Binary, 1)
//{
//	run_body_fulcrum_dynamic();
//}
//using setup_fulcrum_throughput_combined =
//	dynamic_sparse_tunable_throughput_benchmark<
//    kodo_fulcrum::tunable_fulcrum_encoder<fifi::binary8>,
//	kodo_fulcrum::fulcrum_combined_decoder<fifi::binary8>>;
//
//BENCHMARK_F_INLINE(setup_fulcrum_throughput_combined, FulcrumCombined, Binary, 1)
//{
//	run_body_fulcrum_dynamic();
//}
//------------------------------------------------------------------
// Dynamic both - DTEP Fulcrum
//------------------------------------------------------------------
using setup_fulcrum_throughput_inner =
	dynamic_both_throughput_benchmark<
    kodo_fulcrum::tunable_fulcrum_encoder<fifi::binary8>,
	fulcrum_inner_decoder<kodo_fulcrum::fulcrum_inner_decoder<fifi::binary>>>;

BENCHMARK_F_INLINE(setup_fulcrum_throughput_inner, FulcrumInner, Binary, 1)
{
//	run_body_fulcrum_DTEP();
	run_body_fulcrum_DTEP_with_feedback();
}

using setup_fulcrum_throughput_outer =
	dynamic_both_throughput_benchmark<
    kodo_fulcrum::tunable_fulcrum_encoder<fifi::binary8>,
	kodo_fulcrum::fulcrum_outer_decoder<fifi::binary8>>;

BENCHMARK_F_INLINE(setup_fulcrum_throughput_outer, FulcrumOuter, Binary, 3)
{
//	run_body_fulcrum_DTEP();
	run_body_fulcrum_DTEP_with_feedback();
}
//
using setup_fulcrum_throughput_combined =
	dynamic_both_throughput_benchmark<
    kodo_fulcrum::tunable_fulcrum_encoder<fifi::binary8>,
	fulcrum_combined_decoder<kodo_fulcrum::fulcrum_combined_decoder<fifi::binary8>>>;

BENCHMARK_F_INLINE(setup_fulcrum_throughput_combined, FulcrumCombined, Binary, 3)
{
//	run_body_fulcrum_DTEP();
	run_body_fulcrum_DTEP_with_feedback();
}
//------------------------------------------------------------------
// Dynamic both - DTEP Fulcrum - maximum sparity
//------------------------------------------------------------------
//using setup_fulcrum_throughput_inner =
//	dynamic_both_throughput_benchmark<
//    kodo_fulcrum::tunable_fulcrum_encoder<fifi::binary8>,
//	fulcrum_inner_decoder<kodo_fulcrum::fulcrum_inner_decoder<fifi::binary>>>;
//
//BENCHMARK_F_INLINE(setup_fulcrum_throughput_inner, FulcrumInner, Binary, 1)
//{
//	run_body_fulcrum_DTEP_max_sparsity();
//}
//
//using setup_fulcrum_throughput_outer =
//	dynamic_both_throughput_benchmark<
//    kodo_fulcrum::tunable_fulcrum_encoder<fifi::binary8>,
//	kodo_fulcrum::fulcrum_outer_decoder<fifi::binary8>>;
//
//BENCHMARK_F_INLINE(setup_fulcrum_throughput_outer, FulcrumOuter, Binary, 3)
//{
//	run_body_fulcrum_DTEP_max_sparsity();
//}
//
//using setup_fulcrum_throughput_combined =
//	dynamic_both_throughput_benchmark<
//    kodo_fulcrum::tunable_fulcrum_encoder<fifi::binary8>,
//	fulcrum_combined_decoder<kodo_fulcrum::fulcrum_combined_decoder<fifi::binary8>>>;
//
//BENCHMARK_F_INLINE(setup_fulcrum_throughput_combined, FulcrumCombined, Binary, 3)
//{
//	run_body_fulcrum_DTEP_max_sparsity();
//}

//------------------------------------------------------------------
// Dynamic both - STEP Fulcrum: change the expansion at the end the generation, while changing the sparsity in the middle
//------------------------------------------------------------------
//using setup_fulcrum_throughput_inner =
//	dynamic_both_throughput_benchmark<
//    kodo_fulcrum::tunable_fulcrum_encoder<fifi::binary8>,
//	fulcrum_inner_decoder<kodo_fulcrum::fulcrum_inner_decoder<fifi::binary>>>;
//
//BENCHMARK_F_INLINE(setup_fulcrum_throughput_inner, FulcrumInner, Binary, 1)
//{
////	run_body_fulcrum_STEP();
//	run_body_fulcrum_STEP_with_feedback();
//}
//
//using setup_fulcrum_throughput_outer =
//	dynamic_both_throughput_benchmark<
//    kodo_fulcrum::tunable_fulcrum_encoder<fifi::binary8>,
//	kodo_fulcrum::fulcrum_outer_decoder<fifi::binary8>>;
//
//BENCHMARK_F_INLINE(setup_fulcrum_throughput_outer, FulcrumOuter, Binary, 1)
//{
////	run_body_fulcrum_STEP();
//	run_body_fulcrum_STEP_with_feedback();
//}
//
//using setup_fulcrum_throughput_combined =
//	dynamic_both_throughput_benchmark<
//    kodo_fulcrum::tunable_fulcrum_encoder<fifi::binary8>,
//	fulcrum_combined_decoder<kodo_fulcrum::fulcrum_combined_decoder<fifi::binary8>>>;
//
//BENCHMARK_F_INLINE(setup_fulcrum_throughput_combined, FulcrumCombined, Binary, 1)
//{
////	run_body_fulcrum_STEP();
//	run_body_fulcrum_STEP_with_feedback();
//}

int main(int argc, const char* argv[])
{
    srand(static_cast<uint32_t>(time(0)));

    gauge::runner::add_default_printers();
    gauge::runner::run_benchmarks(argc, argv);

    return 0;
}
