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
//#include "operations_benchmark_outer.hpp"
#include "tunable_operations_benchmark_outer.hpp"
//#include "tunable_operations_benchmark_combined.hpp"
//#include "operations_benchmark_combined.hpp"
#include "operations_benchmark_inner.hpp"

#include <tables/table.hpp>

#include <kodo_fulcrum/fulcrum_combined_decoder_counter.hpp>
#include <kodo_fulcrum/fulcrum_codes.hpp>

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

/// Using this macro we may specify options. For specifying options
/// we use the boost program options library. So you may additional
/// details on how to do it in the manual for that library.
BENCHMARK_OPTION(count_operations_options)
{
    gauge::po::options_description options;

    std::vector<uint32_t> symbols;
    //symbols.push_back(16);
    //symbols.push_back(32);
    symbols.push_back(64);
    //symbols.push_back(128);
    //symbols.push_back(256);
    //symbols.push_back(512);
    symbols.push_back(1024);

    auto default_symbols =
        gauge::po::value<std::vector<uint32_t> >()->default_value(
            symbols, "")->multitoken();

    std::vector<uint32_t> symbol_size;
    symbol_size.push_back(1500);

    auto default_symbol_size =
        gauge::po::value<std::vector<uint32_t> >()->default_value(
            symbol_size, "")->multitoken();

    std::vector<std::string> types;
    //types.push_back("encoder");
    //types.push_back("recoder");
    types.push_back("decoder");


    auto default_types =
        gauge::po::value<std::vector<std::string> >()->default_value(
            types, "")->multitoken();


 std::vector<uint32_t> expansion;
    expansion.push_back(1);
//    expansion.push_back(2);
//    expansion.push_back(3);
//    expansion.push_back(4);
//    expansion.push_back(5);
//    expansion.push_back(6);
//    expansion.push_back(7);
//    expansion.push_back(8);
//    expansion.push_back(9);
//    expansion.push_back(10);


    auto default_expansion =
        gauge::po::value<std::vector<uint32_t> >()->default_value(
            expansion, "")->multitoken();


    std::vector<double> erasures;
    erasures.push_back(0);
//    erasures.push_back(0.2);
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

    std::vector<double> density;
//    density.push_back(0.1);
//    density.push_back(0.2);
//    density.push_back(0.3);
//    density.push_back(0.4);
      density.push_back(0.5);

    auto default_density =
        gauge::po::value<std::vector<double> >()->default_value(
            density, "")->multitoken();


    std::vector<uint32_t> threshold;
   threshold.push_back(0);
   auto default_threshold =
		   gauge::po::value<std::vector<uint32_t> >()->default_value(
			   threshold, "")->multitoken();

    options.add_options()
    ("symbols", default_symbols, "Set the number of symbols");

    options.add_options()
    ("symbol_size", default_symbol_size, "Set the symbol size in bytes");

    options.add_options()
    ("type", default_types, "Set type [encoder|decoder]");

   options.add_options()
    ("expansion", default_expansion,
     "Set the expansion of the fulcrum codes");


   options.add_options()
     ("threshold", default_threshold,
      "Set the threshold of the fulcrum codes");

    options.add_options()
    ("erasure", default_erasure,
     "Set the erasure of the fulcrum codes");
    options.add_options()
       ("density", default_density,
        "Set the density of the fulcrum codes");


    gauge::runner::instance().register_options(options);
}


////FULCRUM
//using setup_fulcrum_count_inner =
//    tunable_operations_benchmark_outer<
//    kodo_fulcrum::fulcrum_encoder<fifi::binary8>,
//    kodo_fulcrum::fulcrum_inner_decoder<fifi::binary>,
//	kodo_fulcrum::fulcrum_inner_decoder<fifi::binary>>;
//
//BENCHMARK_F_INLINE(setup_fulcrum_count_inner, FulcrumInner, Binary, 30)
//{
//    run_fulcrum();
//}
//
using setup_fulcrum_count_outer =
    tunable_operations_benchmark_outer<
    kodo_fulcrum::tunable_fulcrum_encoder<fifi::binary8>,
    kodo_fulcrum::fulcrum_inner_decoder<fifi::binary>,
	kodo_fulcrum::fulcrum_outer_decoder<fifi::binary8>>;

BENCHMARK_F_INLINE(setup_fulcrum_count_outer, FulcrumOuter, Binary8, 3)
{
    run_fulcrum();
}

//using setup_fulcrum_count_combined =
//	tunable_operations_benchmark_combined<
//    kodo_fulcrum::fulcrum_encoder<fifi::binary8>,
//    kodo_fulcrum::fulcrum_inner_decoder<fifi::binary>,
//    fulcrum_combined_decoder<kodo_fulcrum::fulcrum_combined_decoder_counter<fifi::binary8>>>;

//BENCHMARK_F_INLINE(setup_fulcrum_count_combined, FulcrumCombined, Binary8, 30)
//{
//    run_fulcrum();
//}


////DYNAMIC TUNABLE FULCRUM
//using setup_tunable_fulcrum_count_outer =
//    tunable_operations_benchmark_outer<
//    kodo_fulcrum::tunable_fulcrum_encoder<fifi::binary8>,
//    kodo_fulcrum::fulcrum_inner_decoder<fifi::binary>,
//	kodo_fulcrum::fulcrum_outer_decoder<fifi::binary8>>;
//
//BENCHMARK_F_INLINE(setup_tunable_fulcrum_count_outer, TunableFulcrumOuter, Binary8, 3)
//{
//	run_dynamic_tep();
//}

//using setup_tunable_fulcrum_count_combined =
//    tunable_operations_benchmark_combined<
//    kodo_fulcrum::tunable_fulcrum_encoder<fifi::binary8>,
//    kodo_fulcrum::fulcrum_inner_decoder<fifi::binary>,
//	fulcrum_combined_decoder<kodo_fulcrum::fulcrum_combined_decoder_counter<fifi::binary8>>>;
//
//BENCHMARK_F_INLINE(setup_tunable_fulcrum_count_combined, TunableFulcrumCombined, Binary8, 3)
//{
//	run_dynamic_tep();
//}

//STATIC TUNABLE FULCRUM
//using setup_tunable_fulcrum_count_outer =
//    tunable_operations_benchmark_outer<
//    kodo_fulcrum::tunable_fulcrum_encoder<fifi::binary8>,
//    kodo_fulcrum::fulcrum_inner_decoder<fifi::binary>,
//	kodo_fulcrum::fulcrum_outer_decoder<fifi::binary8>>;
//
//BENCHMARK_F_INLINE(setup_tunable_fulcrum_count_outer, TunableFulcrumOuter, Binary8, 3)
//{
//    run_static_tep();// cu truoc 14/12 run_region_improvement();
//}

//using setup_tunable_fulcrum_count_combined =
//    tunable_operations_benchmark_combined<
//    kodo_fulcrum::tunable_fulcrum_encoder<fifi::binary8>,
//    kodo_fulcrum::fulcrum_inner_decoder<fifi::binary>,
//	fulcrum_combined_decoder<kodo_fulcrum::fulcrum_combined_decoder_counter<fifi::binary8>>>;
//
//BENCHMARK_F_INLINE(setup_tunable_fulcrum_count_combined, TunableFulcrumCombined, Binary8, 3)
//{
//    run_static_tep();
//}



int main(int argc, const char* argv[])
{
    srand(static_cast<uint32_t>(time(0)));

    gauge::runner::add_default_printers();
    gauge::runner::run_benchmarks(argc, argv);

    return 0;
}
