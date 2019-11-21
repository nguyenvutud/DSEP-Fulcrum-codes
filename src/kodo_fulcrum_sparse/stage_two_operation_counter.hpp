#pragma once

#include <cstdint>

#include <kodo_core/operations_counter.hpp>

namespace kodo_fulcrum
{
    /// @ingroup debug
    ///
    /// This layer "intercepts" all calls to the finite_field_math
    /// layer counting the different operations
    template<class SuperCoder>
    class stage_two_operation_counter : public SuperCoder
    {
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
    	    	 	uint32_t decoder_rank(){
    	    	 		return stage_one_decoder_rank() + stage_two_decoder_rank();
    	    	 	}
    	/// @return The operation counter
    	        kodo_core::operations_counter stage_two_operations_counter()
    	        {
    	            return SuperCoder::stage_two_decoder().get_operations_counter();
    	        }
    	        kodo_core::operations_counter stage_one_operations_counter()
				{
					return SuperCoder::stage_one_decoder().get_operations_counter();
				}

    	        using factory = kodo_core::pool_factory<stage_two_operation_counter>;

    };
}
