#include "okvs.h"
#include "cuda_runtime_api.h"

void trival_test_single_pir_with_gpu(
    size_t num_expansions, 
    size_t further_dims,
    size_t IDX_TARGET
);

void trival_test_single_pir_with_cpu(
    size_t num_expansions, 
    size_t further_dims,
    size_t IDX_TARGET
);

void test_batch_pir_with_pbc_okvs_gpu(
    size_t num_expansions, 
    size_t further_dims,
    size_t query_num
);

void test_batch_pir_with_pbc_okvs_cpu(
    size_t num_expansions, 
    size_t further_dims,
    size_t query_num
);

void test_speed_of_3H_GCT(
    size_t num_expansions, 
    size_t further_dims,
    size_t IDX_TARGET
);

void test_speed_of_without_okvs(
    size_t num_expansions, 
    size_t further_dims,
    size_t IDX_TARGET
);
