#include "parallel.h"
#include "opgpu.cuh"

// int main(int argc, char *argv[]) {
//     // Timer clock;
//     // clock.start_gpu();
//     // print_hello_from_gpu();
//     // clock.stop_gpu();
//     // clock.duration_gpu("Print Hello World from GPU");

//     // return 0;
//     // test_MatAdd();
//     // return 0;

//     size_t num_expansions = 2; // max = 7 //used to be 8
//     size_t further_dims = 2; // v2
//     size_t total_n = (1 << num_expansions) * (1 << further_dims);

//     size_t IDX_TARGET = 14;
//     size_t IDX_DIM0 = IDX_TARGET / (1 << further_dims);


//     // cout << "======================ssssssssss======================\n";
//     #ifndef __EMSCRIPTEN__

//     // omp_set_num_threads(1);

//     build_table();

//     // scratch = (uint64_t *)malloc(crt_count * poly_len * sizeof(uint64_t));

//     ntt_qprime = new intel::hexl::NTT(2048, arb_qprime);

//     bool ubench = false;
//     bool high_rate = false;
//     cout << "argc is:" << argc << endl;

//     if (argc > 1) {
//         num_expansions = strtol(argv[1], NULL, 10); // max = 7 //used to be 8
//         further_dims = strtol(argv[2], NULL, 10);
//         total_n = (1 << num_expansions) * (1 << further_dims);
//         IDX_TARGET = strtol(argv[3], NULL, 10);
//     }

//     // test_GCT_okvs(num_expansions, further_dims);
//     trival_test_single_pir_with_cpu(num_expansions, further_dims, IDX_TARGET);
//     // trival_test_single_pir_with_gpu(num_expansions, further_dims, IDX_TARGET);
//     // test_batch_pir_with_pbc_okvs_gpu(num_expansions, further_dims, IDX_TARGET);
//     // test_batch_pir_with_pbc_okvs_cpu(num_expansions, further_dims, IDX_TARGET);
//     // test_speed_of_3H_GCT(num_expansions, further_dims, IDX_TARGET);
//     // test_speed_of_without_okvs(num_expansions, further_dims, IDX_TARGET);
//     #endif
// }