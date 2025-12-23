// #include "occ.h"

// int main(int argc, char *argv[]) {
//     size_t num_expansions = 2; // max = 7 //used to be 8
//     size_t further_dims = 2; // v2
//     size_t total_n = (1 << num_expansions) * (1 << further_dims);

//     size_t IDX_TARGET = 14;
//     size_t IDX_DIM0 = IDX_TARGET / (1 << further_dims);


//     // cout << "======================ssssssssss======================\n";
//     #ifndef __EMSCRIPTEN__

//     // omp_set_num_threads(1);

//     build_table();

//     scratch = (uint64_t *)malloc(crt_count * poly_len * sizeof(uint64_t));

//     ntt_qprime = new intel::hexl::NTT(2048, arb_qprime);

//     bool ubench = false;
//     bool high_rate = false;
//     cout << "argc is:" << argc << endl;

//     if (argc > 1) {
//         num_expansions = strtol(argv[1], NULL, 10); // max = 7 //used to be 8
//         further_dims = strtol(argv[2], NULL, 10);
//         total_n = (1 << num_expansions) * (1 << further_dims);
//         IDX_TARGET = strtol(argv[3], NULL, 10);
//         IDX_DIM0 = IDX_TARGET / (1 << further_dims);
//     }

//     cout << "Start to test function GenRandVec" << endl;
//     std::vector<size_t> res = GenRandVec(0, 10, 11, 4);
//     for (size_t i = 0 ; i < res.size() ; i++) {
//         cout << i << ' ' << res[i] << endl;
//     }

//     // cout << "Start to test homomorphic operation" << endl;
//     // test_gaussian_elimination();
    
//     // cout << "Start to test sort columns by first non zero" << endl;
//     // test_sort_columns();

//     // cout << "Start to test to solve linear system" << endl;
//     // test_solve_linear_system();

//     // cout << "Start to test plaintext's ntt add" << endl;
//     // test_pt_add_ntt();
    
//     cout << "Start to test to use OCD" << endl;
//     test_mPIR_with_PBC_and_OCD(num_expansions, further_dims);


//     // cout << "Start to test to use OCC" << endl;
//     // test_mPIR_with_PBC_and_OCC(num_expansions, further_dims);

//     // test_mPIR_with_PBC(num_expansions, further_dims);
//     #endif
// }