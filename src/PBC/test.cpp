#include "cuckoo.h"
#include "spiral.h"

// int main(int argc, char *argv[]) {
//     cout << "Hello World\n";
//     cout << "This is a test file for PBC" << endl;
//     return 0;
// }

// int main(int argc, char *argv[]) {
//     size_t num_expansions = 2; // max = 7 //used to be 8
//     size_t further_dims = 2; // v2
//     size_t total_n = (1 << num_expansions) * (1 << further_dims);

//     size_t IDX_TARGET = 14;
//     size_t IDX_DIM0 = IDX_TARGET / (1 << further_dims);


//     cout << "======================ssssssssss======================\n";
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
//         if (argc >= 5) {
//             has_file = true;
//             dbFilename = argv[4];
//             if (dbFilename.length() <= 2) has_file = false;
//             for (size_t i = 5; i < argc; i++) {
//                 if (strcmp(argv[i], "--input") == 0) {
//                     load = true;
//                 }
//                 if (strcmp(argv[i], "--server") == 0) {
//                     server = true;
//                 }
//                 if (strcmp(argv[i], "--client") == 0) {
//                     client = true;
//                 }
//                 if (strcmp(argv[i], "--nonoise") == 0) {
//                     cout << "Using no noise" << endl;
//                     nonoise = true;
//                 }
//                 if (strcmp(argv[i], "--ubench") == 0) {
//                     cout << "Microbenchmarking..." << endl;
//                     ubench = true;
//                 }
//                 if (strcmp(argv[i], "--high-rate") == 0) {
//                     cout << "Using high rate variant..." << endl;
//                     high_rate = true;
//                 }
//                 if (strcmp(argv[i], "--random-data") == 0) {
//                     cout << "Using random data..." << endl;
//                     random_data = true;
//                     dummyWorkingSet = min((1UL << 25) / (total_n), poly_len);
//                     max_trials = (1L << 16) / (total_n);
//                     if (max_trials == 0) {
//                         max_trials = 1;
//                     }
//                 }
//                 if (strcmp(argv[i], "--show-diff") == 0) {
//                     cout << "Showing diff..." << endl;
//                     show_diff = true;
//                 }
//                 if (strcmp(argv[i], "--output-err") == 0) {
//                     cout << "Outputting errors..." << endl;
//                     outputErrFilename = argv[++i];
//                     output_err = true;
//                 }
//                 if (strcmp(argv[i], "--direct-upload") == 0) {
//                     cout << "Direct uploading of query (no compression)" << endl;
//                     direct_upload = true;
//                 }
//                 if (strcmp(argv[i], "--loaddata") == 0) {
//                     has_data = true;
//                     dataFilename = argv[i+1];
//                     i++;
//                 }
//             }
//         }
//     }

//     if (has_file) {
//         if (load) {
//             cout << "loading from file" << endl;
//             fs.open(dbFilename, fstream::in | fstream::binary);
//         } else {
//             cout << "outputting to file" << endl;
//             fs.open(dbFilename, fstream::out | fstream::binary);
//         }
//     }

//     if (has_data) {
//         cout << "loading data from file" << endl;
//         fsData.open(dataFilename, fstream::in | fstream::binary);
//     }

//     do_MatPol_test();

//     setup_constants();

//     generate_gadgets();

//     cout << "Look at this !!!" << endl;
//     if (high_rate) {
//         cout << "target index is:" << IDX_TARGET << endl;
//         testHighRate(num_expansions, further_dims, IDX_TARGET);
//         exit(0);
//     }

//     uint64_t *B;

//     load_db(B, num_expansions, further_dims, IDX_TARGET);
//     // for (size_t i = 0; i < pt_real.rows * pt_real.cols * coeff_count; i++) {
//     //     cout << pt_real.data[i] << endl;
//     // }
//     cout << random_data << endl;
//     cout << "target index is:" << IDX_TARGET << endl;

//     do_test(B, num_expansions, further_dims, IDX_TARGET);
//     #endif
// }








// uint64_t *DB_test;
// uint64_t maskDB = (1ULL << 32) - 1; // 使用1ULL来确保是64位的常量
// MatPoly pts_encd_test(n0, n2);
// MatPoly pt_encd_correct_test(n0, n2);
// MatPoly pt_real_test(n0, n2, false);
// MatPoly pt_test(n0, n0);
// MatPoly pt_correct_test(n0, n0);

// void print_pt_crt(MatPoly pts_encd_test) {
//     // uint64_t *BB = (uint64_t *)malloc(n0 * n2 * crt_count * poly_len * sizeof(uint64_t));
//     for (size_t m = 0; m < n0; m++) {
//         for (size_t c = 0; c < n2; c++) {
//             // memcpy(BB, &pts_encd_test.data[(m * n2 + c) * crt_count * coeff_count], crt_count * coeff_count * sizeof(uint64_t));
//             for (size_t z = 0; z < poly_len; z++) {
//                 cout << pts_encd_test.data[(m * n2 + c) * crt_count * coeff_count + z] << ' ' << pts_encd_test.data[(m * n2 + c) * crt_count * coeff_count + z + poly_len] << ' ';
//             }
//         }
//     }
//     cout << endl;
// }


// void genDB() {
//     size_t dim0 = 1 << 3; // 第一维度的总数
//     size_t num_per = (1 << 5) / dim0; // 数据库总数据量 / 第一维度
    
//     size_t num_bytes_B = sizeof(uint64_t) * dim0 * num_per * n0 * n2 * poly_len; // 2 * poly_len;
//     cout << "num_bytes_B: " << num_bytes_B << endl;
//     DB_test = (uint64_t *)aligned_alloc(64, num_bytes_B);
//     memset(DB_test, 0, num_bytes_B);

//     uint64_t *BB = (uint64_t *)malloc(n0 * n2 * crt_count * poly_len * sizeof(uint64_t));
//     size_t numBytesPlaintextRaw = n0 * n0 * num_bits_q * poly_len / 8;
//     uint64_t *pt_raw = (uint64_t *)malloc(numBytesPlaintextRaw);
//     memset(pt_raw, 0, numBytesPlaintextRaw);
//     uint64_t *pt_buf = (uint64_t *)malloc(n0 * n0 * crt_count * poly_len * sizeof(uint64_t));
//     memset(pt_buf, 0, n0 * n0 * crt_count * poly_len * sizeof(uint64_t));
    
//     cout << "starting generation of db" << endl;
//     MatPoly H_nttd = to_ntt(H_mp);
//     uint64_t *H_encd = H_nttd.data;
//     MatPoly pt_tmp(n0, n0, false);
//     MatPoly pt_encd_raw(n0, n2, false);
//     pts_encd_test = MatPoly(n0, n2);
//     pt_test = MatPoly(n0, n0);
//     for (size_t i = 0 ; i < (1 << 5) ; i++) {
//         generate_random_pt(pt_tmp);

//         pt_encd_raw = pt_tmp;
//         for (size_t pol = 0; pol < n0 * n2 * poly_len; pol++) {
//             int64_t val = (int64_t) pt_encd_raw.data[pol];
//             assert(val >= 0 && val < p_db);
//             if (val >= (p_db / 2)) {
//                 val = val - (int64_t)p_db;
//             }
//             if (val < 0) {
//                 val += Q_i;
//             }
//             assert(val >= 0 && val < Q_i);
//             pt_encd_raw.data[pol] = val;
//         }
//         to_ntt(pts_encd_test, pt_encd_raw);
//         print_pt_crt(pts_encd_test);

//         if (i == 3) {
//             cop(pt_encd_correct_test, pts_encd_test);
//             cop(pt_real_test, pt_tmp);
//             to_ntt(pt_test, pt_tmp);
//             cop(pt_correct_test, pt_test);
//         }

//         // b': i c n z j m
//         size_t ii = i % num_per;
//         size_t j = i / num_per;
//         for (size_t m = 0; m < n0; m++) {
//             for (size_t c = 0; c < n2; c++) {
//                 memcpy(BB, &pts_encd_test.data[(m * n2 + c) * crt_count * coeff_count], crt_count * coeff_count * sizeof(uint64_t));
//                 for (size_t z = 0; z < poly_len; z++) {
//                     size_t idx = z * (num_per * n2 * dim0 * n0) +
//                                  ii * (n2 * dim0 * n0) + 
//                                  c * (dim0 * n0) +
//                                  j * (n0) + m;
//                     DB_test[idx] = BB[z] | (BB[poly_len + z] << 32);
//                 }
//             }
//         }
//     }
//     free(BB);
//     cout << "done loading/generating db." << endl;
// }

// void print_item_of_DB(uint64_t *database, size_t index) {
//     size_t dim0 = 8;
//     size_t num_per = 4;
//     size_t ii = index % num_per;
//     size_t j = index / num_per;
//     for (size_t m = 0 ; m < n0; m++) {
//         for (size_t c = 0 ; c < n2 ; c++) {
//             for (size_t z = 0 ; z < poly_len ; z++) {
//                 size_t idx = z * (num_per * n2 * dim0 * n0) +
//                              ii * (n2 * dim0 * n0) + 
//                              c * (dim0 * n0) +
//                              j * (n0) + m;
//                 cout << (database[idx] & maskDB) << ' ' << ((database[idx] >> 32) & maskDB) << ' ';
//             }
//         }
//     }
//     cout << endl;
// }

// void print_DB(uint64_t *database) {
//     size_t dim0 = 8;
//     size_t num_per = 4;

//     for (size_t i = 0 ; i < (1 << 5) ; i++) {
//         cout << "item " << i;
//         size_t ii = i % num_per;
//         size_t j = i / num_per;
//         cout << " is ";
//         for (size_t m = 0 ; m < n0; m++) {
//             for (size_t c = 0 ; c < n2 ; c++) {
//                 for (size_t z = 0 ; z < poly_len ; z++) {
//                     size_t idx = z * (num_per * n2 * dim0 * n0) +
//                                  ii * (n2 * dim0 * n0) + 
//                                  c * (dim0 * n0) +
//                                  j * (n0) + m;
//                     cout << (database[idx] & maskDB) << ' ' << ((database[idx] >> 32) & maskDB) << ' ';
//                 }
//             }
//         }
//         cout << endl;
//     }
// }

// void print_bucket(std::vector<uint64_t *> collections) {
//     for (size_t i = 0 ; i < collections.size() ; i++) {
//         for (size_t m = 0 ; m < n0; m++) {
//             for (size_t c = 0 ; c < n2 ; c++) {
//                 for (size_t z = 0 ; z < poly_len ; z++) {
//                     size_t idx = m * n2 * poly_len + c * poly_len + z;
//                     cout << (collections[i][idx] & maskDB) << ' ' << ((collections[i][idx] >> 32) & maskDB) << ' ';
//                 }
//             }
//         }
//         cout << endl;
//     }
// }

// void PBC_test(CuckooCode *code) {
//     auto collections = PBC_Encode(code, DB_test);

//     // 生成键
//     std::vector<size_t> keys;
//     for (size_t i = 0; i < code->k; ++i) {
//         keys.push_back(i);
//     }

//     // 得到映射
//     auto schedule = PBC_GetSchedule(code, keys);

//     // 验证调度计划是否有效
//     for (const auto& pair : schedule) {
//         size_t key = pair.first;
//         const auto& buckets = pair.second;
//         assert(buckets.size() == 1);
//         cout << "================== Start Check the " << key << " Buckets =========================" << endl;
//         print_bucket(collections[buckets[0]]);
//         cout << "====================================================================" << endl;
//         print_item_of_DB(DB_test, key);
//         cout << "================== Finish Check the " << key << " Buckets =========================" << endl;
//     }

// }



// int main(int argc, char *argv[]) {
//     cout << "Hello World\n";
//     cout << "This is a test file for PBC" << endl;
//     cout << "=================== Start the Check for DB ===================" << endl;
//     genDB();
//     print_DB(DB_test);
//     cout << "=================== Finish the Check for DB ===================" << endl;
//     size_t test_times = 1;

//     std::mt19937_64 rng(std::random_device{}());
//     for (size_t i = 0 ; i < test_times ; i++) {
//         size_t k = 12 + i + std::uniform_int_distribution<size_t>(0, 7)(rng);
//         CuckooCode code(k, 3, 1.3);
//         PBC_test(&code);
//     }

//     return 0;
// }


