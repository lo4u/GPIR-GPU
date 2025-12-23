// #include "okvs.h"

// void test_speed_of_3H_GCT_and_RB_Matrix(
//     size_t num_expansions, 
//     size_t further_dims,
//     size_t IDX_TARGET
// ) {
//     do_MatPol_test();

//     setup_constants();

//     generate_gadgets();

//     cout << "Look at this !!!" << endl;

//     uint64_t *B;

//     cout << random_data << endl;
//     cout << "the number of target index is:" << IDX_TARGET << endl;

//     std::vector<size_t> IDX_TARGETs(IDX_TARGET);
//     std::iota(IDX_TARGETs.begin(), IDX_TARGETs.end(), 1);

//     // 下面先进行一些准备工作
//     // Both: 新建布谷鸟编码
//     CuckooCode code(IDX_TARGETs.size(), 3, 1.5);

//     OracleTy god = getOracle(&code, num_expansions, further_dims);
//     auto godOracle = std::get<0>(god);
//     auto godSizes = std::get<1>(god);
//     // Client: 得到调度计划
//     auto schedule = PBC_GetSchedule(&code, IDX_TARGETs);
//     // Client: 创建一个新的哈希表，将 schedule 中的每个键映射到一个索引值
//     std::unordered_map<size_t, size_t> indexes;
//     for (const auto& entry : schedule) {
//         size_t key = entry.first;
//         assert(entry.second.size() == 1);
//         size_t bucket = entry.second[0];

//         // 在 oracle 中查找对应的索引
//         auto& oracle_bucket = godOracle[bucket];
//         auto it = std::find_if(oracle_bucket.begin(), oracle_bucket.end(), [key](const size_t& e) {
//             return e == key;
//         });

//         assert(it != oracle_bucket.end());
//         indexes[bucket] = std::distance(oracle_bucket.begin(), it);
//     }

//     // Client: 重构多请求向量
//     std::vector<size_t> real_ind;
//     std::vector<size_t> ind_vec;
//     for (size_t bucket = 0 ; bucket < godOracle.size() ; ++bucket) {
//         if (indexes.find(bucket) != indexes.end()) {
//             ind_vec.push_back(static_cast<uint32_t>(indexes[bucket]));
//             real_ind.push_back(bucket);
//         } else {
//             ind_vec.push_back(static_cast<uint32_t>(0));
//         }
//     }

//     for (size_t i = 0 ; i < ind_vec.size() ; i++) {
//         cout << ind_vec[i] << ' ';
//     }
//     cout << endl;

//     // Client: 创建 Client 类
//     MultiPirClient client_A(godSizes, ind_vec, godOracle, indexes);
//     // 创建密钥
//     client_A.setup();
//     // Client: 生成 packed_poly 向量
//     start_timing();
//     std::vector<MatPoly> packed_pt_vec = client_A.MultiPackedIndexGen();
//     cout << "Time used to do pack is " << end_timing() << endl;

//     // Create OKVS class
//     GCTObvKVStore client_OKVS(real_ind.size(), 3, 1.3, packed_pt_vec.size());
//     cout << "================== Benchmark ==================" << endl;

//     // test: the speed of OKVS
//     start_timing();
//     auto pt_hat = client_OKVS.Encode(packed_pt_vec, real_ind);
//     double time_3H_GCT_Compress = end_timing();
//     cout << "Time used to do OCD.Compress with 3H-GCT is " << time_3H_GCT_Compress << endl;

//     start_timing();
//     auto cvs_hat = client_A.MultiPackedIndexEncrypt(pt_hat); // In NTT
//     double time_encrypt = end_timing();
//     cout << "ALL Time used to do OCD.Compress with 3H-GCT is " << time_3H_GCT_Compress + time_encrypt << endl;
    
//     GCTObvKVStore server_OKVS(real_ind.size(), 3, 1.3, packed_pt_vec.size());
//     start_timing();
//     auto cvs = server_OKVS.Decode(cvs_hat);
//     double time_3H_GCT_Decompress = end_timing();
//     cout << "Time used to do OCD.Decompress with 3H-GCT is " << time_3H_GCT_Decompress << endl;
//     cout << "Total time used to do OCD with 3H-GCT is " << time_3H_GCT_Decompress + time_3H_GCT_Compress + time_encrypt << endl;

//     // // Create OCD class
//     // size_t w = 8 + ceil(log2(real_ind.size()));
//     // LSObvDecompress client_OCD(packed_pt_vec.size(), real_ind.size(), 0.05, w);

//     // // test: the speed of ocd
//     // start_timing();
//     // auto p_hat = client_OCD.Compress(packed_pt_vec, real_ind);
//     // double time_RB_Matrix = end_timing();
//     // cout << "Time used to do OCD with random bond matrix is " << time_RB_Matrix << endl;

//     // start_timing();
//     // auto cvs_hat = client_A.MultiPackedIndexEncrypt(pt_hat);
//     // cout << "Time used to encrypt is " << end_timing() << endl;
// }

// void test_for_noice_of_encrypt(
//     size_t num_expansions,
//     size_t further_dims,
//     size_t IDX_TARGET
// ) {
//     do_MatPol_test();

//     setup_constants();

//     generate_gadgets();

//     cout << "Look at this !!!" << endl;

//     uint64_t *B;

//     cout << random_data << endl;
//     cout << "the number of target index is:" << IDX_TARGET << endl;

//     std::vector<size_t> IDX_TARGETs(IDX_TARGET);
//     std::iota(IDX_TARGETs.begin(), IDX_TARGETs.end(), 1);

//     // 下面先进行一些准备工作
//     // Both: 新建布谷鸟编码
//     CuckooCode code(IDX_TARGETs.size(), 3, 1.5);

//     OracleTy god = getOracle(&code, num_expansions, further_dims);
//     auto godOracle = std::get<0>(god);
//     auto godSizes = std::get<1>(god);
//     // Client: 得到调度计划
//     auto schedule = PBC_GetSchedule(&code, IDX_TARGETs);
//     // Client: 创建一个新的哈希表，将 schedule 中的每个键映射到一个索引值
//     std::unordered_map<size_t, size_t> indexes;
//     for (const auto& entry : schedule) {
//         size_t key = entry.first;
//         assert(entry.second.size() == 1);
//         size_t bucket = entry.second[0];

//         // 在 oracle 中查找对应的索引
//         auto& oracle_bucket = godOracle[bucket];
//         auto it = std::find_if(oracle_bucket.begin(), oracle_bucket.end(), [key](const size_t& e) {
//             return e == key;
//         });

//         assert(it != oracle_bucket.end());
//         indexes[bucket] = std::distance(oracle_bucket.begin(), it);
//     }

//     // Client: 重构多请求向量
//     std::vector<size_t> real_ind;
//     std::vector<size_t> ind_vec;
//     for (size_t bucket = 0 ; bucket < godOracle.size() ; ++bucket) {
//         if (indexes.find(bucket) != indexes.end()) {
//             ind_vec.push_back(static_cast<uint32_t>(indexes[bucket]));
//             real_ind.push_back(bucket);
//         } else {
//             ind_vec.push_back(static_cast<uint32_t>(0));
//         }
//     }

//     for (size_t i = 0 ; i < ind_vec.size() ; i++) {
//         cout << ind_vec[i] << ' ';
//     }
//     cout << endl;

//     // Client: 创建 Client 类
//     MultiPirClient client_A(godSizes, ind_vec, godOracle, indexes);
//     // 创建密钥
//     client_A.setup();

// }

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

//     test_GCT_okvs(num_expansions, further_dims);

//     // cout << "Start to test function GenRandVec" << endl;
//     // std::vector<size_t> res = GenRandVec(0, 10, 11, 4);
//     // for (size_t i = 0 ; i < res.size() ; i++) {
//     //     cout << i << ' ' << res[i] << endl;
//     // }

//     // cout << "Start to test to do OKVS" << endl;
//     // test_speed_of_3H_GCT_and_RB_Matrix(num_expansions, further_dims, IDX_TARGET);

//     // test_mPIR_with_PBC(num_expansions, further_dims);
//     #endif
// }