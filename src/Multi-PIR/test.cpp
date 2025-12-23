// #include <tuple>
// #include "multipir.h"

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
//     }

//     do_MatPol_test();

//     setup_constants();

//     generate_gadgets();

//     cout << "Look at this !!!" << endl;

//     uint64_t *B;

//     cout << random_data << endl;
//     cout << "target index is:" << IDX_TARGET << endl;

//     // std::vector<size_t> IDX_TARGETs;
//     // IDX_TARGETs.push_back(1ULL);

//     // size_t dim0 = 1 << num_expansions; // 第一维度的总数
//     // size_t num_per = 1 << further_dims; // 数据库总数据量 / 第一维度
    
//     // size_t num_bytes_B = sizeof(uint64_t) * dim0 * num_per * n0 * n2 * poly_len; // 2 * poly_len;
//     // cout << "num_bytes_B: " << num_bytes_B << endl;
//     // DB_test = (uint64_t *)aligned_alloc(64, num_bytes_B);
//     // memset(DB_test, 0, num_bytes_B);


//     // genDataBase(B, num_expansions, further_dims, IDX_TARGETs);
    
//     // test_oo(DB_test, num_expansions, further_dims, IDX_TARGET);
//     // test_oo(B, num_expansions, further_dims, IDX_TARGET);
//     // test_oo(DB_test, num_expansions, further_dims, IDX_TARGET);

//     std::vector<size_t> IDX_TARGETs;
//     IDX_TARGETs.push_back(1ULL);
//     IDX_TARGETs.push_back(2ULL);
//     IDX_TARGETs.push_back(3ULL);
//     IDX_TARGETs.push_back(4ULL);
//     IDX_TARGETs.push_back(5ULL);


//     genDataBase(B, num_expansions, further_dims, IDX_TARGETs);



//     // print_DB(B, num_expansions, further_dims);

//     // 下面先进行一些准备工作
//     // Both: 新建布谷鸟编码
//     CuckooCode code(IDX_TARGETs.size(), 3, 1.5);
//     // Server: 对数据库进行 PBC 编码
//     auto collections = PBC_Encode(&code, B, num_expansions, further_dims);
//     // Client: 获取神谕
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
//     // Server: 重建 1.5k 个桶, 变为可用于 Spiral 查询的形式
//     std::vector<uint64_t *> DBs(collections.size());
//     size_t cLen = collections.size();
//     for (size_t i = 0 ; i < cLen ; i++) {
//         cout << "=================== Start to Rebuild the " << i << " Bucket ===================" << endl;
//         DBs[i] = NULL;
//         // print_bucket(collections[i]);
//         reBuildDB(DBs[i], godSizes[i].num_expansions, godSizes[i].further_dims, collections[i]);
//         // cout << "==============================================================" << endl;
//         // print_DB(DBs[i], godSizes[i].num_expansions, godSizes[i].further_dims);
//         cout << "=================== Finish to Rebuild the " << i << " Bucket ===================" << endl;
//     }
//     // print_item_of_DB();
//     // // 验证 PBC 编码的正确性
//     // for (size_t i = 0 ; i < IDX_TARGETs.size() ; i++) {
//     //     size_t bucket = schedule[IDX_TARGETs[i]][0];
//     //     cout << IDX_TARGETs[i] << " \'s" << " bucket is " << bucket << " and its bucket index is " << indexes[bucket] << endl;
//     //     print_item_of_bucket(collections[bucket], indexes[bucket]);
//     //     print_item_of_DB(DBs[bucket], indexes[bucket], godSizes[bucket].num_expansions, godSizes[bucket].further_dims);
//     //     print_item_of_DB(B, IDX_TARGETs[i], num_expansions, further_dims);
//     // }
//     // Client: 重构多请求向量
//     std::vector<size_t> ind_vec;
//     for (size_t bucket = 0 ; bucket < godOracle.size() ; ++bucket) {
//         if (indexes.find(bucket) != indexes.end()) {
//             ind_vec.push_back(static_cast<uint32_t>(indexes[bucket]));
//         } else {
//             ind_vec.push_back(static_cast<uint32_t>(0));
//         }
//     }
//     // free 掉 collections
//     // for (size_t i = 0 ; i < collections.size() ; i++) {
//     //     for (size_t j = 0 ; j < collections[i].size() ; j++) {
//     //         free(collections[i][j]);
//     //     }
//     // }
//     for (size_t i = 0 ; i < ind_vec.size() ; i++) {
//         cout << ind_vec[i] << ' ';
//     }
//     cout << endl;

//     // Client: 创建 Client 类
//     MultiPirClient client_A(godSizes, ind_vec);
//     // 创建密钥
//     client_A.setup();
//     // Client: 发送请求
//     std::vector<MatPoly> cvs = client_A.MultiClientQuery();
//     // Server: 创建 Server 类
//     MultiPirServer server_A(godSizes, DBs, cvs);
//     // Server: 响应请求
//     start_timing();
//     std::vector<MatPoly> rs = server_A.MultiServerAnswer();
//     cout << "======================= Time used with PBC is " << end_timing() << "============================" << endl;
//     // PirClient client_B(ind_vec[0], godSizes[0].num_expansions, godSizes[0].further_dims);
//     // auto cv = client_B.ClientQuery();
//     // PirServer server_B(DBs[0], cv, godSizes[0].num_expansions, godSizes[0].further_dims);
//     // server_B.genQuery();
//     // auto r = server_B.ServerAnswer();
//     // check_corr(server_B.furtherDimsLocals, modswitch_on_server, 4);


//     for (size_t i = 0 ; i < DBs.size() ; i++) {
//         if (indexes.find(i) != indexes.end()) {
//             size_t seq = 0;
//             for (const auto& entry : schedule) {
//                 size_t key = entry.first;
//                 assert(entry.second.size() == 1);
//                 size_t bucket = entry.second[0];
//                 if (bucket == i) {
//                     seq = key;
//                     break;
//                 }
//             }
//             double log_var = check_relation(rs[i], modswitch_on_server, seq-1);
//         }
//     }
//     test_trival_multi_PIR(B, num_expansions, further_dims, IDX_TARGETs);



//     // for (size_t i = 0 ; i < IDX_TARGETs.size() ; i++) {
//     //     cout << "item " << IDX_TARGETs[i] << " is in the " << schedule[IDX_TARGETs[i]][0] << " bucket and its bucket index is " << indexes[schedule[IDX_TARGETs[i]][0]] << endl;
//     // }

//     // for (size_t i = 0 ; i < ind_vec.size(); i++) {
//     //     cout << "indexs[" << i << "] is " << ind_vec[i] << endl;
//     // }


//     // test_trival_multi_PIR(B, num_expansions, further_dims, IDX_TARGETs);

//     // test_oo(B, num_expansions, further_dims, IDX_TARGET);
//     #endif
// }

// // int main(int argc, char *argv[]) {
// //     cout << "Hello World\n";
// //     cout << "This is a test file for PBC" << endl;
// //     cout << "=================== Start the Check for DB ===================" << endl;
// //     genDB(3, 2);
// //     print_DB(DB_test, 3, 2);
// //     cout << "=================== Finish the Check for DB ===================" << endl;
// //     size_t test_times = 1;

// //     std::mt19937_64 rng(std::random_device{}());
// //     for (size_t i = 0 ; i < test_times ; i++) {
// //         size_t k = 12 + i + std::uniform_int_distribution<size_t>(0, 7)(rng);
// //         CuckooCode code(k, 3, 1.3);
// //         PBC_test(&code);
// //         // test_large_db(&code);
// //     }

// //     return 0;
// // }