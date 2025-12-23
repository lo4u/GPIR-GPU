#include "okvs.h"
#include <chrono>
#include <filesystem>

void test_GCT_okvs(
    size_t num_expansions, 
    size_t further_dims
) {
    TestResult testResult;

    do_MatPol_test();

    setup_constants();

    generate_gadgets();

    cout << "Look at this !!!" << endl;

    uint64_t *B;

    cout << random_data << endl;
    cout << "target index is:" << IDX_TARGET << endl;

    std::vector<size_t> IDX_TARGETs(16);
    std::iota(IDX_TARGETs.begin(), IDX_TARGETs.end(), 1);

    if (std::filesystem::exists("./database")) {
        size_t length;
        loadBinary(B, length, "./database");
        assert(length == (1 << num_expansions << further_dims) * n0 * n2 * poly_len);
    } else {
        size_t length = (1 << num_expansions << further_dims) * n0 * n2 * poly_len;
        genDataBase(B, num_expansions, further_dims, IDX_TARGETs);
        saveBinary(B, length, "./database");
    }

    // print_DB(B, num_expansions, further_dims);

    auto startEntire = std::chrono::high_resolution_clock::now();

    // 下面先进行一些准备工作
    // Both: 新建布谷鸟编码
    CuckooCode code(IDX_TARGETs.size(), 3, 1.5);
    // Server: 对数据库进行 PBC 编码
    auto collections = PBC_Encode(&code, B, num_expansions, further_dims);
    // Client: 获取神谕
    OracleTy god = getOracle(&code, num_expansions, further_dims);
    auto godOracle = std::get<0>(god);
    auto godSizes = std::get<1>(god);
    // Client: 得到调度计划
    auto schedule = PBC_GetSchedule(&code, IDX_TARGETs);
    // Client: 创建一个新的哈希表，将 schedule 中的每个键映射到一个索引值
    std::unordered_map<size_t, size_t> indexes;
    for (const auto& entry : schedule) {
        size_t key = entry.first;
        assert(entry.second.size() == 1);
        size_t bucket = entry.second[0];

        // 在 oracle 中查找对应的索引
        auto& oracle_bucket = godOracle[bucket];
        auto it = std::find_if(oracle_bucket.begin(), oracle_bucket.end(), [key](const size_t& e) {
            return e == key;
        });

        assert(it != oracle_bucket.end());
        indexes[bucket] = std::distance(oracle_bucket.begin(), it);
    }
    // Server: 重建 1.5k 个桶, 变为可用于 Spiral 查询的形式
    std::vector<uint64_t *> DBs(collections.size());
    size_t cLen = collections.size();
    for (size_t i = 0 ; i < cLen ; i++) {
        cout << "=================== Start to Rebuild the " << i << " Bucket ===================" << endl;
        DBs[i] = NULL;
        // print_bucket(collections[i]);
        reBuildDB(DBs[i], godSizes[i].num_expansions, godSizes[i].further_dims, collections[i]);
        // cout << "==============================================================" << endl;
        // print_DB(DBs[i], godSizes[i].num_expansions, godSizes[i].further_dims);
        cout << "=================== Finish to Rebuild the " << i << " Bucket ===================" << endl;
    }

    // Client: 重构多请求向量
    std::vector<size_t> real_ind;
    std::vector<size_t> ind_vec;
    for (size_t bucket = 0 ; bucket < godOracle.size() ; ++bucket) {
        if (indexes.find(bucket) != indexes.end()) {
            ind_vec.push_back(static_cast<uint32_t>(indexes[bucket]));
            real_ind.push_back(bucket);
        } else {
            ind_vec.push_back(static_cast<uint32_t>(0));
        }
    }

    for (size_t i = 0 ; i < ind_vec.size() ; i++) {
        cout << ind_vec[i] << ' ';
    }
    cout << endl;

    // Client: 创建 Client 类
    MultiPirClient client_A(godSizes, ind_vec, godOracle, indexes);
    // 创建密钥
    client_A.setup();
    // Client: 生成 packed_poly 向量
    std::vector<MatPoly> packed_pt_vec = client_A.MultiPackedIndexGen();

    // Create OKVS class
    GCTObvKVStore client_OKVS(real_ind.size(), 3, 1.3, packed_pt_vec.size());
    
    auto startCompress = std::chrono::high_resolution_clock::now();
    // test: the correctness of OKVS
    auto pt_hat = client_OKVS.Encode(packed_pt_vec, real_ind);
    auto endCompress = std::chrono::high_resolution_clock::now();
    testResult.compressTime = std::chrono::duration<double, std::milli>(endCompress - startCompress).count();

    cout << "==============================" << endl;

    auto cvs_hat = client_A.MultiPackedIndexEncrypt(pt_hat); // In NTT
    testResult.querySize = calculateMatPolySize(cvs_hat[0]) * cvs_hat.size();

    auto startAnswer = std::chrono::high_resolution_clock::now();

    GCTObvKVStore server_OKVS(real_ind.size(), 3, 1.3, packed_pt_vec.size());

    auto startDecompress = std::chrono::high_resolution_clock::now();
    auto cvs = server_OKVS.Decode(cvs_hat);
    auto endDecompress = std::chrono::high_resolution_clock::now();
    testResult.decompressTime = std::chrono::duration<double, std::milli>(endDecompress - startDecompress).count();

    cout << "==============================" << endl;
    MultiPirServer server_A(godSizes, DBs, cvs);
    // cout << cvs.size() << endl;
    // cout << server_A.handles.size() << endl;
    std::vector<MatPoly> rs = server_A.MultiServerAnswer();

    testResult.responseSize = calculateMatPolySize(rs[0]) * rs.size();

    auto endAnswer = std::chrono::high_resolution_clock::now();
    testResult.answerTime = std::chrono::duration<double, std::milli>(endAnswer - startAnswer).count();

    // cout << "==============================" << endl;
    std::vector<MatPoly> pts = client_A.MultiAnswerDecrypt(rs);

    auto endEntire = std::chrono::high_resolution_clock::now();
    testResult.entireTime = std::chrono::duration<double, std::milli>(endEntire - startEntire).count();


    for (size_t i = 0 ; i < DBs.size() ; i++) {
        if (indexes.find(i) != indexes.end()) {
            size_t seq = 0;
            for (const auto& entry : schedule) {
                size_t key = entry.first;
                assert(entry.second.size() == 1);
                size_t bucket = entry.second[0];
                if (bucket == i) {
                    seq = key;
                    break;
                }
            }
            double log_var = check_relation(rs[i], modswitch_on_server, seq-1);
        }
    }

    testResult.print();

    // test_trival_multi_PIR(B, num_expansions, further_dims, IDX_TARGETs);
}