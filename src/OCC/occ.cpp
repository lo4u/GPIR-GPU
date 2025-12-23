#include "occ.h"

MatPoly genAMatWithV(size_t v) {
    MatPoly A(1, 1, false);
    for (size_t i = 0 ; i < A.rows * A.cols * coeff_count ; i++) {
        A.data[i] = v;
    }
    return A;
}

MatPoly genBMatWithV(size_t v) {
    MatPoly A(n0, n2, false);
    for (size_t i = 0 ; i < A.rows * A.cols * coeff_count ; i++) {
        A.data[i] = v;
    }
    return A;
}

void test_mPIR_with_PBC(
    size_t num_expansions, 
    size_t further_dims
) {
    do_MatPol_test();

    setup_constants();

    generate_gadgets();

    cout << "Look at this !!!" << endl;

    uint64_t *B;

    cout << random_data << endl;
    cout << "target index is:" << IDX_TARGET << endl;

    std::vector<size_t> IDX_TARGETs;
    IDX_TARGETs.push_back(1ULL);
    IDX_TARGETs.push_back(2ULL);
    IDX_TARGETs.push_back(3ULL);
    IDX_TARGETs.push_back(4ULL);
    IDX_TARGETs.push_back(5ULL);
    IDX_TARGETs.push_back(6ULL);
    IDX_TARGETs.push_back(7ULL);
    IDX_TARGETs.push_back(8ULL);
    IDX_TARGETs.push_back(9ULL);
    IDX_TARGETs.push_back(10ULL);
    IDX_TARGETs.push_back(11ULL);
    IDX_TARGETs.push_back(12ULL);


    genDataBase(B, num_expansions, further_dims, IDX_TARGETs);



    // print_DB(B, num_expansions, further_dims);

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
    std::vector<size_t> ind_vec;
    for (size_t bucket = 0 ; bucket < godOracle.size() ; ++bucket) {
        if (indexes.find(bucket) != indexes.end()) {
            ind_vec.push_back(static_cast<uint32_t>(indexes[bucket]));
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
    // // Client: 发送请求
    // std::vector<MatPoly> cvs = client_A.MultiClientQuery();
    // Client: 生成 packed_poly 向量
    std::vector<MatPoly> packed_pt_vec = client_A.MultiPackedIndexGen();
    for (size_t i = 0 ; i < packed_pt_vec.size() ; i++) {
        cout << "====================== Let's check the " << i << " packed_indexes =========================" << endl;
        for (size_t j = 0; j < 10 ; j++) {
            cout << j << " " << packed_pt_vec[i].data[j] << endl;
        }
    }
    // Client: 加密得到 cvs 向量
    std::vector<MatPoly> cvs = client_A.MultiPackedIndexEncrypt(packed_pt_vec);
    // Server: 创建 Server 类
    MultiPirServer server_A(godSizes, DBs, cvs);
    // Server: 响应请求
    start_timing();
    std::vector<MatPoly> rs = server_A.MultiServerAnswer();
    cout << "======================= Time used with PBC is " << end_timing() << "============================" << endl;
    std::vector<MatPoly> pts = client_A.MultiAnswerDecrypt(rs);

    for (size_t i = 0 ; i < pts.size() ; i++) {
        cout << "====================== Let's check the " << i << " pts =========================" << endl;
        for (size_t j = 0; j < 10 ; j++) {
            cout << j << " " << pts[i].data[j] << endl;
        }
    }

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

    // test_trival_multi_PIR(B, num_expansions, further_dims, IDX_TARGETs);
}

void print_const_ntt_matpoly(const MatPoly A) {
    assert(A.cols == 1);
    assert(A.rows == 1);
    assert(A.isNTT);
    MatPoly temp = from_ntt(A);
    for (size_t i = 0 ; i < 10 ; i++) {
        cout << i << ' ' << temp.data[i] << endl;
    }
}

void test_gaussian_elimination() {
    MatPoly C1 = genAMatWithV(10245);
    MatPoly C2 = genAMatWithV(1024);
    MatPoly C3 = genAMatWithV(1025);
    auto D1 = to_ntt(C1);
    auto res = MulScalarAndPt(inv_mod(5, Q_i), MulScalarAndPt(5, D1));
    cout << is_eq(from_ntt(res), C1) << endl;

    cout << "Start to the gaussian elimination" << endl;
    cout << getValue(D1) << endl;
    // print_const_ntt_matpoly(D);
    // std::vector<std::vector<size_t>> A = {
    //     {1, 0, 1, 0},
    //     {0, 0, 1, 0},
    //     {0, 1, 0, 1}
    // };
    std::vector<std::vector<size_t>> A = {
        {1, 1, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 1},
        {0, 1, 0, 0},
        {1, 1, 0, 0}
    };
    std::vector<MatPoly> x(4);
    x[0] = C1;
    x[1] = C2;
    x[2] = C3;
    x[3] = C1;
    auto b = MulRandBondMat(A, x);
    for (size_t i = 0 ; i < b.size() ; i++) {
        b[i] = from_ntt(b[i]);
    }
    auto x_prime = GaussianEliminationWithUnderdeterminedMatrix(A, b);
    for (size_t i = 0 ; i < x_prime.size() ; i++) {
        x_prime[i] = from_ntt(x_prime[i]);
        // cout << i << ' ' << getValue(x_prime[i]) << endl;
    }
    auto b_prime = MulRandBondMat(A, x_prime);
    for (size_t i = 0 ; i < b.size() ; i++) {
        b_prime[i] = from_ntt(b_prime[i]);
    }
    for (size_t i = 0 ; i < b_prime.size() ; i++) {
        cout << i << ' ' << is_eq(b[i], b_prime[i]) << endl;
        // cout << getValue(to_ntt(b[i])) << ' ' << getValue(to_ntt(b_prime[i])) << endl;
    }
}

void test_sort_columns() {
    // 示例矩阵
    std::vector<std::vector<size_t>> matrix = {
        {1, 0, 1, 1},
        {0, 0, 1, 0},
        {0, 1, 0, 1}
    };
    std::vector<std::vector<size_t>> res(matrix.size(), std::vector<size_t>(matrix[0].size(), 0));

    // 对列进行排序并获取原始列的位置
    std::vector<size_t> sortedOrder = sortColumnsByFirstNonZero(res, matrix);

    // 输出排序后的列位置
    for (size_t col : sortedOrder) {
        cout << col << " ";
    }
    cout << endl;

    for (size_t i = 0 ; i < res.size() ; i++) {
        for (size_t j = 0 ; j < res[0].size() ; j++) {
            cout << res[i][j] << ' ';
        }
        cout << endl;
    }

}

void test_solve_linear_system() {
    MatPoly C1 = genAMatWithV(10245);
    MatPoly C2 = genAMatWithV(1024);
    MatPoly C3 = genAMatWithV(1025);

    std::vector<std::vector<size_t>> A = {
        {1, 0, 1, 0},
        {0, 0, 1, 0},
        {0, 1, 0, 1},
        {1, 0, 1, 0},
        {1, 0, 1, 0}
    };

    std::vector<MatPoly> x(3);
    x[0] = C1;
    x[1] = C2;
    x[2] = C3;
    x[3] = C1;

    auto b = MulRandBondMat(A, x);
    for (size_t i = 0 ; i < b.size() ; i++) {
        b[i] = from_ntt(b[i]);
    }
    std::vector<MatPoly> x_prime = SolveLinearSystem(A, b, true);

    auto b_prime = MulRandBondMat(A, x_prime);
    for (size_t i = 0 ; i < b.size() ; i++) {
        b_prime[i] = from_ntt(b_prime[i]);
    }

    for (size_t i = 0 ; i < b_prime.size() ; i++) {
        cout << i << ' ' << is_eq(b[i], b_prime[i]) << endl;
    }
}

void test_mPIR_with_PBC_and_OCD(
    size_t num_expansions, 
    size_t further_dims
) {
    do_MatPol_test();

    setup_constants();

    generate_gadgets();

    cout << "Look at this !!!" << endl;

    uint64_t *B;

    cout << random_data << endl;
    cout << "target index is:" << IDX_TARGET << endl;

    std::vector<size_t> IDX_TARGETs;
    IDX_TARGETs.push_back(1ULL);
    IDX_TARGETs.push_back(2ULL);
    IDX_TARGETs.push_back(3ULL);
    IDX_TARGETs.push_back(4ULL);
    IDX_TARGETs.push_back(5ULL);
    IDX_TARGETs.push_back(6ULL);
    IDX_TARGETs.push_back(7ULL);
    IDX_TARGETs.push_back(8ULL);
    IDX_TARGETs.push_back(9ULL);
    IDX_TARGETs.push_back(10ULL);
    IDX_TARGETs.push_back(11ULL);
    IDX_TARGETs.push_back(12ULL);
    IDX_TARGETs.push_back(13ULL);
    IDX_TARGETs.push_back(14ULL);
    IDX_TARGETs.push_back(15ULL);
    IDX_TARGETs.push_back(16ULL);


    genDataBase(B, num_expansions, further_dims, IDX_TARGETs);



    // print_DB(B, num_expansions, further_dims);

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

    // Client: 创建 OCD 类
    size_t w = 8 + ceil(log2(real_ind.size()));
    LSObvDecompress client_OCD(packed_pt_vec.size(), real_ind.size(), 0.05, w);
    auto p_hat = client_OCD.Compress(packed_pt_vec, real_ind);
    cout << "==============================" << endl;

    auto cvs_hat = client_A.MultiPackedIndexEncrypt(p_hat);

    LSObvDecompress server_OCD(packed_pt_vec.size(), real_ind.size(), 0.05, w);
    auto cvs = server_OCD.ObvDecompress(cvs_hat);
    cout << "==============================" << endl;
    MultiPirServer server_A(godSizes, DBs, cvs);
    cout << cvs.size() << endl;
    cout << server_A.handles.size() << endl;
    std::vector<MatPoly> rs = server_A.MultiServerAnswer();
    // cout << "==============================" << endl;
    std::vector<MatPoly> pts = client_A.MultiAnswerDecrypt(rs);

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

    // test_trival_multi_PIR(B, num_expansions, further_dims, IDX_TARGETs);
}

void test_mPIR_with_PBC_and_OCC(
    size_t num_expansions, 
    size_t further_dims
) {
    do_MatPol_test();

    setup_constants();

    generate_gadgets();

    cout << "Look at this !!!" << endl;

    uint64_t *B;

    cout << random_data << endl;
    cout << "target index is:" << IDX_TARGET << endl;

    std::vector<size_t> IDX_TARGETs;
    IDX_TARGETs.push_back(1ULL);
    IDX_TARGETs.push_back(2ULL);
    IDX_TARGETs.push_back(3ULL);
    IDX_TARGETs.push_back(4ULL);
    IDX_TARGETs.push_back(5ULL);
    IDX_TARGETs.push_back(6ULL);
    IDX_TARGETs.push_back(7ULL);
    IDX_TARGETs.push_back(8ULL);
    IDX_TARGETs.push_back(9ULL);
    IDX_TARGETs.push_back(10ULL);
    IDX_TARGETs.push_back(11ULL);
    IDX_TARGETs.push_back(12ULL);
    IDX_TARGETs.push_back(13ULL);
    IDX_TARGETs.push_back(14ULL);
    IDX_TARGETs.push_back(15ULL);
    IDX_TARGETs.push_back(16ULL);


    genDataBase(B, num_expansions, further_dims, IDX_TARGETs);



    // print_DB(B, num_expansions, further_dims);

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
    // // Client: 发送请求
    // std::vector<MatPoly> cvs = client_A.MultiClientQuery();
    // Client: 生成 packed_poly 向量
    std::vector<MatPoly> packed_pt_vec = client_A.MultiPackedIndexGen();
    for (size_t i = 0 ; i < packed_pt_vec.size() ; i++) {
        cout << "====================== Let's check the " << i << " packed_indexes =========================" << endl;
        for (size_t j = 0; j < 10 ; j++) {
            cout << j << " " << packed_pt_vec[i].data[j] << endl;
        }
    }
    // Client: 加密得到 cvs 向量
    std::vector<MatPoly> cvs = client_A.MultiPackedIndexEncrypt(packed_pt_vec);
    // Server: 创建 Server 类
    MultiPirServer server_A(godSizes, DBs, cvs);
    // Server: 响应请求
    // start_timing();
    // std::vector<MatPoly> rs = server_A.MultiServerAnswer();
    std::vector<MatPoly> cts = server_A.MultiServerRespondWithoutModSwitch();
    size_t w = 8 + ceil(log2(real_ind.size()));
    LSObvCompress server_OCC(cts.size(), real_ind.size(), 0.05, w);
    std::vector<MatPoly> cv = server_OCC.ObvCompress(cts);
    std::vector<MatPoly> r_hats = server_A.MultiServerCompressWithModSwitch(cv);
    // cout << "======================= Time used with PBC is " << end_timing() << "============================" << endl;
    std::vector<MatPoly> pt_hats = client_A.MultiAnswerDecrypt(r_hats);
    LSObvCompress client_OCC(cts.size(), real_ind.size(), 0.05, w);
    auto pts = client_OCC.Decompress(pt_hats, real_ind);

    // for (size_t i = 0 ; i < pts.size() ; i++) {
    //     cout << "====================== Let's check the " << i << " pts =========================" << endl;
    //     for (size_t j = 0; j < 10 ; j++) {
    //         cout << j << " " << pts[i].data[j] << endl;
    //     }
    // }
    
    size_t count = 0;
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
            MatPoly corr = pt_reals[seq-1];
            MatPoly M_result = pts[count++];
            cout << "Is correct?: " << is_eq(corr, M_result) << endl;
            show_diff = true;
            if (show_diff) {
                for (size_t i = 0; i < 10 ; i++) {
                cout << i << " " << corr.data[i] << ", " << M_result.data[i] << endl;
                    if (corr.data[i] != M_result.data[i]) {
                        cout << i << " " << corr.data[i] << ", " << M_result.data[i] << endl;
                        // exit(1);
                    }
                }
            }
        }
    }

    // test_trival_multi_PIR(B, num_expansions, further_dims, IDX_TARGETs);
}

void test_pt_add_ntt() {
    auto A = genBMatWithV(110);
    auto A_prime = genBMatWithV(55);
    auto C = genBMatWithV(0);
    auto a = MulScalarAndPt(inv_mod(2, Q_i), ConvertModToQi(A));
    cout << A.data[0] << endl;
    auto D = ConvertModToPdb(add(MulScalarAndPt(inv_mod(2, Q_i), ConvertModToQi(A)), to_ntt(invert(from_ntt(ConvertModToQi(A_prime))))));
    cout << "Is correct ? " << is_eq(C, D) << ' ' << C.data[0] << ' ' << D.data[0] << endl;
}

