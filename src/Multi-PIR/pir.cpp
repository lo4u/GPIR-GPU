#include "pir.h"

uint64_t maskDB = (1ULL << 32) - 1; // 使用1ULL来确保是64位的常量
MatPoly pts_encd_test(n0, n2);
MatPoly pt_encd_correct_test(n0, n2);
MatPoly pt_test(n0, n0);
MatPoly pt_correct_test(n0, n0);
MatPoly pt_real_test(n0, n2, false);
std::vector<MatPoly> pt_reals;

double check_relation(
    MatPoly total_resp, 
    bool modswitch_on_server,
    size_t seq
) {
    cout << "check " << seq << "'s answer's correctness" << endl;
    MatPoly ct(n1, n2, false);

    MatPoly ct_inp(n1, n2, false);
    // MatPoly total_resp(n1, n2, false);

    MatPoly r_end(S_mp.rows, ct.cols, false);

    MatPoly corr = pt_reals[seq];
    uint64_t q_1 = 4*p_db;
    // MatPoly corr = pt_real_test;
    MatPoly M_result(n0, n0, false);
    MatPoly Sp_mp_nttd_qprime(n0, k_param, false);
    Sp_mp_nttd_qprime = Sp_mp;
    to_ntt_qprime(Sp_mp_nttd_qprime);

    MatPoly first_row_decoded = pick(total_resp, 0, 0, 1, total_resp.cols);
    MatPoly rest_rows_decoded = pick(total_resp, 1, 0, total_resp.rows - 1, total_resp.cols);
    to_ntt_qprime(first_row_decoded); // 方便快速计算乘法
    MatPoly s_prod = mul_over_qprime(Sp_mp_nttd_qprime, first_row_decoded); // 乘法运算，结果矩阵与 rest_rows_decoded 维数相同
    from_ntt_qprime(s_prod);
    for (size_t i = 0; i < s_prod.rows * s_prod.cols * poly_len; i++) { // 循环遍历矩阵中多项式的每个系数
        int64_t val_first = s_prod.data[i];
        if (val_first >= arb_qprime/2) val_first -= arb_qprime; // 对于第一行 Z_q2 中的元素，将其转换到绝对值最小剩余系
        int64_t val_rest = rest_rows_decoded.data[i];
        if (val_rest >= q_1/2) val_rest -= q_1; // 对于剩余行 Z_q1 中的元素，将其转换到绝对值最小剩余系

        uint64_t denom = arb_qprime * (q_1/p_db);

        int64_t r = val_first * q_1;
        r +=  val_rest * arb_qprime;
        // divide r by arb_qprime, rounding
        int64_t sign = r >= 0 ? 1 : -1;
        __int128_t val = r;
        __int128_t result = (r + sign*((int64_t)denom/2)) / (__int128_t)denom;
        result = (result + (denom/p_db)*p_db + 2*p_db) % p_db;

        s_prod.data[i] = (uint64_t)result;
    }
    M_result = s_prod;

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
    return 0;
}

double check_corr(
    FurtherDimsLocals furtherDimsLocals, 
    bool modswitch_on_server, 
    size_t seq
) {
    cout << "check " << seq << "'s answer's correctness" << endl;
    MatPoly ct(n1, n2, false);

    MatPoly ct_inp(n1, n2, false);
    MatPoly total_resp(n1, n2, false);

    MatPoly r_end(S_mp.rows, ct.cols, false);

    MatPoly corr = pt_reals[seq];
    // MatPoly corr = pt_real_test;
    MatPoly M_result(n0, n0, false);

    if (!modswitch_on_server) {
        for (size_t i = 0; i < n1 * n2 * poly_len; i++) {
            ct.data[i] = furtherDimsLocals.cts[i];
        }
        dec_compressed(r_end, to_ntt(ct), to_ntt(S_mp), scale_k);
    } else {

        MatPoly Sp_mp_nttd_qprime(n0, k_param, false);
        Sp_mp_nttd_qprime = Sp_mp;
        to_ntt_qprime(Sp_mp_nttd_qprime);

        for (size_t i = 0; i < n1 * n2 * poly_len; i++) {
            ct_inp.data[i] = furtherDimsLocals.cts[i];
        }
        uint64_t q_1 = 4*p_db;

        // 矩阵的第一行转化为以 q2 为模数，剩余行转化为以 q1 为模数，将其分隔开进行处理
        MatPoly first_row = pick(ct_inp, 0, 0, 1, ct_inp.cols);
        MatPoly first_row_sw = getRescaled(first_row, Q_i, arb_qprime); // Q_i = q -> arb_qprime = q2
        MatPoly rest_rows = pick(ct_inp, 1, 0, ct_inp.rows - 1, ct_inp.cols);
        MatPoly rest_rows_sw = getRescaled(rest_rows, Q_i, q_1); // Q_i = q -> q_1

        // 将模数转换后的两部分重新拼接到一起
        place(total_resp, first_row_sw, 0, 0);
        place(total_resp, rest_rows_sw, 1, 0);

        // ~~transmit~~

        // 重新提取出第一行和后面的部分
        MatPoly first_row_decoded = pick(total_resp, 0, 0, 1, total_resp.cols);
        MatPoly rest_rows_decoded = pick(total_resp, 1, 0, total_resp.rows - 1, total_resp.cols);
        to_ntt_qprime(first_row_decoded); // 方便快速计算乘法
        MatPoly s_prod = mul_over_qprime(Sp_mp_nttd_qprime, first_row_decoded); // 乘法运算，结果矩阵与 rest_rows_decoded 维数相同
        from_ntt_qprime(s_prod);
        for (size_t i = 0; i < s_prod.rows * s_prod.cols * poly_len; i++) { // 循环遍历矩阵中多项式的每个系数
            int64_t val_first = s_prod.data[i];
            if (val_first >= arb_qprime/2) val_first -= arb_qprime; // 对于第一行 Z_q2 中的元素，将其转换到绝对值最小剩余系
            int64_t val_rest = rest_rows_decoded.data[i];
            if (val_rest >= q_1/2) val_rest -= q_1; // 对于剩余行 Z_q1 中的元素，将其转换到绝对值最小剩余系

            uint64_t denom = arb_qprime * (q_1/p_db);

            int64_t r = val_first * q_1;
            r +=  val_rest * arb_qprime;
            // divide r by arb_qprime, rounding
            int64_t sign = r >= 0 ? 1 : -1;
            __int128_t val = r;
            __int128_t result = (r + sign*((int64_t)denom/2)) / (__int128_t)denom;
            result = (result + (denom/p_db)*p_db + 2*p_db) % p_db;

            s_prod.data[i] = (uint64_t)result;
        }
        M_result = s_prod;
    }

    cout << "Is correct?: " << is_eq(corr, M_result) << endl;
    show_diff = true;
    if (show_diff) {
        for (size_t i = 0; i < 10 ; i++) {
            if (corr.data[i] != M_result.data[i]) {
                cout << i << " " << corr.data[i] << ", " << M_result.data[i] << endl;
                // exit(1);
            }
        }
    }
    return 0;
}

void test_oo(
    uint64_t *B, 
    size_t num_expansions, 
    size_t further_dims, 
    size_t IDX_TARGET
) {
    size_t total_n = (1 << num_expansions) * (1 << further_dims);
    size_t dim0 = 1 << num_expansions;
    size_t num_per = total_n / dim0;
    cout << "dim0: " << dim0 << endl;
    cout << "num_per: " << num_per << endl;
    if (total_n % dim0 != 0)
        exit(1);
    bool checking_for_debug = true;
    std::vector<double> vars;
    
    PirClient client_A(IDX_TARGET, num_expansions, further_dims);
    client_A.setup();
    MatPoly cv = client_A.ClientQuery();

    PirServer server_A(B, cv, num_expansions, further_dims);
    server_A.genQuery();
    start_timing();
    MatPoly r = server_A.ServerAnswer();
    cout << "The time used to answer is " << end_timing() << endl;

    if (checking_for_debug) {
        double log_var;
        if (modswitch_on_server) {
            modswitch(server_A.furtherDimsLocals.result, server_A.furtherDimsLocals.cts);
            log_var = check_corr(server_A.furtherDimsLocals, modswitch_on_server, 0ULL);
            // cop(pt_real, pt_real_test);
            // check_final(server_A.furtherDimsLocals, modswitch_on_server);
        } else {
            log_var = check_corr(server_A.furtherDimsLocals, modswitch_on_server, 0ULL);
        }

        vars.push_back(log_var);
    }
}

void genDataBase(
    uint64_t*& DB_test, 
    size_t num_expansions, 
    size_t further_dims, 
    std::vector<size_t> IDX_TARGETs
) {
    size_t total_n = (1 << num_expansions) * (1 << further_dims);
    size_t dim0 = 1 << num_expansions; // 第一维度的总数
    size_t num_per = 1 << further_dims; // 数据库总数据量 / 第一维度
    
    size_t num_bytes_B = sizeof(uint64_t) * dim0 * num_per * n0 * n2 * poly_len; // 2 * poly_len;
    cout << "num_bytes_B: " << num_bytes_B << endl;
    DB_test = (uint64_t *)aligned_alloc(64, num_bytes_B);
    memset(DB_test, 0, num_bytes_B);

    uint64_t *BB = (uint64_t *)malloc(n0 * n2 * crt_count * poly_len * sizeof(uint64_t));
    size_t numBytesPlaintextRaw = n0 * n0 * num_bits_q * poly_len / 8;
    uint64_t *pt_raw = (uint64_t *)malloc(numBytesPlaintextRaw);
    memset(pt_raw, 0, numBytesPlaintextRaw);
    uint64_t *pt_buf = (uint64_t *)malloc(n0 * n0 * crt_count * poly_len * sizeof(uint64_t));
    memset(pt_buf, 0, n0 * n0 * crt_count * poly_len * sizeof(uint64_t));
    
    cout << "starting generation of db" << endl;
    MatPoly H_nttd = to_ntt(H_mp);
    uint64_t *H_encd = H_nttd.data;
    MatPoly pt_tmp(n0, n0, false);
    MatPoly pt_encd_raw(n0, n2, false);
    pts_encd_test = MatPoly(n0, n2);
    pt_test = MatPoly(n0, n0);
    for (size_t i = 0 ; i < total_n ; i++) {
        generate_random_pt(pt_tmp);

        pt_encd_raw = pt_tmp;
        for (size_t pol = 0; pol < n0 * n2 * poly_len; pol++) {
            int64_t val = (int64_t) pt_encd_raw.data[pol];
            assert(val >= 0 && val < p_db);
            if (val >= (p_db / 2)) {
                val = val - (int64_t)p_db;
            }
            if (val < 0) {
                val += Q_i;
            }
            assert(val >= 0 && val < Q_i);
            pt_encd_raw.data[pol] = val;
        }
        to_ntt(pts_encd_test, pt_encd_raw);
        // print_pt_crt(pts_encd_test);

        if (std::find(IDX_TARGETs.begin(), IDX_TARGETs.end(), i) != IDX_TARGETs.end()) {
            // cout << "===============FIND ONE ANSWER" << i << "==================" << endl;
            MatPoly pt_real_test(n0, n2, false);
            cop(pt_encd_correct_test, pts_encd_test);
            cop(pt_real_test, pt_tmp);
            to_ntt(pt_test, pt_tmp);
            cop(pt_correct_test, pt_test);
            pt_reals.push_back(pt_real_test);
            // for (size_t xxx = 0; xxx < 10 ; xxx++) {
            //     cout << xxx << " " << pt_real_test.data[xxx] << endl;
            // }
        }

        // b': i c n z j m
        size_t ii = i % num_per;
        size_t j = i / num_per;
        for (size_t m = 0; m < n0; m++) {
            for (size_t c = 0; c < n2; c++) {
                memcpy(BB, &pts_encd_test.data[(m * n2 + c) * crt_count * coeff_count], crt_count * coeff_count * sizeof(uint64_t));
                for (size_t z = 0; z < poly_len; z++) {
                    size_t idx = z * (num_per * n2 * dim0 * n0) +
                                 ii * (n2 * dim0 * n0) + 
                                 c * (dim0 * n0) +
                                 j * (n0) + m;
                    DB_test[idx] = BB[z] | (BB[poly_len + z] << 32);
                }
            }
        }
    }
    free(BB);
}

void test_trival_multi_PIR(
    uint64_t *B, size_t num_expansions, 
    size_t further_dims, 
    std::vector<size_t> IDX_TARGETs
) {
    size_t test_times = IDX_TARGETs.size();

    size_t total_n = (1 << num_expansions) * (1 << further_dims);
    size_t dim0 = 1 << num_expansions;
    size_t num_per = total_n / dim0;
    cout << "dim0: " << dim0 << endl;
    cout << "num_per: " << num_per << endl;
    if (total_n % dim0 != 0)
        exit(1);
    bool checking_for_debug = true;
    std::vector<double> vars;

    PirClient client_leader(0, num_expansions, further_dims);
    client_leader.setup();
    double time = 0;

    for (size_t i = 0 ; i < test_times ; i++) {
        PirClient client_A(IDX_TARGETs[i], num_expansions, further_dims);
        MatPoly cv = client_A.ClientQuery();

        PirServer server_A(B, cv, num_expansions, further_dims);
        server_A.genQuery();
        start_timing();
        MatPoly r = server_A.ServerAnswer();
        time += end_timing();

        modswitch(server_A.furtherDimsLocals.result, server_A.furtherDimsLocals.cts);
        double log_var = check_corr(server_A.furtherDimsLocals, modswitch_on_server, i);

    }

    cout << "======================= Time used without PBC is " << time << "============================" << endl;
}






