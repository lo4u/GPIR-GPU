#include "parallel.h"
#include "opgpu.cuh"
#include <filesystem>

// client
float time_set_constants = 0.0;
float time_get_schedule = 0.0;
float time_build_original_index = 0.0;
float time_gen_key = 0.0;
float time_build_pt_query = 0.0;
float time_compress_query = 0.0;
float time_encrypt_query = 0.0;
float time_decode_respond = 0.0;
float time_batch_query = 0.0;

// server
float time_encode_db = 0.0;
float time_decompress_query = 0.0;
float time_expand_query = 0.0;
float time_reorient_ct = 0.0;
float time_mul_db_with_query = 0.0;
float time_fold_further_dim = 0.0;
float time_convert_mod = 0.0;
float time_respond = 0.0;
float time_batch_respond = 0.0;

double verifiy_corr(FurtherDimsLocals furtherDimsLocals, bool modswitch_on_server) {
    MatPoly ct(n1, n2, false);

    MatPoly ct_inp(n1, n2, false);
    MatPoly total_resp(n1, n2, false);

    MatPoly r_end(S_mp.rows, ct.cols, false);

    MatPoly corr = pt_real_test;
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
    show_diff = false;
    if (show_diff) {
        for (size_t i = 0; i < M_result.rows * M_result.cols * coeff_count; i++) {
            if (corr.data[i] != M_result.data[i]) {
                cout << i << " " << corr.data[i] << ", " << M_result.data[i] << endl;
                // exit(1);
            }
        }
    }
    return 0;
}


void loadDataBase(
    uint64_t*& DB_test, 
    size_t num_expansions, 
    size_t further_dims, 
    size_t IDX_TARGET
) {
    std::string fileName = "database_gened_" + std::to_string(num_expansions) + "_" + std::to_string(further_dims) + ".dat";

    size_t total_n = (1 << num_expansions) * (1 << further_dims);
    size_t dim0 = 1 << num_expansions; // 第一维度的总数
    size_t num_per = 1 << further_dims; // 数据库总数据量 / 第一维度
    
    size_t num_bytes_B = sizeof(uint64_t) * dim0 * num_per * n0 * n2 * poly_len; // 2 * poly_len;
    cout << "num_bytes_B: " << num_bytes_B << endl;
    DB_test = (uint64_t *)aligned_alloc(64, num_bytes_B);
    memset(DB_test, 0, num_bytes_B);
    if (filesystem::exists(fileName)) {
        cout << "Loading database from file..." << endl;
        size_t len;
        loadBinary(DB_test, len, fileName);
        if (len != num_bytes_B / sizeof(uint64_t)) {
            cout << "Database size mismatch!" << endl;
            exit(1);
        }
        cout << "Database loaded." << endl;
        return;
    }

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

        if (i == IDX_TARGET) {
            cop(pt_encd_correct_test, pts_encd_test);
            cop(pt_real_test, pt_tmp);
            to_ntt(pt_test, pt_tmp);
            cop(pt_correct_test, pt_test);
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
    saveBinary(DB_test, num_bytes_B / sizeof(uint64_t), fileName);
    cout << "done loading/generating db." << endl;
}

void multiplyQueryWithDatabase(
    uint64_t *output,
    const uint64_t *reorientedCiphertexts,
    const uint64_t *database,
    size_t dim0, 
    size_t num_per
) {
    size_t n1_padded = 4; // n1 rounded up to power of 2

    uint32_t low_bits_mask_1 = (1UL<< packed_offset_2) - 1;
    uint32_t low_bits_mask_2 = (1UL<< packed_offset_diff) - 1;

    for (size_t z = 0; z < poly_len; z++) {
        size_t idx_a_base = z * (2*dim0*n1_padded);
        size_t idx_b_base = z * (num_per * n2 * dim0 * n0);
        // if (random_data) idx_b_base = (rand() % dummyWorkingSet) * (num_per * n2 * dim0 * n0);

        for (size_t i = 0; i < num_per; i++) {
            for (size_t c = 0; c < n2; c++) {
                // size_t idx_b = idx_b_base + i * (n2 * dim0 * n0) + c * (dim0 * n0);    
                __uint128_t sums_out_n0_0 = 0, sums_out_n0_1 = 0, sums_out_n0_2 = 0;
                __uint128_t sums_out_n1_0 = 0, sums_out_n1_1 = 0, sums_out_n1_2 = 0;

                // #pragma unroll 16
                for (size_t jm = 0; jm < dim0*2; jm++) {
                    uint64_t b = database[idx_b_base++];
                        
                    const uint64_t *v_a = &reorientedCiphertexts[idx_a_base + jm * n1_padded];
                    uint64_t v_a0 = v_a[0];
                    uint64_t v_a1 = v_a[1];
                    uint64_t v_a2 = v_a[2];

                    uint32_t b_lo = b;
                    uint32_t b_hi = b >> 32L;

                    uint32_t v_a0_lo = v_a0;
                    uint32_t v_a0_hi = v_a0 >> 32L;

                    uint32_t v_a1_lo = v_a1;
                    uint32_t v_a1_hi = v_a1 >> 32L;

                    uint32_t v_a2_lo = v_a2;
                    uint32_t v_a2_hi = v_a2 >> 32L;
                        
                    // do n0
                    sums_out_n0_0 += (v_a0_lo) * (uint64_t)b_lo;
                    sums_out_n0_1 += (v_a1_lo) * (uint64_t)b_lo;
                    sums_out_n0_2 += (v_a2_lo) * (uint64_t)b_lo;

                    // do n1
                    sums_out_n1_0 += ((uint64_t)v_a0_hi) * b_hi;
                    sums_out_n1_1 += ((uint64_t)v_a1_hi) * b_hi;
                    sums_out_n1_2 += ((uint64_t)v_a2_hi) * b_hi;
                }

                // output n0
                size_t n = 0;
                size_t idx_c = i * (n1*n2*crt_count*poly_len) + c * (crt_count*poly_len) + n * (poly_len) + z;
                // cout << z << ',' << i << ',' << c << ',' << n << ':' << idx_c << endl;
                output[idx_c] = sums_out_n0_0 % p_i;
                idx_c += (n2*crt_count*poly_len);
                // cout << z << ',' << i << ',' << c << ',' << n << ':' << idx_c << endl;
                output[idx_c] = sums_out_n0_1 % p_i;
                idx_c += (n2*crt_count*poly_len);
                // cout << z << ',' << i << ',' << c << ',' << n << ':' << idx_c << endl;
                output[idx_c] = sums_out_n0_2 % p_i;

                // output n1
                n = 1;
                idx_c = i * (n1*n2*crt_count*poly_len) + c * (crt_count*poly_len) + n * (poly_len) + z;
                // cout << z << ',' << i << ',' << c << ',' << n << ':' << idx_c << endl;
                output[idx_c] = sums_out_n1_0 % b_i;
                idx_c += (n2*crt_count*poly_len);
                // cout << z << ',' << i << ',' << c << ',' << n << ':' << idx_c << endl;
                output[idx_c] = sums_out_n1_1 % b_i;
                idx_c += (n2*crt_count*poly_len);
                // cout << z << ',' << i << ',' << c << ',' << n << ':' << idx_c << endl;
                output[idx_c] = sums_out_n1_2 % b_i;
            }
        }
    }
}

void processQuery(
    uint64_t *B,
    const uint64_t *expansion_query_ct,      // expand this ct
    const uint64_t *setup_data,
    const uint64_t *further_dims_query_ct,          // further dims query
    const uint64_t *further_dims_query_ct_neg,
    ExpansionLocals expansion_locals,    // must be cleared
    FurtherDimsLocals further_dims_locals, // must be cleared
    size_t num_expansions,
    size_t further_dims
) {
    Timer clock;
    size_t total_n = (1 << num_expansions) * (1 << further_dims);
    size_t dim0 = 1 << num_expansions;
    size_t num_per = total_n / dim0;
    
    clock.start_cpu();
    reorientCiphertexts(
        expansion_locals.reoriented_ciphertexts,
        expansion_locals.cts,
        dim0,
        expansion_locals.n1_padded
    );
    clock.stop_cpu();
    clock.duration_cpu<Timer::ms>("time_reorient_ct");
    time_reorient_ct = clock._timeElasped;

    // clock.start_cpu();
    // clock.start_gpu();
    // cout << "pi, bi is " << p_i << ' ' << b_i << endl;
    multiplyQueryWithDatabaseGPU(
        further_dims_locals.scratch_cts1,
        expansion_locals.reoriented_ciphertexts,
        B,
        dim0,
        num_per
    );
    time_mul_db_with_query = time_mul_DB_with_Q;
    // clock.stop_cpu();
    // clock.duration_cpu<Timer::ms>("Multiply Query with Database");
    // clock.stop_gpu();
    nttInvAndCrtLiftCiphertexts(
        num_per,
        further_dims_locals
    );
    
    size_t cur_dim = 0;
    clock.start_cpu();
    // 迭代处理接下来的 num_per 个维度
    while (num_per >= 2) {
        num_per = num_per / 2;
        foldOneFurtherDimension(cur_dim, num_per, further_dims_query_ct, further_dims_query_ct_neg, further_dims_locals);
        cur_dim++;
    }
    clock.stop_cpu();
    clock.duration_cpu<Timer::ms>("time_fold_further_dim");
    time_fold_further_dim = clock._timeElasped;
}

void processQueryCPU(
    uint64_t *B,
    const uint64_t *expansion_query_ct,      // expand this ct
    const uint64_t *setup_data,
    const uint64_t *further_dims_query_ct,          // further dims query
    const uint64_t *further_dims_query_ct_neg,
    ExpansionLocals expansion_locals,    // must be cleared
    FurtherDimsLocals further_dims_locals, // must be cleared
    size_t num_expansions,
    size_t further_dims
) {
    Timer clock;
    size_t total_n = (1 << num_expansions) * (1 << further_dims);
    size_t dim0 = 1 << num_expansions;
    size_t num_per = total_n / dim0;

    clock.start_cpu();
    reorientCiphertexts(
        expansion_locals.reoriented_ciphertexts,
        expansion_locals.cts,
        dim0,
        expansion_locals.n1_padded
    );
    clock.stop_cpu();
    clock.duration_cpu<Timer::ms>("time_reorient_ct");
    time_reorient_ct = clock._timeElasped;

    clock.start_cpu();
    // clock.start_gpu();
    // cout << "pi, bi is " << p_i << ' ' << b_i << endl;
    multiplyQueryWithDatabase(
        further_dims_locals.scratch_cts1,
        expansion_locals.reoriented_ciphertexts,
        B,
        dim0,
        num_per
    );
    clock.stop_cpu();
    clock.duration_cpu<Timer::ms>("time_mul_db_with_query");
    time_mul_db_with_query = clock._timeElasped;
    // clock.stop_gpu();
    nttInvAndCrtLiftCiphertexts(
        num_per,
        further_dims_locals
    );
    
    size_t cur_dim = 0;
    clock.start_cpu();
    // 迭代处理接下来的 num_per 个维度
    while (num_per >= 2) {
        num_per = num_per / 2;
        foldOneFurtherDimension(cur_dim, num_per, further_dims_query_ct, further_dims_query_ct_neg, further_dims_locals);
        cur_dim++;
    }
    clock.stop_cpu();
    clock.duration_cpu<Timer::ms>("time_fold_further_dim");
    time_fold_further_dim = clock._timeElasped;
    // cout << "done folding" << endl;
}

MatPoly AnswerParallel(
    uint64_t *B,
    MatPoly cv,
    ExpansionLocals expansionLocals,
    FurtherDimsLocals furtherDimsLocals,
    const uint64_t *g_C_fft_crtd,      // expand this ct
    uint64_t *g_Q_crtd,          // further dims query
    const uint64_t *g_Ws_fft,           // use this setup data
    size_t num_expansions,
    size_t further_dims
) {
    size_t setupData_size_bytes = num_expansions * (n1 * mh * poly_len * sizeof(uint64_t));
    size_t query_C_crtd_size_bytes = n1 * m1 * poly_len * sizeof(uint64_t); // 第一个维度
    size_t query_later_inp_cts_crtd_size_bytes = further_dims * n1 * m2 * poly_len * sizeof(uint64_t); // 后续维度
    
    // m2 = m_{GSW} = (n + 1) * t_{GSW}
    // the size of the gadget matrix G_{GSW}
    size_t num_bytes_per_Q = n1 * m2 * crt_count * poly_len * sizeof(uint64_t);
    uint64_t *g_Q = (uint64_t *)malloc(further_dims * num_bytes_per_Q);
    uint64_t *g_Q_neg = (uint64_t *)malloc(further_dims * num_bytes_per_Q);

    uint64_t *g_Q_nttd = (uint64_t *)malloc(crt_count * query_later_inp_cts_crtd_size_bytes);
    
    // the first dim's 2^v1 Matrix Regev ciphertexts for basic vector
    // is restored to expansionLocals.ct
    // the rest dim's v2 GSW Regev ciphertexts for j_{j} in {0, 1}
    // the NTT is restored to g_Q_nttd; the no-NTT is restored to g_Q_crtd
    Timer clock;
    clock.start_cpu();
    cout << num_expansions << ' ' << further_dims << endl;
    ServerExpand(
        cv,
        expansionLocals,
        g_Q_nttd,
        g_Q_crtd,
        num_expansions,
        further_dims
    );
    clock.stop_cpu();
    clock.duration_cpu<Timer::ms>("time_expand_query");
    time_expand_query = clock._timeElasped;

    uint64_t *g_Q_neg_crtd = (uint64_t *)malloc(query_later_inp_cts_crtd_size_bytes);
    #pragma omp parallel for
    for (size_t j = 0; j < further_dims; j++) {
        for (size_t r = 0; r < n1; r++) {
            for (size_t m = 0; m < m2; m++) {
                for (size_t z = 0; z < poly_len; z++) {
                    size_t idx_base = j*(n1 * m2 * poly_len);
                    size_t idx = r*(m2 * poly_len) + m*(poly_len) + z;
                    long val = (long)(G2.data[(r * G2.cols + m) * coeff_count + z]) - (long)g_Q_crtd[idx_base + idx];
                    if (val < 0) val += Q_i_u128;
                    g_Q_neg_crtd[idx_base + idx] = val;
                }
            }
        }
    }
    uint64_t *g_Q_neg_nttd = (uint64_t *)malloc(crt_count * query_later_inp_cts_crtd_size_bytes);
    cpu_crt_to_ucompressed_and_ntt(g_Q_neg_nttd, g_Q_neg_crtd, further_dims * n1 * m2);
    free(g_Q_neg_crtd);
    
    #pragma omp parallel for
    for (size_t j = 0; j < further_dims; j++) {
        size_t idx = j*(n1 * m2 * crt_count * poly_len);
        reorient_Q(&g_Q[idx],       &g_Q_nttd[idx]);
        reorient_Q(&g_Q_neg[idx],   &g_Q_neg_nttd[idx]);
    }

    free(g_Q_nttd);
    free(g_Q_neg_nttd);

    uint64_t *g_C_fft = (uint64_t *)malloc(crt_count * query_C_crtd_size_bytes);

    // start_timing();


    // Timer clock;
    // clock.start_cpu();
    processQuery(
        B,
        g_C_fft,      // expand this ct
        g_Ws_fft,
        g_Q,          // further dims query
        g_Q_neg,
        expansionLocals,
        furtherDimsLocals,
        num_expansions,
        further_dims
    );
    // clock.stop_cpu();
    // clock.duration_cpu<Timer::ms>("Process Query");

    // cout << "Real time used is " << end_timing() << endl;

    // cout << "Done with query processing!" << endl;

    MatPoly ct_inp(n1, n2, false);
    MatPoly total_resp(n1, n2, false);
    
    clock.start_cpu();
    if (modswitch_on_server) {
        modswitch(furtherDimsLocals.result, furtherDimsLocals.cts);
        start_timing();
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
    }
    clock.stop_cpu();
    clock.duration_cpu<Timer::ms>("time_convert_mod");
    time_convert_mod = clock._timeElasped;

    return total_resp;
}

MatPoly AnswerCPU(
    uint64_t *B,
    MatPoly cv,
    ExpansionLocals expansionLocals,
    FurtherDimsLocals furtherDimsLocals,
    const uint64_t *g_C_fft_crtd,      // expand this ct
    uint64_t *g_Q_crtd,          // further dims query
    const uint64_t *g_Ws_fft,           // use this setup data
    size_t num_expansions,
    size_t further_dims
) {
    size_t setupData_size_bytes = num_expansions * (n1 * mh * poly_len * sizeof(uint64_t));
    size_t query_C_crtd_size_bytes = n1 * m1 * poly_len * sizeof(uint64_t); // 第一个维度
    size_t query_later_inp_cts_crtd_size_bytes = further_dims * n1 * m2 * poly_len * sizeof(uint64_t); // 后续维度
    
    // m2 = m_{GSW} = (n + 1) * t_{GSW}
    // the size of the gadget matrix G_{GSW}
    size_t num_bytes_per_Q = n1 * m2 * crt_count * poly_len * sizeof(uint64_t);
    uint64_t *g_Q = (uint64_t *)malloc(further_dims * num_bytes_per_Q);
    uint64_t *g_Q_neg = (uint64_t *)malloc(further_dims * num_bytes_per_Q);

    uint64_t *g_Q_nttd = (uint64_t *)malloc(crt_count * query_later_inp_cts_crtd_size_bytes);
    
    // the first dim's 2^v1 Matrix Regev ciphertexts for basic vector
    // is restored to expansionLocals.ct
    // the rest dim's v2 GSW Regev ciphertexts for j_{j} in {0, 1}
    // the NTT is restored to g_Q_nttd; the no-NTT is restored to g_Q_crtd
    Timer clock;
    clock.start_cpu();
    cout << num_expansions << ' ' << further_dims << endl;
    ServerExpand(
        cv,
        expansionLocals,
        g_Q_nttd,
        g_Q_crtd,
        num_expansions,
        further_dims
    );
    clock.stop_cpu();
    clock.duration_cpu<Timer::ms>("time_expand_query");
    time_expand_query = clock._timeElasped;

    // start_timing();
    uint64_t *g_Q_neg_crtd = (uint64_t *)malloc(query_later_inp_cts_crtd_size_bytes);
    #pragma omp parallel for
    for (size_t j = 0; j < further_dims; j++) {
        for (size_t r = 0; r < n1; r++) {
            for (size_t m = 0; m < m2; m++) {
                for (size_t z = 0; z < poly_len; z++) {
                    size_t idx_base = j*(n1 * m2 * poly_len);
                    size_t idx = r*(m2 * poly_len) + m*(poly_len) + z;
                    long val = (long)(G2.data[(r * G2.cols + m) * coeff_count + z]) - (long)g_Q_crtd[idx_base + idx];
                    if (val < 0) val += Q_i_u128;
                    g_Q_neg_crtd[idx_base + idx] = val;
                }
            }
        }
    }
    uint64_t *g_Q_neg_nttd = (uint64_t *)malloc(crt_count * query_later_inp_cts_crtd_size_bytes);
    cpu_crt_to_ucompressed_and_ntt(g_Q_neg_nttd, g_Q_neg_crtd, further_dims * n1 * m2);
    free(g_Q_neg_crtd);
    
    #pragma omp parallel for
    for (size_t j = 0; j < further_dims; j++) {
        size_t idx = j*(n1 * m2 * crt_count * poly_len);
        reorient_Q(&g_Q[idx],       &g_Q_nttd[idx]);
        reorient_Q(&g_Q_neg[idx],   &g_Q_neg_nttd[idx]);
    }

    free(g_Q_nttd);
    free(g_Q_neg_nttd);

    uint64_t *g_C_fft = (uint64_t *)malloc(crt_count * query_C_crtd_size_bytes);

    // start_timing();


    // Timer clock;
    // clock.start_cpu();
    processQueryCPU(
        B,
        g_C_fft,      // expand this ct
        g_Ws_fft,
        g_Q,          // further dims query
        g_Q_neg,
        expansionLocals,
        furtherDimsLocals,
        num_expansions,
        further_dims
    );
    // clock.stop_cpu();
    // clock.duration_cpu<Timer::ms>("Process Query");

    // cout << "Real time used is " << end_timing() << endl;

    // cout << "Done with query processing!" << endl;

    MatPoly ct_inp(n1, n2, false);
    MatPoly total_resp(n1, n2, false);
     
    clock.start_cpu(); 
    if (modswitch_on_server) {
        modswitch(furtherDimsLocals.result, furtherDimsLocals.cts);
        // start_timing();
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
    }
    clock.stop_cpu();
    clock.duration_cpu<Timer::ms>("time_convert_mod");
    time_convert_mod = clock._timeElasped;

    return total_resp;
}

void trival_test_single_pir_with_gpu(
    size_t num_expansions, 
    size_t further_dims,
    size_t IDX_TARGET
) {
    Timer clock;
    do_MatPol_test();
    clock.start_cpu();
    setup_constants();
    generate_gadgets();
    clock.stop_cpu();
    clock.duration_cpu<Timer::ms>("time_set_constants");
    time_set_constants = clock._timeElasped;

    printf("=================== Test of Single PIR with GPU ===================\n");

    printf("=================== 1. Server: START to Build the DATABASE ===================\n");
    uint64_t *B;

    loadDataBase(B, num_expansions, further_dims, IDX_TARGET);   
    printf("=================== 1. Server: Build the DATABASE ... ok ===================\n");

    size_t total_n = (1 << num_expansions) * (1 << further_dims);
    size_t dim0 = 1 << num_expansions;
    size_t num_per = (1 << further_dims);
    if (total_n % dim0 != 0)
        exit(1);

    size_t num_trials = 1;
    std::vector<double> vars;
    for (size_t trial = 0; trial < num_trials; trial++) {
        printf("=================== 2. Client: START to Gen Key ===================\n");
        PirClient client_A(IDX_TARGET, num_expansions, further_dims);
        clock.start_cpu();
        client_A.setup();
        clock.stop_cpu();
        clock.duration_cpu<Timer::ms>("time_gen_key");
        time_gen_key = clock._timeElasped;
        printf("=================== 2. Client: Gen Key ... ok ===================\n");

        printf("=================== 3. Client: START to Gen Query ===================\n");
        clock.start_cpu();        
        auto cv = client_A.ClientQuery();
        clock.stop_cpu();
        clock.duration_cpu<Timer::ms>("time_encrypt_query");
        time_encrypt_query = clock._timeElasped;
        printf("=================== 3. Client: Gen Query ... ok ===================\n");
        
        printf("=================== 4. Server: START to Respond the Query ===================\n");
        PirServer server_A(B, cv, num_expansions, further_dims);
        server_A.genQuery();
        auto r = AnswerParallel(
            server_A.database,
            server_A.cv,
            server_A.expansionLocals,
            server_A.furtherDimsLocals,
            server_A.g_C_fft_crtd,
            server_A.g_Q_crtd,
            server_A.g_Ws_fft,
            server_A.num_expansions,
            server_A.further_dims
        );
        time_respond = time_expand_query + time_reorient_ct + time_mul_db_with_query + time_fold_further_dim + time_convert_mod;
        printf("time_respond is %f\n", time_respond);
        printf("=================== 4. Server: Respond the Query ... ok ===================\n");
        
        printf("=================== 5. Client: START to Decrypt thr Response ===================\n");
        clock.start_cpu();
        client_A.DecodeAnswer(r, true);
        clock.stop_cpu();
        clock.duration_cpu<Timer::ms>("time_decode_respond");
        time_decode_respond = clock._timeElasped;
        printf("=================== 5. Client: Decrypt thr Response ... ok ===================\n");

        printf("size_of_query is %ld\n", calculateMatPolySize(cv));
        printf("size_of_response is %ld\n", calculateMatPolySize(r));
        
        // double log_var;
        // if (modswitch_on_server) {
        //     modswitch(server_A.furtherDimsLocals.result, server_A.furtherDimsLocals.cts);
        //     log_var = verifiy_corr(server_A.furtherDimsLocals, modswitch_on_server);
        // } else {
        //     log_var = verifiy_corr(server_A.furtherDimsLocals, modswitch_on_server);
        // }

        // vars.push_back(log_var);
    }
    printf("=================== Test of Single PIR with GPU ... ok ===================\n");
}

void trival_test_single_pir_with_cpu(
    size_t num_expansions, 
    size_t further_dims,
    size_t IDX_TARGET
) {
    Timer clock;
    do_MatPol_test();
    clock.start_cpu();
    setup_constants();
    generate_gadgets();
    clock.stop_cpu();
    clock.duration_cpu<Timer::ms>("time_set_constants");
    time_set_constants = clock._timeElasped;

    printf("=================== Test of Single PIR with CPU ===================\n");

    printf("=================== 1. Server: START to Build the DATABASE ===================\n");
    uint64_t *B;

    loadDataBase(B, num_expansions, further_dims, IDX_TARGET);   
    printf("=================== 1. Server: Build the DATABASE ... ok ===================\n");

    size_t total_n = (1 << num_expansions) * (1 << further_dims);
    size_t dim0 = 1 << num_expansions;
    size_t num_per = (1 << further_dims);
    if (total_n % dim0 != 0)
        exit(1);

    size_t num_trials = 1;
    std::vector<double> vars;
    for (size_t trial = 0; trial < num_trials; trial++) {
        printf("=================== 2. Client: START to Gen Key ===================\n");
        PirClient client_A(IDX_TARGET, num_expansions, further_dims);
        clock.start_cpu();
        client_A.setup();
        clock.stop_cpu();
        clock.duration_cpu<Timer::ms>("time_gen_key");
        time_gen_key = clock._timeElasped;
        printf("=================== 2. Client: Gen Key ... ok ===================\n");

        printf("=================== 3. Client: START to Gen Query ===================\n");
        clock.start_cpu();        
        auto cv = client_A.ClientQuery();
        clock.stop_cpu();
        clock.duration_cpu<Timer::ms>("time_encrypt_query");
        time_encrypt_query = clock._timeElasped;
        printf("=================== 3. Client: Gen Query ... ok ===================\n");
        
        printf("=================== 4. Server: START to Respond the Query ===================\n");
        PirServer server_A(B, cv, num_expansions, further_dims);
        server_A.genQuery();
        auto r = AnswerCPU(
            server_A.database,
            server_A.cv,
            server_A.expansionLocals,
            server_A.furtherDimsLocals,
            server_A.g_C_fft_crtd,
            server_A.g_Q_crtd,
            server_A.g_Ws_fft,
            server_A.num_expansions,
            server_A.further_dims
        );
        time_respond = time_expand_query + time_reorient_ct + time_mul_db_with_query + time_fold_further_dim + time_convert_mod;
        printf("time_respond is %f\n", time_respond);
        printf("=================== 4. Server: Respond the Query ... ok ===================\n");
        
        printf("=================== 5. Client: START to Decrypt thr Response ===================\n");
        clock.start_cpu();
        client_A.DecodeAnswer(r, true);
        clock.stop_cpu();
        clock.duration_cpu<Timer::ms>("time_decode_respond");
        time_decode_respond = clock._timeElasped;
        printf("=================== 5. Client: Decrypt thr Response ... ok ===================\n");
        
        // double log_var;
        // if (modswitch_on_server) {
        //     modswitch(server_A.furtherDimsLocals.result, server_A.furtherDimsLocals.cts);
        //     log_var = verifiy_corr(server_A.furtherDimsLocals, modswitch_on_server);
        // } else {
        //     log_var = verifiy_corr(server_A.furtherDimsLocals, modswitch_on_server);
        // }

        // vars.push_back(log_var);
    }
    printf("=================== Test of Single PIR with CPU ... ok ===================\n");
}

// ======================================== benchmark ========================================

void test_batch_pir_with_pbc_okvs_gpu(
    size_t num_expansions, 
    size_t further_dims,
    size_t query_num
) {
    // do some preparation
    Timer clock;
    setGPU();
    do_MatPol_test();
    clock.start_cpu();
    setup_constants();
    generate_gadgets();
    clock.stop_cpu();
    clock.duration_cpu<Timer::ms>("time_set_constants");
    time_set_constants = clock._timeElasped;

    printf("=================== Test of Batch PIR with PBC, OCD and GPU ===================\n");
    
    printf("=================== 1. Server: START to Build the DATABASE ===================\n");
    uint64_t *B;
    
    std::vector<size_t> IDX_TARGETs(query_num);
    std::iota(IDX_TARGETs.begin(), IDX_TARGETs.end(), 1);

    genDataBase(B, num_expansions, further_dims, IDX_TARGETs);
    printf("=================== 1. Server: Build the DATABASE ... ok ===================\n");
    
    CuckooCode code(IDX_TARGETs.size(), 3, 1.5);
    
    printf("=================== 2. Client: START to Get the Schedule to Every Bucket ===================\n");
    OracleTy god = getOracle(&code, num_expansions, further_dims);
    auto godOracle = std::get<0>(god);
    auto godSizes = std::get<1>(god);
    
    clock.start_cpu();
    auto schedule = PBC_GetSchedule(&code, IDX_TARGETs);

    std::unordered_map<size_t, size_t> indexes;
    for (const auto& entry : schedule) {
        size_t key = entry.first;
        assert(entry.second.size() == 1);
        size_t bucket = entry.second[0];
        auto& oracle_bucket = godOracle[bucket];
        auto it = std::find_if(oracle_bucket.begin(), oracle_bucket.end(), [key](const size_t& e) {
            return e == key;
        });

        assert(it != oracle_bucket.end());
        indexes[bucket] = std::distance(oracle_bucket.begin(), it);
    }
    clock.stop_cpu();
    clock.duration_cpu<Timer::ms>("time_get_schedule");
    time_get_schedule = clock._timeElasped;
    printf("=================== 2. Client: Get the Schedule to Every Bucket ... ok ===================\n");
    
    printf("=================== 3. Server: START to PBC.Encode the DATABASE ===================\n");
    clock.start_cpu();
    auto collections = PBC_Encode(&code, B, num_expansions, further_dims);
    clock.stop_cpu();
    clock.duration_cpu<Timer::ms>("time_encode_db");
    time_encode_db = clock._timeElasped;
    free(B);

    std::vector<uint64_t *> DBs(collections.size());
    size_t cLen = collections.size();
    // for (size_t i = 0 ; i < cLen ; i++) {
    //     DBs[i] = NULL;
    //     reBuildDB(DBs[i], godSizes[i].num_expansions, godSizes[i].further_dims, collections[i]);
    // }
    printf("=================== 3. Server: PBC.Encode the DATABASE ... ok ===================\n");
    
    printf("=================== 4. Client: START to Build the Original Index Vector ===================\n");
    clock.start_cpu();
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
    clock.stop_cpu();
    clock.duration_cpu<Timer::ms>("time_build_original_index");
    time_build_original_index = clock._timeElasped;
    printf("=================== 4. Client: Build the Original Index Vector ... ok ===================\n");

    printf("=================== 5. Client: START to Gen Key ===================\n");
    MultiPirClient client_A(godSizes, ind_vec, godOracle, indexes);
    clock.start_cpu();
    client_A.setup();
    clock.stop_cpu();
    clock.duration_cpu<Timer::ms>("time_gen_key");
    time_gen_key = clock._timeElasped;
    printf("=================== 5. Client: Gen Key ... ok ===================\n");
    
    printf("=================== 6. Client: START to Build the Plaintext Query ===================\n");
    clock.start_cpu();
    std::vector<MatPoly> packed_pt_vec = client_A.MultiPackedIndexGen();
    clock.stop_cpu();
    clock.duration_cpu<Timer::ms>("time_build_pt_query");
    time_build_pt_query = clock._timeElasped;
    printf("=================== 6. Client: Build the Plaintext Query ... ok ===================\n");

    printf("=================== 7. Client: START to Compress the Plaintext Query ===================\n");
    GCTObvKVStore client_OKVS(real_ind.size(), 3, 1.3, packed_pt_vec.size());
    clock.start_cpu();
    auto pt_hat = client_OKVS.Encode(packed_pt_vec, real_ind);
    clock.stop_cpu();
    clock.duration_cpu<Timer::ms>("time_compress_query");
    time_compress_query = clock._timeElasped;
    printf("=================== 7. Client: Compress the Plaintext Query ... ok ===================\n");

    printf("=================== 8. Client: START to Encrypt the Plaintext Query ===================\n");
    clock.start_cpu();
    auto cvs_hat = client_A.MultiPackedIndexEncrypt(pt_hat);
    clock.stop_cpu();
    clock.duration_cpu<Timer::ms>("time_encrypt_query");
    time_encrypt_query = clock._timeElasped;
    time_batch_query = time_get_schedule + time_build_original_index + time_gen_key + time_build_pt_query + time_compress_query + time_encrypt_query;
    printf("time_batch_query is %f\n", time_batch_query);
    printf("=================== 8. Client: Encrypt the Plaintext Query ... ok ===================\n");
    
    printf("=================== 9. Server: START to Oblivious Decompress the Query ===================\n");
    GCTObvKVStore server_OKVS(real_ind.size(), 3, 1.3, packed_pt_vec.size());
    clock.start_cpu();
    auto cvs = server_OKVS.Decode(cvs_hat);
    clock.stop_cpu();
    clock.duration_cpu<Timer::ms>("time_decompress_query");
    time_decompress_query = clock._timeElasped;
    printf("=================== 9. Server: Oblivious Decompress the Query ... ok ===================\n");

    printf("=================== 10. Server: START to Respond the Query ===================\n");
    MultiPirServer server_A(godSizes, DBs, cvs);
    // std::vector<MatPoly> rs = server_A.MultiServerAnswer();
    std::vector<MatPoly> rs(server_A.handles.size());
    size_t hLen = server_A.handles.size();
    for (size_t i = 0 ; i < hLen ; i++) {
        reBuildDB(
            server_A.handles[i].database,
            godSizes[i].num_expansions,
            godSizes[i].further_dims,
            collections[i]
        );
        server_A.handles[i].genQuery();
        rs[i] = AnswerParallel(
            server_A.handles[i].database,
            server_A.handles[i].cv,
            server_A.handles[i].expansionLocals,
            server_A.handles[i].furtherDimsLocals,
            server_A.handles[i].g_C_fft_crtd,
            server_A.handles[i].g_Q_crtd,
            server_A.handles[i].g_Ws_fft,
            server_A.handles[i].num_expansions,
            server_A.handles[i].further_dims
        );
        time_respond = time_expand_query + time_reorient_ct + time_mul_db_with_query + time_fold_further_dim + time_convert_mod;
        time_batch_respond += time_respond;
        free(server_A.handles[i].database);
    }
    float time_ave_batch_respond = (time_batch_respond + time_decompress_query) / hLen;
    printf("time_ave_batch_respond is %f\n", time_ave_batch_respond);
    printf("=================== 10. Server: Respond the Query ... ok ===================\n");
    
    printf("=================== 11. Client: START to Decrypt thr Response ===================\n");
    clock.start_cpu();
    std::vector<MatPoly> pts = client_A.MultiAnswerDecrypt(rs);
    clock.stop_cpu();
    clock.duration_cpu<Timer::ms>("time_decode_respond");
    time_decode_respond = clock._timeElasped;
    printf("=================== 11. Client: Decrypt thr Response ... ok ===================\n");
    printf("size_of_query is %ld\n", calculateMatPolySize(cvs_hat[0]) * cvs_hat.size());
    printf("size_of_response is %ld\n", calculateMatPolySize(rs[0]) * rs.size());
    
    finishGPU();
    printf("=================== Test of Batch PIR with PBC, OCD and GPU ... ok ===================\n");   
}

void test_batch_pir_with_pbc_okvs_cpu(
    size_t num_expansions, 
    size_t further_dims,
    size_t query_num
) {
    // do some preparation
    Timer clock;
    setGPU();
    do_MatPol_test();
    clock.start_cpu();
    setup_constants();
    generate_gadgets();
    clock.stop_cpu();
    clock.duration_cpu<Timer::ms>("time_set_constants");
    time_set_constants = clock._timeElasped;

    printf("=================== Test of Batch PIR with PBC, OCD and GPU ===================\n");
    
    printf("=================== 1. Server: START to Build the DATABASE ===================\n");
    uint64_t *B;
    
    std::vector<size_t> IDX_TARGETs(query_num);
    std::iota(IDX_TARGETs.begin(), IDX_TARGETs.end(), 1);

    genDataBase(B, num_expansions, further_dims, IDX_TARGETs);
    printf("=================== 1. Server: Build the DATABASE ... ok ===================\n");
    
    CuckooCode code(IDX_TARGETs.size(), 3, 1.5);
    
    printf("=================== 2. Client: START to Get the Schedule to Every Bucket ===================\n");
    OracleTy god = getOracle(&code, num_expansions, further_dims);
    auto godOracle = std::get<0>(god);
    auto godSizes = std::get<1>(god);
    
    clock.start_cpu();
    auto schedule = PBC_GetSchedule(&code, IDX_TARGETs);

    std::unordered_map<size_t, size_t> indexes;
    for (const auto& entry : schedule) {
        size_t key = entry.first;
        assert(entry.second.size() == 1);
        size_t bucket = entry.second[0];
        auto& oracle_bucket = godOracle[bucket];
        auto it = std::find_if(oracle_bucket.begin(), oracle_bucket.end(), [key](const size_t& e) {
            return e == key;
        });

        assert(it != oracle_bucket.end());
        indexes[bucket] = std::distance(oracle_bucket.begin(), it);
    }
    clock.stop_cpu();
    clock.duration_cpu<Timer::ms>("time_get_schedule");
    time_get_schedule = clock._timeElasped;
    printf("=================== 2. Client: Get the Schedule to Every Bucket ... ok ===================\n");
    
    printf("=================== 3. Server: START to PBC.Encode the DATABASE ===================\n");
    clock.start_cpu();
    auto collections = PBC_Encode(&code, B, num_expansions, further_dims);
    clock.stop_cpu();
    clock.duration_cpu<Timer::ms>("time_encode_db");
    time_encode_db = clock._timeElasped;

    std::vector<uint64_t *> DBs(collections.size());
    size_t cLen = collections.size();
    for (size_t i = 0 ; i < cLen ; i++) {
        DBs[i] = NULL;
        reBuildDB(DBs[i], godSizes[i].num_expansions, godSizes[i].further_dims, collections[i]);
    }
    printf("=================== 3. Server: PBC.Encode the DATABASE ... ok ===================\n");
    
    printf("=================== 4. Client: START to Build the Original Index Vector ===================\n");
    clock.start_cpu();
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
    clock.stop_cpu();
    clock.duration_cpu<Timer::ms>("time_build_original_index");
    time_build_original_index = clock._timeElasped;
    printf("=================== 4. Client: Build the Original Index Vector ... ok ===================\n");

    printf("=================== 5. Client: START to Gen Key ===================\n");
    MultiPirClient client_A(godSizes, ind_vec, godOracle, indexes);
    clock.start_cpu();
    client_A.setup();
    clock.stop_cpu();
    clock.duration_cpu<Timer::ms>("time_gen_key");
    time_gen_key = clock._timeElasped;
    printf("=================== 5. Client: Gen Key ... ok ===================\n");
    
    printf("=================== 6. Client: START to Build the Plaintext Query ===================\n");
    clock.start_cpu();
    std::vector<MatPoly> packed_pt_vec = client_A.MultiPackedIndexGen();
    clock.stop_cpu();
    clock.duration_cpu<Timer::ms>("time_build_pt_query");
    time_build_pt_query = clock._timeElasped;
    printf("=================== 6. Client: Build the Plaintext Query ... ok ===================\n");

    printf("=================== 7. Client: START to Compress the Plaintext Query ===================\n");
    GCTObvKVStore client_OKVS(real_ind.size(), 3, 1.3, packed_pt_vec.size());
    clock.start_cpu();
    auto pt_hat = client_OKVS.Encode(packed_pt_vec, real_ind);
    clock.stop_cpu();
    clock.duration_cpu<Timer::ms>("time_compress_query");
    time_compress_query = clock._timeElasped;
    printf("=================== 7. Client: Compress the Plaintext Query ... ok ===================\n");

    printf("=================== 8. Client: START to Encrypt the Plaintext Query ===================\n");
    clock.start_cpu();
    auto cvs_hat = client_A.MultiPackedIndexEncrypt(pt_hat);
    clock.stop_cpu();
    clock.duration_cpu<Timer::ms>("time_encrypt_query");
    time_encrypt_query = clock._timeElasped;
    time_batch_query = time_get_schedule + time_build_original_index + time_gen_key + time_build_pt_query + time_compress_query + time_encrypt_query;
    printf("time_batch_query is %f\n", time_batch_query);
    printf("=================== 8. Client: Encrypt the Plaintext Query ... ok ===================\n");
    
    printf("=================== 9. Server: START to Oblivious Decompress the Query ===================\n");
    GCTObvKVStore server_OKVS(real_ind.size(), 3, 1.3, packed_pt_vec.size());
    clock.start_cpu();
    auto cvs = server_OKVS.Decode(cvs_hat);
    clock.stop_cpu();
    clock.duration_cpu<Timer::ms>("time_decompress_query");
    time_decompress_query = clock._timeElasped;
    printf("=================== 9. Server: Oblivious Decompress the Query ... ok ===================\n");

    printf("=================== 10. Server: START to Respond the Query ===================\n");
    MultiPirServer server_A(godSizes, DBs, cvs);
    // std::vector<MatPoly> rs = server_A.MultiServerAnswer();
    std::vector<MatPoly> rs(server_A.handles.size());
    size_t hLen = server_A.handles.size();
    start_timing();
    for (size_t i = 0 ; i < hLen ; i++) {
        server_A.handles[i].genQuery();
        rs[i] = AnswerCPU(
            server_A.handles[i].database,
            server_A.handles[i].cv,
            server_A.handles[i].expansionLocals,
            server_A.handles[i].furtherDimsLocals,
            server_A.handles[i].g_C_fft_crtd,
            server_A.handles[i].g_Q_crtd,
            server_A.handles[i].g_Ws_fft,
            server_A.handles[i].num_expansions,
            server_A.handles[i].further_dims
        );
        time_respond = time_decompress_query + time_expand_query + time_reorient_ct + time_mul_db_with_query + time_fold_further_dim + time_convert_mod;
        time_batch_respond += time_respond;
    }
    cout << "======================= Time used with PBC is " << end_timing() << "============================" << endl;
    float time_ave_batch_respond = time_batch_respond / hLen;
    printf("time_ave_batch_respond is %f\n", time_ave_batch_respond);
    printf("=================== 10. Server: Respond the Query ... ok ===================\n");
    
    printf("=================== 11. Client: START to Decrypt thr Response ===================\n");
    clock.start_cpu();
    std::vector<MatPoly> pts = client_A.MultiAnswerDecrypt(rs);
    clock.stop_cpu();
    clock.duration_cpu<Timer::ms>("time_decode_respond");
    time_decode_respond = clock._timeElasped;
    printf("=================== 11. Client: Decrypt thr Response ... ok ===================\n");
    
    finishGPU();
    printf("=================== Test of Batch PIR with PBC, OCD and GPU ... ok ===================\n");   
}

void test_speed_of_3H_GCT(
    size_t num_expansions, 
    size_t further_dims,
    size_t IDX_TARGET
) {
    do_MatPol_test();
    setup_constants();
    generate_gadgets();

    uint64_t *B;
    std::vector<size_t> IDX_TARGETs(IDX_TARGET);
    std::iota(IDX_TARGETs.begin(), IDX_TARGETs.end(), 1);

    CuckooCode code(IDX_TARGETs.size(), 3, 1.5);    
    OracleTy god = getOracle(&code, num_expansions, further_dims);
    auto godOracle = std::get<0>(god);
    auto godSizes = std::get<1>(god);

    start_timing();
    auto schedule = PBC_GetSchedule(&code, IDX_TARGETs);
    std::unordered_map<size_t, size_t> indexes;
    for (const auto& entry : schedule) {
        size_t key = entry.first;
        assert(entry.second.size() == 1);
        size_t bucket = entry.second[0];
        auto& oracle_bucket = godOracle[bucket];
        auto it = std::find_if(oracle_bucket.begin(), oracle_bucket.end(), [key](const size_t& e) {
            return e == key;
        });
        assert(it != oracle_bucket.end());
        indexes[bucket] = std::distance(oracle_bucket.begin(), it);
    }
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
    MultiPirClient client_A(godSizes, ind_vec, godOracle, indexes);
    client_A.setup();
    std::vector<MatPoly> packed_pt_vec = client_A.MultiPackedIndexGen();
    cout << "================== Benchmark ==================" << endl;
    cout << "Time used to do pack is " << end_timing() << endl;
    GCTObvKVStore client_OKVS(real_ind.size(), 3, 1.3, packed_pt_vec.size());

    // test: the speed of OKVS
    start_timing();
    auto pt_hat = client_OKVS.Encode(packed_pt_vec, real_ind);
    double time_3H_GCT_Compress = end_timing();
    cout << "Time used to do OCD.Compress with 3H-GCT is " << time_3H_GCT_Compress << endl;

    start_timing();
    auto cvs_hat = client_A.MultiPackedIndexEncrypt(pt_hat); // In NTT
    double time_encrypt = end_timing();
    cout << "ALL Time used to do OCD.Compress with 3H-GCT is " << time_3H_GCT_Compress + time_encrypt << endl;
    
    GCTObvKVStore server_OKVS(real_ind.size(), 3, 1.3, packed_pt_vec.size());
    start_timing();
    auto cvs = server_OKVS.Decode(cvs_hat);
    double time_3H_GCT_Decompress = end_timing();
    cout << "Time used to do OCD.Decompress with 3H-GCT is " << time_3H_GCT_Decompress << endl;
    cout << "Total time used to do OCD with 3H-GCT is " << time_3H_GCT_Decompress + time_3H_GCT_Compress << endl;
    cout << "Plus time used to do OCD with 3H-GCT is " << time_3H_GCT_Decompress + time_3H_GCT_Compress + time_encrypt<< endl;
}

void test_speed_of_without_okvs(
    size_t num_expansions, 
    size_t further_dims,
    size_t IDX_TARGET
) {
    do_MatPol_test();
    setup_constants();
    generate_gadgets();

    uint64_t *B;
    std::vector<size_t> IDX_TARGETs(IDX_TARGET);
    std::iota(IDX_TARGETs.begin(), IDX_TARGETs.end(), 1);

    CuckooCode code(IDX_TARGETs.size(), 3, 1.5);    
    OracleTy god = getOracle(&code, num_expansions, further_dims);
    auto godOracle = std::get<0>(god);
    auto godSizes = std::get<1>(god);

    start_timing();
    auto schedule = PBC_GetSchedule(&code, IDX_TARGETs);
    std::unordered_map<size_t, size_t> indexes;
    for (const auto& entry : schedule) {
        size_t key = entry.first;
        assert(entry.second.size() == 1);
        size_t bucket = entry.second[0];
        auto& oracle_bucket = godOracle[bucket];
        auto it = std::find_if(oracle_bucket.begin(), oracle_bucket.end(), [key](const size_t& e) {
            return e == key;
        });
        assert(it != oracle_bucket.end());
        indexes[bucket] = std::distance(oracle_bucket.begin(), it);
    }
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
    MultiPirClient client_A(godSizes, ind_vec, godOracle, indexes);
    client_A.setup();
    std::vector<MatPoly> packed_pt_vec = client_A.MultiPackedIndexGen();
    cout << "================== Benchmark ==================" << endl;
    cout << "Time used to do pack is " << end_timing() << endl;

    start_timing();
    auto cvs_hat = client_A.MultiPackedIndexEncrypt(packed_pt_vec); // In NTT
    double time_encrypt = end_timing();
    cout << "ALL Time used to do OCD.Compress with 3H-GCT is " << time_encrypt << endl;
}