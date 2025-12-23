#pragma once

#include <chrono>
#include <fstream>
#include <inttypes.h>
#include <iostream>
#include <random>
#include <sstream>
#include <stdint.h>
#include <stdio.h>
#include <string>
#include <queue>
#include <mutex>
#include <assert.h>

#ifndef __EMSCRIPTEN__
    #include <omp.h>
    // #include <valgrind/callgrind.h>
    #include <unistd.h>
    #include <sys/socket.h>
    #include <sys/un.h>
#endif

using namespace std;

#include "hexl/ntt/ntt.hpp"

#include <bitset>
#include <map>
#include <math.h>

#include "core.h"
#include "poly.h"
#include "values.h"
#include "util.h"
#include "client.h"
#include "testing.h"

#include "common.h"
#include <constants.h>

#include <functional>
#include <iomanip>

#ifdef USE_NFLLIB
#include <nfl.hpp>
#include <nfl/poly.hpp>
#include <nfl/poly_p.hpp>
#include <nfl/core.hpp>
using namespace nfl;
#endif

// 声明 spiral.cpp 中的一些全局变量

extern bool has_file;
extern bool has_data;
extern string dbFilename;
extern string dataFilename;
extern bool load;
extern bool server;
extern bool client;
extern string outputErrFilename;
extern fstream fs;
extern fstream fsData;

extern size_t num_expansions; // max = 7 //used to be 8
extern size_t further_dims; // v2
extern size_t total_n;

extern MatPoly G_hat;

constexpr size_t mh = m_exp;
static_assert(m2 % logQ != 0, "m2 dimension cannot evenly divide logQ");

constexpr size_t m1 = n0;

// constexpr uint64_t switch_factor = 1UL << switch_factor_log2;
// constexpr uint64_t bits_to_hold_q_times_switch_factor = bits_to_hold_q + switch_factor_log2;
// uint64_t factor_to_remove_obliviously = ((uint64_t)round((long double)qq_const / (long double)switch_factor));
// uint64_t factor_init = factor_to_remove_obliviously * switch_factor;

// uint64_t q_i = q_const;
// uint64_t qq_i = qq_const;
extern double F_inv_vals[n2 * n2];

extern MatPoly F_mp;
extern MatPoly H_mp;

// MatPoly S_mp;
// MatPoly Sp_mp;

extern vector<MatPoly> neg1s_mp;
extern MatPoly half_mp;

extern uint64_t *constant_neg1s;
extern uint64_t *constant_half;

extern size_t IDX_TARGET;
extern size_t IDX_DIM0;

class FurtherDimsLocals {
    public:
    uint64_t *result;
    uint64_t *cts;
    uint64_t *scratch_cts1; // 数据库与第一个维度扩展成的密文单位向量点乘后的结果
    uint64_t *scratch_cts2;
    uint64_t *scratch_cts_double1;
    uint64_t *scratch_cts_double2;
    // uint64_t *scratch_cts_double3;
    // uint64_t *scratch_cts_double4;

    size_t num_per;
    size_t num_bytes_C;

    FurtherDimsLocals(size_t np) {
        num_per = np;
        num_bytes_C = sizeof(uint64_t) * num_per * n1 * n2 * 2 * poly_len;
    }

    void allocate() {
        result = (uint64_t*) malloc(2*num_bytes_C);
        cts = (uint64_t*) malloc(2*num_bytes_C);
        scratch_cts1 = (uint64_t*) malloc(num_bytes_C);
        scratch_cts2 = (uint64_t*) malloc(num_bytes_C);
        scratch_cts_double1 = (uint64_t*) malloc(m2 / n1 * num_bytes_C);
        scratch_cts_double2 = (uint64_t*) malloc(m2 / n1 * num_bytes_C);
        // scratch_cts_double3 = (uint64_t*) malloc(m2 / n1 * num_bytes_C);
        // scratch_cts_double4 = (uint64_t*) malloc(m2 / n1 * num_bytes_C);

        clear();
    }

    void clear() {
        memset(result, 0, 2*num_bytes_C);
        memset(cts, 0, 2*num_bytes_C);
        memset(scratch_cts_double1, 0, m2 / n1 * num_bytes_C);
        memset(scratch_cts_double2, 0, m2 / n1 * num_bytes_C);
        // memset(scratch_cts_double3, 0, m2 / n1 * num_bytes_C);
        // memset(scratch_cts_double4, 0, m2 / n1 * num_bytes_C);
        memset(scratch_cts1, 0, num_bytes_C);
        memset(scratch_cts2, 0, num_bytes_C);
    }
};

class ExpansionLocals {
    public:
    uint64_t *cts; // 查询密文
    uint64_t *scratch_cts1;
    uint64_t *scratch_cts2;
    uint64_t *small_coeff_polys;
    uint64_t *reoriented_ciphertexts; // 重新排列后的查询密文

    size_t n1_padded;
    size_t split;

    ExpansionLocals(size_t m = mh) {
        n1_padded = 4; // n1 rounded up to power of 2
        split = m / n1;
    }

    void allocate(size_t num_expansions) {
        size_t num_expanded = 1 << num_expansions;
        size_t C_size_bytes = n1 * m1 * crt_count * poly_len * sizeof(uint64_t);
        size_t L_size_bytes = num_expanded * C_size_bytes;
        size_t smallpols_size_bytes = num_expanded/2 * split * C_size_bytes/2;
        cts = (uint64_t *)malloc(L_size_bytes);
        scratch_cts1 = NULL;//(uint64_t *)malloc(L_size_bytes);
        scratch_cts2 = NULL;//(uint64_t *)malloc(L_size_bytes);
        small_coeff_polys = NULL;//(uint64_t *)malloc(smallpols_size_bytes);

        size_t L_small_size_bytes = num_expanded * n1 * n0 * crt_count * poly_len * sizeof(uint64_t);
        reoriented_ciphertexts = (uint64_t *)aligned_alloc(32, L_small_size_bytes * n1_padded / n1);

        clear(num_expansions);
    }

    void allocate(
        size_t num_expansions,
        size_t m1
    ) {
        size_t num_expanded = 1 << num_expansions;
        size_t C_size_bytes = n1 * m1 * crt_count * poly_len * sizeof(uint64_t);
        size_t L_size_bytes = num_expanded * C_size_bytes;
        size_t smallpols_size_bytes = num_expanded * split * C_size_bytes;
        cts = (uint64_t *)malloc(L_size_bytes);
        scratch_cts1 = NULL;//(uint64_t *)malloc(L_size_bytes);
        scratch_cts2 = NULL;//(uint64_t *)malloc(L_size_bytes);
        small_coeff_polys = NULL;//(uint64_t *)malloc(smallpols_size_bytes);

        size_t L_small_size_bytes = num_expanded * n1 * n0 * crt_count * poly_len * sizeof(uint64_t);
        reoriented_ciphertexts = (uint64_t *)aligned_alloc(64, L_small_size_bytes * n1_padded / n1);

        memset(cts, 0, L_size_bytes);
        // memset(scratch_cts1, 0, L_size_bytes);
        // memset(scratch_cts2, 0, L_size_bytes);
        // memset(small_coeff_polys, 0, smallpols_size_bytes);
        memset(reoriented_ciphertexts, 0, L_small_size_bytes);
    }

    void clear(size_t num_expansions) {
        size_t num_expanded = 1 << num_expansions;
        size_t C_size_bytes = n1 * m1 * crt_count * poly_len * sizeof(uint64_t);
        size_t L_size_bytes = num_expanded * C_size_bytes;
        size_t smallpols_size_bytes = num_expanded/2 * split * C_size_bytes/2;
        size_t L_small_size_bytes = num_expanded * n1 * n0 * crt_count * poly_len * sizeof(uint64_t);

        memset(cts, 0, L_size_bytes);
        // memset(scratch_cts1, 0, L_size_bytes);
        // memset(scratch_cts2, 0, L_size_bytes);
        // memset(small_coeff_polys, 0, smallpols_size_bytes);
        memset(reoriented_ciphertexts, 0, L_small_size_bytes);
    }
};

// namespace uuid {
//     static random_device              rd;
//     static mt19937                    gen(rd());
//     static uniform_int_distribution<> dis(0, 15);
//     static uniform_int_distribution<> dis2(8, 11);

//     string generate_uuid_v4() {
//         stringstream ss;
//         int i;
//         ss << hex;
//         for (i = 0; i < 8; i++) {
//             ss << dis(gen);
//         }
//         ss << "-";
//         for (i = 0; i < 4; i++) {
//             ss << dis(gen);
//         }
//         ss << "-4";
//         for (i = 0; i < 3; i++) {
//             ss << dis(gen);
//         }
//         ss << "-";
//         ss << dis2(gen);
//         for (i = 0; i < 3; i++) {
//             ss << dis(gen);
//         }
//         ss << "-";
//         for (i = 0; i < 12; i++) {
//             ss << dis(gen);
//         };
//         return ss.str();
//     }
// }


extern const char *socket_path;

void cpu_crt(uint64_t *out, const uint64_t *inp, size_t num_polys);
void runConversionTest();
void buildGadget(MatPoly &G);
void build_from_constants(MatPoly &mat, vector<vector<uint64_t> > inp);
uint64_t quickReduceQ(__uint128_t z);
void printScalarMat(MatPoly &A);
void runConversionImprovedTest();
MatPoly getSkVec(size_t idx);
void slow_poly_mul(__uint128_t *res, const __uint128_t *a, const __uint128_t *b, uint64_t modulus);
MatPoly slow_matmul_mod_arbitrary(const MatPoly &a, const MatPoly &b, uint64_t modulus);

void do_MatPol_test();
void setup_constants();
void generate_gadgets();
void load_db(
    uint64_t*& B,
    size_t num_expansions, 
    size_t further_dims, 
    size_t IDX_TARGET
);
MatPoly Query(
    size_t num_expansions, 
    size_t further_dims, 
    size_t IDX_TARGET
);
void ServerExpand(
    MatPoly cv,
    ExpansionLocals expansionLocals,
    uint64_t *g_Q_nttd,
    uint64_t *g_Q_crtd,
    size_t num_expansions,
    size_t further_dims
);
MatPoly Answer(
    uint64_t *B,
    MatPoly cv,
    ExpansionLocals expansionLocals,
    FurtherDimsLocals furtherDimsLocals,
    const uint64_t *g_C_fft_crtd,      // expand this ct
    uint64_t *g_Q_crtd,          // further dims query
    const uint64_t *g_Ws_fft,           // use this setup data
    size_t num_expansions,
    size_t further_dims
); 
void do_test(
    uint64_t *B,
    size_t num_expansions, 
    size_t further_dims, 
    size_t IDX_TARGET
);
void generate_setup_and_query(
    uint64_t **g_C_fft_crtd,      // expand this ct
    uint64_t **g_Q_crtd,          // further dims query
    uint64_t **g_Ws_fft,          // use this setup data
    bool encodeCompressedSetupData,
    size_t num_expansions,
    size_t further_dims
);

void generate_random_pt(MatPoly &M);

void process_crtd_query(
    uint64_t *B,
    ExpansionLocals expansionLocals,
    FurtherDimsLocals furtherDimsLocals,
    const uint64_t *g_C_fft_crtd,      // expand this ct
    const uint64_t *g_Q_crtd,          // further dims query
    const uint64_t *g_Ws_fft,           // use this setup data
    size_t num_expansions,
    size_t further_dims,
    size_t IDX_TARGET
);
void encryptRegev(
    MatPoly &out,           // n1 x m_in
    MatPoly G,               // n0 x m_in
    MatPoly sigma
);

// #ifndef _OPENMP
//     void omp_set_num_threads(int num_threads)
//     {
//     }

//     int omp_get_num_threads(void)
//     {
//     return 1;
//     }

//     int omp_get_max_threads(void)
//     {
//     return 1;
//     }

//     int omp_get_thread_num(void)
//     {
//     return 0;
//     }
// #endif

// 用于面向对象封装
// void generate_setup(
//     uint64_t **g_Ws_fft,          // use this setup data
//     bool encodeCompressedSetupData,
//     bool setNewSecret = true,
//     int num_expansions_h,
//     const MatPoly &G
// );

// void generate_query(
//     uint64_t **g_C_fft_crtd,      // expand this ct
//     uint64_t **g_Q_crtd,           // further dims query
//     size_t further_dims
// );

// void record(const string &s);

// void cpu_crt_to_ucompressed_and_ntt(uint64_t *out, const uint64_t *inp, size_t num_polys);

// void process_query_fast(
//     uint64_t *B,
//     const uint64_t *expansion_query_ct,      // expand this ct
//     const uint64_t *setup_data,
//     const uint64_t *further_dims_query_ct,          // further dims query
//     const uint64_t *further_dims_query_ct_neg,
//     ExpansionLocals expansion_locals,    // must be cleared
//     FurtherDimsLocals further_dims_locals, // must be cleared
//     size_t num_expansions,
//     size_t further_dims
// );

void modswitch(uint64_t *out, const uint64_t *inp);
double check_final(FurtherDimsLocals furtherDimsLocals, bool modswitch_on_server);
void print_summary(size_t num_expansions, size_t further_dims);
void process_crtd_query(
    uint64_t *B,
    ExpansionLocals expansionLocals,
    FurtherDimsLocals furtherDimsLocals,
    const uint64_t *g_C_fft_crtd,      // expand this ct
    uint64_t *g_Q_crtd,          // further dims query
    const uint64_t *g_Ws_fft,           // use this setup data
    size_t num_expansions,
    size_t further_dims,
    size_t IDX_TARGET
);
void generate_setup(
    uint64_t **g_Ws_fft,          // use this setup data
    bool encodeCompressedSetupData,
    bool setNewSecret,
    const MatPoly &G
);
void generate_query(
    uint64_t **g_C_fft_crtd,      // expand this ct
    uint64_t **g_Q_crtd,           // further dims query
    size_t further_dims
);
extern MatPoly pt_real;
void start_timing();
double end_timing();
MatPoly genPtQuery(size_t num_expansions, size_t further_dims, size_t IDX_TARGET);

MatPoly AnswerWithoutModSwitch(
    uint64_t *B,
    MatPoly cv,
    ExpansionLocals expansionLocals,
    FurtherDimsLocals furtherDimsLocals,
    const uint64_t *g_C_fft_crtd,      // expand this ct
    uint64_t *g_Q_crtd,          // further dims query
    const uint64_t *g_Ws_fft,           // use this setup data
    size_t num_expansions,
    size_t further_dims
);

MatPoly CompressWithModSwitch(MatPoly cts);

void cpu_crt_to_ucompressed_and_ntt(uint64_t *out, const uint64_t *inp, size_t num_polys);

void reorient_Q(uint64_t *out, const uint64_t *inp);

void reorientCiphertexts(uint64_t *out, const uint64_t *inp, size_t dim0, size_t n1_padded);

void multiplyQueryByDatabase(
    uint64_t * __restrict__ output,
    const uint64_t * __restrict__ reorientedCiphertexts,
    const uint64_t * __restrict__ database,
    size_t dim0, 
    size_t num_per
);

void nttInvAndCrtLiftCiphertexts(size_t num_per, FurtherDimsLocals furtherDimsLocals);

void foldOneFurtherDimension(
    size_t cur_dim, size_t num_per, 
    const uint64_t *query_ct, const uint64_t *query_ct_neg,
    FurtherDimsLocals locals
);

extern MatPoly G2;
