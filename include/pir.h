#include "spiral.h"
#include "cuckoo.h"
#include <vector>

extern uint64_t maskDB; // 使用1ULL来确保是64位的常量
extern MatPoly pts_encd_test;
extern MatPoly pt_encd_correct_test;
extern MatPoly pt_test;
extern MatPoly pt_correct_test;
extern MatPoly pt_real_test;
extern std::vector<MatPoly> pt_reals;

class PirClient
{
public:

    size_t IDX_TARGET;
    size_t num_expansions;
    size_t further_dims;
    uint64_t **g_Ws_fft;          // use this setup data
    bool encodeCompressedSetupData;

    PirClient(size_t index, size_t num_expansions, size_t further_dims): IDX_TARGET(index), num_expansions(num_expansions), further_dims(further_dims)
    {
        encodeCompressedSetupData = false;
    }

    void setup() {
        generate_setup(
            g_Ws_fft, 
            encodeCompressedSetupData, 
            true,
            G_hat
        );
    }

    MatPoly ClientQuery() {
        MatPoly cv = Query(num_expansions, further_dims, IDX_TARGET);
        return cv;
    }
    
    MatPoly PackIndex() {
        MatPoly sigma = genPtQuery(num_expansions, further_dims, IDX_TARGET);
        return sigma;
    }

    MatPoly DecodeAnswer(MatPoly total_resp, bool modswitch_on_server) {
        MatPoly ct(n1, n2, false);

        MatPoly ct_inp(n1, n2, false);
        // MatPoly total_resp(n1, n2, false);

        MatPoly r_end(S_mp.rows, ct.cols, false);

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
        return M_result;
    }
};

class PirServer
{
public:
     
    uint64_t *database;
    MatPoly cv;
    ExpansionLocals expansionLocals;
    FurtherDimsLocals furtherDimsLocals;
    uint64_t *g_C_fft_crtd;
    uint64_t *g_Q_crtd;
    uint64_t *g_Ws_fft;
    bool encodeCompressedSetupData;
    size_t num_expansions;
    size_t further_dims;

    PirServer(uint64_t *DB, MatPoly cv, size_t num_expansions, size_t further_dims): database(DB), cv(cv), num_expansions(num_expansions), further_dims(further_dims), furtherDimsLocals(1 << further_dims)
    {
        expansionLocals.allocate(num_expansions);

        furtherDimsLocals.allocate();

        encodeCompressedSetupData = false;
    }

    void genQuery() {
        generate_query(
            &g_C_fft_crtd,      // expand this ct
            &g_Q_crtd,           // further dims query
            further_dims
        );
    }

    MatPoly ServerAnswer() {
        MatPoly r = Answer(
            database,
            cv,
            expansionLocals,
            furtherDimsLocals,
            g_C_fft_crtd,
            g_Q_crtd,
            g_Ws_fft,
            num_expansions,
            further_dims
        );
        return r;
    }
    
    MatPoly ServerRespondWithoutModSwitch() {
        return AnswerWithoutModSwitch(
            database,
            cv,
            expansionLocals,
            furtherDimsLocals,
            g_C_fft_crtd,
            g_Q_crtd,
            g_Ws_fft,
            num_expansions,
            further_dims
        );
    }

    MatPoly ServerCompressWithModSwitch(MatPoly ct) {
        return CompressWithModSwitch(ct);
    }
};

double check_corr(
    FurtherDimsLocals furtherDimsLocals, 
    bool modswitch_on_server, 
    size_t seq
);

void test_oo(
    uint64_t *B, 
    size_t num_expansions, 
    size_t further_dims, 
    size_t IDX_TARGET
);

void genDataBase(
    uint64_t*& DB_test, 
    size_t num_expansions, 
    size_t further_dims, 
    std::vector<size_t> IDX_TARGETs
);

void test_trival_multi_PIR(
    uint64_t *B, size_t num_expansions, 
    size_t further_dims, 
    std::vector<size_t> IDX_TARGETs
);

double check_relation(
    MatPoly total_resp, 
    bool modswitch_on_server, 
    size_t seq
);