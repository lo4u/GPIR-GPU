#include "pir.h"

class MultiPirParam
{
public:
    
    size_t num_expansions;
    size_t further_dims;

    MultiPirParam(size_t num_expansions, size_t further_dims): num_expansions(num_expansions), further_dims(further_dims)
    {

    }

};


class MultiPirClient
{
public:
    
    std::vector<PirClient> handles;
    uint64_t **g_Ws_fft;
    bool encodeCompressedSetupData;
    std::vector<std::vector<uint64_t>> godOracle;
    std::unordered_map<size_t, size_t> indexes;

    MultiPirClient(std::vector<MultiPirParam> params, std::vector<size_t> IDX_TARGETs, std::vector<std::vector<uint64_t>> godOracle, std::unordered_map<size_t, size_t> indexes): godOracle(godOracle), indexes(indexes)
    {
        assert(params.size() == IDX_TARGETs.size());
        size_t pLen = params.size();
        for (size_t i = 0 ; i < pLen ; i++) {
            PirClient temp(IDX_TARGETs[i], params[i].num_expansions, params[i].further_dims);
            handles.push_back(temp);
        }
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

    std::vector<MatPoly> MultiClientQuery() {
        std::vector<MatPoly> cvs;
        size_t hLen = handles.size();
        for (size_t i = 0 ; i < hLen ; i++) {
            MatPoly cv = handles[i].ClientQuery();
            cvs.push_back(cv);
        }

        return cvs;
    }

    std::vector<MatPoly> MultiPackedIndexGen() {
        std::vector<MatPoly> packed_ind_pt_vec(godOracle.size());
        #pragma omp parallel for
        for (size_t bucket = 0 ; bucket < godOracle.size() ; ++bucket) {
            if (indexes.find(bucket) != indexes.end()) {
                packed_ind_pt_vec[bucket] = handles[bucket].PackIndex();
            } else {
                MatPoly temp(1, 1, false);
                for (size_t i = 0 ; i < temp.rows * temp.cols * coeff_count ; i++) {
                    temp.data[i] = 0;
                }
                packed_ind_pt_vec[bucket] = temp;
            }
        }
        
        return packed_ind_pt_vec;
    }

    std::vector<MatPoly> MultiPackedIndexEncrypt(std::vector<MatPoly> packed_ind_pt_vec) {
        std::vector<MatPoly> cvs(packed_ind_pt_vec.size());
        #pragma omp parallel for
        for (size_t i = 0 ; i < packed_ind_pt_vec.size() ; i++) {
            assert(packed_ind_pt_vec[i].rows == 1);
            assert(packed_ind_pt_vec[i].cols == 1);
            cvs[i] = encryptSimpleRegev(packed_ind_pt_vec[i]);
        }
        
        return cvs;
    }

    std::vector<MatPoly> MultiAnswerDecrypt(std::vector<MatPoly> rs) {
        std::vector<MatPoly> pts(rs.size());
        #pragma omp parallel for
        for (size_t i = 0 ; i < rs.size() ; i++) {
            pts[i] = handles[0].DecodeAnswer(rs[i], true);
        }
        
        return pts;
    }
};

class MultiPirServer
{
public:

    vector<PirServer> handles;

    MultiPirServer(std::vector<MultiPirParam> params, vector<uint64_t *> DBs, vector<MatPoly> cvs)
    {
        assert(params.size() == DBs.size());
        assert(DBs.size() == cvs.size());
        size_t pLen = params.size();

        for (size_t i = 0 ; i < pLen ; i++) {
            PirServer temp(DBs[i], cvs[i], params[i].num_expansions, params[i].further_dims);
            handles.push_back(temp);
        }
    }

    std::vector<MatPoly> MultiServerAnswer() {
        std::vector<MatPoly> rs(handles.size());
        size_t hLen = handles.size();
        // #pragma omp parallel for
        for (size_t i = 0 ; i < hLen ; i++) {
            // cout << "========" << i << "=========" << endl;
            handles[i].genQuery();
            MatPoly r = handles[i].ServerAnswer();
            rs[i] = r;
        }

        return rs;
    }

    std::vector<MatPoly> MultiServerRespondWithoutModSwitch() {
        std::vector<MatPoly> cts;
        size_t hLen = handles.size();
        for (size_t i = 0 ; i < hLen ; i++) {
            // cout << "========" << i << "=========" << endl;
            handles[i].genQuery();
            MatPoly ct = handles[i].ServerRespondWithoutModSwitch();
            cts.push_back(ct);
        }

        return cts;
    }

    std::vector<MatPoly> MultiServerCompressWithModSwitch(std::vector<MatPoly> cts) {
        std::vector<MatPoly> rs;
        for (size_t i = 0 ; i < cts.size() ; i++) {
            MatPoly r = handles[0].ServerCompressWithModSwitch(cts[i]);
            rs.push_back(r);
        }

        return rs;
    }
};

std::vector<std::vector<size_t>> OracleHelper(
    CuckooCode *self, 
    size_t num_expansions, 
    size_t further_dims
);

// OracleTy类型定义
using OracleTy = std::tuple<std::vector<std::vector<size_t>>, std::vector<MultiPirParam>>;

OracleTy getOracle(
    CuckooCode *code, 
    size_t num_expansions, 
    size_t further_dims
);

void reBuildDB(
    uint64_t*& B, 
    size_t num_expansions, 
    size_t further_dims,
    const std::vector<uint64_t *> bucket
);

extern uint64_t *DB_test;

void print_pt_crt(MatPoly pts_encd_test);

void genDB(size_t num_expansions, size_t further_dims);

void print_item_of_DB(uint64_t *database, size_t index, size_t num_expansions, size_t further_dims);

void print_bucket(std::vector<uint64_t *> collections);

void print_item_of_oracle(std::vector<std::vector<size_t>> oracle, size_t bucket);

void print_DB(uint64_t *database, size_t num_expansions, size_t further_dims);

void PBC_test(CuckooCode *code);

void print_item_of_oracle(std::vector<std::vector<size_t>> oracle, size_t bucket);

void test_large_db(CuckooCode *code);

void print_item_of_bucket(std::vector<uint64_t *> collections, size_t index);