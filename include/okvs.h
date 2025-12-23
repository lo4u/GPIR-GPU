#include "3H_GCT.h"

class GCTObvKVStore
{
public:
    
    CuckooCode code;
    size_t num_all_pts;

    GCTObvKVStore(size_t k, size_t d, double r, size_t n): code(k, d, r), num_all_pts(n) 
    {
        
    }

    std::vector<MatPoly> Encode(
        std::vector<MatPoly> pts, 
        std::vector<size_t> realIndexes
    ) {
        assert(pts.size() == num_all_pts);
        assert(realIndexes.size() == code.k);
        size_t m = ceil(code.r * (double)code.k);
        // cout << "m is " << m << endl;

        // Initialize empty vectors
        std::vector<MatPoly> L(m);
        // Initialize empty stack
        std::stack<VertexOriented> Stack;
        
        // Get real plaintexts
        std::vector<MatPoly> real_pts;
        for (size_t i = 0 ; i < realIndexes.size() ; i++) {
            real_pts.push_back(pts[realIndexes[i]]);
        }

        std::unordered_map<size_t, std::vector<size_t>> positions; // map containing idx -> [position 1, ..., position d]
        std::vector<Vertex> HyperGraph(m);
        HyperGraphSetup(HyperGraph);
        std::vector<std::vector<size_t>> NeiborMatrix(real_pts.size(), std::vector<size_t>(m, 0));
        size_t count = 0;

        for (const auto& idx : realIndexes) {
            const size_t size = sizeof(size_t);
            std::vector<unsigned char> bytes(size, 0);
            bytes = serialize(idx);
            std::vector<size_t> position_choices;

            // Map entry's key to d positions
            for (size_t id = 0; id < code.d; ++id) {
                size_t nonce = 0;
                // The following computes position = sha_d(idx) % m;
                size_t position = hash_and_mod(id, nonce, bytes, m);
            
                // Ensure each key maps to *different* buckets
                while (std::find(position_choices.begin(), position_choices.end(), position) != position_choices.end()) {
                    nonce++;
                    position = hash_and_mod(id, nonce, bytes, m);
                }

                position_choices.push_back(position);
                NeiborMatrix[count][position] = 1;
            }
              
            for (size_t jj = 0; jj < code.d; jj++) {
                HyperGraph[position_choices[jj]].vertex_id = position_choices[jj];
                HyperGraph[position_choices[jj]].edge_table[count] = Edge(count, position_choices[(jj+1)%3], position_choices[(jj+2)%3]);
            }

            count++;
            positions[idx] = position_choices;
        }

        std::vector<bool> flagList = peel_the_hypergraph(Stack, HyperGraph, NeiborMatrix);

        std::vector<std::vector<uint64_t>> core2 = choice_rows_in_neiborMatirx(flagList, NeiborMatrix);

        std::vector<MatPoly> core2_value = choice_position_in_real_pt(flagList, real_pts);

        // print_some_vars(
        //     NeiborMatrix,
        //     HyperGraph,
        //     flagList,
        //     Stack,
        //     core2
        // );
         
        std::vector<MatPoly> solution(m, MatPoly(1, 1, false));
        if (core2_value.size() > 0) {
            solution = SolveLinearSystem(core2, core2_value, true);
            // print_solutions(solution);
        }

        auto res = get_whole_solution(solution, Stack, real_pts, positions, realIndexes);
        // print_solutions(res);

        // for (size_t ii = 0; ii < realIndexes.size(); ii++) {
        //     auto position_choice = positions[realIndexes[ii]];
        //     auto temp = add(add(to_ntt(res[position_choice[0]]), to_ntt(res[position_choice[1]])), to_ntt(res[position_choice[2]]));
        //     temp = from_ntt(temp);
        //     cout << ii << ' ' << is_eq(temp, real_pts[ii]) << endl;
        // }

        return res;
    }

    std::vector<MatPoly> Decode(std::vector<MatPoly> cts_hat) {
        size_t m = ceil(code.r * (double)code.k);
        std::vector<MatPoly> res;

        // cout << "m is " << m << endl;

        for(size_t i = 0 ; i < num_all_pts ; i++) {
            const size_t size = sizeof(size_t);
            std::vector<unsigned char> bytes(size, 0);
            bytes = serialize(i);
            std::vector<size_t> position_choices;

            // Map entry's key to d positions
            for (size_t id = 0; id < code.d; ++id) {
                size_t nonce = 0;
                // The following computes position = sha_d(idx) % m;
                size_t position = hash_and_mod(id, nonce, bytes, m);
            
                // Ensure each key maps to *different* buckets
                while (std::find(position_choices.begin(), position_choices.end(), position) != position_choices.end()) {
                    nonce++;
                    position = hash_and_mod(id, nonce, bytes, m);
                }

                position_choices.push_back(position);
            }

            auto temp = add(add(cts_hat[position_choices[0]], cts_hat[position_choices[1]]), cts_hat[position_choices[2]]);
            res.push_back(temp);
        }

        return res;
    }
};

void test_GCT_okvs(
    size_t num_expansions, 
    size_t further_dims
);
