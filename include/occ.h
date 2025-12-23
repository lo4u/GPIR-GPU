#include "rbmatrices.h"

class LSObvCompress
{
public:

    size_t num_all_cts; // 压缩前全部密文的数量（包括零的加密）
    size_t num_non_zero_pts; // 非零明文的数量
    double epsilon; // m <- (1 + epsilon) * t
    size_t bond_width; // 随机带的宽度

    LSObvCompress(size_t n, size_t t, double epsilon, size_t w): num_all_cts(n), num_non_zero_pts(t), epsilon(epsilon), bond_width(w)
    {

    }

    std::vector<MatPoly> ObvCompress(std::vector<MatPoly> cts) {
        assert(cts.size() == num_all_cts);
        size_t m = ceil((1 + epsilon) * (double)num_non_zero_pts);
        // 创建一个 m * n 的二维向量，初始值为 0
        std::vector<std::vector<size_t>> M(m, std::vector<size_t>(num_all_cts, 0));
        // 对矩阵的列进行迭代
        for (size_t i = 0 ; i < num_all_cts ; i++) {
            size_t R = HASH1(i, serialize(i), maskDB); // 获取随机性
            auto vi = GenRandVec(i, m, R, bond_width); // 生成随机向量
            for (size_t j = 0 ; j < m ; j++) {
                M[j][i] = vi[j];
            }
        }

        std::vector<MatPoly> ans = MulRBMatWithCts(M, cts);
        return ans;
    }

    std::vector<MatPoly> Decompress(
        std::vector<MatPoly> pts, 
        std::vector<size_t> realIndexes
    ) {
        size_t m = ceil((1 + epsilon) * (double)num_non_zero_pts);
        assert(pts.size() == m);
        assert(realIndexes.size() == num_non_zero_pts);
        // 创建一个 m * n 的二维向量，初始值为 0
        std::vector<std::vector<size_t>> M(m, std::vector<size_t>(num_non_zero_pts, 0));
        // 对矩阵的列进行迭代
        for (size_t i = 0 ; i < num_non_zero_pts ; i++) {
            size_t R = HASH1(realIndexes[i], serialize(realIndexes[i]), maskDB); // 获取随机性
            auto vi = GenRandVec(realIndexes[i], m, R, bond_width); // 生成随机向量
            for (size_t j = 0 ; j < m ; j++) {
                M[j][i] = vi[j];
            }
        }

        // solve linear system
        auto p_I = SolveLinearSystem(M, pts, false);

        return p_I;
    }
};


class LSObvDecompress
{
public:

    size_t num_all_pts; // 压缩前全部明文的数量（包括零的加密）
    size_t num_non_zero_pts; // 非零明文的数量
    double epsilon; // m <- (1 + epsilon) * t
    size_t bond_width; // 随机带的宽度

    LSObvDecompress(size_t n, size_t t, double epsilon, size_t w): num_all_pts(n), num_non_zero_pts(t), epsilon(epsilon), bond_width(w)
    {

    }

    std::vector<MatPoly> Compress(
        std::vector<MatPoly> pts, 
        std::vector<size_t> realIndexes
    ) {
        assert(pts.size() == num_all_pts);
        assert(realIndexes.size() == num_non_zero_pts);
        size_t m = ceil((1 + epsilon) * (double)num_non_zero_pts);
        cout << "m is " << m << endl;

        // 创建一个 t * m 的二维向量，初始值为 0
        std::vector<std::vector<size_t>> M(num_non_zero_pts, std::vector<size_t>(m, 0));
        for(size_t i = 0 ; i < M.size() ; i++) {
            size_t R = HASH1(realIndexes[i], serialize(realIndexes[i]), maskDB); // 获取随机性
            M[i] = GenRandVec(realIndexes[i], m, R, bond_width); // 生成随机向量
            // cout << realIndexes[i] << ": ";
            // for (size_t k = 0 ; k < M[i].size() ; k++) {
            //     cout << M[i][k] << ' ';
            // }
            // cout << endl;
        }

        std::vector<MatPoly> real_pts;
        for (size_t i = 0 ; i < realIndexes.size() ; i++) {
            real_pts.push_back(pts[realIndexes[i]]);
        }

        // solve linear system
        auto p_hat = SolveLinearSystem(M, real_pts, true);

        // test:
        // test_correctness_of_OCD(M, p_hat, real_pts);

        return p_hat;
    }

    std::vector<MatPoly> ObvDecompress(std::vector<MatPoly> cts_hat) {
        size_t m = ceil((1 + epsilon) * (double)num_non_zero_pts);
        std::vector<MatPoly> res;

        cout << "m is " << m << endl;

        // 创建一个 n * m 的二维向量，初始值为 0
        for(size_t i = 0 ; i < num_all_pts ; i++) {
            size_t R = HASH1(i, serialize(i), maskDB); // 获取随机性
            auto a = GenRandVec(i, m, R, bond_width); // 生成随机向量
            res.push_back(InnerVecMul(a, cts_hat));
            cout << i << ": ";
            for (size_t k = 0 ; k < a.size() ; k++) {
                cout << a[k] << ' ';
            }
            cout << endl;
        }

        return res;
    }
};

MatPoly genAMatWithV(size_t v);

MatPoly genBMatWithV(size_t v);

void test_mPIR_with_PBC(
    size_t num_expansions, 
    size_t further_dims
);

void print_const_ntt_matpoly(const MatPoly A);

void test_gaussian_elimination();

void test_sort_columns();

void test_solve_linear_system();

void test_mPIR_with_PBC_and_OCD(
    size_t num_expansions, 
    size_t further_dims
);

void test_mPIR_with_PBC_and_OCC(
    size_t num_expansions, 
    size_t further_dims
);

void test_pt_add_ntt();
