#include "rbmatrices.h"

// 生成条带开始位置的哈希函数
size_t HASH1(
    size_t R, 
    const std::vector<unsigned char> data, 
    size_t modulus
) {
    // Use SHA-256 for hashing
    SHA256_CTX sha256;
    unsigned char hash[SHA256_DIGEST_LENGTH];

    // Initialize SHA-256 context
    SHA256_Init(&sha256);

    // Hash the custom string (id and nonce)
    std::stringstream ss;
    ss << R;
    SHA256_Update(&sha256, ss.str().c_str(), ss.str().size());

    // Hash the data
    SHA256_Update(&sha256, data.data(), data.size());

    // Finalize the hash
    SHA256_Final(hash, &sha256);

    // Convert hash to a big integer
    boost::multiprecision::cpp_int int_value;
    for (unsigned int i = 0; i < SHA256_DIGEST_LENGTH; ++i) {
        int_value = (int_value << 8) | hash[i];
    }

    // Perform modulo operation
    return static_cast<size_t>(int_value % modulus);
}

// 基于SHA256的哈希函数
std::string Sha256Hash(const std::string& input) {
    unsigned char hash[SHA256_DIGEST_LENGTH];
    SHA256_CTX sha256;
    SHA256_Init(&sha256);
    SHA256_Update(&sha256, input.c_str(), input.size());
    SHA256_Final(hash, &sha256);

    std::stringstream ss;
    for (int i = 0; i < SHA256_DIGEST_LENGTH; i++) {
        ss << std::hex << std::setw(2) << std::setfill('0') << (int)hash[i];
    }
    return ss.str();
}

// 生成随机向量的哈希函数
std::vector<size_t> HASH2(size_t vLen, size_t R, size_t i) {
    std::stringstream ss;
    ss << R << "||" << i;
    std::string input = ss.str();
    std::string hash = Sha256Hash(input);

    std::vector<size_t> random_vec(vLen, 0);

    // 将哈希值转换为0和1的向量
    for (size_t j = 0 ; j < vLen ; ++j) {
        // 从哈希值中取位
        size_t bit_pos = j % (SHA256_DIGEST_LENGTH * 8);
        size_t byte_pos = bit_pos / 8;
        size_t bit = ((hash[byte_pos] >> (bit_pos % 8)) & 1);
        // cout << bit << ' ';
        random_vec[j] = bit;
        // cout << random_vec[j] << ' ';
    }

    return random_vec;
}

// 生成随机条带向量
std::vector<size_t> GenRandVec(
    size_t colIndex, 
    size_t vLen, 
    size_t R, 
    size_t bondWidth
) {
    auto s = HASH1(R, serialize(colIndex), vLen - bondWidth + 1);
    auto u = HASH2(bondWidth, R, colIndex);
    std::vector<size_t> vi(vLen, 0);

    for (size_t j = 0 ; j < bondWidth ; j++) {
        vi[s + j] = u[j];
    }
    return vi;
}

std::vector<MatPoly> MulRBMatWithCts(
    std::vector<std::vector<size_t>> M, 
    std::vector<MatPoly> cts
) {
    std::vector<MatPoly> res;
    size_t rows = M.size();
    size_t cols = M[0].size();
    assert(cols == cts.size());
    // 将 cts 转换到 NTT 域上
    for (size_t i = 0 ; i < cols ; i++) {
        cts[i] = to_ntt(cts[i]);
    }
    // 运行矩阵乘法
    for (size_t i = 0 ; i < rows ; i++) {
        MatPoly temp(n1, n2);
        for (size_t j = 0 ; j < cols ; j++) {
            if (M[i][j] == 1) {
                temp = add(temp, cts[j]);
            }
        }
        res.push_back(from_ntt(temp));
    }

    return res;
}

MatPoly MulScalarAndPt(size_t scalar, MatPoly pt) {
    if (!pt.isNTT) {
        std::runtime_error("MulScalarAndPt Failed");
    }
    MatPoly scalar_prime(1, 1, false);
    scalar_prime.data[0] = scalar;
    return mul_by_const(to_ntt(scalar_prime), pt);
}

// result is in NTT
MatPoly ConvertModToQi(const MatPoly pt_tmp) {
    MatPoly pt_encd_raw(n0, n2, false);
    for (size_t pol = 0; pol < n0 * n2 * poly_len; pol++) {
        int64_t val = (int64_t) pt_tmp.data[pol];
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

    return to_ntt(pt_encd_raw);
}

// result is not in NTT
MatPoly ConvertModToPdb(MatPoly pt_encd_raw) {
    MatPoly A = from_ntt(pt_encd_raw);
    MatPoly res(pt_encd_raw.rows, pt_encd_raw.cols, false);
    for (size_t i = 0; i < A.rows * A.cols * poly_len; i++) {
        int64_t a = (int64_t) (A.data[i] % Q_i);
        if (a >= Q_i / 2) a -= Q_i;
        while (a < 0)
        {
            a += p_db;
        }
        res.data[i] = a % p_db;
    }
    return res;
}

const double EPSILON = 1e-10;

// 检查解是否有效
bool isValid(double value) {
    return std::abs(value) > EPSILON;
}

void print_matrix_of_scalar(std::vector<std::vector<MatPoly>> Matrix) {
    cout << "================================" << endl;
    for (size_t i = 0 ; i < Matrix.size() ; i++) {
        for (size_t j = 0 ; j < Matrix[0].size() ; j++) {
            cout << getValue(Matrix[i][j]) << ' ';
        }
        cout << endl;
    }
    cout << "================================" << endl;
}

// 高斯消元法求解线性方程组
std::vector<double> gaussianElimination(
    const std::vector<std::vector<double>>& A, 
    const std::vector<double>& b
) {
    int m = A.size(); // 行数
    int n = A[0].size(); // 列数

    // 构建增广矩阵
    std::vector<std::vector<double>> augmentedMatrix(m, std::vector<double>(n + 1, 0.0));
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            augmentedMatrix[i][j] = A[i][j];
        }
        augmentedMatrix[i][n] = b[i];
    }

    // 前向消元
    for (int i = 0; i < std::min(m, n); ++i) {
        // 找到最大的主元
        int maxRow = i;
        for (int k = i + 1; k < m; ++k) {
            if (std::abs(augmentedMatrix[k][i]) > std::abs(augmentedMatrix[maxRow][i])) {
                maxRow = k;
            }
        }

        // 交换行
        std::swap(augmentedMatrix[maxRow], augmentedMatrix[i]);

        // 检查是否有解
        if (!isValid(augmentedMatrix[i][i])) {
            throw std::runtime_error("No solution or infinite solutions");
        }

        // 消元
        for (int k = i + 1; k < m; ++k) {
            double factor = augmentedMatrix[k][i] / augmentedMatrix[i][i];
            for (int j = i; j <= n; ++j) {
                augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
            }
        }
    }

    // 回代
    std::vector<double> x(n, 0.0);
    for (int i = std::min(m, n) - 1; i >= 0; --i) {
        x[i] = augmentedMatrix[i][n] / augmentedMatrix[i][i];
        for (int k = i - 1; k >= 0; --k) {
            augmentedMatrix[k][n] -= augmentedMatrix[k][i] * x[i];
        }
    }

    // 如果行数大于列数，处理自由变量
    if (m > n) {
        std::vector<double> solution(n, 0.0);
        for (int i = 0; i < n; ++i) {
            solution[i] = x[i];
        }
        for (int i = n; i < m; ++i) {
            double sum = augmentedMatrix[i][n];
            for (int j = 0; j < n; ++j) {
                sum -= augmentedMatrix[i][j] * solution[j];
            }
            if (isValid(augmentedMatrix[i][i])) {
                throw std::runtime_error("No solution");
            }
        }
        return solution;
    }

    return x;
}

// 获取 value 值
uint64_t getValue(MatPoly input) {
    if (!input.isNTT) {
        std::runtime_error("Get Value Failed");
    }
    assert(input.rows == 1);
    assert(input.cols == 1);
    assert(input.isNTT);
    auto temp = from_ntt(input);
    return temp.data[0];
}

// 检查解是否有效
bool isZero(MatPoly A) {
    if (!A.isNTT) {
        std::runtime_error("Is zero Failed");
    }
    auto temp = from_ntt(A);
    auto helper = MatPoly(1, 1, false);
    return is_eq(temp, helper);
}

// return -a*b in NTT
MatPoly MulNeg(MatPoly a, MatPoly b) {
    if (!a.isNTT || !b.isNTT) {
        std::runtime_error("Mul Neg Failed");
    }
    return to_ntt(invert(from_ntt(mul_by_const(a, b))));
}

std::vector<MatPoly> GaussianElimination(
    const std::vector<std::vector<size_t>>& A, 
    const std::vector<MatPoly>& b
) {
    size_t rows = A.size();
    size_t cols = A[0].size();
    assert(rows == b.size());
    // print_matrix_of_scalar(A);

    // 构建增广矩阵
    std::vector<std::vector<MatPoly>> augmentedMatrix(rows, std::vector<MatPoly>(cols + 1));
    for (size_t i = 0 ; i < rows ; ++i) {
        for (size_t j = 0; j < cols ; ++j) {
            augmentedMatrix[i][j] = MatPoly(1, 1, false);
            augmentedMatrix[i][j].data[0] = A[i][j];
            augmentedMatrix[i][j] = to_ntt(augmentedMatrix[i][j]);
        }
        augmentedMatrix[i][cols] = ConvertModToQi(b[i]);
    }
    print_matrix_of_scalar(augmentedMatrix);

    size_t r = 0;
    for (size_t c = 0 ; c < cols ; c++) {
        if (r >= rows) {
            break;
        }
        // 找到最大的主元
        size_t maxRow = r;
        for (size_t k = r + 1 ; k < rows ; ++k) {
            if (getValue(augmentedMatrix[k][c]) > getValue(augmentedMatrix[maxRow][c])) {
                maxRow = k;
            }
        }

        if (isZero(augmentedMatrix[maxRow][c])) {
            continue;
        }
        // 交换行
        std::swap(augmentedMatrix[maxRow], augmentedMatrix[r]);

        for (size_t i = 0 ; i < rows ; i++) {
            if (i == r) {
                continue;
            }
            MatPoly factor = MulScalarAndPt(
                          inv_mod(getValue(augmentedMatrix[r][c]), Q_i), 
                          augmentedMatrix[i][c]
                          );
            for (size_t j = c ; j <= cols ; j++) {
                augmentedMatrix[i][j] = add(
                                    augmentedMatrix[i][j], 
                                    MulNeg(factor, augmentedMatrix[r][j])
                                    );
            }
        }
        r++;
        print_matrix_of_scalar(augmentedMatrix);
    }
    // print_matrix_of_scalar(augmentedMatrix);

    if (r < rows) {
        for (size_t i = r ; i < rows ; i++) {
            if (isZero(augmentedMatrix[i][cols])) {
                std::runtime_error("No solution");
            }
        }
    }

    std::vector<MatPoly> x(cols, MatPoly(1, 1, true));
    for (size_t i = 0 ; i < rows ; i++) {
        size_t j;
        for (j = i ; j < cols ; j++) {
            if (!isZero(augmentedMatrix[i][j])) {
                x[j] = MulScalarAndPt(
                    inv_mod(getValue(augmentedMatrix[i][j]), Q_i),
                    augmentedMatrix[i][cols]
                );
                break;
            }
        }
        if (j == cols) {
            break;
        }
    }

    return x;
}

std::vector<MatPoly> GaussianEliminationWithUnderdeterminedMatrix(
    const std::vector<std::vector<size_t>>& A, 
    const std::vector<MatPoly>& b
) {
    size_t rows = A.size();
    size_t cols = A[0].size();
    assert(rows == b.size());
    // print_matrix_of_scalar(A);

    // 构建增广矩阵
    std::vector<std::vector<MatPoly>> augmentedMatrix(rows, std::vector<MatPoly>(cols + 1));
    for (size_t i = 0 ; i < rows ; ++i) {
        for (size_t j = 0; j < cols ; ++j) {
            augmentedMatrix[i][j] = MatPoly(1, 1, false);
            augmentedMatrix[i][j].data[0] = A[i][j];
            augmentedMatrix[i][j] = to_ntt(augmentedMatrix[i][j]);
        }
        augmentedMatrix[i][cols] = to_ntt(b[i]);
    }
    // print_matrix_of_scalar(augmentedMatrix);

    size_t r = 0;
    for (size_t c = 0 ; c < cols ; c++) {
        if (r >= rows) {
            break;
        }
        // 找到最大的主元
        size_t maxRow = r;
        for (size_t k = r + 1 ; k < rows ; ++k) {
            if (getValue(augmentedMatrix[k][c]) > getValue(augmentedMatrix[maxRow][c])) {
                maxRow = k;
            }
        }

        if (isZero(augmentedMatrix[maxRow][c])) {
            continue;
        }
        // 交换行
        std::swap(augmentedMatrix[maxRow], augmentedMatrix[r]);

        for (size_t i = 0 ; i < rows ; i++) {
            if (i == r) {
                continue;
            }
            MatPoly factor = MulScalarAndPt(
                          inv_mod(getValue(augmentedMatrix[r][c]), Q_i), 
                          augmentedMatrix[i][c]
                          );
            for (size_t j = c ; j <= cols ; j++) {
                augmentedMatrix[i][j] = add(
                                    augmentedMatrix[i][j], 
                                    MulNeg(factor, augmentedMatrix[r][j])
                                    );
            }
        }
        r++;
        // print_matrix_of_scalar(augmentedMatrix);
    }
    // print_matrix_of_scalar(augmentedMatrix);

    if (r < rows) {
        for (size_t i = r ; i < rows ; i++) {
            if (isZero(augmentedMatrix[i][cols])) {
                std::runtime_error("No solution");
            }
        }
    }

    std::vector<MatPoly> x(cols, MatPoly(1, 1, true));
    for (size_t i = 0 ; i < rows ; i++) {
        size_t j;
        for (j = i ; j < cols ; j++) {
            if (!isZero(augmentedMatrix[i][j])) {
                x[j] = MulScalarAndPt(
                    inv_mod(getValue(augmentedMatrix[i][j]), Q_i),
                    augmentedMatrix[i][cols]
                );
                break;
            }
        }
        if (j == cols) {
            break;
        }
    }

    return x;
}

std::vector<MatPoly> MulRandBondMat_1(
    const std::vector<std::vector<size_t>>& A, 
    const std::vector<MatPoly>& x
) {
    size_t rows = A.size();
    size_t cols = A[0].size();
    assert(cols == x.size());

    std::vector<MatPoly> x_temp(cols);
    for (size_t i = 0 ; i < cols ; i++) {
        x_temp[i] = to_ntt(x[i]);
    }

    std::vector<MatPoly> res(rows);
    for (size_t i = 0 ; i < rows ; i++) {
        res[i] = MatPoly(1, 1, false);
        res[i] = to_ntt(res[i]);
        for (size_t j = 0 ; j < cols ; j++) {
            if (A[i][j] == 1) {
                res[i] = add(res[i], x_temp[j]);
            }
        }
    }

    return res;
}

std::vector<MatPoly> MulRandBondMat(
    const std::vector<std::vector<size_t>>& A, 
    const std::vector<MatPoly>& x
) {
    size_t rows = A.size();
    size_t cols = A[0].size();
    assert(cols == x.size());

    std::vector<std::vector<MatPoly>> A_temp(rows, std::vector<MatPoly>(cols));
    for (size_t i = 0 ; i < rows ; ++i) {
        for (size_t j = 0; j < cols ; ++j) {
            A_temp[i][j] = MatPoly(1, 1, false);
            A_temp[i][j].data[0] = A[i][j];
            A_temp[i][j] = to_ntt(A_temp[i][j]);
        }
    }

    std::vector<MatPoly> x_temp(cols);
    for (size_t i = 0 ; i < cols ; i++) {
        x_temp[i] = to_ntt(x[i]);
    }

    std::vector<MatPoly> res(rows);
    for (size_t i = 0 ; i < rows ; i++) {
        res[i] = MatPoly(1, 1, false);
        res[i] = to_ntt(res[i]);
        for (size_t j = 0 ; j < cols ; j++) {
            res[i] = add(res[i], mul_by_const(A_temp[i][j], x_temp[j]));
        }
    }

    return res;
}

// 函数用于找到每列中第一个非零元素的行索引
std::vector<size_t> findFirstNonZeroPerColumn(const std::vector<std::vector<size_t>>& matrix) {
    std::vector<size_t> firstNonZero(matrix[0].size(), matrix.size()); // 初始化为-1，表示未找到非零元素
    for (size_t j = 0 ; j < matrix[0].size() ; j++) {
        for (size_t i = 0 ; i < matrix.size() ; i++) {
            if (matrix[i][j] != 0 && firstNonZero[j] == matrix.size()) {
                firstNonZero[j] = i;
                break;
            }
        }
    }
    return firstNonZero;
}

// 比较函数，用于根据第一个非零元素的行索引对列进行排序
bool compareColumns(const size_t& col1, const size_t& col2, const std::vector<size_t>& firstNonZero) {
    return firstNonZero[col1] < firstNonZero[col2];
}

// 主函数，用于对矩阵的列进行排序并返回原始列的位置
std::vector<size_t> sortColumnsByFirstNonZero(
    std::vector<std::vector<size_t>>& out, 
    const std::vector<std::vector<size_t>>& matrix
) {
    std::vector<size_t> firstNonZero = findFirstNonZeroPerColumn(matrix);
    std::vector<size_t> originalOrder(matrix[0].size()); // 保存原始列的顺序
    std::iota(originalOrder.begin(), originalOrder.end(), 0); // 填充为0, 1, 2, ..., n-1

    // 根据第一个非零元素的行索引对列进行排序
    std::sort(originalOrder.begin(), originalOrder.end(), [&firstNonZero](size_t col1, size_t col2) {
        return compareColumns(col1, col2, firstNonZero);
    });

    for (size_t j = 0 ; j < out[0].size() ; j++) {
        for (size_t i = 0 ; i < out.size() ; i++) {
            out[i][j] = matrix[i][originalOrder[j]];
        }
    }

    return originalOrder;
}

// 返回值不在 NTT 域中
std::vector<MatPoly> SolveLinearSystem(
    const std::vector<std::vector<size_t>>& A, 
    const std::vector<MatPoly>& b,
    bool mode = true
) {
    std::vector<std::vector<size_t>> res(A.size(), std::vector<size_t>(A[0].size(), 0));
    std::vector<size_t> sortedOrder = sortColumnsByFirstNonZero(res, A);
    std::vector<MatPoly> p_pi;
    // print_matrix_of_scalar(A);

    if (mode) {
        p_pi = GaussianEliminationWithUnderdeterminedMatrix(res, b);
    } else {
        p_pi = GaussianElimination(res, b);
    }
    std::vector<MatPoly> p(p_pi.size());
    for (size_t i = 0 ; i < p.size() ; i++) {
        p[sortedOrder[i]] = from_ntt(p_pi[i]);
    }

    return p;
}

void test_correctness_of_OCD(
    std::vector<std::vector<size_t>> M, 
    std::vector<MatPoly> p_hat,
    std::vector<MatPoly> real_pts
) {
    auto b_prime = MulRandBondMat_1(M, p_hat);
    for (size_t i = 0 ; i < b_prime.size() ; i++) {
        // cout << b_prime[i].isNTT << endl;
        // cout << getValue(b_prime[i]) << ' ';
        b_prime[i] = from_ntt(b_prime[i]);
    }
    // cout << endl;
    for (size_t i = 0 ; i < b_prime.size() ; i++) {
        cout << i << ' ' << is_eq(real_pts[i], b_prime[i]) << endl;
        // for (size_t k = 0; k < real_pts[i].rows * b_prime[i].cols * coeff_count; k++) {
        //     if (real_pts[i].data[k] != b_prime[i].data[k]) {
        //         cout << k << " " << real_pts[i].data[k] << ", " << b_prime[i].data[k] << endl;
        //     }
        // }
    }
}

MatPoly InnerVecMul(
    std::vector<size_t> a,
    std::vector<MatPoly> b
) {
    assert(a.size() == b.size());
    MatPoly res(n0, 1, true);
    
    for(size_t i = 0 ; i < a.size() ; i++) {
        if (a[i] == 1) {
            res = add(res, b[i]);
        }
    }
    return res;
}

