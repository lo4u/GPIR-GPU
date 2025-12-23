#include <tuple>
#include "multipir.h"

size_t HASH1(
    size_t R,
    const std::vector<unsigned char> data,
    size_t modulus
);

std::vector<size_t> GenRandVec(
    size_t colIndex, 
    size_t vLen, 
    size_t R, 
    size_t bondWidth
);

std::vector<MatPoly> MulRBMatWithCts(
    std::vector<std::vector<size_t>> M, 
    std::vector<MatPoly> cts
);

MatPoly MulScalarAndPt(size_t scalar, MatPoly pt);

MatPoly ConvertModToQi(MatPoly pt_tmp);

MatPoly ConvertModToPdb(MatPoly pt_encd_raw);

std::vector<double> gaussianElimination(
    const std::vector<std::vector<double>>& A, 
    const std::vector<double>& b
);

std::vector<MatPoly> GaussianElimination(
    const std::vector<std::vector<size_t>>& A, 
    const std::vector<MatPoly>& b
);

uint64_t getValue(MatPoly input);

std::vector<MatPoly> MulRandBondMat(
    const std::vector<std::vector<size_t>>& A, 
    const std::vector<MatPoly>& x
);

std::vector<size_t> sortColumnsByFirstNonZero(
    std::vector<std::vector<size_t>>& out, 
    const std::vector<std::vector<size_t>>& matrix
);

std::vector<MatPoly> SolveLinearSystem(
    const std::vector<std::vector<size_t>>& A, 
    const std::vector<MatPoly>& b,
    bool mode
);

std::vector<MatPoly> GaussianEliminationWithUnderdeterminedMatrix(
    const std::vector<std::vector<size_t>>& A, 
    const std::vector<MatPoly>& b
);

void test_correctness_of_OCD(
    std::vector<std::vector<size_t>> M, 
    std::vector<MatPoly> p_hat,
    std::vector<MatPoly> real_pts
);

MatPoly InnerVecMul(
    std::vector<size_t> a,
    std::vector<MatPoly> b
);

