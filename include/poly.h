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
#include <stdexcept>

#include <assert.h>

using namespace std;

#include <constants.h>
#include "core.h"
#include "values.h"

// extern uint64_t *scratch;


struct MatPoly {
    size_t rows; // 矩阵的行数
    size_t cols; // 矩阵的列数
    uint64_t *data; // 矩阵中的数据
    bool isNTT; // 矩阵是否处于快速傅里叶变换形式

    MatPoly(size_t rs, size_t cs) : rows(rs), cols(cs), isNTT(true) // 新建函数 1
    {
        data = (uint64_t *) calloc(rows * cols * crt_count * coeff_count, sizeof(uint64_t));
    }

    MatPoly(size_t rs, size_t cs, bool isNTTFlag) : rows(rs), cols(cs), isNTT(isNTTFlag) // 新建函数 2
    {
        size_t factor = isNTT ? crt_count : 1;
        data = (uint64_t *) calloc(rows * cols * factor * coeff_count, sizeof(uint64_t));
    }

    MatPoly() : rows(0), cols(0), data(NULL), isNTT(true) // 新建函数 3
    {
        
    }

    uint64_t *operator[](size_t index) const { // 使用数组索引的方式访问矩阵的行。它返回指向当前行数据的指针
        return &data[index * cols * coeff_count];
    }

    MatPoly& operator=(const MatPoly& other) // 赋值操作符，用于复制一个 MatPoly 对象到另一个
    {
        if (this != &other) // not a self-assignment
        {
            this->rows = other.rows;
            this->cols = other.cols;
            this->isNTT = other.isNTT;
            size_t factor = (this->isNTT) ? crt_count : 1;
            data = (uint64_t *) calloc(rows * cols * factor * coeff_count, sizeof(uint64_t));
            memcpy(this->data, other.data, rows * cols * factor * coeff_count * sizeof(uint64_t));
        }
        return *this;
    }

};

// template<size_t s_rows, size_t s_cols>
// struct MatPolyNttd : MatPoly {
//     uint64_t data_backing[s_rows * s_cols * crt_count * poly_len * sizeof(uint64_t)];

//     MatPolyNttd() {
//         this->rows = s_rows;
//         this->cols = s_cols;
//         this->isNTT = true;
//         this->data = data_backing;
//     }

//     MatPoly& operator=(const MatPoly& other) = delete;
// };

void multiply(MatPoly& c, const MatPoly& a, const MatPoly& b);
void multiply_no_reduce(MatPoly& out, const MatPoly& a, const MatPoly& b);
MatPoly multiply(const MatPoly& a, const MatPoly& b); // return a * b

void to_ntt(MatPoly& out, const MatPoly& a);
void to_ntt_no_reduce(MatPoly& out, const MatPoly& a);
MatPoly to_ntt(const MatPoly& a);
void from_ntt(MatPoly& out, const MatPoly& a);
MatPoly from_ntt(const MatPoly& a);

void add(MatPoly& c, const MatPoly& a, const MatPoly& b);
MatPoly add(const MatPoly& a, const MatPoly& b); // return a + b
void add_into(MatPoly& c, const MatPoly& a, const MatPoly& b, size_t t_row, size_t t_col);
void invert(MatPoly& out, const MatPoly& a);
MatPoly invert(const MatPoly& a);
void vertical_merge(MatPoly& c, const MatPoly& a, const MatPoly& b);
void automorph(MatPoly& out, const MatPoly& a, uint64_t t);
MatPoly automorph(const MatPoly& a, uint64_t t);

MatPoly single_poly(uint64_t value);
void mul_by_const(MatPoly& out, const MatPoly& singlePoly, const MatPoly& a);
MatPoly mul_by_const(const MatPoly& singlePoly, const MatPoly& a);

uint64_t mul_and_reduce_qq(uint64_t a, uint64_t b);
uint64_t reduce_qq(uint64_t z);
uint64_t reduce_qq_from_u128(__uint128_t z);

void to_simple_crtd(uint64_t *out, const MatPoly& a);

void reduce_mod(MatPoly &a, uint64_t mod);
void cop(MatPoly& out, const MatPoly& a);
void cop(
    MatPoly& out, const MatPoly& a, 
    size_t s_row, size_t s_col, 
    size_t t_row, size_t t_col,
    size_t num_row, size_t num_col
);
void place(MatPoly& out, const MatPoly& a, size_t t_row, size_t t_col);
void pick(MatPoly& out, const MatPoly& a, size_t t_row, size_t t_col);
MatPoly pick(const MatPoly& a, size_t t_row, size_t t_col, size_t num_rows, size_t num_cols);

bool is_eq(const MatPoly &A, const MatPoly &B);

uint64_t crt_compose_pa(uint64_t x, uint64_t y);
uint64_t crt_compose(uint64_t x, uint64_t y, uint64_t z);

// void serialize(const MatPoly &mat, fstream &f);
// void deserialize(MatPoly &mat, size_t rs, size_t cs, fstream &f);

void decode_mat(MatPoly &mat, uint64_t *inp);

uint64_t rescale(uint64_t a, uint64_t inp_mod, uint64_t out_mod, bool nosign = false);
MatPoly getRescaled(const MatPoly &a, uint64_t inp_mod, uint64_t out_mod);

void divide_by_const(MatPoly& out, const MatPoly& a, uint64_t divisor, uint64_t mod);
void divide_by_const(MatPoly& out, const MatPoly& a, uint64_t numer, uint64_t denom, uint64_t mod);

inline uint64_t barrett_raw_u64(uint64_t input, uint64_t const_ratio_1, uint64_t modulus) {
    unsigned long long tmp[2];
    tmp[1] = static_cast<unsigned long long>(((static_cast<__uint128_t>(input) * static_cast<__uint128_t>(const_ratio_1)) >> 64));

    // Barrett subtraction
    tmp[0] = input - tmp[1] * modulus;

    // One more subtraction is enough
    return (tmp[0] >= modulus ? (tmp[0] - modulus) : (tmp[0]));
}

// return val % moduli[n];
inline uint64_t barrett_coeff(uint64_t val, size_t n) {
    // return val % moduli[n];
    return (n == 0) ? barrett_raw_u64(val, cr1_p, p_i)
                    : barrett_raw_u64(val, cr1_b, b_i);
    return val;
}

// 优化大数模运算
#ifdef MODULUS_56_BIT
inline uint64_t reduction_u128_qq_rough(uint64_t val) {
    return val;
}
inline uint64_t reduction_u128_qq(uint64_t val) {
    return val % b_i;
}
#else
inline uint64_t reduction_u128_qq_rough(uint64_t val) {
    uint64_t x = val & ((1UL << 31) - 1);
    uint64_t y = val >> 31;
    return (x + (y << 17) - y);
}
inline uint64_t reduction_u128_qq(uint64_t val) {
    uint64_t x = val & ((1UL << 31) - 1);
    uint64_t y = val >> 31;
    return (x + (y << 17) - y) % b_i;
}
#endif

// ---------------- 写入函数 ----------------
inline void saveMatPolyVector(const std::vector<MatPoly>& vec, const std::string& filename) {
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        throw std::runtime_error("Failed to open file for writing");
    }

    size_t vec_size = vec.size();
    out.write(reinterpret_cast<const char*>(&vec_size), sizeof(vec_size));

    for (const auto& mat : vec) {
        out.write(reinterpret_cast<const char*>(&mat.rows), sizeof(mat.rows));
        out.write(reinterpret_cast<const char*>(&mat.cols), sizeof(mat.cols));
        out.write(reinterpret_cast<const char*>(&mat.isNTT), sizeof(mat.isNTT));

        size_t factor = mat.isNTT ? crt_count : 1;
        size_t data_len = mat.rows * mat.cols * poly_len * factor;

        out.write(reinterpret_cast<const char*>(mat.data), data_len * sizeof(uint64_t));
    }

    out.close();
}

// ---------------- 读取函数 ----------------
inline std::vector<MatPoly> loadMatPolyVector(const std::string& filename) {
    std::ifstream in(filename, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Failed to open file for reading");
    }

    size_t vec_size;
    in.read(reinterpret_cast<char*>(&vec_size), sizeof(vec_size));

    std::vector<MatPoly> vec(vec_size);

    for (auto& mat : vec) {
        in.read(reinterpret_cast<char*>(&mat.rows), sizeof(mat.rows));
        in.read(reinterpret_cast<char*>(&mat.cols), sizeof(mat.cols));
        in.read(reinterpret_cast<char*>(&mat.isNTT), sizeof(mat.isNTT));

        size_t factor = mat.isNTT ? crt_count : 1;
        size_t data_len = mat.rows * mat.cols * poly_len * factor;

        mat.data = new uint64_t[data_len];
        in.read(reinterpret_cast<char*>(mat.data), data_len * sizeof(uint64_t));
    }

    in.close();
    return vec;
}

// ---------------- 写入二进制函数 ----------------
// ptr: 指向数组的指针
// length: 元素个数
// filename: 输出文件名
inline void saveBinary(const uint64_t* ptr, size_t length, const std::string& filename) {
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        throw std::runtime_error("Failed to open file for writing");
    }

    // 写入元素个数（可选，用于读取时知道长度）
    out.write(reinterpret_cast<const char*>(&length), sizeof(length));

    // 写入数据
    out.write(reinterpret_cast<const char*>(ptr), length * sizeof(uint64_t));
    out.close();
}

// ---------------- 读取二进制函数 ----------------
// 返回新分配的数组指针，length 会被设置为元素个数
inline void loadBinary(uint64_t *ptr, size_t& length, const std::string& filename) {
    std::ifstream in(filename, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Failed to open file for reading");
    }

    // 先读取元素个数
    in.read(reinterpret_cast<char*>(&length), sizeof(length));
    in.read(reinterpret_cast<char*>(ptr), length * sizeof(uint64_t));
    in.close();
}

// 计算一个多项式矩阵的大小,单位是字节
inline size_t calculateMatPolySize(const MatPoly& mat) {
    size_t factor = mat.isNTT ? crt_count : 1;
    size_t data_len = mat.rows * mat.cols * poly_len * factor;
    return sizeof(size_t) * 2 + data_len * sizeof(uint64_t) + sizeof(bool); // rows, cols, isNTT + data
}