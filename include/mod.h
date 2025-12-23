#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <openssl/sha.h>
#include <openssl/evp.h>
#include <boost/multiprecision/cpp_int.hpp>
#include "spiral.h"

struct Tuple {
    size_t K; // 键（即为 DB 中的索引）
    uint64_t *V; // 值（即为所存储明文矩阵 MatPoly 的 data 字段）

    Tuple(size_t key) : K(key)
    {
        V = (uint64_t *) calloc(n0 * n2 * poly_len, sizeof(uint64_t));
    }
};

size_t hash_and_mod(
    size_t id, 
    size_t nonce, 
    const std::vector<unsigned char> data, 
    size_t modulus
);