#include <iostream>
#include <unordered_map>
#include <vector>
#include <functional>
#include <limits>
#include <mod.h>

struct CuckooCode {
    size_t k;
    size_t d; // d choices
    double r; // total buckets = ceil(k * r)

    CuckooCode(size_t k, size_t d, double r) : k(k), d(d), r(r)
    {

    }
};

// 插入元素
bool CuckooInsert(
    std::unordered_map<size_t, size_t>& elements, // 输出结果：保存 bucket -> [current key] 
    std::unordered_map<size_t, std::vector<size_t>>& buckets, // 输入：保存 K -> [bucket 1, ..., bucket d] 的映射
    const size_t& key, // 输入：待插入的键（数据库的索引）
    size_t attempt, // 尝试次数
    std::mt19937_64 rng // 随机性
);

std::vector<std::vector<uint64_t *>> PBC_Encode(
    CuckooCode *self, 
    uint64_t *DB, 
    size_t num_expansions, 
    size_t further_dims
);

std::unordered_map<size_t, std::vector<size_t>> PBC_GetSchedule(
    CuckooCode *self, 
    const std::vector<size_t>& keys
);

MatPoly PBC_Decode(std::vector<MatPoly> respond);
std::vector<unsigned char> serialize(size_t value);
