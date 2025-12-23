#include "cuckoo.h"
// issue: The buckets need to be padded to power of 2

void print_item_of_buckets(std::vector<uint64_t *> collections, size_t index) {
    uint64_t maskDB = (1ULL << 32) - 1; // 使用1ULL来确保是64位的常量
    size_t count = 0;
    for (size_t m = 0 ; m < n0; m++) {
        for (size_t c = 0 ; c < n2 ; c++) {
            for (size_t z = 0 ; z < poly_len ; z++) {
                size_t idx = m * n2 * poly_len + c * poly_len + z;
                cout << (collections[index][idx] & maskDB) << ' ' << ((collections[index][idx] >> 32) & maskDB) << ' ';
                if (count++ > 5) {
                    cout << endl;
                    return;
                }
            }
        }
    }
    cout << endl;
}

const size_t MAX_ATTEMPTS = 1000;
uint64_t maskDB_prime = (1ULL << 32) - 1; // 使用1ULL来确保是64位的常量
void print_value(uint64_t *value) {
    for (size_t i = 0 ; i < n0 * n2 * poly_len ; i++) {
        cout << (value[i] & maskDB_prime) << " " << ((value[i] >> 32) & maskDB_prime) << ' ';
    }
    cout << endl;
}

std::vector<unsigned char> serialize(size_t value) {
    const size_t size = sizeof(size_t);
    std::vector<unsigned char> bytes(size, 0);

    // 将 size_t 值的每个字节复制到 vector 中
    for (size_t i = 0; i < size; ++i) {
        bytes[i] = static_cast<unsigned char>((value >> (i * 8)) & 0xFF);
    }

    return bytes;
}

// 插入元素
bool CuckooInsert(
    std::unordered_map<size_t, size_t>& elements, // 输出结果：保存 bucket -> [current key] 
    std::unordered_map<size_t, std::vector<size_t>>& buckets, // 输入：保存 K -> [bucket 1, ..., bucket d] 的映射
    const size_t& key, // 输入：待插入的键（数据库的索引）
    size_t attempt = 0, // 尝试次数
    std::mt19937_64 rng = std::mt19937_64() // 随机性
) {
    if (attempt >= MAX_ATTEMPTS) {
        return false;
    }

    // Case 1: 检查是否有空桶可以插入
    for (auto bucket : buckets[key]) {
        if (elements.find(bucket) == elements.end()) {
            elements[bucket] = key;
            return true;
        }
    }

    // Case 2: 所有可能的桶都满了，需要重新定位一个现有元素
    size_t index = std::uniform_int_distribution<size_t>(0, buckets[key].size() - 1)(rng);
    size_t chosen_bucket = buckets[key][index];

    // 插入新键，并获取之前插入的键
    auto old_key = elements[chosen_bucket];
    elements[chosen_bucket] = key;

    // 重新插入被移动的键
    return CuckooInsert(elements, buckets, old_key, attempt + 1, rng);
}

std::vector<std::vector<uint64_t *>> PBC_Encode(
    CuckooCode *self, 
    uint64_t *DB, 
    size_t num_expansions, 
    size_t further_dims
) {
    size_t dim0 = 1 << num_expansions; // 第一维度的总数
    size_t num_per = 1 << further_dims; // 数据库总数据量 / 第一维度
    size_t total_n = dim0 * num_per;
    // cout << "==================== Start to Encode the Database ====================" << endl;

    size_t total_buckets = (size_t)ceil((double)self->k * self->r);
    std::vector<std::vector<uint64_t *>> collections(total_buckets);
    // cout << "total_buckets is " << total_buckets << endl;

    for (size_t i = 0 ; i < total_n ; i++) {
        // b': i c n z j m
        size_t ii = i % num_per;
        size_t j = i / num_per;
        size_t count = 0;
        uint64_t *value = (uint64_t *)malloc(n0 * n2 * poly_len * sizeof(uint64_t));
        // 得到索引 i 处数据库中的元素
        for (size_t m = 0; m < n0; m++) {
            for (size_t c = 0; c < n2; c++) {
                for (size_t z = 0; z < poly_len; z++) {
                    size_t idx = z * (num_per * n2 * dim0 * n0) +
                                 ii * (n2 * dim0 * n0) + 
                                 c * (dim0 * n0) +
                                 j * (n0) + m;
                    value[count++] = DB[idx];
                    // free(&DB[idx]);
                }
            }
        }
        // 对键（索引 i）进行序列化
        const size_t size = sizeof(size_t);
        std::vector<unsigned char> bytes(size, 0);
        bytes = serialize(i);

        std::vector<size_t> bucket_choices;

        // cout << "item " << i << " is ";       
        // print_value(value);
        // cout << "It's in the ";
        // 将键映射到 d 个桶（不重复）
        for (size_t id = 0 ; id < self->d ; id++) {
            size_t nonce = 0;

            // 计算 bucket = bucket = sha_d(key) % k;
            auto bucket = hash_and_mod(id, nonce, bytes, total_buckets);

            // 确保每个键被映射到 *不同的* 桶
            while (std::find(bucket_choices.begin(), bucket_choices.end(), bucket) != bucket_choices.end()) {
                nonce++;
                bucket = hash_and_mod(id, nonce, bytes, total_buckets);
            }

            // cout << bucket << ", ";
            bucket_choices.push_back(bucket);
            collections[bucket].push_back(value);
            // if (i == 5) {
            //     cout << "5 \'s" << " bucket is " << bucket << endl;
            //     print_item_of_buckets(collections[bucket], collections[bucket].size() - 1);
            // }
        }
        
        // cout << "buckets" << endl;
    }

    // cout << "==================== Finish to Encode the Database ====================" << endl;

    return collections;
}

std::unordered_map<size_t, std::vector<size_t>> PBC_GetSchedule(CuckooCode *self, const std::vector<size_t>& keys) {
    assert(keys.size() <= self->k);

    size_t total_buckets = static_cast<size_t>(std::ceil(self->k * self->r));

    std::unordered_map<size_t, std::vector<size_t>> buckets; // map containing K -> [bucket 1, ..., bucket d]

    for (const auto& key : keys) {
        const size_t size = sizeof(size_t);
        std::vector<unsigned char> bytes(size, 0);
        bytes = serialize(key);
        std::vector<size_t> bucket_choices;

        // Map entry's key to d buckets (no repeats)
        for (size_t id = 0; id < self->d; ++id) {
            size_t nonce = 0;
            // The following computes bucket = sha_d(key) % k;
            size_t bucket = hash_and_mod(id, nonce, bytes, total_buckets);
            
            // Ensure each key maps to *different* buckets
            while (std::find(bucket_choices.begin(), bucket_choices.end(), bucket) != bucket_choices.end()) {
                nonce++;
                bucket = hash_and_mod(id, nonce, bytes, total_buckets);
            }

            bucket_choices.push_back(bucket);
            // if (key == 5) {
            //     cout << "5 \'s bucket is " << bucket << endl;
            // }
        }

        buckets[key] = bucket_choices;
    }

    // This is a variant of the Insert algorithm in cuckoo hashing (Pagh and Rodler).
    // One difference is we only have 1 table and d hash functions.
    // The difference is that we are doing this for retrieval rather than insertion! (the keys
    // have already been inserted). What this means is that cuckoo hashing is being applied
    // with respect to the client's keys (not the keys that the storage server received).
    // This is a crucial but subtle difference.

    std::mt19937_64 rng(std::random_device{}());
    std::unordered_map<size_t, size_t> elements; // map containing bucket -> [current key]

    for (const auto& key : keys) {
        if (!CuckooInsert(elements, buckets, key, 0, rng)) {
            throw std::runtime_error("无法生成有效的映射");
        }
    }

    std::unordered_map<size_t, std::vector<size_t>> schedule;

    for (const auto& pair : elements) {
        schedule[pair.second].push_back(pair.first);
    }

    assert(schedule.size() == keys.size());

    return schedule;
}

MatPoly PBC_Decode(std::vector<MatPoly> respond) {
    assert(respond.size() == 1);
    return respond[0];
}

