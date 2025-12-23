#include "multipir.h"

std::vector<std::vector<size_t>> OracleHelper(
    CuckooCode *self,
    size_t num_expansions, 
    size_t further_dims
) {
    size_t dim0 = 1 << num_expansions; // 第一维度的总数
    size_t num_per = 1 << further_dims; // 数据库总数据量 / 第一维度
    size_t total_n = dim0 * num_per;

    size_t total_buckets = (size_t)ceil((double)self->k * self->r);
    std::vector<std::vector<size_t>> collections(total_buckets);

    for (size_t i = 0 ; i < total_n ; i++) {
        // 对键（索引 i）进行序列化
        const size_t size = sizeof(size_t);
        std::vector<unsigned char> bytes(size, 0);
        bytes = serialize(i);

        std::vector<size_t> bucket_choices;

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
            collections[bucket].push_back(i);
        }
        
        // cout << "buckets" << endl;
    }

    return collections;
}

size_t nextPowerOfTwo(size_t n) {
    // 如果n已经是2的幂，直接返回n
    if ((n & (n - 1)) == 0) {
        return n;
    }
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n |= n >> 32;
    return n + 1;
}

size_t log2OfPowerOfTwo(size_t n) {
    size_t log = 0;
    while (n > 1) {
        n >>= 1;
        log++;
    }
    return log;
}

OracleTy getOracle(
    CuckooCode *code, 
    size_t num_expansions, 
    size_t further_dims
) {
    auto oracle = OracleHelper(code, num_expansions, further_dims);
    std::vector<MultiPirParam> sizes;
    for (const auto& vec : oracle) {
        size_t bLen = vec.size();
        // cout << "bLen is " << bLen << endl;
        size_t bLen_padded = nextPowerOfTwo(bLen);
        // cout << "bLen_padded is " << bLen_padded << endl;
        size_t log_bLen_padded = log2OfPowerOfTwo(bLen_padded);
        size_t num_expansions_prime = (size_t)ceil(((double)log_bLen_padded) / 2);
        if (log_bLen_padded % 2 == 0) {
            num_expansions_prime++;
        } 
        // cout << "log_bLen_padded is " << log_bLen_padded << endl;
        // cout << "num_expansions_prime " << num_expansions_prime << endl;
        size_t further_dims_prime = log_bLen_padded - num_expansions_prime;
        MultiPirParam temp(num_expansions_prime, further_dims_prime);
        sizes.push_back(temp);
    }

    return std::make_tuple(oracle, sizes);
}

void print_pt_crt(MatPoly pts_encd_test) {
    // uint64_t *BB = (uint64_t *)malloc(n0 * n2 * crt_count * poly_len * sizeof(uint64_t));
    for (size_t m = 0; m < n0; m++) {
        for (size_t c = 0; c < n2; c++) {
            // memcpy(BB, &pts_encd_test.data[(m * n2 + c) * crt_count * coeff_count], crt_count * coeff_count * sizeof(uint64_t));
            for (size_t z = 0; z < poly_len; z++) {
                cout << pts_encd_test.data[(m * n2 + c) * crt_count * coeff_count + z] << ' ' << pts_encd_test.data[(m * n2 + c) * crt_count * coeff_count + z + poly_len] << ' ';
            }
        }
    }
    cout << endl;
}

void reBuildDB(
    uint64_t*& B, 
    size_t num_expansions, 
    size_t further_dims,
    const std::vector<uint64_t *> bucket
) {
    size_t dim0 = 1 << num_expansions; // 第一维度的总数
    size_t num_per = 1 << further_dims; // 数据库总数据量 / 第一维度
    
    size_t num_bytes_B = sizeof(uint64_t) * dim0 * num_per * n0 * n2 * poly_len; // 2 * poly_len;
    B = (uint64_t *)aligned_alloc(64, num_bytes_B);
    memset(B, 0, num_bytes_B);

    // pts_encd_test = MatPoly(n0, n2);
    size_t bLen = bucket.size();

    uint64_t *BB = (uint64_t *)malloc(n0 * n2 * poly_len * sizeof(uint64_t));
    for (size_t i = 0 ; i < bLen ; i++) {
        // cop(pts_encd_test, bucket[i]);
        // print_pt_crt(pts_encd_test);
        size_t ii = i % num_per;
        size_t j = i / num_per;
        for (size_t m = 0 ; m < n0; m++) {
            for (size_t c = 0 ; c < n2 ; c++) {
                memcpy(BB, &bucket[i][(m * n2 + c) * coeff_count], coeff_count * sizeof(uint64_t));
                for (size_t z = 0 ; z < poly_len ; z++) {
                    size_t idx = z * (num_per * n2 * dim0 * n0) +
                                 ii * (n2 * dim0 * n0) + 
                                 c * (dim0 * n0) +
                                 j * (n0) + m;
                    B[idx] = BB[z];
                }
                // free(&bucket[i][(m * n2 + c) * coeff_count]);
            }
        }
    }
    free(BB);
}

uint64_t *DB_test;

void genDB(size_t num_expansions, size_t further_dims) {
    size_t dim0 = 1 << num_expansions; // 第一维度的总数
    size_t num_per = 1 << further_dims; // 数据库总数据量 / 第一维度
    size_t total_n = dim0 * num_per;
    
    size_t num_bytes_B = sizeof(uint64_t) * dim0 * num_per * n0 * n2 * poly_len; // 2 * poly_len;
    cout << "num_bytes_B: " << num_bytes_B << endl;
    DB_test = (uint64_t *)aligned_alloc(64, num_bytes_B);
    memset(DB_test, 0, num_bytes_B);

    uint64_t *BB = (uint64_t *)malloc(n0 * n2 * crt_count * poly_len * sizeof(uint64_t));
    size_t numBytesPlaintextRaw = n0 * n0 * num_bits_q * poly_len / 8;
    uint64_t *pt_raw = (uint64_t *)malloc(numBytesPlaintextRaw);
    memset(pt_raw, 0, numBytesPlaintextRaw);
    uint64_t *pt_buf = (uint64_t *)malloc(n0 * n0 * crt_count * poly_len * sizeof(uint64_t));
    memset(pt_buf, 0, n0 * n0 * crt_count * poly_len * sizeof(uint64_t));
    
    cout << "starting generation of db" << endl;
    MatPoly H_nttd = to_ntt(H_mp);
    uint64_t *H_encd = H_nttd.data;
    MatPoly pt_tmp(n0, n0, false);
    MatPoly pt_encd_raw(n0, n2, false);
    pts_encd_test = MatPoly(n0, n2);
    pt_test = MatPoly(n0, n0);
    for (size_t i = 0 ; i < total_n ; i++) {
        generate_random_pt(pt_tmp);

        pt_encd_raw = pt_tmp;
        for (size_t pol = 0; pol < n0 * n2 * poly_len; pol++) {
            int64_t val = (int64_t) pt_encd_raw.data[pol];
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
        to_ntt(pts_encd_test, pt_encd_raw);
        print_pt_crt(pts_encd_test);

        // b': i c n z j m
        size_t ii = i % num_per;
        size_t j = i / num_per;
        for (size_t m = 0; m < n0; m++) {
            for (size_t c = 0; c < n2; c++) {
                memcpy(BB, &pts_encd_test.data[(m * n2 + c) * crt_count * coeff_count], crt_count * coeff_count * sizeof(uint64_t));
                for (size_t z = 0; z < poly_len; z++) {
                    size_t idx = z * (num_per * n2 * dim0 * n0) +
                                 ii * (n2 * dim0 * n0) + 
                                 c * (dim0 * n0) +
                                 j * (n0) + m;
                    DB_test[idx] = BB[z] | (BB[poly_len + z] << 32);
                }
            }
        }
    }
    free(BB);
    cout << "done loading/generating db." << endl;
}

void print_item_of_DB(uint64_t *database, size_t index, size_t num_expansions, size_t further_dims) {
    size_t dim0 = 1 << num_expansions;
    size_t num_per = 1 << further_dims;
    size_t ii = index % num_per;
    size_t j = index / num_per;
    size_t count = 0;
    for (size_t m = 0 ; m < n0; m++) {
        for (size_t c = 0 ; c < n2 ; c++) {
            for (size_t z = 0 ; z < poly_len ; z++) {
                size_t idx = z * (num_per * n2 * dim0 * n0) +
                             ii * (n2 * dim0 * n0) + 
                             c * (dim0 * n0) +
                             j * (n0) + m;
                cout << (database[idx] & maskDB) << ' ' << ((database[idx] >> 32) & maskDB) << ' ';
                if (count++ > 5) {
                    cout << endl;
                    return;
                }
            }
        }
    }
    cout << endl;
}

void print_DB(uint64_t *database, size_t num_expansions, size_t further_dims) {
    size_t dim0 = 1 << num_expansions;
    size_t num_per = 1 << further_dims;
    size_t total_n = dim0 * num_per;

    for (size_t i = 0 ; i < total_n ; i++) {
        cout << "item " << i;
        size_t ii = i % num_per;
        size_t j = i / num_per;
        cout << " is ";
        for (size_t m = 0 ; m < n0; m++) {
            for (size_t c = 0 ; c < n2 ; c++) {
                for (size_t z = 0 ; z < poly_len ; z++) {
                    size_t idx = z * (num_per * n2 * dim0 * n0) +
                                 ii * (n2 * dim0 * n0) + 
                                 c * (dim0 * n0) +
                                 j * (n0) + m;
                    cout << (database[idx] & maskDB) << ' ' << ((database[idx] >> 32) & maskDB) << ' ';
                }
            }
        }
        cout << endl;
    }
}

void print_bucket(std::vector<uint64_t *> collections) {
    for (size_t i = 0 ; i < collections.size() ; i++) {
        for (size_t m = 0 ; m < n0; m++) {
            for (size_t c = 0 ; c < n2 ; c++) {
                for (size_t z = 0 ; z < poly_len ; z++) {
                    size_t idx = m * n2 * poly_len + c * poly_len + z;
                    cout << (collections[i][idx] & maskDB) << ' ' << ((collections[i][idx] >> 32) & maskDB) << ' ';
                }
            }
        }
        cout << endl;
    }
}

void print_item_of_oracle(std::vector<std::vector<size_t>> oracle, size_t bucket) {
    size_t bLen = oracle[bucket].size();
    for (size_t i = 0 ; i < bLen ; i++) {
        cout << oracle[bucket][i] << ' ';
    }
    cout << endl;
}

void PBC_test(CuckooCode *code) {
    auto collections = PBC_Encode(code, DB_test, 3, 2);
    OracleTy god = getOracle(code, 3, 2);
    auto godness = std::get<0>(god);
    auto godson = std::get<1>(god);

    // 生成键
    std::vector<size_t> keys;
    for (size_t i = 0; i < code->k; ++i) {
        keys.push_back(i);
    }

    // 得到映射
    auto schedule = PBC_GetSchedule(code, keys);

    // 验证调度计划是否有效
    for (const auto& pair : schedule) {
        size_t key = pair.first;
        const auto& buckets = pair.second;
        assert(buckets.size() == 1);
        cout << "================== Start Check Whether the " << key << " Element in " << buckets[0] << " Bucket =========================" << endl;
        print_bucket(collections[buckets[0]]);
        cout << "====================================================================" << endl;
        print_item_of_DB(DB_test, key, num_expansions, further_dims);
        print_item_of_oracle(godness, buckets[0]);
        cout << "================== Finish Check Whether the " << key << " Element in " << buckets[0] << " Bucket =========================" << endl;
    }

    std::vector<uint64_t *> DBs(collections.size());
    size_t cLen = collections.size();
    for (size_t i = 0 ; i < cLen ; i++) {
        cout << "=================== Start to Rebuild the " << i << " Bucket ===================" << endl;
        DBs[i] = NULL;
        print_bucket(collections[i]);
        reBuildDB(DBs[i], godson[i].num_expansions, godson[i].further_dims, collections[i]);
        cout << "==============================================================" << endl;
        print_DB(DBs[i], godson[i].num_expansions, godson[i].further_dims);
        cout << "=================== Finish to Rebuild the " << i << " Bucket ===================" << endl;
    }
}

void test_large_db(CuckooCode *code) {
    auto collections = PBC_Encode(code, DB_test, 7, 5);
    OracleTy god = getOracle(code, 7, 5);
    auto godness = std::get<0>(god);
    auto godson = std::get<1>(god);

    // 生成键
    std::vector<size_t> keys;
    for (size_t i = 0 ; i < code->k ; ++i) {
        keys.push_back(i);
    }

    // 得到映射
    auto schedule = PBC_GetSchedule(code, keys);

    // 验证调度计划是否有效
    for (const auto& pair : schedule) {
        size_t key = pair.first;
        const auto& buckets = pair.second;
        assert(buckets.size() == 1);
        // cout << "================== Start Check Whether the " << key << " Element in " << buckets[0] << " Bucket =========================" << endl;
        // print_bucket(collections[buckets[0]]);
        // cout << "====================================================================" << endl;
        // print_item_of_DB(DB_test, key);
        // print_item_of_oracle(godness, buckets[0]);
        // cout << "================== Finish Check Whether the " << key << " Element in " << buckets[0] << " Bucket =========================" << endl;
    }

    std::vector<uint64_t *> DBs(collections.size());
    size_t cLen = collections.size();
    for (size_t i = 0 ; i < cLen ; i++) {
        cout << "=================== Start to Rebuild the " << i << " Bucket ===================" << endl;
        DBs[i] = NULL;
        // print_bucket(collections[i]);
        reBuildDB(DBs[i], godson[i].num_expansions, godson[i].further_dims, collections[i]);
        // cout << "==============================================================" << endl;
        // print_DB(DBs[i], godson[i].num_expansions, godson[i].further_dims);
        cout << "=================== Finish to Rebuild the " << i << " Bucket ===================" << endl;
    }

}

void print_item_of_bucket(std::vector<uint64_t *> collections, size_t index) {
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



