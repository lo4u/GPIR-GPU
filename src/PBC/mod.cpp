#include "mod.h"
const int SM3_DIGEST_LENGTH = 32;
// Utility function that computes a hash of the key and mods it by the given modulus
size_t hash_and_mod(
    size_t id, 
    size_t nonce, 
    const std::vector<unsigned char> data, 
    size_t modulus
) {
    // Use SM3 for hashing
    EVP_MD_CTX* mdctx;
    unsigned char hash[SM3_DIGEST_LENGTH]; // SM3 hash length is 32 bytes

    // Initialize SM3 context
    mdctx = EVP_MD_CTX_new();
    const EVP_MD* sm3 = EVP_sm3(); // Get SM3 algorithm from OpenSSL
    EVP_DigestInit_ex(mdctx, sm3, NULL);

    // Hash the custom string (id and nonce)
    std::stringstream ss;
    ss << id << nonce;
    EVP_DigestUpdate(mdctx, ss.str().c_str(), ss.str().size());

    // Hash the data
    EVP_DigestUpdate(mdctx, data.data(), data.size());

    // Finalize the hash
    EVP_DigestFinal_ex(mdctx, hash, NULL);

    // Clean up context
    EVP_MD_CTX_free(mdctx);

    // Convert hash to a big integer
    boost::multiprecision::cpp_int int_value = 0;
    for (unsigned int i = 0; i < SM3_DIGEST_LENGTH; ++i) {
        int_value = (int_value << 8) | hash[i];
    }

    // Perform modulo operation
    return static_cast<size_t>(int_value % modulus);
}

