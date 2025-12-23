#include "opgpu.cuh"

float time_mul_DB_with_Q = 0.0;

__global__ void hello_from_gpu() {
    const int bid = blockIdx.x;
    const int tid = threadIdx.x;
    const int id = blockIdx.x * blockDim.x + threadIdx.x;

    printf("Hello World from block %d and the thread %d, global id %d the GPU\n", bid, tid, id);
}

void print_hello_from_gpu() {
    hello_from_gpu<<<2, 4>>>();
    cudaDeviceSynchronize();
}

cudaError_t ErrorCheck(
    cudaError_t error_code,
    const char* filename,
    int lineNumber
) {
    if (error_code != cudaSuccess) {
        printf("CUDA error:\r\ncode=%d, name=%s, description=%s\r\nfile=%s, line%d\r\n",
                error_code, cudaGetErrorName(error_code), cudaGetErrorString(error_code), filename, lineNumber);
        return error_code;
    }
    return error_code;
}

void setGPU() {
    int iDeviceCount = 0;
    cudaError_t error = cudaGetDeviceCount(&iDeviceCount);

    if (error != cudaSuccess || iDeviceCount == 0) {
        printf("No CUDA campatable GPU found!\n");
        exit(-1);
    } else {
        printf("The count of GPUs is %d.\n", iDeviceCount);
    }

    int iDev = 0;
    error = cudaSetDevice(iDev);
    if (error != cudaSuccess) {
        printf("Fail to set GPU 0 for computing.\n");
        exit(-1);
    } else {
        printf("Set GPU 0 for computing.\n");
    }
}

Timer::Timer(){
    _timeElasped = 0;
    _cStart = std::chrono::high_resolution_clock::now();
    _cStop = std::chrono::high_resolution_clock::now();
    cudaEventCreate(&_gStart);
    cudaEventCreate(&_gStop);
}

Timer::~Timer(){
    cudaFree(_gStart);
    cudaFree(_gStop);
}

void Timer::start_gpu() {
    cudaEventRecord(_gStart, 0);
}

void Timer::stop_gpu() {
    cudaEventRecord(_gStop, 0);
}

void Timer::start_cpu() {
    _cStart = std::chrono::high_resolution_clock::now();
}

void Timer::stop_cpu() {
    _cStop = std::chrono::high_resolution_clock::now();
}

void Timer::duration_gpu(std::string msg){
    CUDA_CHECK(cudaEventSynchronize(_gStart));
    CUDA_CHECK(cudaEventSynchronize(_gStop));
    cudaEventElapsedTime(&_timeElasped, _gStart, _gStop);

    LOG("%-60s uses %.6lf ms", msg.c_str(), _timeElasped);
}

__global__ void addMatrix(int *out, int *in1, int *in2, const int nx, const int ny) {
    int ix = threadIdx.x + blockIdx.x * blockDim.x;
    int iy = threadIdx.y + blockIdx.y * blockDim.y;
    unsigned int idx = iy * nx + ix;
    if (ix < nx && iy < ny) {
        out[idx] = in1[idx] + in2[idx];
    }
}

void test_MatAdd() {
    setGPU();

    int nx = 16;
    int ny = 8;
    int nxy = nx * ny;
    size_t stBytesCount = nxy * sizeof(int);

    int *ipHost_A, *ipHost_B, *ipHost_C;
    ipHost_A = (int *)malloc(stBytesCount);
    ipHost_B = (int *)malloc(stBytesCount);
    ipHost_C = (int *)malloc(stBytesCount);
    if (ipHost_A != NULL && ipHost_B != NULL && ipHost_C != NULL) {
        for (size_t i = 0; i < nxy; i++) {
            ipHost_A[i] = i;
            ipHost_B[i] = i + 1;
        }
        memset(ipHost_C, 0, stBytesCount);
    } else {
        printf("Fail to allocate host memory!\n");
        exit(-1);
    }

    int *ipDevice_A, *ipDevice_B, *ipDevice_C;
    ErrorCheck(cudaMalloc((int **)&ipDevice_A, stBytesCount), __FILE__, __LINE__);
    ErrorCheck(cudaMalloc((int **)&ipDevice_B, stBytesCount), __FILE__, __LINE__);
    ErrorCheck(cudaMalloc((int **)&ipDevice_C, stBytesCount), __FILE__, __LINE__);
    if (ipDevice_A != NULL && ipDevice_B != NULL && ipDevice_C != NULL) {
        ErrorCheck(cudaMemcpy(ipDevice_A, ipHost_A, stBytesCount, cudaMemcpyHostToDevice), __FILE__, __LINE__);
        ErrorCheck(cudaMemcpy(ipDevice_B, ipHost_B, stBytesCount, cudaMemcpyHostToDevice), __FILE__, __LINE__);
        ErrorCheck(cudaMemcpy(ipDevice_C, ipHost_C, stBytesCount, cudaMemcpyHostToDevice), __FILE__, __LINE__);
    } else {
        printf("Fail to allocate device memory\n");
        free(ipHost_A);
        free(ipHost_B);
        free(ipHost_C);
        exit(1);
    }

    dim3 block(4, 4);
    dim3 grid((nx + block.x - 1) / block.x, (ny + block.y - 1) / block.y);
    printf("thread config:grid:<%d, %d>, block:<%d, %d>\n", grid.x, grid.y, block.x, block.y);

    addMatrix<<<grid, block>>>(ipDevice_C, ipDevice_A, ipDevice_B, nx, ny);
    ErrorCheck(cudaMemcpy(ipHost_C, ipDevice_C, stBytesCount, cudaMemcpyDeviceToHost), __FILE__, __LINE__);
    for (int i = 0; i < 10; i++) {
        printf("id=%d, matrix_A=%d, matrix_B=%d, result=%d\n", i + 1, ipHost_A[i], ipHost_B[i], ipHost_C[i]);
    }

    free(ipHost_A);
    free(ipHost_B);
    free(ipHost_C);

    ErrorCheck(cudaFree(ipDevice_A), __FILE__, __LINE__);
    ErrorCheck(cudaFree(ipDevice_B), __FILE__, __LINE__);
    ErrorCheck(cudaFree(ipDevice_C), __FILE__, __LINE__);
    
    ErrorCheck(cudaDeviceReset(), __FILE__, __LINE__);
}

__global__ void MulQueryDBKernal(
    uint64_t *output,
    const uint64_t *reorientedCiphertexts,
    const uint64_t *database,
    size_t dim0, 
    size_t num_per,
    const size_t poly_len = 2048,
    const size_t n0 = 2,
    const size_t n1 = 3,
    const size_t n2 = 2,
    const size_t crt_count = 2
) {
    size_t n1_padded = 4; // n1 rounded up to power of 2

    size_t coeff_id = threadIdx.x + blockIdx.x * blockDim.x;
    size_t num_per_id = threadIdx.y + blockIdx.y * blockDim.y;
    if (coeff_id < 2048 && num_per_id < num_per) {
        size_t idx_a_base = coeff_id * (2 * dim0 * n1_padded);
        size_t idx_b_base = coeff_id * (num_per * n2 * dim0 * n0);

        const uint32_t p_i = 268369921;
        const uint32_t b_i = 249561089;

        for (size_t c = 0; c < n2; c++) {
            // size_t idx_b = idx_b_base + i * (n2 * dim0 * n0) + c * (dim0 * n0);    
            __uint128_t sums_out_n0_0 = 0, sums_out_n0_1 = 0, sums_out_n0_2 = 0;
            __uint128_t sums_out_n1_0 = 0, sums_out_n1_1 = 0, sums_out_n1_2 = 0;

            // #pragma unroll 16
            for (size_t jm = 0; jm < dim0*2; jm++) {
                uint64_t b = database[idx_b_base++];
                        
                const uint64_t *v_a = &reorientedCiphertexts[idx_a_base + jm * n1_padded];
                uint64_t v_a0 = v_a[0];
                uint64_t v_a1 = v_a[1];
                uint64_t v_a2 = v_a[2];

                uint32_t b_lo = b;
                uint32_t b_hi = b >> 32L;

                uint32_t v_a0_lo = v_a0;
                uint32_t v_a0_hi = v_a0 >> 32L;

                uint32_t v_a1_lo = v_a1;
                uint32_t v_a1_hi = v_a1 >> 32L;

                uint32_t v_a2_lo = v_a2;
                uint32_t v_a2_hi = v_a2 >> 32L;
                        
                // do n0
                sums_out_n0_0 += (v_a0_lo) * (uint64_t)b_lo;
                sums_out_n0_1 += (v_a1_lo) * (uint64_t)b_lo;
                sums_out_n0_2 += (v_a2_lo) * (uint64_t)b_lo;

                // do n1
                sums_out_n1_0 += ((uint64_t)v_a0_hi) * b_hi;
                sums_out_n1_1 += ((uint64_t)v_a1_hi) * b_hi;
                sums_out_n1_2 += ((uint64_t)v_a2_hi) * b_hi;
            }

            // output n0
            size_t n = 0;
            size_t idx_c = num_per_id * (n1*n2*crt_count*poly_len) + c * (crt_count*poly_len) + n * (poly_len) + coeff_id;
            // printf("%lu,%lu,%lu,%lu:%lu\n", coeff_id, num_per_id, c, n, idx_c);
            output[idx_c] = sums_out_n0_0 % p_i;
            idx_c += (n2*crt_count*poly_len);
            // printf("%lu,%lu,%lu,%lu:%lu\n", coeff_id, num_per_id, c, n, idx_c);
            output[idx_c] = sums_out_n0_1 % p_i;
            idx_c += (n2*crt_count*poly_len);
            // printf("%lu,%lu,%lu,%lu:%lu\n", coeff_id, num_per_id, c, n, idx_c);
            output[idx_c] = sums_out_n0_2 % p_i;

            // output n1
            n = 1;
            idx_c = num_per_id * (n1*n2*crt_count*poly_len) + c * (crt_count*poly_len) + n * (poly_len) + coeff_id;
            // printf("%lu,%lu,%lu,%lu:%lu\n", coeff_id, num_per_id, c, n, idx_c);
            output[idx_c] = sums_out_n1_0 % b_i;
            idx_c += (n2*crt_count*poly_len);
            // printf("%lu,%lu,%lu,%lu:%lu\n", coeff_id, num_per_id, c, n, idx_c);
            output[idx_c] = sums_out_n1_1 % b_i;
            idx_c += (n2*crt_count*poly_len);
            // printf("%lu,%lu,%lu,%lu:%lu\n", coeff_id, num_per_id, c, n, idx_c);
            output[idx_c] = sums_out_n1_2 % b_i;
        }
    }
    // printf("%lu,%lu\n", coeff_id, num_per_id);
    // printf("%lu,%lu\n", coeff_id, num_per_id);
}

// // --- Kernel 定义 (修改版) ---
// __global__ void MulQueryDBKernal(
//     uint64_t *output,
//     const uint64_t *reorientedCiphertexts,
//     const uint64_t *database, // 注意：这里传入的是当前分片的数据库指针
//     size_t dim0, 
//     size_t num_per,
//     size_t coeff_offset,      // [新增] 当前分片在全局系数中的起始偏移量
//     const size_t poly_len = 2048,
//     const size_t n0 = 2,
//     const size_t n1 = 3,
//     const size_t n2 = 2,
//     const size_t crt_count = 2
// ) {
//     size_t n1_padded = 4; // n1 rounded up to power of 2

//     // relative_coeff_id: 在当前 GPU 显存 Buffer 中的相对索引 (0 ~ chunk_size-1)
//     size_t relative_coeff_id = threadIdx.x + blockIdx.x * blockDim.x;
    
//     // global_coeff_id: 在整个多项式中的真实逻辑索引 (0 ~ 2047)
//     // 用于定位 Output 和 ReorientedCiphertexts
//     size_t global_coeff_id = relative_coeff_id + coeff_offset; 

//     size_t num_per_id = threadIdx.y + blockIdx.y * blockDim.y;

//     // 边界检查：确保没有越界
//     // 注意这里我们检查的是 global_coeff_id 是否超过 poly_len
//     if (global_coeff_id < poly_len && num_per_id < num_per) {
        
//         // 1. 计算 ReorientedCiphertexts 的索引 (使用全局索引，因为它是全量传进去的)
//         size_t idx_a_base = global_coeff_id * (2 * dim0 * n1_padded);
        
//         // 2. 计算 Database 的索引 (关键修改点!)
//         // 显存里的 'database' 指针指向的是当前这一小块(Shard)的起始位置
//         // 所以这里必须使用 relative_coeff_id，否则会越界访问到未分配的显存
//         size_t idx_b_base = relative_coeff_id * (num_per * n2 * dim0 * n0);

//         const uint32_t p_i = 268369921;
//         const uint32_t b_i = 249561089;

//         for (size_t c = 0; c < n2; c++) {
//             __uint128_t sums_out_n0_0 = 0, sums_out_n0_1 = 0, sums_out_n0_2 = 0;
//             __uint128_t sums_out_n1_0 = 0, sums_out_n1_1 = 0, sums_out_n1_2 = 0;

//             for (size_t jm = 0; jm < dim0*2; jm++) {
//                 // 读取当前分片内的数据库数据
//                 uint64_t b = database[idx_b_base++]; 
                        
//                 // 读取全量的密文数据
//                 const uint64_t *v_a = &reorientedCiphertexts[idx_a_base + jm * n1_padded];
//                 uint64_t v_a0 = v_a[0];
//                 uint64_t v_a1 = v_a[1];
//                 uint64_t v_a2 = v_a[2];

//                 // --- 保持原有的数学运算逻辑不变 ---
//                 uint32_t b_lo = b;
//                 uint32_t b_hi = b >> 32L;

//                 uint32_t v_a0_lo = v_a0; uint32_t v_a0_hi = v_a0 >> 32L;
//                 uint32_t v_a1_lo = v_a1; uint32_t v_a1_hi = v_a1 >> 32L;
//                 uint32_t v_a2_lo = v_a2; uint32_t v_a2_hi = v_a2 >> 32L;
                        
//                 // do n0
//                 sums_out_n0_0 += (v_a0_lo) * (uint64_t)b_lo;
//                 sums_out_n0_1 += (v_a1_lo) * (uint64_t)b_lo;
//                 sums_out_n0_2 += (v_a2_lo) * (uint64_t)b_lo;

//                 // do n1
//                 sums_out_n1_0 += ((uint64_t)v_a0_hi) * b_hi;
//                 sums_out_n1_1 += ((uint64_t)v_a1_hi) * b_hi;
//                 sums_out_n1_2 += ((uint64_t)v_a2_hi) * b_hi;
//             }

//             // --- 输出结果 ---
//             // 写入 Output 时，必须使用 global_coeff_id，因为 Output 也是全量分配的
//             size_t n = 0;
//             size_t idx_c = num_per_id * (n1*n2*crt_count*poly_len) + c * (crt_count*poly_len) + n * (poly_len) + global_coeff_id;
            
//             output[idx_c] = sums_out_n0_0 % p_i;
//             idx_c += (n2*crt_count*poly_len);
//             output[idx_c] = sums_out_n0_1 % p_i;
//             idx_c += (n2*crt_count*poly_len);
//             output[idx_c] = sums_out_n0_2 % p_i;

//             n = 1;
//             idx_c = num_per_id * (n1*n2*crt_count*poly_len) + c * (crt_count*poly_len) + n * (poly_len) + global_coeff_id;
            
//             output[idx_c] = sums_out_n1_0 % b_i;
//             idx_c += (n2*crt_count*poly_len);
//             output[idx_c] = sums_out_n1_1 % b_i;
//             idx_c += (n2*crt_count*poly_len);
//             output[idx_c] = sums_out_n1_2 % b_i;
//         }
//     }
// }

void multiplyQueryWithDatabaseGPU(
    uint64_t *output,
    const uint64_t *reorientedCiphertexts,
    const uint64_t *database,
    size_t dim0, 
    size_t num_per,
    const size_t poly_len,
    const size_t n0,
    const size_t n1,
    const size_t n2,
    const size_t crt_count
) {
    size_t n1_padded = 4; // n1 rounded up to power of 2
    size_t nx = poly_len;
    size_t ny = num_per;
    size_t L_small_size_bytes = dim0 * n1 * n0 * crt_count * poly_len * sizeof(uint64_t);
    size_t output_size = sizeof(uint64_t) * num_per * n1 * n2 * 2 * poly_len;
    size_t reC_size = L_small_size_bytes * n1_padded / n1;
    size_t DB_size = sizeof(uint64_t) * dim0 * num_per * n0 * n2 * poly_len;
    
    uint64_t *ipDeviceOutput;
    uint64_t *ipDeviceReorientedCiphertexts;
    uint64_t *ipDeviceDB;
    CUDA_CHECK(cudaMalloc((uint64_t **)&ipDeviceOutput, output_size));
    CUDA_CHECK(cudaMalloc((uint64_t **)&ipDeviceReorientedCiphertexts, reC_size));
    CUDA_CHECK(cudaMalloc((uint64_t **)&ipDeviceDB, DB_size));
    if (ipDeviceOutput != NULL && ipDeviceReorientedCiphertexts != NULL && ipDeviceDB != NULL) {
        CUDA_CHECK(cudaMemcpy(ipDeviceOutput, output, output_size, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(ipDeviceReorientedCiphertexts, reorientedCiphertexts, reC_size, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(ipDeviceDB, database, DB_size, cudaMemcpyHostToDevice));
    } else {
        printf("Fail to allocate device memory\n");
        exit(1);
    }

    dim3 block(8, 4);
    dim3 grid((nx + block.x - 1) / block.x, (ny + block.y - 1) / block.y);
    // printf("thread config:grid:<%d, %d>, block:<%d, %d>\n", grid.x, grid.y, block.x, block.y);
    Timer clock;
    clock.start_gpu();
    MulQueryDBKernal<<<grid, block>>>(
        ipDeviceOutput,
        ipDeviceReorientedCiphertexts, 
        ipDeviceDB,
        dim0, num_per
    );
    CUDA_CHECK(cudaMemcpy(output, ipDeviceOutput, output_size, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaDeviceSynchronize());
    // CUDA_KERNEL_CHECK();
    
    clock.stop_gpu();
    clock.duration_gpu("time_mul_db_with_query");
    time_mul_DB_with_Q = clock._timeElasped;

    CUDA_CHECK(cudaFree(ipDeviceDB));
    CUDA_CHECK(cudaFree(ipDeviceOutput));
    CUDA_CHECK(cudaFree(ipDeviceReorientedCiphertexts));

    // CUDA_CHECK(cudaDeviceReset());
}

// // --- Host 函数 (修改版) ---
// void multiplyQueryWithDatabaseGPU(
//     uint64_t *output,
//     const uint64_t *reorientedCiphertexts,
//     const uint64_t *database,
//     size_t dim0, 
//     size_t num_per,
//     const size_t poly_len,
//     const size_t n0,
//     const size_t n1,
//     const size_t n2,
//     const size_t crt_count
// ) {
//     size_t n1_padded = 4;
    
//     // 1. 计算总内存需求
//     size_t L_small_size_bytes = dim0 * n1 * n0 * crt_count * poly_len * sizeof(uint64_t);
//     size_t output_size = sizeof(uint64_t) * num_per * n1 * n2 * 2 * poly_len;
//     size_t reC_size = L_small_size_bytes * n1_padded / n1;
//     size_t DB_total_size = sizeof(uint64_t) * dim0 * num_per * n0 * n2 * poly_len;

//     printf("\n====== [GPU Memory Analysis] ======\n");
//     printf("Expected Total DB Size:  %10.4f GB\n", (double)DB_total_size / (1024.0*1024.0*1024.0));
//     printf("Reoriented Ciphertexts:  %10.4f GB\n", (double)reC_size / (1024.0*1024.0*1024.0));
//     printf("Output Buffer Size:      %10.4f MB\n", (double)output_size / (1024.0*1024.0));
//     printf("===================================\n");

//     // 2. 分段策略配置 (解决 OOM 的核心)
//     // RTX 3090 (24GB) 无法一次放下 64GB 数据库。
//     // 我们将 2048 个系数切分为 4 份，每份 512 个系数。
//     // 每份 DB 大小约为 16GB，加上其他变量，正好能塞进显存。
//     int num_shards = 4; 
//     size_t coeffs_per_shard = poly_len / num_shards; // 例如 512
//     size_t db_shard_size = DB_total_size / num_shards; // 单个分片的大小

//     printf("[GPU Setup] Total DB: %.2f GB. Splitting into %d shards (%.2f GB each).\n", 
//            (double)DB_total_size/(1024*1024*1024), num_shards, (double)db_shard_size/(1024*1024*1024));

//     // 3. 显存分配
//     uint64_t *ipDeviceOutput;
//     uint64_t *ipDeviceReorientedCiphertexts;
//     uint64_t *ipDeviceDB_Shard; // 注意：这里只需要申请一份分片的大小

//     CUDA_CHECK(cudaMalloc((uint64_t **)&ipDeviceOutput, output_size));
//     CUDA_CHECK(cudaMalloc((uint64_t **)&ipDeviceReorientedCiphertexts, reC_size));
    
//     // 尝试分配 1/4 大小的数据库缓存区
//     cudaError_t err = cudaMalloc((uint64_t **)&ipDeviceDB_Shard, db_shard_size);
//     if (err != cudaSuccess) {
//         printf("ERROR: Failed to allocate DB Shard memory (%.2f GB)!\n", (double)db_shard_size/(1024*1024*1024));
//         exit(1);
//     }

//     // 4. 拷贝不变的数据 (Output 和 Ciphertexts 比较小，直接全量拷贝)
//     CUDA_CHECK(cudaMemcpy(ipDeviceOutput, output, output_size, cudaMemcpyHostToDevice));
//     CUDA_CHECK(cudaMemcpy(ipDeviceReorientedCiphertexts, reorientedCiphertexts, reC_size, cudaMemcpyHostToDevice));

//     // 5. 配置 Kernel 线程块
//     dim3 block(8, 4);
//     // Grid.x 只需要覆盖 coeffs_per_shard (512)，而不是 poly_len (2048)
//     dim3 grid((coeffs_per_shard + block.x - 1) / block.x, (num_per + block.y - 1) / block.y);

//     Timer clock;
//     clock.start_gpu();

//     // 6. 循环分片处理
//     for (int s = 0; s < num_shards; s++) {
//         size_t current_coeff_offset = s * coeffs_per_shard;
        
//         // 计算当前分片在 Host 端数据库数组中的字节偏移量
//         // 公式逻辑：每个系数 (coeff) 对应的数据块大小 * 之前的系数个数
//         // 原公式: idx_b_base = coeff_id * (num_per * n2 * dim0 * n0);
//         // 所以 Host 偏移就是: current_coeff_offset * (num_per * n2 * dim0 * n0) * sizeof(uint64_t)
//         size_t host_byte_offset = s * db_shard_size;

//         // A. 拷贝当前这部分数据库 (Host -> Device)
//         // 注意：这里 database 指针做了偏移
//         CUDA_CHECK(cudaMemcpy(ipDeviceDB_Shard, 
//                               (uint8_t*)database + host_byte_offset, 
//                               db_shard_size, 
//                               cudaMemcpyHostToDevice));

//         // B. 启动 Kernel
//         // 传入 ipDeviceDB_Shard (它只包含当前的 1/4 数据)
//         // 传入 current_coeff_offset (告诉 Kernel 这是第几段)
//         MulQueryDBKernal<<<grid, block>>>(
//             ipDeviceOutput,
//             ipDeviceReorientedCiphertexts, 
//             ipDeviceDB_Shard,
//             dim0, num_per,
//             current_coeff_offset // [Critical] 偏移量
//         );
        
//         // 可选：在这里加同步，方便调试。生产环境去掉以允许 CPU 准备下一块数据。
//         // CUDA_CHECK(cudaDeviceSynchronize());
//     }

//     // 7. 拷贝结果回 Host
//     CUDA_CHECK(cudaMemcpy(output, ipDeviceOutput, output_size, cudaMemcpyDeviceToHost));
    
//     // 确保所有任务完成
//     CUDA_CHECK(cudaDeviceSynchronize());
    
//     clock.stop_gpu();
//     clock.duration_gpu("time_mul_db_with_query");
//     time_mul_DB_with_Q = clock._timeElasped;

//     // 8. 释放显存
//     CUDA_CHECK(cudaFree(ipDeviceDB_Shard));
//     CUDA_CHECK(cudaFree(ipDeviceOutput));
//     CUDA_CHECK(cudaFree(ipDeviceReorientedCiphertexts));
// }

void finishGPU() {
    CUDA_CHECK(cudaDeviceReset());
}


