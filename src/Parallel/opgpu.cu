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

void finishGPU() {
    CUDA_CHECK(cudaDeviceReset());
}


