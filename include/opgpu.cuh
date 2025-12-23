#include <stdio.h>
#include <cuda_runtime.h>
#include <system_error>
#include <stdarg.h>
#include <chrono>
extern float time_mul_DB_with_Q;

void print_hello_from_gpu();

void setGPU();

cudaError_t ErrorCheck(
    cudaError_t error_code,
    const char* filename,
    int lineNumber
);

#define CUDA_CHECK(call)             __cudaCheck(call, __FILE__, __LINE__)
#define CUDA_KERNEL_CHECK()          __kernelCheck(__FILE__, __LINE__)
#define LOG(...)                     __log_info(__VA_ARGS__)

inline static void __cudaCheck(cudaError_t err, const char* file, const int line) 
{
    if (err != cudaSuccess) 
    {
        printf("ERROR: %s:%d, ", file, line);
        printf("CODE:%s, DETAIL:%s\n", cudaGetErrorName(err), cudaGetErrorString(err));
        exit(1);
    }
}

inline static void __kernelCheck(const char* file, const int line) 
{
    cudaError_t err = cudaPeekAtLastError();
    if (err != cudaSuccess) 
    {
        printf("ERROR: %s:%d, ", file, line);
        printf("CODE:%s, DETAIL:%s\n", cudaGetErrorName(err), cudaGetErrorString(err));
        exit(1);
    }
}

static void __log_info(const char* format, ...) 
{
    char msg[1000];
    va_list args;
    va_start(args, format);

    vsnprintf(msg, sizeof(msg), format, args);

    fprintf(stdout, "%s\n", msg);
    va_end(args);
}

class Timer {
public:
    using s  = std::ratio<1, 1>;
    using ms = std::ratio<1, 1000>;
    using us = std::ratio<1, 1000000>;
    using ns = std::ratio<1, 1000000000>;

public:
    Timer();
    ~Timer();

public:
    void start_cpu();
    void start_gpu();
    void stop_cpu();
    void stop_gpu();

    template <typename span>
    void duration_cpu(std::string msg);

    void duration_gpu(std::string msg);

private:
    std::chrono::time_point<std::chrono::high_resolution_clock> _cStart;
    std::chrono::time_point<std::chrono::high_resolution_clock> _cStop;
    cudaEvent_t _gStart;
    cudaEvent_t _gStop;
public:
    float _timeElasped;
};

template <typename span>
void Timer::duration_cpu(std::string msg){
    std::string str;

    if(std::is_same<span, s>::value) { str = "s"; }
    else if(std::is_same<span, ms>::value) { str = "ms"; }
    else if(std::is_same<span, us>::value) { str = "us"; }
    else if(std::is_same<span, ns>::value) { str = "ns"; }

    std::chrono::duration<double, span> time = _cStop - _cStart;
    _timeElasped = time.count();
    LOG("%-40s uses %.6lf %s", msg.c_str(), time.count(), str.c_str());
}

__global__ void addMatrix(int *out, int *in1, int *in2, const int nx, const int ny);

void test_MatAdd();

void multiplyQueryWithDatabaseGPU(
    uint64_t *output,
    const uint64_t *reorientedCiphertexts,
    const uint64_t *database,
    size_t dim0, 
    size_t num_per,
    const size_t poly_len = 2048,
    const size_t n0 = 2,
    const size_t n1 = 3,
    const size_t n2 = 2,
    const size_t crt_count =2
);

void finishGPU();