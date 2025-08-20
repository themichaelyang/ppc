/* From https://ppc.cs.aalto.fi/ch4/v0/ */
#include <cstdlib>
#include <iostream>
#include <cuda_runtime.h>

static inline void check(cudaError_t err, const char* context) {
    if (err != cudaSuccess) {
        std::cerr << "CUDA error: " << context << ": "
            << cudaGetErrorString(err) << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

#define CHECK(x) check(x, #x)

static inline int divup(int a, int b) {
    return (a + b - 1)/b;
}

/* static inline int roundup(int a, int b) { */
/*     return divup(a, b) * b; */
/* } */

__global__ void correlation_kernel(float *result, const float *data, int row_width, int rows);
__device__ double get_data(int x, int y, const float *data, int nx);
__device__ double get_deno_term(int n, double sq_sum, double sum);

/*
This is the function you need to implement. Quick reference:
- input rows: 0 <= y < ny
- input columns (row_width): 0 <= x < nx
- element at row y and column x is stored in data[x + y*nx]
- correlation between rows i and row j has to be stored in result[i + j*ny]
- only parts with 0 <= j <= i < ny need to be filled
*/
void correlate(int rows, int row_width, const float *data, float *result) {
  // Allocate memory & copy data to GPU
  float* data_GPU = NULL;
  float* result_GPU = NULL;

  CHECK(cudaMalloc((void**)&data_GPU, rows * row_width * sizeof(float)));
  CHECK(cudaMalloc((void**)&result_GPU, rows * rows * sizeof(float)));
  CHECK(cudaMemcpy(data_GPU, data, rows * row_width * sizeof(float), cudaMemcpyHostToDevice));

  // Run kernel
  // block x, thread x
  correlation_kernel<<<rows, rows>>>(result_GPU, data_GPU, row_width, rows);
  CHECK(cudaGetLastError());

  // Copy data back to CPU & release memory
  CHECK(cudaMemcpy(result, result_GPU, rows * rows * sizeof(float), cudaMemcpyDeviceToHost));
  CHECK(cudaFree(data_GPU));
  CHECK(cudaFree(result_GPU));
}

__global__ void correlation_kernel(float *result, const float *data, int row_width, int rows) {
  int row_i = blockIdx.x;
  int row_j = threadIdx.x;

  if (row_i >= rows || row_j >= rows || row_i < row_j)
    return;

  float row_dot = 0;
  float row_i_sum = 0;
  float row_j_sum = 0;

  float row_i_sq_sum = 0;
  float row_j_sq_sum = 0;

  for (int col=0; col < row_width; col++) {
    float row_i_col = get_data(col, row_i, data, row_width);
    float row_j_col = get_data(col, row_j, data, row_width);

    row_dot += row_i_col * row_j_col;

    row_i_sum += row_i_col;
    row_j_sum += row_j_col;

    row_i_sq_sum += row_i_col * row_i_col;
    row_j_sq_sum += row_j_col * row_j_col;
  }

  float numerator = (row_width * row_dot) - (row_i_sum * row_j_sum);
  float deno_row_i_term = get_deno_term(row_width, row_i_sq_sum, row_i_sum);
  float deno_row_j_term = get_deno_term(row_width, row_j_sq_sum, row_j_sum);

  result[row_i + row_j*rows] = numerator / (deno_row_i_term * deno_row_j_term);
}

__device__ double get_data(int x, int y, const float *data, int nx) {
  return data[x + y*nx];
}

__device__ double get_deno_term(int n, double sq_sum, double sum) {
  double sum_squared = sum * sum;
  return std::sqrt((n * sq_sum) - sum_squared);
}
