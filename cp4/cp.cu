/* From https://ppc.cs.aalto.fi/ch4/v0/ */

// A reasonable way to calculate all pairwise correlations is the following:

// First normalize the input rows so that each row has the arithmetic mean of 0 — 
// be careful to do the normalization so that you do not change pairwise correlations.
//
// Then normalize the input rows so that for each row the sum of the squares of the elements is 1 — 
// again, be careful to do the normalization so that you do not change pairwise correlations.
//
// Let X be the normalized input matrix.
// Calculate the (upper triangle of the) matrix product Y = XX^T.
//
// Now matrix Y contains all pairwise correlations. The only computationally-intensive part is the
// computation of the matrix product; the normalizations can be done in linear time in the input size.

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
__host__ __device__ float get_data(int x, int y, const float *data, int nx);
__host__ __device__ int get_index(int x, int y, int nx);

/*
This is the function you need to implement. Quick reference:
- input rows: 0 <= y < ny
- input columns (row_width): 0 <= x < nx
- element at row y and column x is stored in data[x + y*nx]
- correlation between rows i and row j has to be stored in result[i + j*ny]
- only parts with 0 <= j <= i < ny need to be filled
*/

// TODO: look into Godbolt
void correlate(int ny, int nx, const float *data, float *result) {
  int rows = ny;
  int row_width = nx;

  // float, to match the size and data expected by GPU. 
  // calloc instead of malloc to initialize memory to pass compute-sanitizer.
  float *normalized = (float *)calloc(nx * ny, sizeof(float));

  for (int row=0; row < rows; row++) {
    float sum = 0;

    // get the mean
    for (int element=0; element < row_width; element++) {
      sum += get_data(element, row, data, row_width);
    }
    double mean = sum / row_width;

    float squared_sum = 0;
    // zero the mean
    for (int element=0; element < row_width; element++) {
      double zero_meaned = get_data(element, row, data, row_width) - mean;
      normalized[get_index(element, row, row_width)] = zero_meaned;
      squared_sum += (zero_meaned * zero_meaned);
    }

    squared_sum = std::sqrt(squared_sum);
    for (int element=0; element < row_width; element++) {
      normalized[get_index(element, row, row_width)] /= squared_sum;
    }
  }

  // Allocate memory & copy data to GPU
  float* data_GPU = NULL;
  float* result_GPU = NULL;

  CHECK(cudaMalloc((void**)&data_GPU, rows * row_width * sizeof(float)));
  CHECK(cudaMalloc((void**)&result_GPU, rows * rows * sizeof(float)));

  // Need to initialize memory to avoid angering compute-sanitizer
  CHECK(cudaMemset(result_GPU, 0, ny * ny * sizeof(float)));
  CHECK(cudaMemcpy(data_GPU, normalized, rows * row_width * sizeof(float), cudaMemcpyHostToDevice));

  // Run kernel
  dim3 dimBlock(16, 16); // dimensions of a block
  dim3 dimGrid(divup(ny, dimBlock.x), divup(ny, dimBlock.y));

  correlation_kernel<<<dimGrid, dimBlock>>>(result_GPU, data_GPU, row_width, rows);
  CHECK(cudaGetLastError());

  // Copy data back to CPU & release memory
  CHECK(cudaMemcpy(result, result_GPU, rows * rows * sizeof(float), cudaMemcpyDeviceToHost));
  CHECK(cudaFree(data_GPU));
  CHECK(cudaFree(result_GPU));

  free(normalized);
}

__global__ void correlation_kernel(float *result, const float *normalized_data, int row_width, int rows) {
  int row_i = threadIdx.x + blockIdx.x * blockDim.x;
  int row_j = threadIdx.y + blockIdx.y * blockDim.y;

  if (row_i >= rows || row_j >= rows || row_i < row_j)
    return;

  float correlation = 0;
  for (int element=0; element < row_width; element++) {
    correlation += (get_data(element, row_i, normalized_data, row_width) * get_data(element, row_j, normalized_data, row_width));
  }

  // result is sized row x row, since row-wise correlations
  result[get_index(row_i, row_j, rows)] = correlation;
}

__host__ __device__ float get_data(int x, int y, const float *data, int nx) {
  return data[x + y*nx];
}

__host__ __device__ int get_index(int x, int y, int nx) {
  return x + y*nx;
}
