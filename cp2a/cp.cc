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
#include <cmath>

double correlate_x(int x, int row_i, int row_j, double *normalized, int row_width);
double get_data(int x, int y, double *data, int nx);
float get_data(int x, int y, const float *data, int nx);
int get_index(int x, int y, int nx);

const int parallel_instr = 4;

/*
This is the function you need to implement. Quick reference:
- input rows: 0 <= y < ny
- input columns (row_width): 0 <= x < nx
- element at row y and column x is stored in data[x + y*nx]
- correlation between rows i and row j has to be stored in result[i + j*ny]
- only parts with 0 <= j <= i < ny need to be filled
*/
void correlate(int ny, int nx, const float *data, float *result) {
  int rows = ny;
  int row_width = nx;
  int normalized_row_width = row_width;

  // make row width divisible by # parallel instructions wanted, reduce branching in inner loop
  if (row_width % parallel_instr != 0) {
    normalized_row_width = row_width + (parallel_instr - row_width % parallel_instr);
  }

  double* normalized = (double*) calloc(normalized_row_width * ny, sizeof(double));

  for (int row=0; row < rows; row++) {
    double sum = 0;

    // get the mean
    for (int element=0; element < row_width; element++) {
      sum += get_data(element, row, data, row_width);
    }
    double mean = sum / row_width;

    // zero the mean, calculate squared sum
    double squared_sum = 0;
    for (int element=0; element < row_width; element++) {
      double zero_meaned = get_data(element, row, data, row_width) - mean;
      normalized[get_index(element, row, normalized_row_width)] = zero_meaned;
      squared_sum += (zero_meaned * zero_meaned);
    }

    // one the row magnitude
    double magnitude = std::sqrt(squared_sum);
    for (int element=0; element < row_width; element++) {
      normalized[get_index(element, row, normalized_row_width)] /= magnitude;
    }
  }

  double* partials = (double*) calloc(normalized_row_width, sizeof(double));

  // only populate upper triangle
  for (int row_j=0; row_j < rows; row_j++) {
    for (int row_i=row_j; row_i < rows; row_i++) {
      double correlation = 0;
      asm("# inner loop start");
      for (int offset=0; offset < parallel_instr; offset++) {
        partials[offset] = 0;
      }

      // parallelize instructions by reducing data dependencies
      for (int s=0; s < normalized_row_width / parallel_instr; s++) {
        for (int offset=0; offset < parallel_instr; offset++) {
          partials[offset] += correlate_x(s * parallel_instr + offset, row_i, row_j, normalized, normalized_row_width);
        }
      }

      for (int offset=0; offset < parallel_instr; offset++) {
        correlation += partials[offset];
      }
      result[get_index(row_i, row_j, rows)] = correlation;
      asm("# inner loop end");
    }
  }

  free(normalized);
  free(partials);
}

double correlate_x(int x, int row_i, int row_j, double *normalized, int row_width) {
  return (get_data(x, row_i, normalized, row_width) * get_data(x, row_j, normalized, row_width));
}

double get_data(int x, int y, double *data, int nx) {
  return data[x + y*nx];
}

float get_data(int x, int y, const float *data, int nx) {
  return data[x + y*nx];
}

int get_index(int x, int y, int nx) {
  return x + y*nx;
}
