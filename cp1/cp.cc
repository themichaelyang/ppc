#include <cmath>

float correlation(int row_i, int row_j, const float *data, int row_width);
float get_data(int x, int y, const float *data, int nx);
float get_deno_term(int n, float sq_sum, float sum);

/*
This is the function you need to implement. Quick reference:
- input rows: 0 <= y < ny
- input columns: 0 <= x < nx
- element at row y and column x is stored in data[x + y*nx]
- correlation between rows i and row j has to be stored in result[i + j*ny]
- only parts with 0 <= j <= i < ny need to be filled
*/
void correlate(int ny, int nx, const float *data, float *result) {
  int rows = ny;
  int row_width = nx;

  for (int i=0; i < rows; i++) {
    for (int j=0; j < rows; j++) {
      result[i + j*rows] = correlation(i, j, data, row_width);
    }
  }
}

/* https://en.wikipedia.org/wiki/Pearson_correlation_coefficient#For_a_sample
 *
 * Note that Wikipedia uses x and y as the two samples, whereas x and y in our
 * code relates to the row/col coordinates. Our samples are the rows.
 * */
float correlation(int row_i, int row_j, const float *data, int row_width) {
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

  return numerator / (deno_row_i_term * deno_row_j_term);
}

float get_data(int x, int y, const float *data, int nx) {
  return data[x + y*nx];
}

float get_deno_term(int n, float sq_sum, float sum) {
  float sum_squared = sum * sum;
  return std::sqrt((n * sq_sum) - sum_squared);
}
