#include <vector>
#include <cmath>

double* make_summed_area(int height, int width, const float *data);
double get_summed_area(double* areas, int x, int y, int c, int width, int height, int nx);
int index(int c, int x, int y, int nx);
double sum(double* arr, int len);

struct Result {
  int y0;
  int x0;
  int y1;
  int x1;
  float outer[3];
  float inner[3];
};

/*
This is the function you need to implement. Quick reference:
- x coordinates: 0 <= x < nx
- y coordinates: 0 <= y < ny
- color components: 0 <= c < 3
- input: data[c + 3 * x + 3 * nx * y]
*/

// To min the cost of the inner and outer individually, set the averages.
// Then, to find the inner and outer pair with the lowest cost:
//
Result segment(int ny, int nx, const float *data) {
  Result result{0, 0, 0, 0, {0, 0, 0}, {0, 0, 0}};

  double* areas = make_summed_area(ny, nx, data);

  int total_area = nx * ny;
  double total_sum[3] = {0.0, 0.0, 0.0};
  for (int c = 0; c < 3; c++) {
    total_sum[c] = get_summed_area(areas, 0, 0, c, nx - 1, ny - 1, nx);
  }

  double max = 0;

  for (int width = 1; width < nx - x; width++) {
    for (int height = 1; height < ny - y; height++) {
      for (int y = 0; y < ny; y++) {
        for (int x = 0; x < nx; x++) {

          int inner_area = width * height;
          double inner_sum[3] = {0.0, 0.0, 0.0};
          double inner_avg[3] = {0.0, 0.0, 0.0};

          int outer_area = total_area - inner_area;
          double outer_sum[3] = {0.0, 0.0, 0.0};
          double outer_avg[3] = {0.0, 0.0, 0.0};

          for (int c = 0; c < 3; c++) {
            inner_sum[c] = get_summed_area(areas, x, y, c, width, height, nx);
            outer_sum[c] = total_sum[c] - inner_sum[c];
          }

          for (int c = 0; c < 3; c++) {
            inner_avg[c] = inner_sum[c] / inner_area;
            outer_avg[c] = outer_sum[c] / outer_area;
          }

          double diff = std::abs(inner_avg - outer_avg);

          if (diff > max) {
            int x1 = x + width;
            int y1 = y + height;
            result.x0 = x;
            result.y0 = y;
            result.x1 = x1;
            result.y1 = y1;

            for (int c = 0; c < 3; c++) {
              result.outer[c] = outer_avg[c];
              result.inner[c] = inner_avg[c];
            }
          }
        }
      }
    }
  }

  free(areas);

  return result;
}

// essentially a 2d prefix sum
double* make_summed_area(int ny, int nx, const float *data) {
  double* areas = (double*) calloc(nx * ny * 3, sizeof(double));

  for (int y = 0; y < ny; y++) {
    for (int x = 0; x < nx; x++) {
      for (int c = 0; c < 3; c++) {
        int i = index(c, x, y, nx);
        areas[i] = data[index(c, x, y, nx)];

        if (x - 1 >= 0) {
          areas[i] += areas[index(c, x-1, y, nx)];
        }

        if (y - 1 >= 0) {
          areas[i] += areas[index(c, x, y-1, nx)];
        }

        if (x - 1 >= 0 && y - 1 >= 0) {
          areas[i] += areas[index(c, x-1, y-1, nx)];
        }
      }
    }
  }

  return areas;
}

double get_summed_area(double* areas, int x, int y, int c, int width, int height, int nx) {
  double sum = areas[index(c, x + width, y + height, nx)];
  sum -= areas[index(c, x, y, nx)];
  sum -= areas[index(c, x, y + height, nx)];
  sum += areas[index(c, x, y, nx)];
  return sum;
}

int index(int c, int x, int y, int nx) {
  return c + 3*x + 3*y*nx;
}

double sum(double* arr, int len) {
  double ret = 0;

  for (int i = 0; i < len; i++) {
    ret += arr[i];
  }

  return ret;
}

