// #include <cmath>
#include <iomanip>
#include <iostream>
#include <ostream>

// precompute color sums for all rectangles anchored at (0, 0) using a single-pass
// summed area table algorithm, essentially a 2d prefix sum

// precompute color squared sum in summed area table: (a_i)^2

// iterate through all combinations of rectangles for each color channel and keep track of minimum SSE
// SSE = sum_i((p - a_i)^2) = n*p^2 - 2p*sum_i(a_i) + sum_i(a_i^2)
// where p is predicted and a_i is color at pixel i

// do this for both inner and outer rectangles since they have different p
// then add together (can break up sum because of associativity)

// use precomputed summed area tables to O(1) calculate terms of SSE for every x, y, width, height
// note: inclusive upper left, exclusive bottom right. in general be careful of boundaries.

struct precomp_values {
  double *areas;
  double *component_squared_areas;
};

size_t index(int c, int x, int y, int nx) { return c + 3 * x + 3 * nx * y; }

precomp_values *make_summed_area(int ny, int nx, const float *data) {
  double *areas = (double *)calloc(nx * ny * 3, sizeof(double));
  double *component_squared_areas =
      (double *)calloc(nx * ny * 3, sizeof(double));

  precomp_values *pv = (precomp_values *)calloc(1, sizeof(precomp_values));

  for (int y = 0; y < ny; y++) {
    for (int x = 0; x < nx; x++) {
      for (int c = 0; c < 3; c++) {
        int i = index(c, x, y, nx);
        float pixel_color = data[index(c, x, y, nx)];
        areas[i] = pixel_color;
        component_squared_areas[i] = pixel_color * pixel_color;

        if (x - 1 >= 0) {
          areas[i] += areas[index(c, x - 1, y, nx)];
          component_squared_areas[i] +=
              component_squared_areas[index(c, x - 1, y, nx)];
        }

        if (y - 1 >= 0) {
          areas[i] += areas[index(c, x, y - 1, nx)];
          component_squared_areas[i] +=
              component_squared_areas[index(c, x, y - 1, nx)];
        }

        if (x - 1 >= 0 && y - 1 >= 0) {
          areas[i] -= areas[index(c, x - 1, y - 1, nx)];
          component_squared_areas[i] -=
              component_squared_areas[index(c, x - 1, y - 1, nx)];
        }
      }
    }
  }

  pv->areas = areas;
  pv->component_squared_areas = component_squared_areas;
  return pv;
}

struct Result {
  // upper left corner
  int y0;
  int x0;
  // lower right corner (exclusive)
  int y1;
  int x1;
  // color of background
  float outer[3];
  // color of rectangle
  float inner[3];
};

double get_summed_area(double *areas, int x, int y, int c, int width, int height, int nx) {
  int x1 = x + width - 1;
  int y1 = y + height - 1;
  double sum = areas[index(c, x1, y1, nx)];
  if (y > 0)
    sum -= areas[index(c, x1, y - 1, nx)];
  if (x > 0)
    sum -= areas[index(c, x - 1, y1, nx)];
  if (x > 0 && y > 0)
    sum += areas[index(c, x - 1, y - 1, nx)];
  return sum;
}


/*
This is the function you need to implement. Quick reference:
- x coordinates: 0 <= x < nx
- y coordinates: 0 <= y < ny
- color components: 0 <= c < 3
- input: data[c + 3 * x + 3 * nx * y]
*/
Result segment(int ny, int nx, const float *data) {
  Result result{0, 0, 0, 0, {0, 0, 0}, {0, 0, 0}};

  // precomputation
  precomp_values *pv = make_summed_area(ny, nx, data);

  double min_error = 1e10;

  for (int y0 = 0; y0 < ny; ++y0) {
    for (int x0 = 0; x0 < nx; ++x0) {
      for (int y1 = y0 + 1; y1 <= ny; ++y1) {
        for (int x1 = x0 + 1; x1 <= nx; ++x1) {
          int width = x1 - x0;
          int height = y1 - y0;

          Result candidate = {
              .y0 = y0,
              .x0 = x0,
              .y1 = y1,
              .x1 = x1,
              .outer = {0, 0, 0},
              .inner = {0, 0, 0},
          };
          double err = 0.0;

          for (int c = 0; c < 3; ++c) {
            double total = get_summed_area(pv->areas, 0, 0, c, nx, ny, nx);
            double summed_area =
                get_summed_area(pv->areas, x0, y0, c, x1 - x0, y1 - y0, nx);

            double total_squared_components = get_summed_area(
                pv->component_squared_areas, 0, 0, c, nx, ny, nx);
            double summed_component_squared_area = get_summed_area(
                pv->component_squared_areas, x0, y0, c, x1 - x0, y1 - y0, nx);

            double inner_color = summed_area / (width * height);
            double outer_color =
                (total - summed_area) / ((nx * ny) - (width * height));

            candidate.inner[c] = inner_color;
            candidate.outer[c] = outer_color;

            double err_inner = (width * height) * (inner_color * inner_color) -
                               (2 * summed_area * inner_color) +
                               summed_component_squared_area;
            double err_outer =
                ((nx * ny) - (width * height)) * (outer_color * outer_color) -
                (2 * (total - summed_area) * outer_color) +
                (total_squared_components - summed_component_squared_area);

            err += err_inner + err_outer;
          }

          if (err < min_error) {
            result = candidate;
            min_error = err;
          }
        }
      }
    }
  }

  // appease address sanitizer
  free(pv->areas);
  free(pv->component_squared_areas);
  free(pv);

  return result;
}
