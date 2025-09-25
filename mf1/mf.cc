#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <iostream>

float median_pixel(int x, int y, int nx, int ny, int hx, int hy, const float *in);

/*
This is the function you need to implement. Quick reference:
- input rows: 0 <= y < ny
- input columns: 0 <= x < nx
- element at row y and column x is stored in in[x + y*nx]
- for each pixel (x, y), store the median of the pixels (a, b) which satisfy
  max(x-hx, 0) <= a < min(x+hx+1, nx), max(y-hy, 0) <= b < min(y+hy+1, ny)
  in out[x + y*nx].
*/
void mf(int ny, int nx, int hy, int hx, const float *in, float *out) {
  for (int y = 0; y < ny; y++) {
    for (int x = 0; x < nx; x++) {
      out[x + y * nx] = median_pixel(x, y, nx, ny, hx, hy, in);
    }
  }
}

float median_pixel(int x, int y, int nx, int ny, int hx, int hy, const float *in) {
  int x_low = std::max(x - hx, 0);
  int x_hi = std::min(x + hx, nx - 1);

  int y_low = std::max(y - hy, 0);
  int y_hi = std::min(y + hy, ny - 1);

  int length = (y_hi - y_low + 1) * (x_hi - x_low + 1);
  float* pixels = (float*) calloc(length, sizeof(float));
  int i = 0;

  for (int yi = y_low; yi <= y_hi; yi++) {
    int cols = yi * nx;
    for (int xi = x_low; xi <= x_hi; xi++) {
      pixels[i] = in[xi + cols];
      i++;
    }
  }

  int halfway = length/2;
  std::nth_element(pixels, pixels + halfway, pixels + length);
  float median;

  if (length % 2 == 0) {
    std::nth_element(pixels, pixels + halfway - 1, pixels + length);
    median = (pixels[halfway] + pixels[halfway - 1]) / 2;
  }
  else {
    median = pixels[halfway];
  }
  free(pixels);

  return median;
}
