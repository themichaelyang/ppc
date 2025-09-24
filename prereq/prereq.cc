struct Result {
    float avg[3];
};

/*
This is the function you need to implement. Quick reference:
- x coordinates: 0 <= x < nx
- y coordinates: 0 <= y < ny
- horizontal position: 0 <= x0 < x1 <= nx
- vertical position: 0 <= y0 < y1 <= ny
- color components: 0 <= c < 3
- input: data[c + 3 * x + 3 * nx * y]
- output: avg[c]
*/
Result calculate(int ny, int nx, const float *data, int y0, int x0, int y1, int x1) {
  double totals[3] = {0.0d, 0.0d, 0.0d};
  Result result{{0.0f, 0.0f, 0.0f}};
  int size = (x1 - x0) * (y1 - y0);

  for (int x = x0; x < x1; x++) {
    int row_offset = 3*x;
    for (int y = y0; y < y1; y++) {
      int column_offset = 3*nx*y;
      for (int c = 0; c < 3; c++) {
        totals[c] += data[c + row_offset + column_offset];
      }
    }
  }

  for (int c = 0; c < 3; c++) {
    result.avg[c] = totals[c] / size;
  }

  return result;
}
