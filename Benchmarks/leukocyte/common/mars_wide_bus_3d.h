/************************************************************************************
 *  (c) Copyright 2014-2015 Falcon Computing Solutions, Inc. All rights reserved.
 *
 *  This file contains confidential and proprietary information
 *  of Falcon Computing Solutions, Inc. and is protected under U.S. and
 *  international copyright and other intellectual property laws.
 *
 ************************************************************************************/

/************************************************************************************
 *  (c) Copyright 2014-2015 Falcon Computing Solutions, Inc. All rights reserved.
 *
 *  This file contains confidential and proprietary information
 *  of Falcon Computing Solutions, Inc. and is protected under U.S. and
 *  international copyright and other intellectual property laws.
 *
 ************************************************************************************/

#include<assert.h>
#include<stdlib.h>
#define BUS_WIDTH LARGE_BUS / 8
#define cpp_get_range(tmp, x, y) tmp(x, y)
#define c_get_range(tmp, x, y) apint_get_range(tmp, x, y)
#define cpp_set_range(tmp, x, y, val) tmp(x, y) = val
#define c_set_range(tmp, x, y, val) tmp = apint_set_range(tmp, x, y, val)
#ifdef __cplusplus
#define tmp2(x, y) cpp_get_range(tmp, x, y)
#define tmp3(x, y, val) cpp_set_range(tmp, x, y, val)
#else
#define tmp2(x, y) c_get_range(tmp, x, y)
#define tmp3(x, y, val) c_set_range(tmp, x, y, val)
#endif

#define memcpy_wide_bus(access, type, size1, size2) instant_size(access, type, size1, size2)
#define instant_size(access, type, size1, size2) memcpy_wide_bus_##access##_##type##_3d_##size1##_##size2


static void
memcpy_wide_bus(read, char, SIZE_1, SIZE_2)(char a_buf[][SIZE_1][SIZE_2], 
                            size_t index3_offset, size_t index2_offset,
                            size_t index1_offset, MARS_WIDE_BUS_TYPE *a,
                            size_t offset_byte, size_t size_byte) {
#pragma HLS inline self
  const size_t data_width = sizeof(char);
  const size_t bus_width = BUS_WIDTH;
  const size_t num_elements = bus_width / data_width;
  size_t buf_size = size_byte / data_width;
  size_t offset = offset_byte / data_width;
  size_t head_align = offset & (num_elements - 1);
  size_t new_offset = offset + buf_size;
  size_t tail_align = (new_offset - 1) & (num_elements - 1);
  size_t start = offset / num_elements;
  size_t end = (offset + buf_size + num_elements - 1) / num_elements;
  //MARS_WIDE_BUS_TYPE *a_offset = a + start;
  size_t i, j;
  size_t index3 = index3_offset, index2 = index2_offset, index1 = index1_offset;
  const size_t bound1 = SIZE_2 / num_elements;
  const size_t SIZE = SIZE_1 * SIZE_2;
  const size_t bound2 = SIZE / num_elements;
  const size_t bound0 = num_elements / SIZE;
  int aligned = 0;
  if ((SIZE_2 % num_elements) == 0 && (index1_offset % num_elements) == 0) {
    aligned = 1;
    index1 = index1_offset / num_elements;
  }
  else if ((SIZE % num_elements) == 0 && 
           (num_elements % SIZE_2) == 0 &&
      ((index2_offset * SIZE_2) % num_elements) == 0 &&
      !index1_offset) {
    aligned = 2;
    index2 = index2_offset * SIZE_2 / num_elements;
  }
  else if ((num_elements % SIZE) == 0 && 
           (index3_offset * SIZE % num_elements) == 0 &&
           !index2_offset && !index1_offset) {
    index3 = index3_offset * SIZE  / num_elements;
    aligned = -1;
  }
  const size_t index_offset = index3_offset * SIZE + index2_offset * 
                              SIZE_2 + index1_offset;
  int len = end - start;
  assert(len <= buf_size / num_elements + 2);
  assert(len >= buf_size / num_elements);
  if (1 == len) {
#ifdef __cplusplus
    MARS_WIDE_BUS_TYPE tmp(a[start]);
#else
    MARS_WIDE_BUS_TYPE tmp = a[start];
#endif
    for (j = 0; j < num_elements; ++j) {
       if (j < head_align || j > tail_align)
         continue;
       size_t buffer_index = j - head_align + index_offset;
       a_buf[buffer_index / SIZE][(buffer_index % SIZE) / SIZE_2][(buffer_index % SIZE) % SIZE_2] =
           tmp2(((j + 1) * data_width * 8 - 1), (j * data_width * 8));
    }
    return;
  }
  for (i = 0; i < len; ++i) {
#pragma HLS pipeline
#ifdef __cplusplus
    MARS_WIDE_BUS_TYPE tmp(a[i + start]);
#else
    MARS_WIDE_BUS_TYPE tmp = a[i + start];
#endif
    if (head_align == 0 && aligned) {
      if (aligned == 1) {
        for (j = 0; j < num_elements; ++j) {
          if (i == end - start - 1 && j > tail_align)
            continue;
          size_t buffer_index = i * num_elements + j - 0;
          a_buf[index3][index2][index1 * num_elements + j] =
                tmp2(((j + 1) * data_width * 8 - 1), (j * data_width * 8));
        }
        index1++;
        if (index1 == bound1) {
          index1 = 0;
          ++index2;
        }
        if (index2 == SIZE_1) {
          index2 = 0;
          ++index3;
        }
      } else if (aligned == 2) {
        int j0, j1;
        for (j0 = 0; j0 < num_elements / SIZE_2; ++j0) {
          for (j1 = 0; j1 < SIZE_2; ++j1) {
            j = j0 * SIZE_2 + j1;
            if (i == end - start - 1 && j > tail_align)
              continue;
            a_buf[index3][index2 * num_elements / SIZE_2  + j0][j1] =
                tmp2(((j + 1) * data_width * 8 - 1), (j * data_width * 8));
          }
        }
        ++index2;
        if (index2 == bound2) {
          index2 = 0;
           ++index3;
        }
      } else if (aligned == -1) {
        int j0, j1, j2;
        for (j0 = 0; j0 < num_elements / SIZE; ++j0) {
          for (j1 = 0; j1 < SIZE_1; ++j1) {
            for (j2 = 0; j2 < SIZE_2; ++j2) {
              j = j0 * SIZE + j1 * SIZE_2 + j2;
              if (i == end - start - 1 && j > tail_align)
                continue;
              a_buf[index3 * bound0 + j0][j1][j2] =
                  tmp2(((j + 1) * data_width * 8 - 1), (j * data_width * 8));
            }
          }
        }
        ++index3;
      }
    } else {
      for (j = 0; j < num_elements; ++j) {
        if (i == 0 && j < head_align)
          continue;
        if (i == end - start - 1 && j > tail_align)
          continue;
        size_t buffer_index = i * num_elements + j - head_align + index_offset;
        a_buf[buffer_index / SIZE][(buffer_index % SIZE) / SIZE_2][(buffer_index % SIZE) % SIZE_2] =
            tmp2(((j + 1) * data_width * 8 - 1), (j * data_width * 8));
      }
    }
  }
}

static void
memcpy_wide_bus(write, char, SIZE_1, SIZE_2)(MARS_WIDE_BUS_TYPE *c, char c_buf[][SIZE_1][SIZE_2],
                             size_t index3_offset, size_t index2_offset, 
                             size_t index1_offset,
                             size_t offset_byte, size_t size_byte) {
#pragma HLS inline self
  const size_t data_width = sizeof(char);
  const size_t bus_width = BUS_WIDTH;
  const size_t num_elements = bus_width / data_width;
  size_t buf_size = size_byte / data_width;
  size_t offset = offset_byte / data_width;
  size_t head_align = offset & (num_elements - 1);
  size_t new_offset = offset + buf_size;
  size_t tail_align = (new_offset - 1) & (num_elements - 1);
  size_t start = offset / num_elements;
  size_t end = (offset + buf_size + num_elements - 1) / num_elements;
  size_t len = end - start;
  size_t i, j;
  const size_t SIZE = SIZE_1 * SIZE_2;
  size_t index_offset = SIZE * index3_offset + SIZE_2 * index2_offset +
                        index1_offset;
  if (head_align == 0)
    len = (buf_size + num_elements - 1) / num_elements;
  size_t align = 0;
  if (len == 1) {
    MARS_WIDE_BUS_TYPE tmp;
    if (head_align != 0 || tail_align != (num_elements - 1))
      tmp = c[start];
    for (j = 0; j < num_elements; ++j) {
      if (j < head_align)
        continue;
      if (j > tail_align)
        continue;
      size_t buffer_index = j - head_align + index_offset;
      tmp3(((j + 1) * data_width * 8 - 1), (j * data_width * 8),
           c_buf[buffer_index / SIZE][(buffer_index % SIZE) / SIZE_2][(buffer_index % SIZE) %
                                                           SIZE_2]);
    }
    c[start] = tmp;
    return;
  }
  if (head_align != 0) {
    MARS_WIDE_BUS_TYPE tmp = c[start];
    for (j = 0; j < num_elements; ++j) {
      if (j < head_align)
        continue;
      size_t buffer_index = j - head_align + index_offset;
      tmp3(((j + 1) * data_width * 8 - 1), (j * data_width * 8),
           c_buf[buffer_index / SIZE][(buffer_index % SIZE) / SIZE_2][
                 (buffer_index % SIZE) % SIZE_2]);
    }
    c[start] = tmp;
    start++;
    align++;
  }
  if (tail_align != (num_elements - 1))
    align++;
  size_t index3 = index3_offset, index2 = index2_offset, index1 = index1_offset;
  const size_t bound1 = SIZE_2 / num_elements;
  const size_t bound2 = SIZE / num_elements;
  const size_t bound0 = num_elements / SIZE;
  int aligned = 0;
  if ((SIZE_2 % num_elements) == 0 && (index1_offset % num_elements) == 0) {
    aligned = 1;
    index1 = index1_offset / num_elements;
  }
  else if ((SIZE % num_elements) == 0 && 
           (num_elements % SIZE_2) == 0 &&
      ((index2_offset * SIZE_2) % num_elements) == 0 &&
      !index1_offset) {
    aligned = 2;
    index2 = index2_offset * SIZE_2 / num_elements;
  }
  else if ((num_elements % SIZE) == 0 && 
           (index3_offset * SIZE % num_elements) == 0 &&
           !index2_offset && !index1_offset) {
    index3 = index3_offset * SIZE  / num_elements;
    aligned = -1;
  }
  int burst_len = len - align;
  assert(burst_len <= buf_size / num_elements);
  for (i = 0; i < burst_len; ++i) {
#pragma HLS pipeline
    MARS_WIDE_BUS_TYPE tmp;
    if (head_align == 0 && aligned) {
      if (aligned == 1) {
        for (j = 0; j < num_elements; ++j) {
          int val =  c_buf[index3][index2][index1 * num_elements + j];
          tmp3(((j + 1) * data_width * 8 - 1), (j * data_width * 8), val);
        }
        index1++;
        if (index1 == bound1) {
          index1 = 0;
          ++index2;
        }
        if (index2 == SIZE_1) {
          index2 = 0;
          ++index3;
        }
      } else if (aligned == 2) {
        int j0, j1;
        for (j0 = 0; j0 < num_elements / SIZE_2; ++j0) {
          for (j1 = 0; j1 < SIZE_2; ++j1) {
            j = j0 * SIZE_2 + j1;
            int val = c_buf[index3][index2 * num_elements / SIZE_2 + j0][j1];
            tmp3(((j + 1) * data_width * 8 - 1), (j * data_width * 8), val);
          }
        }
        ++index2;
        if (index2 == bound2) {
          index2 = 0;
          ++index3;
        }
      } else if (aligned == -1) {
        int j0, j1, j2;
        for (j0 = 0; j0 < bound0; ++j0) {
          for (j1 = 0; j1 < SIZE_1; ++j1) {
            for (j2 = 0; j2 < SIZE_2; ++j2) {
              j = j0 * SIZE + j1 * SIZE_2 + j2;
              int val = c_buf[index3 *bound0 + j0][j1][j2];
              tmp3(((j + 1) * data_width * 8 - 1), (j * data_width * 8), val);
            }
          }
        }
        ++index3;
      }
    } else {
      for (j = 0; j < num_elements; ++j) {
        size_t buffer_index = i * num_elements + j + num_elements - head_align + index_offset;
        tmp3(((j + 1) * data_width * 8 - 1), (j * data_width * 8),
             c_buf[buffer_index / SIZE][(buffer_index % SIZE) / SIZE_2][(buffer_index % SIZE) %
                                                             SIZE_2]);
      }
    }
    c[i + start] = tmp;
  }
  if (tail_align != num_elements - 1) {
    MARS_WIDE_BUS_TYPE tmp = c[end - 1];
    size_t pos = (len - align) * num_elements + index_offset;
    pos += (num_elements - head_align) % num_elements;
    for (j = 0; j < num_elements; ++j) {
      if (j > tail_align)
        continue;
      size_t buffer_index = pos + j;
      tmp3(((j + 1) * data_width * 8 - 1), (j * data_width * 8),
           c_buf[buffer_index / SIZE][(buffer_index % SIZE) / SIZE_2][(buffer_index % SIZE) %
                                                           SIZE_2]);
    }
    c[end - 1] = tmp;
  }
}

static void
memcpy_wide_bus(read, int, SIZE_1, SIZE_2)(int a_buf[][SIZE_1][SIZE_2], 
                                                   size_t index3_offset,
                                                   size_t index2_offset,
                                                   size_t index1_offset,
                            MARS_WIDE_BUS_TYPE *a,
                            size_t offset_byte, size_t size_byte) {
#pragma HLS inline self
  const size_t data_width = sizeof(int);
  const size_t bus_width = BUS_WIDTH;
  const size_t num_elements = bus_width / data_width;
  size_t buf_size = size_byte / data_width;
  size_t offset = offset_byte / data_width;
  size_t head_align = offset & (num_elements - 1);
  size_t new_offset = offset + buf_size;
  size_t tail_align = (new_offset - 1) & (num_elements - 1);
  size_t start = offset / num_elements;
  size_t end = (offset + buf_size + num_elements - 1) / num_elements;
  //MARS_WIDE_BUS_TYPE *a_offset = a + start;
  size_t i, j;
  size_t index3 = index3_offset, index2 = index2_offset, index1 = index1_offset;
  const size_t bound1 = SIZE_2 / num_elements;
  const size_t SIZE = SIZE_1 * SIZE_2;
  const size_t bound2 = SIZE / num_elements;
  const size_t bound0 = num_elements / SIZE;
  int aligned = 0;
  if ((SIZE_2 % num_elements) == 0 && (index1_offset % num_elements) == 0) {
    aligned = 1;
    index1 = index1_offset / num_elements;
  }
  else if ((SIZE % num_elements) == 0 && 
           (num_elements % SIZE_2) == 0 &&
      ((index2_offset * SIZE_2) % num_elements) == 0 &&
      !index1_offset) {
    aligned = 2;
    index2 = index2_offset * SIZE_2 / num_elements;
  }
  else if ((num_elements % SIZE) == 0 && 
           (index3_offset * SIZE % num_elements) == 0 &&
           !index2_offset && !index1_offset) {
    index3 = index3_offset * SIZE  / num_elements;
    aligned = -1;
  }

  const size_t index_offset = index3_offset * SIZE + index2_offset * SIZE_2 +
                              index1_offset;
  int len = end - start;
  assert(len <= buf_size / num_elements + 2);
  assert(len >= buf_size / num_elements);
  if (1 == len) {
#ifdef __cplusplus
    MARS_WIDE_BUS_TYPE tmp(a[start]);
#else
    MARS_WIDE_BUS_TYPE tmp = a[start];
#endif
    for (j = 0; j < num_elements; ++j) {
       if (j < head_align || j > tail_align)
         continue;
       size_t buffer_index = j - head_align + index_offset;
       a_buf[buffer_index / SIZE][(buffer_index % SIZE) / SIZE_2][(buffer_index % SIZE) % SIZE_2] =
           tmp2(((j + 1) * data_width * 8 - 1), (j * data_width * 8));
    }
    return;
  }
  for (i = 0; i < len; ++i) {
#pragma HLS pipeline
#ifdef __cplusplus
    MARS_WIDE_BUS_TYPE tmp(a[i + start]);
#else
    MARS_WIDE_BUS_TYPE tmp = a[i + start];
#endif
    if (head_align == 0 && aligned) {
      if (aligned == 1) {
        for (j = 0; j < num_elements; ++j) {
          if (i == end - start - 1 && j > tail_align)
            continue;
          size_t buffer_index = i * num_elements + j - 0;
          a_buf[index3][index2][index1 * num_elements + j] =
                tmp2(((j + 1) * data_width * 8 - 1), (j * data_width * 8));
        }
        index1++;
        if (index1 == bound1) {
          index1 = 0;
          ++index2;
        }
        if (index2 == SIZE_1) {
          index2 = 0;
          ++index3;
        }
      } else if (aligned == 2) {
        int j0, j1;
        for (j0 = 0; j0 < num_elements / SIZE_2; ++j0) {
          for (j1 = 0; j1 < SIZE_2; ++j1) {
            j = j0 * SIZE_2 + j1;
            if (i == end - start - 1 && j > tail_align)
              continue;
            a_buf[index3][index2 * num_elements / SIZE_2  + j0][j1] =
                tmp2(((j + 1) * data_width * 8 - 1), (j * data_width * 8));
          }
        }
        ++index2;
        if (index2 == bound2) {
          index2 = 0;
           ++index3;
        }
      } else if (aligned == -1) {
        int j0, j1, j2;
        for (j0 = 0; j0 < bound0; ++j0) {
          for (j1 = 0; j1 < SIZE_1; ++j1) {
            for (j2 = 0; j2 < SIZE_2; ++j2) {
              j = j0 * SIZE + j1 * SIZE_2 + j2;
              if (i == end - start - 1 && j > tail_align)
                continue;
              a_buf[index3 * bound0 + j0][j1][j2] =
                  tmp2(((j + 1) * data_width * 8 - 1), (j * data_width * 8));
            }
          }
        }
        ++index3;
      }
    } else {
      for (j = 0; j < num_elements; ++j) {
        if (i == 0 && j < head_align)
          continue;
        if (i == end - start - 1 && j > tail_align)
          continue;
        size_t buffer_index = i * num_elements + j - head_align + index_offset;
        a_buf[buffer_index / SIZE][(buffer_index % SIZE) / SIZE_2][(buffer_index % SIZE) % SIZE_2] =
            tmp2(((j + 1) * data_width * 8 - 1), (j * data_width * 8));
      }
    }
  }
}

static void
memcpy_wide_bus(write, int, SIZE_1, SIZE_2)(MARS_WIDE_BUS_TYPE *c, int c_buf[][SIZE_1][SIZE_2],
                                                   size_t index3_offset,
                                                   size_t index2_offset,
                                                   size_t index1_offset,
                             size_t offset_byte, size_t size_byte) {
#pragma HLS inline self
  const size_t data_width = sizeof(int);
  const size_t bus_width = BUS_WIDTH;
  const size_t num_elements = bus_width / data_width;
  size_t buf_size = size_byte / data_width;
  size_t offset = offset_byte / data_width;
  size_t head_align = offset & (num_elements - 1);
  size_t new_offset = offset + buf_size;
  size_t tail_align = (new_offset - 1) & (num_elements - 1);
  size_t start = offset / num_elements;
  size_t end = (offset + buf_size + num_elements - 1) / num_elements;
  size_t len = end - start;
  size_t i, j;
  const size_t SIZE = SIZE_1 * SIZE_2;
  const size_t index_offset = index3_offset * SIZE + index2_offset * SIZE_2 +
                              index1_offset;
  if (head_align == 0)
    len = (buf_size + num_elements - 1) / num_elements;
  size_t align = 0;
  if (len == 1) {
    MARS_WIDE_BUS_TYPE tmp;
    if (head_align != 0 || tail_align != (num_elements - 1))
      tmp = c[start];
    for (j = 0; j < num_elements; ++j) {
      if (j < head_align)
        continue;
      if (j > tail_align)
        continue;
      size_t buffer_index = j - head_align + index_offset ;
      tmp3(((j + 1) * data_width * 8 - 1), (j * data_width * 8),
           c_buf[buffer_index / SIZE][(buffer_index % SIZE) / SIZE_2][(buffer_index % SIZE) %
                                                           SIZE_2]);
    }
    c[start] = tmp;
    return;
  }
  if (head_align != 0) {
    MARS_WIDE_BUS_TYPE tmp = c[start];
    for (j = 0; j < num_elements; ++j) {
      if (j < head_align)
        continue;
      size_t buffer_index = j - head_align + index_offset;
      tmp3(((j + 1) * data_width * 8 - 1), (j * data_width * 8),
           c_buf[buffer_index / SIZE][(buffer_index % SIZE) / SIZE_2][
                 (buffer_index % SIZE) % SIZE_2]);
    }
    c[start] = tmp;
    start++;
    align++;
  }
  if (tail_align != (num_elements - 1))
    align++;
  size_t index3 = index3_offset, index2 = index2_offset, index1 = index1_offset;
  const size_t bound1 = SIZE_2 / num_elements;
  const size_t bound2 = SIZE / num_elements;
  const size_t bound0 = num_elements / SIZE;
  int aligned = 0;
  if ((SIZE_2 % num_elements) == 0 && (index1_offset % num_elements) == 0) {
    aligned = 1;
    index1 = index1_offset / num_elements;
  }
  else if ((SIZE % num_elements) == 0 && 
           (num_elements % SIZE_2) == 0 &&
      ((index2_offset * SIZE_2) % num_elements) == 0 &&
      !index1_offset) {
    aligned = 2;
    index2 = index2_offset * SIZE_2 / num_elements;
  }
  else if ((num_elements % SIZE) == 0 && 
           (index3_offset * SIZE % num_elements) == 0 &&
           !index2_offset && !index1_offset) {
    index3 = index3_offset * SIZE  / num_elements;
    aligned = -1;
  }
  int burst_len = len - align;
  assert(burst_len <= buf_size / num_elements);
  for (i = 0; i < burst_len; ++i) {
#pragma HLS pipeline
    MARS_WIDE_BUS_TYPE tmp;
    if (head_align == 0 && aligned) {
      if (aligned == 1) {
        for (j = 0; j < num_elements; ++j) {
          int val =  c_buf[index3][index2][index1 * num_elements + j];
          tmp3(((j + 1) * data_width * 8 - 1), (j * data_width * 8), val);
        }
        index1++;
        if (index1 == bound1) {
          index1 = 0;
          ++index2;
        }
        if (index2 == SIZE_1) {
          index2 = 0;
          ++index3;
        }
      } else if (aligned == 2) {
        int j0, j1;
        for (j0 = 0; j0 < num_elements / SIZE_2; ++j0) {
          for (j1 = 0; j1 < SIZE_2; ++j1) {
            j = j0 * SIZE_2 + j1;
            int val = c_buf[index3][index2 * num_elements / SIZE_2 + j0][j1];
            tmp3(((j + 1) * data_width * 8 - 1), (j * data_width * 8), val);
          }
        }
        ++index2;
        if (index2 == bound2) {
          index2 = 0;
          ++index3;
        }
      } else if (aligned == -1) {
        int j0, j1, j2;
        for (j0 = 0; j0 < bound0; ++j0) {
          for (j1 = 0; j1 < SIZE_1; ++j1) {
            for (j2 = 0; j2 < SIZE_2; ++j2) {
              j = j0 * SIZE + j1 * SIZE_2 + j2;
              int val = c_buf[index3 * bound0 + j0][j1][j2];
              tmp3(((j + 1) * data_width * 8 - 1), (j * data_width * 8), val);
            }
          }
        }
        ++index3;
      }
    } else {
      for (j = 0; j < num_elements; ++j) {
        size_t buffer_index = i * num_elements + j + num_elements - head_align + index_offset;
        tmp3(((j + 1) * data_width * 8 - 1), (j * data_width * 8),
             c_buf[buffer_index / SIZE][(buffer_index % SIZE) / SIZE_2][(buffer_index % SIZE) %
                                                             SIZE_2]);
      }
    }
    c[i + start] = tmp;
  }
  if (tail_align != num_elements - 1) {
    MARS_WIDE_BUS_TYPE tmp = c[end - 1];
    size_t pos = (len - align) * num_elements+ index_offset;
    pos += (num_elements - head_align) % num_elements;
    for (j = 0; j < num_elements; ++j) {
      if (j > tail_align)
        continue;
      size_t buffer_index = pos + j;
      tmp3(((j + 1) * data_width * 8 - 1), (j * data_width * 8),
           c_buf[buffer_index / SIZE][(buffer_index % SIZE) / SIZE_2][(buffer_index % SIZE) %
                                                           SIZE_2]);
    }
    c[end - 1] = tmp;
  }
}

static void memcpy_wide_bus(read, float, SIZE_1, SIZE_2)(float a_buf[][SIZE_1][SIZE_2],
                                                   size_t index3_offset,
                                                   size_t index2_offset,
                                                   size_t index1_offset,
                                                 MARS_WIDE_BUS_TYPE *a,
                                                 size_t offset_byte,
                                                 size_t size_byte) {
#pragma HLS inline self
  const size_t data_width = sizeof(float);
  const size_t bus_width = BUS_WIDTH;
  const size_t num_elements = bus_width / data_width;
  size_t buf_size = size_byte / data_width;
  size_t offset = offset_byte / data_width;
  size_t head_align = offset & (num_elements - 1);
  size_t new_offset = offset + buf_size;
  size_t tail_align = (new_offset - 1) & (num_elements - 1);
  size_t start = offset / num_elements;
  size_t end = (offset + buf_size + num_elements - 1) / num_elements;
  //MARS_WIDE_BUS_TYPE *a_offset = a + start;
  size_t i, j;
  size_t index3 = index3_offset, index2 = index2_offset, index1 = index1_offset;
  const size_t bound1 = SIZE_2 / num_elements;
  const size_t SIZE = SIZE_1 * SIZE_2;
  const size_t bound2 = SIZE / num_elements;
  const size_t bound0 = num_elements / SIZE;
  int aligned = 0;
  if ((SIZE_2 % num_elements) == 0 && (index1_offset % num_elements) == 0) {
    aligned = 1;
    index1 = index1_offset / num_elements;
  }
  else if ((SIZE % num_elements) == 0 && 
           (num_elements % SIZE_2) == 0 &&
      ((index2_offset * SIZE_2) % num_elements) == 0 &&
      !index1_offset) {
    aligned = 2;
    index2 = index2_offset * SIZE_2 / num_elements;
  }
  else if ((num_elements % SIZE) == 0 && 
           (index3_offset * SIZE % num_elements) == 0 &&
           !index2_offset && !index1_offset) {
    index3 = index3_offset * SIZE  / num_elements;
    aligned = -1;
  }
  const size_t index_offset = index3_offset * SIZE + index2_offset * SIZE_2 +
                              index1_offset;
  int len = end - start;
  assert(len <= buf_size / num_elements + 2);
  assert(len >= buf_size / num_elements);
  if (1 == len) {
#ifdef __cplusplus
    MARS_WIDE_BUS_TYPE tmp(a[start]);
#else
    MARS_WIDE_BUS_TYPE tmp = a[start];
#endif
    for (j = 0; j < num_elements; ++j) {
      if (j < head_align || j > tail_align)
         continue;
      size_t buffer_index = j - head_align + index_offset;
      int raw_bits =
          tmp2(((j + 1) * data_width * 8 - 1), (j * data_width * 8));
      a_buf[buffer_index / SIZE][(buffer_index % SIZE) / SIZE_2][(buffer_index % SIZE) % SIZE_2] =
          *(float *)(&raw_bits);
    }
    return;
  }
  for (i = 0; i < len; ++i) {
#pragma HLS pipeline
#ifdef __cplusplus
    MARS_WIDE_BUS_TYPE tmp(a[i + start]);
#else
    MARS_WIDE_BUS_TYPE tmp = a[i + start];
#endif

    if (head_align == 0 && aligned) {
      if (aligned == 1) {
        for (j = 0; j < num_elements; ++j) {
          if (i == end - start - 1 && j > tail_align)
            continue;
          int raw_bits =
              tmp2(((j + 1) * data_width * 8 - 1), (j * data_width * 8));
          a_buf[index3][index2][index1 * num_elements + j] =
                *(float *)(&raw_bits);
        }
        index1++;
        if (index1 == bound1) {
          index1 = 0;
          ++index2;
        }
        if (index2 == SIZE_1) {
          index2 = 0;
          ++index3;
        }
      } else if (aligned == 2) {
        int j0, j1;
        for (j0 = 0; j0 < num_elements / SIZE_2; ++j0) {
          for (j1 = 0; j1 < SIZE_2; ++j1) {
            j = j0 * SIZE_2  + j1;
            if (i == end - start - 1 && j > tail_align)
              continue;
            int raw_bits =
                tmp2(((j + 1) * data_width * 8 - 1), (j * data_width * 8));
            a_buf[index3][index2 * num_elements / SIZE_2 + j0][j1] =
                  *(float *)(&raw_bits);
          }
        }
        ++index2;
        if (index2 == bound2) {
          index2 = 0;
          ++index3;
        }
      } else if (aligned == -1) {
        int j0, j1, j2;
        for (j0 = 0; j0 < bound0; ++j0) {
          for (j1 = 0; j1 < SIZE_1; ++j1) {
            for (j2 = 0; j2 < SIZE_2; ++j2) {
              j = j0 * SIZE + j1 * SIZE_2 + j2;
              if (i == end - start - 1 && j > tail_align)
                continue;
              int raw_bits =
                  tmp2(((j + 1) * data_width * 8 - 1), (j * data_width * 8));
              a_buf[index3 * bound0 + j0][j1][j2] =
                    *(float *)(&raw_bits);
            }
          }
        }
        ++index3;
      }
    } else {
      for (j = 0; j < num_elements; ++j) {
        if (i == 0 && j < head_align)
          continue;
        if (i == end - start - 1 && j > tail_align)
          continue;
        size_t buffer_index = i * num_elements + j - head_align + index_offset;
        int raw_bits =
            tmp2(((j + 1) * data_width * 8 - 1), (j * data_width * 8));
        a_buf[buffer_index / SIZE][(buffer_index % SIZE) / SIZE_2][(buffer_index % SIZE) % SIZE_2] =
            *(float *)(&raw_bits);
      }
    }
  }
}

static void memcpy_wide_bus(write, float, SIZE_1, SIZE_2)(MARS_WIDE_BUS_TYPE *c,
                                                  float c_buf[][SIZE_1][SIZE_2],
                                                   size_t index3_offset,
                                                   size_t index2_offset,
                                                   size_t index1_offset,
                                                  size_t offset_byte,
                                                  size_t size_byte) {
#pragma HLS inline self
  const size_t data_width = sizeof(float);
  const size_t bus_width = BUS_WIDTH;
  const size_t num_elements = bus_width / data_width;
  size_t buf_size = size_byte / data_width;
  size_t offset = offset_byte / data_width;
  size_t head_align = offset & (num_elements - 1);
  size_t new_offset = offset + buf_size;
  size_t tail_align = (new_offset - 1) & (num_elements - 1);
  size_t start = offset / num_elements;
  size_t end = (offset + buf_size + num_elements - 1) / num_elements;
  size_t len = end - start;
  size_t i, j;
  const size_t SIZE = SIZE_1 * SIZE_2;
  const size_t index_offset = index3_offset * SIZE + index2_offset * SIZE_2 +
                              index1_offset;
  if (head_align == 0)
    len = (buf_size + num_elements - 1) / num_elements;
  if (len == 1) {
    MARS_WIDE_BUS_TYPE tmp;
    if (head_align != 0 || tail_align != (num_elements - 1))
      tmp = c[start];
    for (j = 0; j < num_elements; ++j) {
      if (j < head_align)
        continue;
      if (j > tail_align)
        continue;
      size_t buffer_index = j - head_align + index_offset;
      float buf_tmp = c_buf[buffer_index / SIZE][
        (buffer_index % SIZE) / SIZE_2][(buffer_index % SIZE) % SIZE_2];
      int raw_bits = *(int *)&buf_tmp;
      tmp3(((j + 1) * data_width * 8 - 1), (j * data_width * 8), raw_bits);
    }
    c[start] = tmp;
    return;
  }
  unsigned align = 0;
  if (head_align != 0) {
    MARS_WIDE_BUS_TYPE tmp = c[start];
    for (j = 0; j < num_elements; ++j) {
      if (j < head_align)
        continue;
      size_t buffer_index = j - head_align + index_offset;
      float buf_tmp = c_buf[buffer_index / SIZE][
        (buffer_index % SIZE) / SIZE_2][(buffer_index % SIZE) % SIZE_2];
      int raw_bits = *(int *)&buf_tmp;
      tmp3(((j + 1) * data_width * 8 - 1), (j * data_width * 8), raw_bits);
    }
    c[start] = tmp;
    start++;
    align++;
  }
  if (tail_align != (num_elements - 1))
    align++;
  size_t index3 = index3_offset, index2 = index2_offset, index1 = index1_offset;
  const size_t bound1 = SIZE_2 / num_elements;
  const size_t bound2 = SIZE / num_elements;
  const size_t bound0 = num_elements / SIZE;
  int aligned = 0;
  if ((SIZE_2 % num_elements) == 0 && (index1_offset % num_elements) == 0) {
    aligned = 1;
    index1 = index1_offset / num_elements;
  }
  else if ((SIZE % num_elements) == 0 && 
           (num_elements % SIZE_2) == 0 &&
      ((index2_offset * SIZE_2) % num_elements) == 0 &&
      !index1_offset) {
    aligned = 2;
    index2 = index2_offset * SIZE_2 / num_elements;
  }
  else if ((num_elements % SIZE) == 0 && 
           (index3_offset * SIZE % num_elements) == 0 &&
           !index2_offset && !index1_offset) {
    index3 = index3_offset * SIZE  / num_elements;
    aligned = -1;
  }
  int burst_len = len - align;
  assert(burst_len <= buf_size / num_elements);
  for (i = 0; i < burst_len; ++i) {
#pragma HLS pipeline
    MARS_WIDE_BUS_TYPE tmp;
    if (head_align == 0 && aligned) {
      if (aligned == 1) {
        for (j = 0; j < num_elements; ++j) {

          float buf_tmp =
              c_buf[index3][index2][index1 * num_elements + j];
          int raw_bits = *(int *)&buf_tmp;
          tmp3(((j + 1) * data_width * 8 - 1), (j * data_width * 8), raw_bits);
        }
        index1++;
        if (index1 == bound1) {
          index1 = 0;
          ++index2;
        }
        if (index2 == SIZE_1) {
          index2 = 0;
          ++index3;
        }
      } else if (aligned == 2) {
        int j0, j1;
        for (j0 = 0; j0 < num_elements / SIZE_2; ++j0) {
          for (j1 = 0; j1 < SIZE_2; ++j1) {
            j = j0 * SIZE_2 + j1;
            float buf_tmp =
                c_buf[index3][index2 * num_elements / SIZE_2 + j0][j1];
            int raw_bits = *(int *)&buf_tmp;
            tmp3(((j + 1) * data_width * 8 - 1), (j * data_width * 8), raw_bits);
          }
        }
        ++index2;
        if (index2 == bound2) {
          index2 = 0;
          ++index3;
        }
      } else if (aligned == -1) {
        int j0, j1, j2;
        for (j0 = 0; j0 < bound0; ++j0) {
          for (j1 = 0; j1 < SIZE_1; ++j1) {
            for (j2 = 0; j2 < SIZE_2; ++j2) {
              j = j0 * SIZE + j1 * SIZE_2 + j2;
              float buf_tmp =
                  c_buf[index3 *bound0 + j0][j1][j2];
              int raw_bits = *(int *)&buf_tmp;
              tmp3(((j + 1) * data_width * 8 - 1), (j * data_width * 8), raw_bits);
            }
          }
        }
        ++index3;
      }
    }
    else {
      for (j = 0; j < num_elements; ++j) {
        size_t buffer_index = i * num_elements + j + num_elements - head_align + index_offset;
        float buf_tmp = c_buf[buffer_index / SIZE][(buffer_index % SIZE) /SIZE_2][(buffer_index % SIZE) % SIZE_2];
        int raw_bits = *(int *)&buf_tmp;
        tmp3(((j + 1) * data_width * 8 - 1), (j * data_width * 8), raw_bits);
      }
    }
    c[i + start] = tmp;
  }
  if (tail_align != num_elements - 1) {
    MARS_WIDE_BUS_TYPE tmp = c[end - 1];
    size_t pos = (len - align) * num_elements + index_offset;
    pos += (num_elements - head_align) % num_elements;
    for (j = 0; j < num_elements; ++j) {
      if (j > tail_align)
        continue;
      size_t buffer_index = pos + j;
      float buf_tmp = c_buf[buffer_index / SIZE][(buffer_index % SIZE) /SIZE_2][(buffer_index % SIZE) % SIZE_2];
      int raw_bits = *(int *)&buf_tmp;
      tmp3(((j + 1) * data_width * 8 - 1), (j * data_width * 8), raw_bits);
    }
    c[end - 1] = tmp;
  }
}

static void
memcpy_wide_bus(read, double, SIZE_1, SIZE_2)(double a_buf[][SIZE_1][SIZE_2],
                                                   size_t index3_offset,
                                                   size_t index2_offset,
                                                   size_t index1_offset,
                               MARS_WIDE_BUS_TYPE *a, size_t offset_byte,
                               size_t size_byte) {
#pragma HLS inline self
  const size_t data_width = sizeof(double);
  const size_t bus_width = BUS_WIDTH;
  const size_t num_elements = bus_width / data_width;
  size_t buf_size = size_byte / data_width;
  size_t offset = offset_byte / data_width;
  size_t head_align = offset & (num_elements - 1);
  size_t new_offset = offset + buf_size;
  size_t tail_align = (new_offset - 1) & (num_elements - 1);
  size_t start = offset / num_elements;
  size_t end = (offset + buf_size + num_elements - 1) / num_elements;
  //MARS_WIDE_BUS_TYPE *a_offset = a + start;
  size_t i, j;
  size_t index3 = index3_offset, index2 = index2_offset, index1 = index1_offset;
  const size_t SIZE = SIZE_1 * SIZE_2;
  const size_t index_offset = index3_offset * SIZE + index2_offset * SIZE_2 +
                              index1_offset;
  const size_t bound1 = SIZE_2 / num_elements;
  const size_t bound2 = SIZE / num_elements;
  const size_t bound0 = num_elements / SIZE;
  int aligned = 0;
  if ((SIZE_2 % num_elements) == 0 && (index1_offset % num_elements) == 0) {
    aligned = 1;
    index1 = index1_offset / num_elements;
  }
  else if ((SIZE % num_elements) == 0 && 
           (num_elements % SIZE_2) == 0 &&
      ((index2_offset * SIZE_2) % num_elements) == 0 &&
      !index1_offset) {
    aligned = 2;
    index2 = index2_offset * SIZE_2 / num_elements;
  }
  else if ((num_elements % SIZE) == 0 && 
           (index3_offset * SIZE % num_elements) == 0 &&
           !index2_offset && !index1_offset) {
    index3 = index3_offset * SIZE  / num_elements;
    aligned = -1;
  }
  int len = end - start;
  assert(len <= buf_size / num_elements + 2);
  assert(len >= buf_size / num_elements);
  if (1 == len) {
#ifdef __cplusplus
    MARS_WIDE_BUS_TYPE tmp(a[start]);
#else
    MARS_WIDE_BUS_TYPE tmp = a[start];
#endif
    for (j = 0; j < num_elements; ++j) {
       if (j < head_align || j > tail_align)
         continue;

        size_t buffer_index = j - head_align + index_offset;
        long long raw_bits =
            tmp2(((j + 1) * data_width * 8 - 1), (j * data_width * 8));
        a_buf[buffer_index / SIZE][(buffer_index % SIZE) / SIZE_2][(buffer_index % SIZE) % SIZE_2] =
            *(double *)(&raw_bits);
    }
    return;
  }
  for (i = 0; i < len; ++i) {
#pragma HLS pipeline
#ifdef __cplusplus
    MARS_WIDE_BUS_TYPE tmp(a[i + start]);
#else
    MARS_WIDE_BUS_TYPE tmp = a[i + start];
#endif
    if (head_align == 0 && aligned) {
      if (aligned == 1) {
        for (j = 0; j < num_elements; ++j) {
          if (i == end - start - 1 && j > tail_align)
            continue;
          long long raw_bits =
              tmp2(((j + 1) * data_width * 8 - 1), (j * data_width * 8));
          a_buf[index3][index2][index1 * num_elements + j] =
                *(double *)(&raw_bits);
        }
        index1++;
        if (index1 == bound1) {
          index1 = 0;
          ++index2;
        }
        if (index2 == SIZE_1) {
          index2 = 0;
          ++index3;
        }
      } else if (aligned == 2) {
        int j0, j1;
        for (j0 = 0; j0 < num_elements / SIZE_2; ++j0) {
          for (j1 = 0; j1 < SIZE_2; ++j1) {
            j = j0 * SIZE_2  + j1;
            if (i == end - start - 1 && j > tail_align)
              continue;
            long long raw_bits =
                tmp2(((j + 1) * data_width * 8 - 1), (j * data_width * 8));
            a_buf[index3][index2 * num_elements / SIZE_2 + j0][j1] =
                  *(double *)(&raw_bits);
          }
        }
        ++index2;
        if (index2 == bound2) {
          index2 = 0;
          ++index3;
        }
      } else if (aligned == -1) {
        int j0, j1, j2;
        for (j0 = 0; j0 < bound0; ++j0) {
          for (j1 = 0; j1 < SIZE_1; ++j1) {
            for (j2 = 0; j2 < SIZE_2; ++j2) {
              j = j0 * SIZE + j1 * SIZE_2 + j2;
              if (i == end - start - 1 && j > tail_align)
                continue;
              long long raw_bits =
                  tmp2(((j + 1) * data_width * 8 - 1), (j * data_width * 8));
              a_buf[index3 * bound0 + j0][j1][j2] =
                    *(double *)(&raw_bits);
            }
          }
        }
        ++index3;
      }
    } else {
      for (j = 0; j < num_elements; ++j) {
        if (i == 0 && j < head_align)
          continue;
        if (i == end - start - 1 && j > tail_align)
          continue;
        size_t buffer_index = i * num_elements + j - head_align + index_offset;
        long long raw_bits =
            tmp2(((j + 1) * data_width * 8 - 1), (j * data_width * 8));
        a_buf[buffer_index / SIZE][(buffer_index % SIZE) / SIZE_2][(buffer_index % SIZE) % SIZE_2] =
            *(double *)(&raw_bits);
      }
    }
  }
}

static void memcpy_wide_bus(write, double, SIZE_1, SIZE_2)(
    MARS_WIDE_BUS_TYPE *c, double c_buf[][SIZE_1][SIZE_2],
                                                   size_t index3_offset,
                                                   size_t index2_offset,
                                                   size_t index1_offset,
    size_t offset_byte, size_t size_byte) {
#pragma HLS inline self
  const size_t data_width = sizeof(double);
  const size_t bus_width = BUS_WIDTH;
  const size_t num_elements = bus_width / data_width;
  size_t buf_size = size_byte / data_width;
  size_t offset = offset_byte / data_width;
  size_t head_align = offset & (num_elements - 1);
  size_t new_offset = offset + buf_size;
  size_t tail_align = (new_offset - 1) & (num_elements - 1);
  size_t start = offset / num_elements;
  size_t end = (offset + buf_size + num_elements - 1) / num_elements;
  size_t len = end - start;
  size_t i, j;
  const size_t bound1 = SIZE_2 / num_elements;
  const size_t SIZE = SIZE_1 * SIZE_2;
  const size_t bound2 = SIZE / num_elements;
  const size_t bound0 = num_elements / SIZE;
  int aligned = 0;
  size_t index3 = index3_offset, index2 = index2_offset, index1 = index1_offset;
  if ((SIZE_2 % num_elements) == 0 && (index1_offset % num_elements) == 0) {
    aligned = 1;
    index1 = index1_offset / num_elements;
  }
  else if ((SIZE % num_elements) == 0 &&
           (num_elements % SIZE_2) == 0 &&
      ((index2_offset * SIZE_2) % num_elements) == 0 &&
      !index1_offset) {
    aligned = 2;
    index2 = index2_offset * SIZE_2 / num_elements;
  }
  else if ((num_elements % SIZE) == 0 && 
           (index3_offset * SIZE % num_elements) == 0 &&
           !index2_offset && !index1_offset) {
    index3 = index3_offset * SIZE  / num_elements;
    aligned = -1;
  }
  const size_t index_offset = index3_offset * SIZE + index2_offset * SIZE_2 +
                              index1_offset;
  if (head_align == 0)
    len = (buf_size + num_elements - 1) / num_elements;
  if (len == 1) {
    MARS_WIDE_BUS_TYPE tmp;
    if (head_align != 0 || tail_align != (num_elements - 1))
      tmp = c[start];
    for (j = 0; j < num_elements; ++j) {
      if (j < head_align)
        continue;
      if (j > tail_align)
        continue;
      size_t buffer_index = j - head_align + index_offset;
      double buf_tmp = c_buf[buffer_index / SIZE][(buffer_index % SIZE) /SIZE_2][(buffer_index % SIZE) % SIZE_2];
      long long raw_bits = *(long long *)&buf_tmp;
      tmp3(((j + 1) * data_width * 8 - 1), (j * data_width * 8), raw_bits);
    }
    c[start] = tmp;
    return;
  }
  size_t align = 0;
  if (head_align != 0) {
    MARS_WIDE_BUS_TYPE tmp = c[start];
    for (j = 0; j < num_elements; ++j) {
      if (j < head_align)
        continue;
      size_t buffer_index = j - head_align + index_offset;
      double buf_tmp = c_buf[buffer_index / SIZE][(buffer_index % SIZE) /SIZE_2][(buffer_index % SIZE) % SIZE_2];
      long long raw_bits = *(long long *)&buf_tmp;
      tmp3(((j + 1) * data_width * 8 - 1), (j * data_width * 8), raw_bits);
    }
    c[start] = tmp;
    start++;
    align++;
  }
  if (tail_align != (num_elements - 1))
    align++;
  int burst_len = len - align;
  assert(burst_len <= buf_size / num_elements);
  for (i = 0; i < burst_len; ++i) {
#pragma HLS pipeline
    MARS_WIDE_BUS_TYPE tmp;
    if (head_align == 0 && aligned) {
      if (aligned == 1) {
        for (j = 0; j < num_elements; ++j) {

          double buf_tmp =
              c_buf[index3][index2][index1 * num_elements + j];
          long long raw_bits = *(long long *)&buf_tmp;
          tmp3(((j + 1) * data_width * 8 - 1), (j * data_width * 8), raw_bits);
        }
        index1++;
        if (index1 == bound1) {
          index1 = 0;
          ++index2;
        }
        if (index2 == SIZE_1) {
          index2 = 0;
          ++index3;
        }
      } else if (aligned == 2) {
        int j0, j1;
        for (j0 = 0; j0 < num_elements / SIZE_2; ++j0) {
          for (j1 = 0; j1 < SIZE_2; ++j1) {
            j = j0 * SIZE_2 + j1;
            double buf_tmp =
                c_buf[index3][index2 * num_elements / SIZE_2 + j0][j1];
            long long raw_bits = *(long long *)&buf_tmp;
            tmp3(((j + 1) * data_width * 8 - 1), (j * data_width * 8), raw_bits);
          }
        }
        ++index2;
        if (index2 == bound2) {
          index2 = 0;
          ++index3;
        }
      } else if (aligned == -1) {
        int j0, j1, j2;
        for (j0 = 0; j0 < bound0; ++j0) {
          for (j1 = 0; j1 < SIZE_1; ++j1) {
            for (j2 = 0; j2 < SIZE_2; ++j2) {
              j = j0 * SIZE + j1 * SIZE_2 + j2;
              double buf_tmp =
                  c_buf[index3 *bound0 + j0][j1][j2];
              long long raw_bits = *(long long *)&buf_tmp;
              tmp3(((j + 1) * data_width * 8 - 1), (j * data_width * 8), raw_bits);
            }
          }
        }
        ++index3;
      }
    }
    else {
      for (j = 0; j < num_elements; ++j) {
        size_t buffer_index = i * num_elements + j + num_elements - head_align + index_offset;
        double buf_tmp =
            c_buf[buffer_index / SIZE][(buffer_index % SIZE) / SIZE_2][(buffer_index % SIZE) %SIZE_2];
        long long raw_bits = *(long long *)&buf_tmp;
        tmp3(((j + 1) * data_width * 8 - 1), (j * data_width * 8), raw_bits);
      }
    }
    c[i + start] = tmp;
  }
  if (tail_align != num_elements - 1) {
    MARS_WIDE_BUS_TYPE tmp = c[end - 1];
    size_t pos = (len - align) * num_elements + index_offset;
    pos += (num_elements - head_align) % num_elements;
    for (j = 0; j < num_elements; ++j) {
      if (j > tail_align)
        continue;
      size_t buffer_index = pos + j;
      double buf_tmp = c_buf[buffer_index / SIZE][(buffer_index % SIZE) /SIZE_2][(buffer_index % SIZE) % SIZE_2];
      long long raw_bits = *(long long *)&buf_tmp;
      tmp3(((j + 1) * data_width * 8 - 1), (j * data_width * 8), raw_bits);
    }
    c[end - 1] = tmp;
  }
}

//#undef memcpy_wide_bus
//#undef instant_size
//
//#undef LARGE_BUS
//#undef BUS_WIDTH
//#undef cpp_get_range
//#undef c_get_range
//#undef cpp_set_range
//#undef c_set_range
//#undef tmp2
//#undef tmp3
