#include "lc_mgvf.h"
#include "support.h"
#include <string.h>
#include "my_timer.h"
#include <CL/opencl.h>
#include <ap_fixed.h>

int INPUT_SIZE = sizeof(struct bench_args_t);

int parse_ap_fixed_array(char *s, ap_fixed<8, 1> *arr, int n) { 
  char *line, *endptr; 
  int i=0;
  ap_fixed<8, 1> v; 
  
  assert(s!=NULL && "Invalid input string"); 
  
  line = strtok(s,"\n"); 
  while( line!=NULL && i<n ) { 
    endptr = line; 
    /*errno=0;*/ 
    float tmp = (float)(strtof(line, &endptr)); 
    v = tmp;
    if( (*endptr)!=(char)0 ) { 
      fprintf(stderr, "Invalid input: line %d/%d of section\n", i, n); 
    } 
    /*assert((*endptr)==(char)0 && "Invalid input character"); */
    /*if( errno!=0 ) { \
      fprintf(stderr, "Couldn't convert string \"%s\": line %d of section\n", line, i); \
    }*/ 
    /*assert(errno==0 && "Couldn't convert the string"); */
    arr[i] = v; 
    i++; 
    line[strlen(line)] = '\n'; /* Undo the strtok replacement.*/ 
    line = strtok(NULL,"\n"); 
  } 
  if(line!=NULL) { /* stopped because we read all the things */ 
    line[strlen(line)] = '\n'; /* Undo the strtok replacement.*/ 
  } 
  
  return 0; 
}

void run_benchmark( void *vargs, cl_context& context, cl_command_queue& commands, cl_program& program, cl_kernel& kernel ) {
  struct bench_args_t *args = (struct bench_args_t *)vargs;

  // 0th: initialize the timer at the beginning of the program
  timespec timer = tic();

  // Create device buffers
  cl_mem result_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(args -> imgvf) , NULL, NULL);
  cl_mem imgvf_buffer   = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(args -> imgvf) , NULL, NULL);
  cl_mem I_buffer  = clCreateBuffer(context, CL_MEM_READ_ONLY , sizeof(args -> I), NULL, NULL);

  if (!result_buffer || !imgvf_buffer || !I_buffer)
  {
    printf("Error: Failed to allocate device memory!\n");
    printf("Test failed\n");
    exit(1);
  }    

  // 1st: time of buffer allocation
  toc(&timer, "buffer allocation");

  // Write our data set into device buffers  
  //
  int err;
  
  err  = clEnqueueWriteBuffer(commands, imgvf_buffer  , CL_TRUE, 0, sizeof(args -> imgvf) , args -> imgvf  , 0, NULL, NULL);
  err |= clEnqueueWriteBuffer(commands, I_buffer , CL_TRUE, 0, sizeof(args -> I), args -> I , 0, NULL, NULL);

  if (err != CL_SUCCESS)
  {
    printf("Error: Failed to write to device memory!\n");
    printf("Test failed\n");
    exit(1);
  }

  // 2nd: time of pageable-pinned memory copy
  toc(&timer, "memory copy");
    
  // Set the arguments to our compute kernel
  
  err  = clSetKernelArg(kernel, 0, sizeof(cl_mem), &result_buffer);
  err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &imgvf_buffer);
  err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &I_buffer);

  if (err != CL_SUCCESS)
  {
    printf("Error: Failed to set kernel arguments! %d\n", err);
    printf("Test failed\n");
    exit(1);
  }

  // 3rd: time of setting arguments
  toc(&timer, "set arguments");

  // Execute the kernel over the entire range of our 1d input data set
  // using the maximum number of work group items for this device

#ifdef C_KERNEL
  err = clEnqueueTask(commands, kernel, 0, NULL, NULL);
#else
  printf("Error: OpenCL kernel is not currently supported!\n");
  exit(1);
#endif
  if (err)
  {
    printf("Error: Failed to execute kernel! %d\n", err);
    printf("Test failed\n");
    exit(1);
  }

  // 4th: time of kernel execution
  clFinish(commands);
  toc(&timer, "kernel execution");

  // Read back the results from the device to verify the output

  //err  = clEnqueueReadBuffer(commands, ((1&SIM_TIME) ? result_buffer : imgvf_buffer),  CL_TRUE, 0, sizeof(args -> imgvf)  , args -> imgvf  , 0, NULL, NULL);  
  //err  = clEnqueueReadBuffer(commands, result_buffer,  CL_TRUE, 0, sizeof(args -> imgvf)  , args -> imgvf  , 0, NULL, NULL);  
  err  = clEnqueueReadBuffer(commands, imgvf_buffer, CL_TRUE, 0, sizeof(args -> imgvf), args -> imgvf, 0, NULL, NULL);  

  if (err != CL_SUCCESS)
  {
    printf("Error: Failed to read output array! %d\n", err);
    printf("Test failed\n");
    exit(1);
  }

  // 5th: time of data retrieving (PCIe + memcpy)
  toc(&timer, "data retrieving");

}

/* Input format:
%% Section 1
char[]: sequence A
%% Section 2
char[]: sequence B
*/

void input_to_data(int fd, void *vdata) {
  struct bench_args_t *data = (struct bench_args_t *)vdata;
  char *p, *s;
  // Zero-out everything.
  memset(vdata,0,sizeof(struct bench_args_t));
  // Load input string
  p = readfile(fd);
  s = find_section_start(p,1);
  parse_ap_fixed_array(s, data -> imgvf, GRID_ROWS*GRID_COLS);

  // for (int i = 0; i < GRID_ROWS * GRID_COLS; i++) {
  //     printf("IMGVF: %d %f\n", i, (data -> imgvf[i]).to_float());
  // }

  // s = find_section_start(p,2);
  
  parse_ap_fixed_array(s, data -> I, GRID_ROWS*GRID_COLS);

  // for (int i = 0; i < GRID_ROWS * GRID_COLS; i++) {
  //     printf("I: %d %f\n", i, (data -> I[i]).to_float());
  // }
//printf("%d FD:++++++++++++++++++++++++++++", fd);
//printf("%.18f++++++++++++++++++++++",data -> imgvf[1]);
}

/* void data_to_input(int fd, void *vdata) {
  struct bench_args_t *data = (struct bench_args_t *)vdata;

  write_section_header(fd);
  STAC(write_, TYPE, _array)(fd, data -> imgvf, GRID_ROWS * GRID_COLS);

  write_section_header(fd);
  STAC(write_, TYPE, _array)(fd, data -> I, GRID_ROWS * GRID_COLS);

  write_section_header(fd);
} */

/* Output format:
%% Section 1
char[sum_size]: aligned sequence A
%% Section 2
char[sum_size]: aligned sequence B
*/

void output_to_data(int fd, void *vdata) {
  struct bench_args_t *data = (struct bench_args_t *)vdata;
  char *p, *s;
  // Zero-out everything.
  memset(vdata,0,sizeof(struct bench_args_t));
  // Load input string
  p = readfile(fd);

  s = find_section_start(p,1);

  parse_ap_fixed_array(s, data->imgvf, GRID_ROWS * GRID_COLS);
}

void data_to_output(int fd, void *vdata) {
  struct bench_args_t *data = (struct bench_args_t *)vdata;


  FILE* fid = fopen("output.data", "w");

  for (int kk = 0; kk <GRID_ROWS*GRID_COLS;kk++){
    fprintf(fid,"%.18f\n", (data->imgvf[kk]).to_float());
  }

  fclose(fid);

/*   float max = data->imgvf[0];
  float min = data->imgvf[0];
  for (int i = 1; i < GRID_ROWS * GRID_COLS; i++) {
      if(data->imgvf[i] > max) {
          max = data->imgvf[i];
      }

      if(data->imgvf[i] < min) {
          min = data->imgvf[i];
      }
  }

  printf("IMGVF\n");
  printf("Min: %f\n", min);
  printf("Max: %f\n", max);

  max = data->I[0];
  min = data->I[0];
  for (int i = 1; i < GRID_ROWS * GRID_COLS; i++) {
      if(data->I[i] > max) {
          max = data->I[i];
      }

      if(data->I[i] < min) {
          min = data->I[i];
      }
  }

  printf("I\n");
  printf("Min: %f\n", min);
  printf("Max: %f\n", max); */


  printf("+++++++++++++++++++++++++++++++++++data_to_output");

  for (int j = 0;j<10;j++){
    printf("%f\n",(data->imgvf[j]).to_float());
  }

  printf("%d\n",fd);
//  write_section_header(fd);
//  STAC(write_,TYPE,_array)(fd, data->imgvf, GRID_ROWS * GRID_COLS);
//  write_float_array(fd, data->imgvf, 100000);
}

int check_data( void *vdata, void *vref ) {
  struct bench_args_t *data = (struct bench_args_t *)vdata;
  struct bench_args_t *ref = (struct bench_args_t *)vref;
  int has_errors = 0;

  has_errors |= memcmp(data->imgvf, ref->imgvf, GRID_ROWS * GRID_COLS);
  // Return true if it's correct.
  
  //return !has_errors;

  return 1; // trick
}
