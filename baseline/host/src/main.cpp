/**********************************************************************
 * Nadesh Ramanathan, Imperial College London
 *
 * File: main.cpp
 *
 * Additional Comments: distributed under a BSD license, see LICENSE.txt
 *
 **********************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstring>
#include "CL/opencl.h"
#include "AOCLUtils/aocl_utils.h"

#include "filtering_algorithm_top.h"
#include "filtering_algorithm_util.h"
#include "build_kdTree.h"

using namespace aocl_utils;

#define STRING_BUFFER_LEN 1024

// OpenCL runtime configuration
static cl_platform_id platform = NULL;
static cl_device_id device = NULL;
static cl_context context = NULL;
static cl_command_queue queue = NULL;
static cl_kernel kernel = NULL;
static cl_program program = NULL;

// Function prototypes
bool init();
void cleanup();
static void device_info_ulong( cl_device_id device, cl_device_info param, const char* name);
static void device_info_uint( cl_device_id device, cl_device_info param, const char* name);
static void device_info_bool( cl_device_id device, cl_device_info param, const char* name);
static void device_info_string( cl_device_id device, cl_device_info param, const char* name);
static void display_device_info( cl_device_id device );


data_type initial_centre_positions[K];

kdTree_type tree_image_out[P][TREE_HEAP_SIZE];

void make_clusters_out_file_name(char *result, uint n, uint k, uint d, double std_dev, bool hex)
{
    if (!hex)
        sprintf(result,"golden_ref/clusters_out_N%d_K%d_D%d_s%.2f.mat",n,k,d,std_dev);
    else
        sprintf(result,"golden_ref/clusters_out_N%d_K%d_D%d_s%.2f.hex",n,k,d,std_dev);
}

void make_distortion_out_file_name(char *result, uint n, uint k, uint d, double std_dev, bool hex)
{
    if (!hex)
        sprintf(result,"golden_ref/distortion_out_N%d_K%d_D%d_s%.2f.mat",n,k,d,std_dev);
    else
        sprintf(result,"golden_ref/distortion_out_N%d_K%d_D%d_s%.2f.hex",n,k,d,std_dev);
}

void init_tree_node_memory( kdTree_type *in_r[P],
        uint * size,
        kdTree_type (*out)[TREE_HEAP_SIZE]
        ){

    for(int p =0; p < P; p++){
        kdTree_type *in = in_r[p];

        for (node_pointer i=0; i<size[p]; i++) {
            kdTree_type tmp_node;
            tmp_node = in[i];

            uint idx;
            data_type_ext dummy_wgtCent;
            data_type dummy_midPoint;
            data_type dummy_bnd_hi, dummy_bnd_lo;
            coord_type_ext dummy_sum_sq;
            coord_type_ext dummy_count;
            node_pointer dummy_left, dummy_right;
            get_kd_tree_type_items(	tmp_node,
                    &idx,
                    &dummy_wgtCent,
                    &dummy_midPoint,
                    &dummy_bnd_hi,
                    &dummy_bnd_lo,
                    &dummy_sum_sq,
                    &dummy_count,
                    &dummy_left,
                    &dummy_right);

            uint node_address = idx;


            out[p][node_address] = tmp_node;
        }
        printf("\n");
    }
}




// Entry point.
int main() {
    cl_int status;

    if(!init()) {
        return -1;
    }

    const uint k = 128;

    data_type clusters_out[K];
    coord_type_ext distortion_out[K];
    uint visited_nodes[L][P];
    uint nc_pairs[L][P];
    uint alloc_set[L][P];

    node_pointer root[P];
    uint size[P];

    posix_memalign ((void**)(&tree_image_out), 64, sizeof(tree_image_out));
    posix_memalign ((void**)(&initial_centre_positions), 64, sizeof(initial_centre_positions));
    posix_memalign ((void**)(&clusters_out), 64, sizeof(clusters_out));
    posix_memalign ((void**)(&distortion_out), 64, sizeof(distortion_out));
    posix_memalign ((void**)(&visited_nodes), 64, sizeof(visited_nodes));
    posix_memalign ((void**)(&visited_nodes), 64, sizeof(nc_pairs));
    posix_memalign ((void**)(&visited_nodes), 64, sizeof(alloc_set));

    cl_int err;

    cl_mem kdTree_buffer; 
    kdTree_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY, 
            sizeof(tree_image_out), NULL , &err);
    if(err < 0) {
        perror("Couldn't create a buffer");
        return -1;   
    }

    cl_mem initial_centre_buffer; 
    initial_centre_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY, 
            sizeof(initial_centre_positions), NULL , &err);
    if(err < 0) {
        perror("Couldn't create a buffer");
        return -1;   
    }

    cl_mem root_buffer; 
    root_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY, 
            sizeof(root), NULL , &err);
    if(err < 0) {
        perror("Couldn't create a buffer");
        return -1;   
    }

    const uint k_value = k;

    cl_mem k_buffer; 
    k_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY, 
            sizeof(k_value), NULL , &err);
    if(err < 0) {
        perror("Couldn't create a buffer");
        return -1;   
    }

    centre_heap_type heap[P][CENTRESET_HEAP_SIZE];
    cl_mem heap_buffer; 
    heap_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY, 
            sizeof(heap), NULL , &err);
    if(err < 0) {
        perror("Couldn't create a buffer");
        return -1;   
    }

    centre_list_pointer centre_freelist[P][CENTRESET_HEAP_SIZE];
    cl_mem list_buffer; 
    list_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY, 
            sizeof(centre_freelist), NULL , &err);
    if(err < 0) {
        perror("Couldn't create a buffer");
        return -1;   
    }


    cl_mem distortion_out_buffer; 
    distortion_out_buffer = clCreateBuffer(context, CL_MEM_WRITE_ONLY, 
            sizeof(distortion_out), NULL , &err);
    if(err < 0) {
        perror("Couldn't create a buffer");
        return -1;   
    }

    cl_mem clusters_out_buffer; 
    clusters_out_buffer = clCreateBuffer(context, CL_MEM_WRITE_ONLY, 
            sizeof(clusters_out), NULL , &err);
    if(err < 0) {
        perror("Couldn't create a buffer");
        return -1;   
    }

    cl_mem visited_nodes_buffer; 
    visited_nodes_buffer = clCreateBuffer(context, CL_MEM_WRITE_ONLY, 
            sizeof(visited_nodes), NULL , &err);
    if(err < 0) {
        perror("Couldn't create a buffer");
        return -1;   
    }

    cl_mem nc_pairs_buffer; 
    nc_pairs_buffer = clCreateBuffer(context, CL_MEM_WRITE_ONLY, 
            sizeof(nc_pairs), NULL , &err);
    if(err < 0) {
        perror("Couldn't create a buffer");
        return -1;   
    }


    cl_mem alloc_set_buffer; 
    alloc_set_buffer = clCreateBuffer(context, CL_MEM_WRITE_ONLY, 
            sizeof(alloc_set), NULL , &err);
    if(err < 0) {
        perror("Couldn't create a buffer");
        return -1;   
    }

    uint n = 1024*1024;
    double std_dev = 0.10;

    kdTree_type* tree_image_r[P];

    for (uint runs=0; runs<2; runs++) {

        // initialize random seed
        srand (time(NULL));

        if(runs==0){
            printf("Building balanced tree\n");
            if ( !build_subtrees(P, n, root, size, tree_image_r, initial_centre_positions, n, k, std_dev))
                return 1;
        }else {
            printf("Building unbalanced tree\n");
            if ( !build_degenerate_subtrees(P, n, root, size, tree_image_r, initial_centre_positions, n, k, std_dev))
                return 1;
        }

        init_tree_node_memory(tree_image_r, size, tree_image_out);

        printf("\nKernel initialization is complete.\n");
        printf("Launching the kernel...\n\n");


        // Configure work set over which the kernel will execute
        size_t wgSize[3] = {P, 1, 1};
        size_t gSize[3] = {P, 1, 1};


        err = clEnqueueWriteBuffer(queue, kdTree_buffer, CL_TRUE, 0,sizeof(tree_image_out), tree_image_out, 0, 0, 0);
        if (err != CL_SUCCESS) {
            return -1;
        }


        err = clEnqueueWriteBuffer(queue, initial_centre_buffer, CL_TRUE, 0,sizeof(initial_centre_positions), initial_centre_positions, 0, 0, 0);
        if (err != CL_SUCCESS) {
            return -1;
        }


        err = clEnqueueWriteBuffer(queue, root_buffer, CL_TRUE, 0,sizeof(root), root, 0, 0, 0);
        if (err != CL_SUCCESS) {
            return -1;
        }


        err = clEnqueueWriteBuffer(queue, k_buffer, CL_TRUE, 0,sizeof(k_value), &k_value, 0, 0, 0);
        if (err != CL_SUCCESS) {
            return -1;
        }

        int arg = 0;

        err = clSetKernelArg(kernel, arg++, sizeof(cl_mem), &kdTree_buffer);
        if (err != CL_SUCCESS) {
            return -1;
        };

        err = clSetKernelArg(kernel, arg++, sizeof(cl_mem), &initial_centre_buffer);
        if (err != CL_SUCCESS) {
            return -1;
        };

        err = clSetKernelArg(kernel, arg++, sizeof(cl_mem), &root_buffer);
        if (err != CL_SUCCESS) {
            return -1;
        };

        err = clSetKernelArg(kernel, arg++, sizeof(cl_mem), &k_buffer);
        if (err != CL_SUCCESS) {
            return -1;
        };

        err = clSetKernelArg(kernel, arg++, sizeof(cl_mem), &heap_buffer);
        if (err != CL_SUCCESS) {
            return -1;
        };

        err = clSetKernelArg(kernel, arg++, sizeof(cl_mem), &distortion_out_buffer);
        if (err != CL_SUCCESS) {
            return -1;
        };

        err = clSetKernelArg(kernel, arg++, sizeof(cl_mem), &clusters_out_buffer);
        if (err != CL_SUCCESS) {
            return -1;
        };

        err = clSetKernelArg(kernel, arg++, sizeof(cl_mem), &visited_nodes_buffer);
        if (err != CL_SUCCESS) {
            return -1;
        };

        err = clSetKernelArg(kernel, arg++, sizeof(cl_mem), &nc_pairs_buffer);
        if (err != CL_SUCCESS) {
            return -1;
        };


        err = clSetKernelArg(kernel, arg++, sizeof(cl_mem), &alloc_set_buffer);
        if (err != CL_SUCCESS) {
            return -1;
        };

        cl_event event;

        // Launch the kernel
        status = clEnqueueNDRangeKernel(queue, kernel, 1, NULL, gSize, wgSize, 0, NULL, &event);
        checkError(status, "Failed to launch kernel");

        clWaitForEvents(1,&event);
        // Wait for command queue to complete pending events
        status = clFinish(queue);
        checkError(status, "Failed to finish");

        cl_ulong time_start, time_end;

        double total_time;

        clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
        clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);
        total_time = time_end-time_start;

        printf("\nKernel execution is complete.\n");

        err = clEnqueueReadBuffer(queue, clusters_out_buffer, CL_TRUE, 0,sizeof(clusters_out), clusters_out, 0, 0, 0);
        if (err != CL_SUCCESS) {
            return -1;
        }

        err = clEnqueueReadBuffer(queue, distortion_out_buffer, CL_TRUE, 0,sizeof(distortion_out), distortion_out, 0, 0, 0);
        if (err != CL_SUCCESS) {
            return -1;
        }

        err = clEnqueueReadBuffer(queue, visited_nodes_buffer, CL_TRUE, 0,sizeof(visited_nodes), visited_nodes, 0, 0, 0);
        if (err != CL_SUCCESS) {
            return -1;
        }

        err = clEnqueueReadBuffer(queue, nc_pairs_buffer, CL_TRUE, 0,sizeof(nc_pairs), nc_pairs, 0, 0, 0);
        if (err != CL_SUCCESS) {
            return -1;
        }

        err = clEnqueueReadBuffer(queue, alloc_set_buffer, CL_TRUE, 0,sizeof(alloc_set), alloc_set, 0, 0, 0);
        if (err != CL_SUCCESS) {
            return -1;
        }

        printf("New centres after clustering\n");

        for (uint i=0; i<k; i++) {
            printf("%d: ",i);
            for (uint d=0; d<D-1; d++) {
                printf("%d ",(int)get_coord_type_vector_item(clusters_out[i].value, d));
            }
            printf("%d\n",(int)get_coord_type_vector_item(clusters_out[i].value, D-1));

        }

        if(runs==0){        
            printf("\nBalanced tree stats for baseline k-means\n");
        }else{
            printf("\nUnbalanced tree stats for baseline k-means\n");
        }

        printf("Visited Node Statistics\n");
        printf("Work-item:\t");
        for(uint p=0; p<P; p++){
            printf("%d\t",p);
        }
        printf("\n");
        for(uint l=0; l < L; l++){
            printf("%d:\t\t",l);
            uint total = 0;
            for(uint p=0; p<P; p++){
                printf("%d\t",visited_nodes[l][p] );
                total+=visited_nodes[l][p];
            }
            printf("\t = %d\n",total);
        } 

        printf("Node-Center Pairs Statistics\n");
        printf("Work-item:\t");
        for(uint p=0; p<P; p++){
            printf("%d\t",p);
        }
        printf("\n");
        for(uint l=0; l < L; l++){
            printf("%d:\t\t",l);
            uint total = 0;
            for(uint p=0; p<P; p++){
                printf("%d\t",nc_pairs[l][p] );
                total+=nc_pairs[l][p];
            }
            printf("\t = %d\n",total);
        }


        printf("OpenCl Execution time is: %0.3f milliseconds \n",total_time/1000000.0);
    }

    for (uint p=0; p<P; p++) {
        delete[] tree_image_r[p];
    }
    // free memory allocated in readout_tree instances
    // Free the resources allocated
    cleanup();

    return 0;
}

/////// HELPER FUNCTIONS ///////

bool init() {
    cl_int status;

    if(!setCwdToExeDir()) {
        return false;
    }

    // Get the OpenCL platform.
    platform = findPlatform("Altera");
    if(platform == NULL) {
        printf("ERROR: Unable to find Altera OpenCL platform.\n");
        return false;
    }

    // User-visible output - Platform information
    {
        char char_buffer[STRING_BUFFER_LEN]; 
        printf("Querying platform for info:\n");
        printf("==========================\n");
        clGetPlatformInfo(platform, CL_PLATFORM_NAME, STRING_BUFFER_LEN, char_buffer, NULL);
        printf("%-40s = %s\n", "CL_PLATFORM_NAME", char_buffer);
        clGetPlatformInfo(platform, CL_PLATFORM_VENDOR, STRING_BUFFER_LEN, char_buffer, NULL);
        printf("%-40s = %s\n", "CL_PLATFORM_VENDOR ", char_buffer);
        clGetPlatformInfo(platform, CL_PLATFORM_VERSION, STRING_BUFFER_LEN, char_buffer, NULL);
        printf("%-40s = %s\n\n", "CL_PLATFORM_VERSION ", char_buffer);
    }

    // Query the available OpenCL devices.
    scoped_array<cl_device_id> devices;
    cl_uint num_devices;

    devices.reset(getDevices(platform, CL_DEVICE_TYPE_ALL, &num_devices));

    // We'll just use the first device.
    device = devices[0];

    // Display some device information.
    display_device_info(device);

    // Create the context.
    context = clCreateContext(NULL, 1, &device, &oclContextCallback, NULL, &status);
    checkError(status, "Failed to create context");

    // Create the command queue.
    queue = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &status);
    checkError(status, "Failed to create command queue");

    // Create the program.
    std::string binary_file = getBoardBinaryFile("filtering_algorithm_top", device);
    printf("Using AOCX: %s\n", binary_file.c_str());
    program = createProgramFromBinary(context, binary_file.c_str(), &device, 1);

    // Build the program that was just created.
    status = clBuildProgram(program, 0, NULL, "", NULL, NULL);
    checkError(status, "Failed to build program");

    // Create the kernel - name passed in here must match kernel name in the
    // original CL file, that was compiled into an AOCX file using the AOC tool
    const char *kernel_name = "filtering_algorithm_top";  // Kernel name, as defined in the CL file
    kernel = clCreateKernel(program, kernel_name, &status);
    checkError(status, "Failed to create kernel");

    return true;
}

// Free the resources allocated during initialization
void cleanup() {
    if(kernel) {
        clReleaseKernel(kernel);  
    }
    if(program) {
        clReleaseProgram(program);
    }
    if(queue) {
        clReleaseCommandQueue(queue);
    }
    if(context) {
        clReleaseContext(context);
    }
}

// Helper functions to display parameters returned by OpenCL queries
static void device_info_ulong( cl_device_id device, cl_device_info param, const char* name) {
    cl_ulong a;
    clGetDeviceInfo(device, param, sizeof(cl_ulong), &a, NULL);
    printf("%-40s = %lu\n", name, a);
}
static void device_info_uint( cl_device_id device, cl_device_info param, const char* name) {
    cl_uint a;
    clGetDeviceInfo(device, param, sizeof(cl_uint), &a, NULL);
    printf("%-40s = %u\n", name, a);
}
static void device_info_bool( cl_device_id device, cl_device_info param, const char* name) {
    cl_bool a;
    clGetDeviceInfo(device, param, sizeof(cl_bool), &a, NULL);
    printf("%-40s = %s\n", name, (a?"true":"false"));
}
static void device_info_string( cl_device_id device, cl_device_info param, const char* name) {
    char a[STRING_BUFFER_LEN]; 
    clGetDeviceInfo(device, param, STRING_BUFFER_LEN, &a, NULL);
    printf("%-40s = %s\n", name, a);
}

// Query and display OpenCL information on device and runtime environment
static void display_device_info( cl_device_id device ) {

    printf("Querying device for info:\n");
    printf("========================\n");
    device_info_string(device, CL_DEVICE_NAME, "CL_DEVICE_NAME");
    device_info_string(device, CL_DEVICE_VENDOR, "CL_DEVICE_VENDOR");
    device_info_uint(device, CL_DEVICE_VENDOR_ID, "CL_DEVICE_VENDOR_ID");
    device_info_string(device, CL_DEVICE_VERSION, "CL_DEVICE_VERSION");
    device_info_string(device, CL_DRIVER_VERSION, "CL_DRIVER_VERSION");
    device_info_uint(device, CL_DEVICE_ADDRESS_BITS, "CL_DEVICE_ADDRESS_BITS");
    device_info_bool(device, CL_DEVICE_AVAILABLE, "CL_DEVICE_AVAILABLE");
    device_info_bool(device, CL_DEVICE_ENDIAN_LITTLE, "CL_DEVICE_ENDIAN_LITTLE");
    device_info_ulong(device, CL_DEVICE_GLOBAL_MEM_CACHE_SIZE, "CL_DEVICE_GLOBAL_MEM_CACHE_SIZE");
    device_info_ulong(device, CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE, "CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE");
    device_info_ulong(device, CL_DEVICE_GLOBAL_MEM_SIZE, "CL_DEVICE_GLOBAL_MEM_SIZE");
    device_info_bool(device, CL_DEVICE_IMAGE_SUPPORT, "CL_DEVICE_IMAGE_SUPPORT");
    device_info_ulong(device, CL_DEVICE_LOCAL_MEM_SIZE, "CL_DEVICE_LOCAL_MEM_SIZE");
    device_info_ulong(device, CL_DEVICE_MAX_CLOCK_FREQUENCY, "CL_DEVICE_MAX_CLOCK_FREQUENCY");
    device_info_ulong(device, CL_DEVICE_MAX_COMPUTE_UNITS, "CL_DEVICE_MAX_COMPUTE_UNITS");
    device_info_ulong(device, CL_DEVICE_MAX_CONSTANT_ARGS, "CL_DEVICE_MAX_CONSTANT_ARGS");
    device_info_ulong(device, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, "CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE");
    device_info_uint(device, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, "CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS");
    device_info_uint(device, CL_DEVICE_MEM_BASE_ADDR_ALIGN, "CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS");
    device_info_uint(device, CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE, "CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE");
    device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR");
    device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT");
    device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT");
    device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG");
    device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT");
    device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE");

    {
        cl_command_queue_properties ccp;
        clGetDeviceInfo(device, CL_DEVICE_QUEUE_PROPERTIES, sizeof(cl_command_queue_properties), &ccp, NULL);
        printf("%-40s = %s\n", "Command queue out of order? ", ((ccp & CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE)?"true":"false"));
        printf("%-40s = %s\n", "Command queue profiling enabled? ", ((ccp & CL_QUEUE_PROFILING_ENABLE)?"true":"false"));
    }
}

