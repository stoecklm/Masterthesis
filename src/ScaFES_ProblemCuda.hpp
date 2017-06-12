/* ScaFES
 * Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file ScaFES_ProblemCuda.hpp
 *  @brief Contains methods related to CUDA and the class template Problem.
 */

#ifndef SCAFES_PROBLEMCUDA_HPP_
#define SCAFES_PROBLEMCUDA_HPP_

#include <vector>

#ifdef SCAFES_HAVE_CUDA
#include <cuda.h>
#include<cuda_runtime_api.h>
#endif

#define BLOCK_SIZE 32

/** Get the error code from cuda */
void check_cuda_errors(const char* filename, const int line_number);

/** Computes method \c updateInner using CUDA. */
void updateGlobalInnerUsingCudaVect(short int* const& kindElemData,
                                    std::vector<double*>& dfDepElemData,
                                    std::vector<double*> const& dfIndepElemData,
                                    int const& timestep,
                                    std::vector<int> const& nodeType,
                                    int const& nNodesTotalPartition);

int host_kernelUpdateInnerVect(
           dim3 blocksPerGrid,
           dim3 threadsPerBlock,
    short int* const kindElemData, std::vector<double*> dfDepElemData,
    std::vector<double*> const dfIndepElemData, int const timeIter,
    std::vector<int> const setGridNodes, int const nNodesTotalPartition);


/*******************************************************************************
*******************************************************************************/
inline void check_cuda_errors(const char* filename, const int line_number)
{
    cudaError_t error = cudaGetLastError();

    if (error != cudaSuccess)
    {
        printf("CUDA error at %s:%i: %s\n", filename, line_number,
               cudaGetErrorString(error));
        exit(-1);
    }
}

/*******************************************************************************
*******************************************************************************/
inline void getCudaInfo()
{
    int devicenr = 0;

    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, devicenr);
    printf("NVIDIA -Card Properties\n");
    printf(" name: %s\n", prop.name);
    printf(" global memory: %d MByte(s)\n", prop.totalGlobalMem / 1024 / 1024);
    printf(" shared memory: %d kByte(s)\n", prop.sharedMemPerBlock / 1024);
    printf(" register per multiprocessor: %d\n", prop.regsPerBlock);
    printf(" warp size: %d threads\n", prop.warpSize);
    printf(" max. threads per block: %d\n", prop.maxThreadsPerBlock);
    printf(" max. threadblock size: %dx%dx%d\n", prop.maxThreadsDim[0],
           prop.maxThreadsDim[1], prop.maxThreadsDim[2]);
    printf(" max. blockgrid size: %dx%dx%d\n", prop.maxGridSize[0],
           prop.maxGridSize[1], prop.maxGridSize[2]);
    printf(" constant memory: %d Byte(s)\n", prop.totalConstMem);
    printf(" compute capability: %d.%d\n", prop.major, prop.minor);
    printf(" clock rate: %1.0f MHz\n", (float)prop.clockRate / 1000.0);
    printf(" multiprocessor count: %d\n\n", prop.multiProcessorCount);
}

/*******************************************************************************
*******************************************************************************/
void updateGlobalInnerUsingCudaVect(short int* const& kindElemData,
                                    std::vector<double*>& dfDepElemData,
                                    std::vector<double*> const& dfIndepElemData,
                                    std::vector<int> const& setGridNodes,
                                    int const& timeIter,
                                    int const& nNodesTotalPartition)
{
    /*---------------------------------------------------------------------------
    | dimBlock = number of threads per block
    | dimGrid = number of blocks
    ---------------------------------------------------------------------------*/
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE, 1);
    dim3 dimGrid(ceil((nNodesTotalPartition + dimBlock.x - 1) / dimBlock.x),
                 ceil((nNodesTotalPartition + dimBlock.y - 1) / dimBlock.y));

    host_kernelUpdateInnerVect(dimGrid, dimBlock, kindElemData,
                              dfDepElemData, dfIndepElemData, timeIter,
                              setGridNodes, nNodesTotalPartition);

    cudaThreadSynchronize();
    check_cuda_errors(__FILE__, __LINE__);
}

// /*******************************************************************************
// *******************************************************************************/
// TODO: Call of updateGlobalInnerUsingCuda() within
// Problem<OWNPRBLM, CT, DIM>::evalTypeTwoAt():
// [...]
// std::vector<double*> dfDepElemData;
// std::vector<double*> dfIndepElemData;
// for (std::size_t ii = 0; ii < dfIndep.size(); ++ii)
// {    dfIndepElemData.push_back(dfIndep.data());
// }
// for (std::size_t ii = 0; ii < dfDep.size(); ++ii)
// {    dfDepElemData.push_back(dfDep.data());
// }
// int nodeType = 0;
// int nNodesTotalPartition = setGridNodes.size();
// updateGlobalInnerUsingCudaVect(this->mKind.data(),
//                                dfDepElemData,
//                                dfIndepElemData,
//                                timeIter,
//                                setGridNodes,
//                                nNodesTotalPartition);
// [...]
//

#endif
