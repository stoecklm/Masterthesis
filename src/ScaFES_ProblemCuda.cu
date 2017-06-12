#include <vector>

/*******************************************************************************
*******************************************************************************/
__global__ void kernelUpdateInnerVect(
    short int* const kindElemData, std::vector<double*> dfDepElemData,
    std::vector<double*> const dfIndepElemData, int const timeIter,
    std::vector<int> const setGridNodes, int const nNodesTotalPartition)
{
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int ii = col * nNodesTotalPartition + row;
    // TODO: Index globNode = ...(ii);

    if (ii < nNodesTotalPartition)
    {
        // Call of OwnProblem::updateInner(...);
    }
}

/*******************************************************************************
*******************************************************************************/
int host_kernelUpdateInnerVect(
           dim3 blocksPerGrid,
           dim3 threadsPerBlock,
    short int* const kindElemData, std::vector<double*> dfDepElemData,
    std::vector<double*> const dfIndepElemData, int const timeIter,
    std::vector<int> const setGridNodes, int const nNodesTotalPartition)
{
    kernelUpdateInnerVect<<<blocksPerGrid,threadsPerBlock>>>(kindElemData,
                                                             dfDepElemData,
                                                             dfIndepElemData,
                                                             timeIter,
                                                             setGridNodes,
                                                             nNodesTotalPartition);
    return EXIT_SUCCESS;
}

// /*******************************************************************************
// *******************************************************************************/
// template<class OWNPRBLM, std::size_t DIM, typename CT>
// __global__ void kernelUpdateInnerCuda(DataField<short
// int,DIM> const& mKind,
//                                   std::vector< DataField<CT,DIM> >& dfDep,
//                                   std::vector< DataField<CT,DIM> > const& dfIndep,
//                                   int const& timeIter,
//                                   std::vector<int> const& setGridNodes,
//                                   int const& nNodesTotalPartition)
// {
//    int col = blockIdx.x * blockDim.x + threadIdx.x;
//    int row = blockIdx.y * blockDim.y + threadIdx.y;
//    int ii = col * nNodesTotalPartition + row;
//
//    if (ii < nNodesTotalPartition) {
//      Ntuple<int, DIM> idxNode = mKind.memoryPos2IdxNode(ii);
//      static_cast<OWNPRBLM*>(this)->template updateInner<CT>(dfDep, dfIndep,
//                                                             idxNode, timeIter);
//    }
// }
//
// /*******************************************************************************
// *******************************************************************************/
// template<class OWNPRBLM, std::size_t DIM>
// template<typename CT>
// inline void Problem<OWNPRBLM,
// DIM>::updateGlobalInnerUsingCuda(DataField<short int,DIM> const& mKind,
//                                   std::vector< DataField<CT,DIM> >& dfDep,
//                                   std::vector< DataField<CT,DIM> > const&
// dfIndep,
//                                   int const& timeIter,
//                                   std::vector<int> const& setGridNodes,
//                                   int const& nNodesTotalPartition)
// {
//    /*---------------------------------------------------------------------------
//    | dimBlock = number of threads per block
//    | dimGrid = number of blocks
//    ---------------------------------------------------------------------------*/
//    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE, 1);
//    dim3 dimGrid(ceil((nNodesTotalPartition + dimBlock.x - 1) / dimBlock.x),
//                 ceil((nNodesTotalPartition + dimBlock.y - 1) / dimBlock.y));
//
//    kernelUpdateInner<<<dimGrid,dimBlock>>>(kind, dfDep, dfIndep,
//                                            timeIter,
//                                            nodeType,
//                                            nNodesTotalPartition);
//    cudaThreadSynchronize();
//    check_cuda_errors(__FILE__, __LINE__);
// }
//
