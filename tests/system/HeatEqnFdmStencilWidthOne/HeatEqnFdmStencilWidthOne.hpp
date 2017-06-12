/*  Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
 *  For details, see the files COPYING and LICENSE in the base directory
 *  of the package.
 */

/**
 *  @file HeatEqnFdmStencilWidthOne.hpp
 *
 *  @brief System test: Contains the class HeatEqnFdmStencilWidthOne.
 */

#include <vector>
#include <string>

#include "ScaFES.hpp"

/**
 * \class HeatEqnFdmStencilWidthOne
 *  @brief System test: Scalar field with stencil width of 1.
 *
 * System test.
 *
 * Initialisation:
 * \code$f_{i,j,k}^0 = i+j+k \quad \forall \ i, j, k$\endcode
 *
 * Computation:
 *  \code
    \begin{align*}
       f_{i,j,k}^{n+1} = \frac{2}{7} \left(  f_{i,j,k}^n
        + f_{i-1,j,k}^n + f_{i+1,j,k}^n + f_{i,j+1,k}^n + f_{i,j-1,k}^n
                                + f_{i,j,k-1}^n + f_{i,j,k+1}^n \right)  \\
                                \forall \ i = 1, 2, \dots, N_0-2 \nonumber \\
                                        \ j = 1, 2, \dots, N_1-2 \nonumber \\
                                        \ k = 1, 2, \dots, N_2-2 \nonumber
    \end{align*}
    \endcode
 *
 * Verification:
 * \code$f_{i,j,k}^n = 2^n(i+j+k)\quad \forall \ n= 0, 1, \dots, T$\endcode
 *
 * Function
 * \code
   \begin{align}
   fana(x,y,z,t) &:= x/hx + y/hy + z/hz + 2 t/tau \\
   ==> fana(x_i,y_j,z_k,t_l) &= x_i/h_0 + y_j/h_1 + z_k/h_2 + 2 t_l/tau \\
                           &= i h_0/h_0 + j h_1/h_1 + k h_2/h_2 + 2 l tau/tau\\
                           &= i + j + k + 2 l.
   \end{align}
 * \endcode
 *
 * \code $f(locNode)$\endcode is calculated as the arithmetic average of
 * the values of the six direct neighbour nodes and itself.
 * If the exchange works fine, $f(locNode)$ is constant for all time steps.
 *
 * Simplest case at the boundary: Each node needs the values at one neighbour
 * node:
 * \code t = 1: f[i,j,k] = (i+j+k) + (i+1+j+k) - 1 = 2*(i,j,k)\endcode
 * More complicated case at the boundary:
 * Each node needs the values at all direct neighbour nodes,
 * corner nodes will need the values at three neighbour nodes.
 */
template<typename CT, std::size_t DIM>
class HeatEqnFdmStencilWidthOne
    : public ScaFES::Problem<HeatEqnFdmStencilWidthOne<CT, DIM>, CT, DIM>
{
    public:
        HeatEqnFdmStencilWidthOne(ScaFES::Parameters const& cl,
                                  ScaFES::GridGlobal<DIM> const& gg,
            bool useLeapfrog,
            std::vector<std::string> const& nameDatafield,
            std::vector<int> const& stencilWidth,
            std::vector<bool> const& isKnownDf,
            std::vector<int> const& nLayers = std::vector<int>(),
            std::vector<CT> const& defaultValue = std::vector<CT>(),
            std::vector<ScaFES::WriteHowOften> const& writeToFile
                = std::vector<ScaFES::WriteHowOften>(),
            std::vector<bool> const& computeError = std::vector<bool>(),
            std::vector<CT> const& geomparamsInit = std::vector<CT>() )
            : ScaFES::Problem<HeatEqnFdmStencilWidthOne<CT, DIM>, CT, DIM>(
                     cl, gg, useLeapfrog, nameDatafield, stencilWidth,
                     isKnownDf, nLayers, defaultValue, writeToFile,
                     computeError, geomparamsInit) { }

        void evalInner(std::vector< ScaFES::DataField<CT, DIM> >& vNew,
                      ScaFES::Ntuple<int,DIM> const& idxNode,
                      int const& timestep) {
            for (std::size_t ii = 0; ii < vNew.size(); ++ii) {
                vNew[ii](idxNode) = 0.0;
                for (std::size_t pp = 0; pp < DIM; ++pp) {
                    vNew[ii](idxNode) += idxNode[pp];
                }
                vNew[ii](idxNode) *= pow(2.0, timestep);
            }
        }
        void evalBorder(std::vector< ScaFES::DataField<CT, DIM> >& vNew,
                      ScaFES::Ntuple<int,DIM> const& idxNode,
                      int const& timestep) {
            this->evalInner(vNew, idxNode, timestep);
        }

        template<typename TT>
        void initInner(std::vector< ScaFES::DataField<TT, DIM> >& vNew,
                       std::vector<TT> const& /*vOld*/,
                       ScaFES::Ntuple<int,DIM> const& idxNode,
                       int const& /*timestep*/) {
            for (std::size_t ii = 0; ii < vNew.size(); ++ii) {
                vNew[ii](idxNode) = 0.0;
                for (std::size_t pp = 0; pp < DIM; ++pp) {
                    vNew[ii](idxNode) += idxNode[pp];
                }
            }
        }
        template<typename TT>
        void initBorder(std::vector< ScaFES::DataField<TT, DIM> >& vNew,
                        std::vector<TT> const& vOld,
                        ScaFES::Ntuple<int,DIM> const& idxNode,
                        int const& timestep) {
            this->template initInner<TT>(vNew, vOld, idxNode, timestep);
        }

        template<typename TT>
        void updateInner(std::vector< ScaFES::DataField<TT, DIM> >& vNew,
                         std::vector< ScaFES::DataField<TT, DIM> > const& vOld,
                         ScaFES::Ntuple<int,DIM> const& idxNode,
                         int const& /*timestep*/) {
            for (std::size_t ii = 0; ii < vNew.size(); ++ii) {
                vNew[ii](idxNode) = vOld[ii](idxNode);
                for (std::size_t pp = 0; pp < 2*DIM; ++pp) {
                    vNew[ii](idxNode) += vOld[ii](this->connect(idxNode,pp));
                }
                vNew[ii](idxNode) *= ( 2.0 / (2 * DIM + 1) );
            }
        }

        template<typename TT>
        void updateBorder(std::vector< ScaFES::DataField<TT, DIM> >& vNew,
                          std::vector< ScaFES::DataField<TT, DIM> > const& vOld,
                          ScaFES::Ntuple<int,DIM> const& idxNode,
                          int const& timestep) {
            for (std::size_t ii = 0; ii < vNew.size(); ++ii) {
                // ATTENTION: timestep must be traced, too!
                TT correction = pow(2.0, timestep - 1);
                for (std::size_t pp = 0; pp < 2 * DIM; ++pp) {
                    // ATTENTION: The following if condition has to be rewritten!
                    // Otherwise, ADOL-C mode will not work!
                    // At boundary?
                    if (ScaFES::Ntuple<int,DIM>(-(pp+1))
                        == this->connect(idxNode,pp)) {
                        // LEFT, FRONT, BOTTOM...
//                  if (0 == pp % 2) {
//                      vNew[ii](idxNode) = vOld[ii](idxNode)
//                              + vOld[ii](this->connect(idxNode,pp+1))
//                              - correction;
//                 }
//                 // RIGHT, BACK, TOP...
//                 else {
//                      vNew[ii](idxNode) = vOld[ii](idxNode)
//                              + vOld[ii](this->connect(idxNode,pp-1))
//                              + correction;
//                 }
                        vNew[ii](idxNode) = vOld[ii](idxNode)
                         + vOld[ii](this->connect(idxNode,pp+pow(-1.0,pp)))
                          - pow(-1.0,pp) * correction;
                    }
                }
            }
        }

       template<typename TT>
       void updateInner2(std::vector<ScaFES::DataField<TT,DIM>>&,
                         std::vector<ScaFES::DataField<TT,DIM>> const&,
                         ScaFES::Ntuple<int,DIM> const&,
                         int const&) {}

       template<typename TT>
       void updateBorder2(std::vector<ScaFES::DataField<TT,DIM>>&,
                          std::vector<ScaFES::DataField<TT,DIM>>const&,
                          ScaFES::Ntuple<int,DIM> const&,
                          int const&) {}
};
