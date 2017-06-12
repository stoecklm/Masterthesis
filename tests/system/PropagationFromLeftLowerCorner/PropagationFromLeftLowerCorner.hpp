/*  Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
 *  For details, see the files COPYING and LICENSE in the base directory
 *  of the package.
 */

/**
 *  @file PropagationFromLeftLowerCorner.hpp
 *
 *  @brief System test: Contains the class PropagationFromLeftLowerCorner.
 */

#include <vector>
#include <string>

#include "ScaFES.hpp"

const double EPS = 2.2e-2;
const double DEFAULT_VALUE = -25.0;

/**
 * \class PropagationFromLeftLowerCorner
 *  @brief System test: Propagation from left lower corner through grid.
 *
 * System test.
 *
 *  Scalar field with border width of 1.
 *  \n Initialisation:
    \code\begin{flalign*}
      d &= 15,
      u_{i,j,k}^0  &= d & \forall \ i = 1, 2, \dots, N_0-1, &\\
                    &                            & j = 1, 2, \dots, N_1-1, \\
                    &                            & k = 1, 2, \dots, N_2-1  \\
      \text{und } &u_{0,0,0}^0  = 1.0
    \end{flalign*}\endcode
 *
 *  Computation:
 *  \code
    \begin{flalign*}
      u_{i,j,k}^{n+1} = \begin{cases}
                    n+1, & \mbox{if } u_{i-1,j,k} \ne d, \\
                    n+1, & \mbox{if } u_{i,j-1,k} \ne d, \\
                    n+1, & \mbox{if } u_{i,j,k-1} \ne d, \\
                    u_{i,j,k}^{n}  & \mbox{else.}
                  \end{cases}
    \end{flalign*}
    \endcode
 *
 *  Verification:
 *  \code$u_{i,j,k}^n = i+j+k\quad \mbox{for } \ i+j+k \le n$\endcode
 */
template<typename CT, std::size_t DIM>
class  PropagationFromLeftLowerCorner
    : public ScaFES::Problem< PropagationFromLeftLowerCorner<CT, DIM>, CT, DIM>
{
    public:
        PropagationFromLeftLowerCorner(ScaFES::Parameters const& params,
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
            : ScaFES::Problem<PropagationFromLeftLowerCorner<CT, DIM>, CT, DIM>(
                     params, gg, useLeapfrog, nameDatafield, stencilWidth,
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
                int sumIndex = 0.0;
                for (std::size_t pp = 0; pp < DIM; ++pp) {
                    sumIndex += idxNode[pp];
                }
                if (sumIndex > timestep) {
                    vNew[ii](idxNode) = DEFAULT_VALUE;
                }
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
                vNew[ii](idxNode) = DEFAULT_VALUE;
                ScaFES::Ntuple<double,DIM> x = this->coordinates(idxNode);
                if (x < EPS) {
                    vNew[ii](idxNode) = 0.0;
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
                         int const& timestep) {
            for (std::size_t ii = 0; ii < vNew.size(); ++ii) {
                CT defvalCmpx = DEFAULT_VALUE;
                bool updateUnew = false;
                for (std::size_t pp = 0; pp < DIM; ++pp) {
                    if ( (fabs(defvalCmpx
                          - vOld[ii](this->connect(idxNode,2*pp))) >= EPS) ) {
                        updateUnew = true;
                        break;
                    }
                }
                if (updateUnew && (fabs(defvalCmpx
                      - vOld[ii](idxNode)) < EPS) ) {
                    vNew[ii](idxNode) = static_cast<double>(timestep);
                }
                else {
                    vNew[ii](idxNode) = vOld[ii](idxNode);
                }
            }
        }

        template<typename TT>
        void updateBorder(std::vector< ScaFES::DataField<TT, DIM> >& vNew,
                          std::vector< ScaFES::DataField<TT, DIM> > const& vOld,
                          ScaFES::Ntuple<int,DIM> const& idxNode,
                          int const& timestep) {
            for (std::size_t ii = 0; ii < vNew.size(); ++ii) {
                TT defvalCmpx = DEFAULT_VALUE;
                TT cmplx = defvalCmpx;
                bool updateUnew = false;
                for (std::size_t pp = 0; pp < DIM; ++pp) {
                    if (vOld[ii].gridGlobal().discreteDomain().idxNodeFirst(pp)
                          < idxNode.elem(pp)) {
                        cmplx = vOld[ii](this->connect(idxNode, 2*pp));
                    }
                    else {
                        cmplx = defvalCmpx;
                    }
                    if ( (fabs(defvalCmpx - cmplx) >= EPS) ) {
                        updateUnew = true;
                        break;
                    }
                }
                if (updateUnew
                    && (fabs(defvalCmpx - vOld[ii](idxNode)) < EPS) ) {
                    vNew[ii](idxNode) = static_cast<double>(timestep);
                }
                else {
                    vNew[ii](idxNode) = vOld[ii](idxNode);
                }
            }
        }

       template<typename TT>
       void updateInner2(std::vector<ScaFES::DataField<TT,DIM>>& ,
                        std::vector<ScaFES::DataField<TT,DIM>> const& ,
                        ScaFES::Ntuple<int,DIM> const& ,
                       int const& ) {}

       template<typename TT>
       void updateBorder2(std::vector<ScaFES::DataField<TT,DIM>>&,
                         std::vector<ScaFES::DataField<TT,DIM>>const&,
                         ScaFES::Ntuple<int,DIM> const&,
                         int const&) {}

};
