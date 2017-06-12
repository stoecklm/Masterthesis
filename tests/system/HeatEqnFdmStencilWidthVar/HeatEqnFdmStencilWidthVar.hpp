/*  Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
 *  For details, see the files COPYING and LICENSE in the base directory
 *  of the package.
 */

/**
 *  @file HeatEqnFdmStencilWidthVar.hpp
 *
 *  @brief System test: Contains the class HeatEqnFdmStencilWidthVar.
 */

#include <vector>
#include <string>

#include "ScaFES.hpp"

/**
 * \class HeatEqnFdmStencilWidthVar
 *  @brief System test: Scalar field with variable stencil width.
 *
 * System test.
 *
 * Initialisation:
 * \code $\varepsilon_{i,j,k}^0 = i+j+k \quad\forall \ i, j, k$\endcode
 *
 * Computation: Star with witdh "bw" or less if near to boundary.
 *
 * Verification:
 * \code$\varepsilon_{i,j,k}^n = i+j+k+2^n \quad\forall \ n= 0,1,\dots$\endcode
 *
 * Time step \code $t = 1, 2, \dots$\endcode:
 * \code
    \begin{flalign}
    v &= \sum\limits_{r=1}^{m_i} \varepsilon_{i-r,j,k}^n
      + \sum\limits_{r=1}^{n_i} \varepsilon_{i+r,j,k}^n
      + \sum\limits_{r=1}^{m_i} \varepsilon_{i,j-r,k}^n
      + \sum\limits_{r=1}^{n_i} \varepsilon_{i,j+r,k}^n
      + \sum\limits_{r=1}^{m_i} \varepsilon_{i,j,k-r}^n
      + \sum\limits_{r=1}^{n_i} \varepsilon_{i,j,k+r}^n \nonumber \\
     \varepsilon_{i,j,k}^{n+1} &= \underbrace{
     \frac{\varepsilon_{i,j,k}^n + v - cc}{num} }_{= \varepsilon_{i,j,k}^n} +
       \ 2.0
   \end{flalign}
   \endcode
 * \code $m_s$\endcode and \code$n_s$\endcode are the numbers of neighboured
 * grid nodes (in one given direction). The values at these neighboured grid
 * nodes will be added.
 * \code
   \begin{flalign*}
     m_s &:= \min(s, epsBw) \mbox{ f\"ur } s =i,j,k \\
     \mbox{and }  n_s &:= min(N_{ss} - s - 1, epsBw)
    \mbox{ for } s = i (ss=0), j (ss=1), k (ss=2)
    \end{flalign*}
    \endcode
 *
 *  The value \code$\varepsilon_{i,j,k}^{n+1}$\endcode is the sum of the
 *  values at all neighboured nodes which are in a circle with a radius
 *  smaller than of equal to the given stencil width.
 *  If the distance of a grid node to the global boundary is less than the
 *  given stencil width, then the sum must be corrected by \code$cc$\endcode.
 */
template<typename CT, std::size_t DIM>
class HeatEqnFdmStencilWidthVar
    : public ScaFES::Problem<HeatEqnFdmStencilWidthVar<CT, DIM>, CT, DIM>
{
    public:
        HeatEqnFdmStencilWidthVar(ScaFES::Parameters const& params,
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
            : ScaFES::Problem<HeatEqnFdmStencilWidthVar<CT, DIM>, CT, DIM>(
                     params, gg, useLeapfrog, nameDatafield, stencilWidth,
                     isKnownDf, nLayers, defaultValue, writeToFile,
                     computeError, geomparamsInit) { }
        void evalInner(std::vector< ScaFES::DataField<CT, DIM> >& vNew,
                       ScaFES::Ntuple<int,DIM> const& idxNode,
                       int const& timestep) {
            for (std::size_t ii = 0; ii < vNew.size(); ++ii) {
                vNew[ii](idxNode) = 2.0 * timestep;
                for (std::size_t pp = 0; pp < DIM; ++pp) {
                    vNew[ii](idxNode) += idxNode[pp];
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
                // How many grid nodes are between the current global node
                // and the boundary?
                ScaFES::Ntuple<int,DIM> diff2Start
                    = idxNode - vOld[ii].gridGlobal().discreteDomain().idxNodeFirst();
                ScaFES::Ntuple<int,DIM> diff2End
                    = vOld[ii].gridGlobal().discreteDomain().idxNodeLast() - idxNode;

                int epsBw = vOld[ii].stencilWidth();
                TT sumValuesAdded = 0.0;
                // Correcture term: +1+2...+bw and -1-2-...
                // are not equal at boundary nodes.
                int cc = 0;
                int nValuesAdded = 1;
                int minStartBw;
                int minEndBw;
                ScaFES::Ntuple<int,DIM> tmp;
                ScaFES::Ntuple<int,DIM> tmp2;
                for (std::size_t kk = 0; kk < DIM; ++kk) {
                    // To the left, front, bottom...
                    minStartBw = std::min(diff2Start[kk], epsBw);
                    for (int rr = 1; rr <= minStartBw; ++rr) {
                        tmp = idxNode;
                        tmp[kk] -= rr;
                        sumValuesAdded += vOld[ii](tmp);
                        cc -= rr;
                    }
                    // To the right, back, top...
                    minEndBw = std::min(diff2End[kk], epsBw);
                    for (int rr = 1; rr <= minEndBw; ++rr) {
                        tmp = idxNode;
                        tmp[kk] += rr;
                        sumValuesAdded += vOld[ii](tmp);
                        cc += rr;
                    }
                    nValuesAdded += (minStartBw + minEndBw);
                }
                vNew[ii](idxNode) =  2.0
                              + (sumValuesAdded + vOld[ii](idxNode) - cc)
                                  / nValuesAdded;
            }
        }

        template<typename TT>
        void updateBorder(std::vector< ScaFES::DataField<TT, DIM> >& vNew,
                          std::vector< ScaFES::DataField<TT, DIM> > const& vOld,
                          ScaFES::Ntuple<int,DIM> const& idxNode,
                          int const& timestep) {
            this->updateInner(vNew, vOld, idxNode, timestep);
        }

        template<typename TT>
        void updateInner2(std::vector< ScaFES::DataField<TT,DIM> >& /*vNew*/,
                          std::vector< ScaFES::DataField<TT,DIM> > const& /*vOld*/,
                          ScaFES::Ntuple<int,DIM> const& /*idxNode*/,
                          int const& /*timestep*/) { }

        template<typename TT>
        void updateBorder2(std::vector< ScaFES::DataField<TT,DIM> >& /*vNew*/,
                           std::vector< ScaFES::DataField<TT,DIM> > const& /*vOld*/,
                           ScaFES::Ntuple<int,DIM> const& /*idxNode*/,
                           int const& /*timestep*/) { }
};
