/* ScaFES
 * Copyright (c) 2017, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file analyticalSolutions.hpp
 *  @brief Contains functions which can be used as analytical solutions.
 */

#ifndef ANALYTICAL_SOLUTIONS_HPP_
#define ANALYTICAL_SOLUTIONS_HPP_

/* Functions for constant analytical solutions. */
template<typename TT>
TT consFuncArg(void) {
    return 1.0;
}

template<typename TT>
TT consFuncTimeDerivative(void) {
    return consFuncArg<TT>();
}

template<typename TT>
TT consFuncSumOfSpaceDerivatives2ndOrder(void) {
    return 0.0;
}

template<typename TT>
TT consFunc(TT const& t) {
    return t*consFuncArg<TT>();
}

/* Functions for linear analytical solutions. */
template<typename TT, std::size_t DIM>
TT linFuncArg(ScaFES::Ntuple<TT,DIM> const& x) {
    TT argument = 1.0;
    for (std::size_t pp = 0; pp < DIM; ++pp) {
        argument += x[pp];
    }
    return argument;
}

template<typename TT, std::size_t DIM>
TT linFuncTimeDerivative(ScaFES::Ntuple<TT,DIM> const& x) {
    return linFuncArg<TT,DIM>(x);
}

template<typename TT>
TT linFuncSumOfSpaceDerivatives2ndOrder(void) {
    return 0.0;
}

template<typename TT, std::size_t DIM>
TT linFunc(ScaFES::Ntuple<TT,DIM> const& x, TT const& t) {
    return t*linFuncArg<TT,DIM>(x);
}

/* Functions for quadratic analytical solutions. */
template<typename TT, std::size_t DIM>
TT quadFuncArg(ScaFES::Ntuple<TT,DIM> const& x) {
    TT argument = 1.0;
    for (std::size_t pp = 0; pp < DIM; ++pp) {
        argument += x[pp];
        for (std::size_t ii = pp; ii < DIM; ++ii) {
            argument += x[pp]*x[ii];
        }
    }
    return argument;
}

template<typename TT, std::size_t DIM>
TT quadFuncTimeDerivative(ScaFES::Ntuple<TT,DIM> const& x) {
    return quadFuncArg<TT,DIM>(x);
}

template<typename TT, std::size_t DIM>
TT quadFuncSumOfSpaceDerivatives2ndOrder(ScaFES::Ntuple<TT,DIM> const& x,
                                         TT const& t) {
    TT sumSpaceDerivative = 0.0;
    TT spaceDerivative = 0.0;
    for (std::size_t pp = 0; pp < DIM; pp++) {
        spaceDerivative = 2.0;
        for (std::size_t ii = 0; ii < DIM; ++ii) {
            if (pp != ii) {
                spaceDerivative += 2.0*x[ii];
            } else {
                spaceDerivative += 6.0*x[ii];
            }
        }
        sumSpaceDerivative += t * spaceDerivative;
    }
    return sumSpaceDerivative;
}

template<typename TT, std::size_t DIM>
TT quadFunc(ScaFES::Ntuple<TT,DIM> const& x, TT const& t) {
    return t*quadFuncArg<TT,DIM>(x);
}

/* Functions for cubic analytical solutions. */
template<typename TT, std::size_t DIM>
TT cubicFuncArg(ScaFES::Ntuple<TT,DIM> const& x) {
    TT argument = 1.0;
    for (std::size_t pp = 0; pp < DIM; ++pp) {
        argument += x[pp];
        for (std::size_t ii = pp; ii < DIM; ++ii) {
            argument += x[pp]*x[ii];
        }
        for (std::size_t ii = 0; ii < DIM; ++ii) {
            argument += x[pp]*x[pp]*x[ii];
        }
        for (std::size_t ii = pp+1; ii < DIM; ++ii) {
            for (std::size_t jj = ii+1; jj < DIM; ++jj) {
                argument += x[pp]*x[ii]*x[jj];
            }
        }
    }
    return argument;
}

template<typename TT, std::size_t DIM>
TT cubicFuncTimeDerivative(ScaFES::Ntuple<TT,DIM> const& x) {
    return cubicFuncArg<TT,DIM>(x);
}

template<typename TT, std::size_t DIM>
TT cubicFunc(ScaFES::Ntuple<TT,DIM> const& x, TT const& t) {
    return t*cubicFuncArg<TT,DIM>(x);
}

#endif
