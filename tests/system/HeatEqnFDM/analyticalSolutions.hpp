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
TT constantFunctionArgument(void) {
    return 1.0;
}

template<typename TT>
TT constantFunctionTimeDerivative(void) {
    return constantFunctionArgument<TT>();
}

template<typename TT>
TT constantFunctionSpaceDerivatives(void) {
    return 0.0;
}

template<typename TT>
TT constantFunction(TT const& t) {
    return t*constantFunctionArgument<TT>();
}

/* Functions for linear analytical solutions. */
template<typename TT, std::size_t DIM>
TT linearFunctionArgument(ScaFES::Ntuple<TT,DIM> const& x) {
    TT argument = 1.0;
    for (std::size_t pp = 0; pp < DIM; ++pp) {
        argument += x[pp];
    }
    return argument;
}

template<typename TT, std::size_t DIM>
TT linearFunctionTimeDerivative(ScaFES::Ntuple<TT,DIM> const& x) {
    return linearFunctionArgument<TT,DIM>(x);
}

template<typename TT>
TT linearFunctionSpaceDerivatives(void) {
    return 0.0;
}

template<typename TT, std::size_t DIM>
TT linearFunction(ScaFES::Ntuple<TT,DIM> const& x, TT const& t) {
    return t*linearFunctionArgument<TT,DIM>(x);
}

/* Functions for quadratic analytical solutions. */
template<typename TT, std::size_t DIM>
TT quadraticFunctionArgument(ScaFES::Ntuple<TT,DIM> const& x) {
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
TT quadraticFunctionTimeDerivative(ScaFES::Ntuple<TT,DIM> const& x) {
    return quadraticFunctionArgument<TT,DIM>(x);
}

template<typename TT, std::size_t DIM>
TT quadraticFunctionSpaceDerivatives(ScaFES::Ntuple<TT,DIM> const& x, TT const& t) {
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
TT quadraticFunction(ScaFES::Ntuple<TT,DIM> const& x, TT const& t) {
    return t*quadraticFunctionArgument<TT,DIM>(x);
}

/* Functions for cubic analytical solutions. */
template<typename TT, std::size_t DIM>
TT cubicFunctionArgument(ScaFES::Ntuple<TT,DIM> const& x) {
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
TT cubicFunctionTimeDerivative(ScaFES::Ntuple<TT,DIM> const& x) {
    return cubicFunctionArgument<TT,DIM>(x);
}

template<typename TT, std::size_t DIM>
TT cubicFunction(ScaFES::Ntuple<TT,DIM> const& x, TT const& t) {
    return t*cubicFunctionArgument<TT,DIM>(x);
}

#endif
