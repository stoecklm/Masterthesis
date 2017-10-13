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

/** Calculates value of argument (value inside parenthesis)
 *  of a constant function for given grid point.
 *  @param x Coordinates for this grid point.
 */
template<typename TT, std::size_t DIM>
TT consFuncArg(ScaFES::Ntuple<TT,DIM> const& /*x*/ ) {
    return 1.0;
}

/** Calculates value of time derivative of a constant function
 *  for given grid point.
 *  @param x Coordinates for this grid point.
 */
template<typename TT, std::size_t DIM>
TT consFuncTimeDerivative(ScaFES::Ntuple<TT,DIM> const& /*x*/) {
    return consFuncArg<TT>();
}

/** Calculates value of first order space derivative of a constant function
 *  for given grid point, time and orientation.
 *  @param x Coordinates for this grid point.
 *  @param t Absolute time.
 *  @param currDIM Current dimension (orientation of points).
 */
template<typename TT, std::size_t DIM>
TT consFuncSpaceDerivative1stOrder(ScaFES::Ntuple<TT,DIM> const& /*x*/,
                                   TT const& /*t*/, std::size_t /*currDIM*/) {
    return 0.0;
}

/** Calculates sum of second order space derivatives for all dimensions
 *  of a constant function for given grid point and time.
 *  @param x Coordinates for this grid point.
 *  @param t Absolute time.
 */
template<typename TT, std::size_t DIM>
TT consFuncSumOfSpaceDerivatives2ndOrder(ScaFES::Ntuple<TT,DIM> const& /*x*/,
                                         TT const& /*t*/) {
    return 0.0;
}

/** Calculates value of a constant function for given grid point and time.
 *  @param x Coordinates for this grid point.
 *  @param t Absolute time.
 */
template<typename TT, std::size_t DIM>
TT consFunc(ScaFES::Ntuple<TT,DIM> const& /*x*/, TT const& t) {
    return t*consFuncArg<TT>();
}

/** Calculates value of time derivative of a linear function
 *  for given grid point.
 *  @param x Coordinates for this grid point.
 */
template<typename TT, std::size_t DIM>
TT linFuncArg(ScaFES::Ntuple<TT,DIM> const& x) {
    TT argument = 1.0;
    for (std::size_t pp = 0; pp < DIM; ++pp) {
        argument += x[pp];
    }
    return argument;
}

/** Calculates value of time derivative of a linear function
 *  for given grid point.
 *  @param x Coordinates for this grid point.
 */
template<typename TT, std::size_t DIM>
TT linFuncTimeDerivative(ScaFES::Ntuple<TT,DIM> const& x) {
    return linFuncArg<TT,DIM>(x);
}

/** Calculates value of first order space derivative of a linear function
 *  for given grid point, time and orientation.
 *  @param x Coordinates for this grid point.
 *  @param t Absolute time.
 *  @param currDIM Current dimension (orientation of points).
 */
template<typename TT, std::size_t DIM>
TT linFuncSpaceDerivative1stOrder(ScaFES::Ntuple<TT,DIM> const& /*x*/,
                                  TT const& t, std::size_t /*currDIM*/) {
    return t;
}

/** Calculates sum of second order space derivatives for all dimensions
 *  of a linear function for given grid point and time.
 *  @param x Coordinates for this grid point.
 *  @param t Absolute time.
 */
template<typename TT, std::size_t DIM>
TT linFuncSumOfSpaceDerivatives2ndOrder(ScaFES::Ntuple<TT,DIM> const& /*x*/,
                                        TT const& /*t*/) {
    return 0.0;
}

/** Calculates value of a linear function for given grid point and time.
 *  @param x Coordinates for this grid point.
 *  @param t Absolute time.
 */
template<typename TT, std::size_t DIM>
TT linFunc(ScaFES::Ntuple<TT,DIM> const& x, TT const& t) {
    return t*linFuncArg<TT,DIM>(x);
}

/** Calculates value of time derivative of a quadratic function
 *  for given grid point.
 *  @param x Coordinates for this grid point.
 */
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

/** Calculates value of time derivative of a quadratic function
 *  for given grid point.
 *  @param x Coordinates for this grid point.
 */
template<typename TT, std::size_t DIM>
TT quadFuncTimeDerivative(ScaFES::Ntuple<TT,DIM> const& x) {
    return quadFuncArg<TT,DIM>(x);
}

/** Calculates value of first order space derivative of a quadratic function
 *  for given grid point, time and orientation.
 *  @param x Coordinates for this grid point.
 *  @param t Absolute time.
 *  @param currDIM Current dimension (orientation of points).
 */
template<typename TT, std::size_t DIM>
TT quadFuncSpaceDerivative1stOrder(ScaFES::Ntuple<TT,DIM> const& x,
                                   TT const& t, std::size_t currDIM) {
    TT spaceDerivative = 1.0;
    for (std::size_t pp = 0; pp < DIM; pp++) {
        if (pp == currDIM) {
            spaceDerivative += 2.0*x[pp];
        } else {
            spaceDerivative += x[pp];
        }
    }
    return t*spaceDerivative;
}

/** Calculates sum of second order space derivatives for all dimensions
 *  of a quadratic function for given grid point and time.
 *  @param x Coordinates for this grid point.
 *  @param t Absolute time.
 */
template<typename TT, std::size_t DIM>
TT quadFuncSumOfSpaceDerivatives2ndOrder(ScaFES::Ntuple<TT,DIM> const& /*x*/,
                                         TT const& t) {
    return DIM*2.0*t;
}

/** Calculates value of a quadratic function for given grid point and time.
 *  @param x Coordinates for this grid point.
 *  @param t Absolute time.
 */
template<typename TT, std::size_t DIM>
TT quadFunc(ScaFES::Ntuple<TT,DIM> const& x, TT const& t) {
    return t*quadFuncArg<TT,DIM>(x);
}

/** Calculates value of time derivative of a cubic function
 *  for given grid point.
 *  @param x Coordinates for this grid point.
 */
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

/** Calculates value of time derivative of a cubic function
 *  for given grid point.
 *  @param x Coordinates for this grid point.
 */
template<typename TT, std::size_t DIM>
TT cubicFuncTimeDerivative(ScaFES::Ntuple<TT,DIM> const& x) {
    return cubicFuncArg<TT,DIM>(x);
}

/** Calculates value of first order space derivative of a cubic function
 *  for given grid point, time and orientation.
 *  @param x Coordinates for this grid point.
 *  @param t Absolute time.
 *  @param currDIM Current dimension (orientation of points).
 */
template<typename TT, std::size_t DIM>
TT cubicFuncSpaceDerivative1stOrder(ScaFES::Ntuple<TT,DIM> const& x,
                                    TT const& t, std::size_t currDIM) {
    TT spaceDerivative = 1.0;
    for (std::size_t pp = 0; pp < DIM; pp++) {
        if (pp == currDIM) {
            spaceDerivative += 2.0*x[pp];
        } else {
            spaceDerivative += x[pp];
        }
    }
    for (std::size_t pp = 0; pp < DIM; pp++) {
        if (pp == currDIM) {
            spaceDerivative += 3.0*x[pp]*x[pp];
        } else {
            spaceDerivative += x[pp]*x[pp];
        }
    }
    for (std::size_t pp = 0; pp < DIM; pp++) {
        if (pp != currDIM) {
            spaceDerivative += 2.0*x[currDIM]*x[pp];
        }
    }
    if (DIM == 3) {
        if (currDIM == 0) {
            spaceDerivative += x[1]*x[2];
        } else if (currDIM == 1) {
            spaceDerivative += x[0]*x[2];
        } else /* (currDIM == 2) */ {
            spaceDerivative += x[0]*x[1];
        }
    }
    return t*spaceDerivative;
}

/** Calculates sum of second order space derivatives for all dimensions
 *  of a cubic function for given grid point and time.
 *  @param x Coordinates for this grid point.
 *  @param t Absolute time.
 */
template<typename TT, std::size_t DIM>
TT cubicFuncSumOfSpaceDerivatives2ndOrder(ScaFES::Ntuple<TT,DIM> const& x,
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

/** Calculates value of a cubic function for given grid point and time.
 *  @param x Coordinates for this grid point.
 *  @param t Absolute time.
 */
template<typename TT, std::size_t DIM>
TT cubicFunc(ScaFES::Ntuple<TT,DIM> const& x, TT const& t) {
    return t*cubicFuncArg<TT,DIM>(x);
}

#endif
