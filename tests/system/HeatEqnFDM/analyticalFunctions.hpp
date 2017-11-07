/* ScaFES
 * Copyright (c) 2017, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file analyticalFunctions.hpp
 *  @brief Contains functions which can be used as analytical solutions.
 */

#ifndef ANALYTICAL_FUNCTIONS_HPP_
#define ANALYTICAL_FUNCTIONS_HPP_

/** Calculates value of a linear time function for given time.
 *  @param t Value of absolute time.
 */
template<typename TT>
TT timeLin(TT const& t) {
    return (1.0 + t);
}

/** Calculates value of a constant space function for given grid point.
 *  @param x Coordinates of given grid point.
 */
template<typename TT, std::size_t DIM>
TT spaceConst(ScaFES::Ntuple<TT,DIM> const& /*x*/) {
    return 1.0;
}

/** Calculates value of a linear space function for given grid point.
 *  @param x Coordinates of given grid point.
 */
template<typename TT, std::size_t DIM>
TT spaceLin(ScaFES::Ntuple<TT,DIM> const& x) {
    TT argument = 1.0;
    for (std::size_t pp = 0; pp < DIM; ++pp) {
        argument += x[pp];
    }
    return argument;
}

/** Calculates value of time derivative of a linear time and
 *  a constant space function for given grid point.
 *  @param x Coordinates of given grid point.
 */
template<typename TT, std::size_t DIM>
TT timeLinSpaceConstdTime(ScaFES::Ntuple<TT,DIM> const& x) {
    return spaceConst<TT,DIM>(x);
}

/** Calculates value of time derivative of a linear time and
 *  a linear space function for given grid point.
 *  @param x Coordinates of given grid point.
 */
template<typename TT, std::size_t DIM>
TT timeLinSpaceLindTime(ScaFES::Ntuple<TT,DIM> const& x) {
    return spaceLin<TT,DIM>(x);
}

/** Calculates value of first order space derivative of a linear time
 *  and a constant space function for given grid point, time and orientation.
 *  @param x Coordinates of given grid point.
 *  @param t Value of absolute time.
 *  @param currDIM Current dimension (orientation of points).
 */
template<typename TT, std::size_t DIM>
TT timeLinSpaceConstdSpace1stOrder(ScaFES::Ntuple<TT,DIM> const& /*x*/,
                                   TT const& /*t*/, std::size_t /*currDIM*/) {
    return 0.0;
}

/** Calculates value of first order space derivative of a linear time
 *  and a linear space function for given grid point, time and orientation.
 *  @param x Coordinates of given grid point.
 *  @param t Value of absolute time.
 *  @param currDIM Current dimension (orientation of points).
 */
template<typename TT, std::size_t DIM>
TT timeLinSpaceLindSpace1stOrder(ScaFES::Ntuple<TT,DIM> const& /*x*/,
                                 TT const& t, std::size_t /*currDIM*/) {
    return timeLin<TT>(t);
}

/** Calculates sum of second order space derivatives for all dimensions of
 *  a linear time and a constant space function for given grid point and time.
 *  @param x Coordinates of given grid point.
 *  @param t Value of absolute time.
 */
template<typename TT, std::size_t DIM>
TT timeLinSpaceConstSumOfdSpace2ndOrder(ScaFES::Ntuple<TT,DIM> const& /*x*/,
                                        TT const& /*t*/) {
    return 0.0;
}

/** Calculates sum of second order space derivatives for all dimensions
 *  of a linear time and a linear space function for given grid point and time.
 *  @param x Coordinates of given grid point.
 *  @param t Value of absolute time.
 */
template<typename TT, std::size_t DIM>
TT timeLinSpaceLinSumOfdSpace2ndOrder(ScaFES::Ntuple<TT,DIM> const& /*x*/,
                                      TT const& /*t*/) {
    return 0.0;
}

/** Calculates value of a linear time and a constant space function
 *  for given grid point and time.
 *  @param x Coordinates of given grid point.
 *  @param t Value of absolute time.
 */
template<typename TT, std::size_t DIM>
TT timeLinSpaceConstFunc(ScaFES::Ntuple<TT,DIM> const& x, TT const& t) {
    return (timeLin<TT>(t) * spaceConst<TT,DIM>(x));
}

/** Calculates value of a linear time and a linear space function
 *  for given grid point and time.
 *  @param x Coordinates of given grid point.
 *  @param t Value of absolute time.
 */
template<typename TT, std::size_t DIM>
TT timeLinSpaceLinFunc(ScaFES::Ntuple<TT,DIM> const& x, TT const& t) {
    return (timeLin<TT>(t) * spaceLin<TT,DIM>(x));
}

#endif
