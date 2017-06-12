/* ScaFES
 * Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

#ifndef INDICATORELLIPSOID_HPP
#define    INDICATORELLIPSOID_HPP

#include <ScaFES.hpp>

template <typename T>
class IndicatorEllipsoid {
    // steepness param h
    const T h_=1.e-8;

public:
    typedef ScaFES::Ntuple<T, 3> T3;

    IndicatorEllipsoid(const IndicatorEllipsoid& orig)=delete;
    ~IndicatorEllipsoid()=default;

    /**
     * Construct an ellipsoid with indicator function.
     * @param x center
     * @param y center
     * @param z center
     * @param a radius
     * @param b radius
     * @param c radius
     * @param value_inside
     * @param value_outside
     */
    IndicatorEllipsoid(
            const T& x, const T& y, const T& z,
            const T& a, const T& b, const T& c,
            const T& value_inside, const T& value_outside)
    : center_(x, y, z)
    , radius(a, b, c)
    , inside_(value_inside)
    , outside_(value_outside) { };

    /**
     * Compute parameter value for provided coordinate.
     * @param addr
     * @return value of material parameter at given coordinate.
     */
    T operator()(const T3& addr) const
    {
        T3 faddr=addr-center_;

        // the same hints as in class EMFDTD apply: do not use auto, keep to the prescribed order
        T tmp((faddr.elem(0)*faddr.elem(0))+(faddr.elem(1)*faddr.elem(1)));
        T theta(atan2(sqrt(tmp), faddr.elem(2)));
        T phi(atan2(faddr.elem(1), faddr.elem(0)));
        T ro(sqrt(tmp+(faddr.elem(2)*faddr.elem(2))));

        T x1(radius.elem(0)*cos(theta)*cos(phi));
        T y1(radius.elem(1)*cos(theta)*sin(phi));
        T z1(radius.elem(2)*sin(theta));

        T rn=sqrt((x1*x1)+(y1*y1)+(z1*z1));

        T ret=0.5*inside_*(((ro+rn)/(sqrt((ro+rn)*(ro+rn)+h_*h_))) -
            ((ro-rn)/(sqrt((ro-rn)*(ro-rn)+h_*h_))));

        return ret;
    };

private:
    const T3 center_, radius;
    const T& inside_, outside_;
};

#endif    /* INDICATORELLIPSOID_HPP */
