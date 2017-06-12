/* ScaFES
 * Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

#ifndef MODULATEDGAUSSIANSOURCE_HPP_
#define    MODULATEDGAUSSIANSOURCE_HPP_

#include <boost/property_tree/ptree.hpp>

template <typename T>
class ModulatedGaussianSource {
public:
    ModulatedGaussianSource(const ModulatedGaussianSource&)=delete;

    ~ModulatedGaussianSource() { };

    ModulatedGaussianSource(const T& dt, const boost::property_tree::ptree& p)
    : dt_(dt)
    , f_mod_(p.get<T>("source.f_mod"))
    , t0_(p.get<T>("source.t0"))
    , spread_(p.get<T>("source.spread")) { };

    template <typename TT> TT operator()(const TT & n) const
    {
        TT v=(t0_-n)/spread_;
        v*=v;
        v*=-.5;
        v=exp(v);
        TT mod=cos(2.*M_PI*f_mod_*n*dt_);
        TT ret=v*mod;
        return ret;
    }

private:
    const T dt_;
    const T f_mod_;
    const T t0_;
    const T spread_;
};

#endif    /* MODULATEDGAUSSIANSOURCE_HPP_ */

