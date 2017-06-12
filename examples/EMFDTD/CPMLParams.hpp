/* ScaFES
 * Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

#ifndef CPMLPARAMS_HPP
#define CPMLPARAMS_HPP

#include <cmath>
#include <algorithm>
#include <vector>

#include "physical_constants.hpp"

template <typename T>
struct CPMLParams {
    CPMLParams()=delete;
    CPMLParams(const CPMLParams&)=delete;
    CPMLParams& operator=(const CPMLParams&)=delete;
    ~CPMLParams()=default;

    CPMLParams(const T dt,
        std::size_t D_,
        const T m_a, const T m,
        const T k_max, const T s_max, const T a_max)
    : b(D_+1,T(0.)), c(D_+1,T(0.)), ook(D_+1,T(1.)), D(D_)
    {
        const T d=T(D);

        for (std::size_t w=0; w<D_; ++w) {
            const T wd=std::pow<T>(w/(d-1.), m);
            const T dw_a=std::pow<T>((d-1.-w)/(d-1.), m_a);

            const T s=wd*s_max;
            const T k=1.+((k_max-1.)*wd);

            const T a_=dw_a*a_max;

            b[w]=std::exp(-dt*(((s/k)+a_)/EPSILON_0));
            c[w]=(s*(b[w]-1.))/((s*k)+(k*k*a_));
            if (b[w]==1.)
                c[w]=0.;
            ook[w]=1./k;
        }
        std::reverse(b.begin(),--b.end());
        std::reverse(c.begin(),--c.end());
        std::reverse(ook.begin(),--ook.end());
    };

    void show() const
    {
        std::cout<<"# w b c 1/k"<<std::endl;
        for (std::size_t w=0; w<=D; ++w)
            std::cout<<w<<" "<<b[w]<<" "<<c[w]<<" "<<ook[w]<<std::endl;
    };

    std::vector<T> b, c, ook;
    std::size_t D;
};

#endif    /* CPMLPARAMS_HPP */
