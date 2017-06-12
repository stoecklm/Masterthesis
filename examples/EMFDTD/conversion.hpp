/* ScaFES
 * Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

#ifndef CONVERSION_HPP_
#define CONVERSION_HPP_

#include <boost/mpl/if.hpp>
#include <adolc/adolc.h>

template <typename T>
struct pass_through {
    T operator()(const T& s) const
    {
        return s;
    }
};

template <typename D, typename S>
struct filter { };

template <typename D>
struct filter<D,adouble> {
    D operator()(const adouble& s) const
    {
        return s.value();
    };
};

template <typename D, typename S>
struct convert {
    typedef typename
    boost::mpl::if_<boost::is_same<D, S>, pass_through<D>, filter<D, S> >::type executor;

    D operator()(const S& s) const
    {
        return executor()(s);
    };
};

#endif // CONVERSION_HPP_

