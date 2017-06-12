/* ScaFES
 * Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

#ifndef COMPLETECPMLCURLSTENCIL_HPP
#define    COMPLETECPMLCURLSTENCIL_HPP

template <typename T, typename I,
typename I::value_type off_l0,
typename I::value_type off_l1,
typename I::value_type off_l2,
typename I::value_type off_r0,
typename I::value_type off_r1,
typename I::value_type off_r2>
class CompleteCurlCPMLStencil {
private:
    typedef T value_type;
    typedef I Index;
    typedef typename I::value_type index_type;

public:

    CompleteCurlCPMLStencil(const CPMLParams<value_type>& p,
                const Index& start, const Index& end,
                const Index& last_node,
                const value_type constant)
    : p_(p)
    , constant_(constant)
    , I_(last_node[0]+end[0])
    , J_(last_node[1]+end[1])
    , K_(last_node[2]+end[2])
    , start_(start)
    , end_(I_, J_, K_) { };

    template <typename DF>
    typename DF::value_type operator()(const Index& idx,
        const DF& dest_old,
        const DF& from_l, const DF& psi_l_old, DF& psi_l_new,
        const DF& from_r, const DF& psi_r_old, DF& psi_r_new
        ) const
    {
        // skip invalid elements
        for (unsigned long int i=0; i<Index::DIM; ++i) {
            if ((start_[i]>idx[i])||idx[i]>end_[i]) {
                return 0.;
            }
        }

        typedef typename DF::value_type TT;

        const index_type i=idx[0];
        const index_type j=idx[1];
        const index_type k=idx[2];

        const index_type dlo=
            std::min<index_type>(i, I_-i)*std::abs(off_l0)+
            std::min<index_type>(j, J_-j)*std::abs(off_l1)+
            std::min<index_type>(k, K_-k)*std::abs(off_l2);
        const index_type dl=std::min<index_type>(dlo, p_.D);

        TT dfdl=from_l(Index(i+off_l0, j+off_l1, k+off_l2))-from_l(idx);
        TT psi_l_n=(p_.b[dl]*psi_l_old(idx)) + (p_.c[dl]*dfdl);
        dfdl=(p_.ook[dl]*dfdl)+psi_l_n;
        psi_l_new(idx)=psi_l_n;

        const index_type dro=std::min<index_type>(i, I_-i)*std::abs(off_r0)+
            std::min<index_type>(j, J_-j)*std::abs(off_r1)+
            std::min<index_type>(k, K_-k)*std::abs(off_r2);
        const index_type dr=std::min<index_type>(dro, p_.D);

        TT dfdr=from_r(Index(i+off_r0, j+off_r1, k+off_r2))-from_r(idx);
        TT psi_r_n=(p_.b[dr]*psi_r_old(idx)) + (p_.c[dr]*dfdr);
        dfdr=(p_.ook[dr]*dfdr)+psi_r_n;
        psi_r_new(idx)=psi_r_n;

        TT ret=dest_old(idx)+(constant_*(dfdl-dfdr));
        return ret;
    }

private:
    const CPMLParams<value_type>& p_;
    const value_type constant_;
    const index_type I_, J_, K_;
    const Index start_, end_;
};


#endif    /* COMPLETECPMLCURLSTENCIL_HPP */

