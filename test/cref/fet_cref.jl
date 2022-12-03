#!/usr/bin/env julia

import fet.hgacc_t

CREF_IMAGE = string(@__DIR__, "/libkfunc.so")

function kf_lgamma_cref(z::Cdouble) ::Cdouble
    return ccall((:kf_lgamma, CREF_IMAGE),
               Cdouble, (Cdouble,),
               z)
end # kf_lgamma_cref


function kf_erfc_cref(x::Cdouble) ::Cdouble
    return ccall((:kf_erfc, CREF_IMAGE),
                 Cdouble, (Cdouble,),
                x)
end # kf_erfc_cref


function _kf_gammap_cref(s::Cdouble, z::Cdouble) ::Cdouble
    return ccall((:_kf_gammap, CREF_IMAGE),
                 Cdouble, (Cdouble, Cdouble),
                 s, z)
end # _kf_gammap_cref


function _kf_gammaq_cref(s::Cdouble, z::Cdouble) ::Cdouble
    return ccall((:_kf_gammaq, CREF_IMAGE),
                 Cdouble, (Cdouble, Cdouble),
                 s, z)
end # _kf_gammaq_cref


function kf_gammap_cref(s::Cdouble, z::Cdouble) ::Cdouble
    return ccall((:kf_gammap, CREF_IMAGE),
                 Cdouble, (Cdouble, Cdouble),
                 s, z)
end # kf_gammap_cref


function kf_gammaq_cref(s::Cdouble, z::Cdouble) ::Cdouble
    return ccall((:kf_gammaq, CREF_IMAGE),
                 Cdouble, (Cdouble, Cdouble),
                 s, z)
end # kf_gammaq_cref


function kf_betai_aux_cref(a::Cdouble, b::Cdouble, x::Cdouble) ::Cdouble
    return ccall((:kf_betai_aux, CREF_IMAGE),
                Cdouble, (Cdouble, Cdouble, Cdouble),
                a, b, x)
end # kf_betai_cref


function kf_betai_cref(a::Cdouble, b::Cdouble, x::Cdouble) ::Cdouble
    return ccall((:kf_betai, CREF_IMAGE),
                Cdouble, (Cdouble, Cdouble, Cdouble),
                a, b, x)
end # kf_betai_cref


function lbinom_cref(n::Clonglong, k::Clonglong) ::Cdouble
    return ccall((:lbinom, CREF_IMAGE),
                 Cdouble, (Clonglong, Clonglong),
                 n, k)
end # lbinom_cref


function hypergeo_cref(n11::Clonglong, n1_::Clonglong, n_1::Clonglong, n::Clonglong) ::Cdouble
    return ccall((:hypergeo, CREF_IMAGE),
                 Cdouble, (Clonglong, Clonglong, Clonglong, Clonglong),
                 n11, n1_, n_1, n)
end # hypergeo_cref


function hypergeo_acc_cref(n11::Clonglong, n1_::Clonglong, n_1::Clonglong, n::Clonglong, aux::Ptr{hgacc_t}) ::Cdouble
    return ccall((:hypergeo_acc, CREF_IMAGE),
                 Cdouble, (Clonglong, Clonglong, Clonglong, Clonglong, Ptr{hgacc_t}),
                 n11, n1_, n_1, n, aux)
end # hypergeo_acc_cref


function kt_fisher_exact_cref(n11::Clonglong, n12::Clonglong, n21::Clonglong, n22::Clonglong) ::Tuple{Cdouble, Cdouble, Cdouble, Cdouble}
    p_l = Ref(Cdouble(0.0));
    p_r = Ref(Cdouble(0.0));
    p_t = Ref(Cdouble(0.0));
    p = ccall((:kt_fisher_exact, CREF_IMAGE),
              Cdouble, (Clonglong, Clonglong, Clonglong, Clonglong, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), 
              n11, n12, n21, n22, p_l, p_r, p_t);

    return (p, p_l.x, p_r.x, p_t.x)
end
