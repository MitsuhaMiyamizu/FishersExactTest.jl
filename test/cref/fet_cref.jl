#!/usr/bin/env julia

function fisher_exact_cref(n11::Clonglong, n12::Clonglong, n21::Clonglong, n22::Clonglong) ::Tuple{Cdouble, Cdouble, Cdouble, Cdouble}
    p_l = Ref(Cdouble(0.0));
    p_r = Ref(Cdouble(0.0));
    p_t = Ref(Cdouble(0.0));
    p = ccall((:kt_fisher_exact, string(@__DIR__, "/libkfunc.so")),
              Cdouble, (Clonglong, Clonglong, Clonglong, Clonglong, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), 
              n11, n12, n21, n22, p_l, p_r, p_t);

    return (p, p_l.x, p_r.x, p_t.x)
end
