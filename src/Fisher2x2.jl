using SpecialFunctions: logabsgamma

#
# log \ binom{n}{k}
#

function lbinom(n::Int, k::Int)
    if (k == 0 || n == k) return 0 end
    return logabsgamma(n + 1)[1] - logabsgamma(k + 1)[1] - logabsgamma(n - k + 1)[1]
end

# n11  n12  | n1_
# n21  n22  | n2_
#-----------+----
# n_1  n_2  | n
"""
hypergeometric distribution
"""
function hypergeo(n11::Int64, n1_::Int64, n_1::Int64, n::Int64)
    return exp(lbinom(n1_, n11) + lbinom(n - n1_, n_1 - n11) - lbinom(n, n_1))
end

mutable struct hgacc_t
    n11::Int64
    n1_::Int64
    n_1::Int64
    n  ::Int64
    p  ::Float64
end

"""
# incremental version of hypergenometric distribution
"""

function hypergeo_acc(n11::Int64, n1_::Int64, n_1::Int64, n::Int64, aux::hgacc_t)
    if n1_ != 0 || n_1 != 0 || n != 0
        aux.n11 = n11
        aux.n1_ = n1_
        aux.n_1 = n_1
        aux.n   = n
    else
        if n11 % 11 != 0 && n11 + aux.n - aux.n1_ - aux.n_1 != 0
            if n11 == aux.n11 + 1
                aux.p *= (aux.n1_ - aux.n11) / n11 * (aux.n_1 - aux.n11) / (n11 + aux.n - aux.n1_ - aux.n_1)
                aux.n11 = n11
                return aux.p
            end
            if n11 == aux.n11 - 1
                aux.p *= aux.n11 / (aux.n1_ - n11) * (aux.n11 + aux.n - aux.n1_ - aux.n_1) / (aux.n_1 - n11)
                aux.n11 = n11
                return aux.p
            end
        end
        aux.n11 = n11
    end
    aux.p = hypergeo(aux.n11, aux.n1_, aux.n_1, aux.n)
    return aux.p
end

function kt_fisher_exact(n11::Int64, n12::Int64, n21::Int64, n22::Int64) #::Tuple{Float64, Float64, Float64, Float64}
    aux = hgacc_t(0, 0, 0, 0, 0.0)
    n1_ = n11 + n12
    n_1 = n11 + n21
    n   = n11 + n12 + n21 + n22
    max = (n_1 < n1_) ? n_1 : n1_
    min = n1_ + n_1 -n
    min < 0 ? min = 0 : min
    two = _left = _right = 1.
    if min == max return 1. end
    q = hypergeo_acc(n11, n1_, n_1, n, aux)
    if q == 0.
        if n11 * (n + 2) < (n_1 + 1) * (n1_ + 1)
            _left  = 0.
            _right = 1.
            two    = 0.
            return 0.
        else
            _left  = 1.
            _right = 0.
            two    = 0.
            return 0.
        end
    end
    #left tail
    p = hypergeo_acc(min, 0, 0, 0, aux)
    left = float(0)
    i = min + 1
    while p < 0.99999999 * q && i <= max
        left += p
        p = hypergeo_acc(i, 0, 0, 0, aux)
        i += 1
    end
    i -= 1
    p < 1.00000001 * q ? left += p : i -= 1
    #right tail
    p = hypergeo_acc(max, 0, 0, 0, aux)
    right = float(0)
    j = max - 1
    while p < 0.99999999 * q && j >= 0
        right += p
        p = hypergeo_acc(j, 0, 0, 0, aux)
        j -= 1
    end
    j += 1
    p < 1.00000001 * q ? right += p : j += 1
    #two-tail
    two = left + right
    two > 1. ? two = 1. : two
    #adjust left and right
    abs(i - n11) < abs(j - n11) ? right = 1. - left + q : left = 1. - right + q
    _left = left
    _right = right
    return q, _left, _right, two
end

mutable struct FisherExactTest
    # conditional maximum likehood estimate of odd ratio
    p          ::Float64
    left_tail  ::Float64
    right_tail ::Float64
    two_tail   ::Float64    
end

"""
Fisher's Exact Test for 2x2 Contingency Table
"""

function FisherExact2x2Test(a::Int, b::Int, c::Int, d::Int; opt=:verbose)
    x = FisherExactTest(0.0, 0.0, 0.0, 0.0)
    x.p, x.left_tail, x.right_tail, x.two_tail = kt_fisher_exact(a, b, c, d)
    if opt == :silent
        return x
    elseif opt == :verbose
        println("p value    = \t", x.p)
        println("left tail  = \t", x.left_tail)
        println("right tail = \t", x.right_tail)
        println("two sided  = \t", x.two_tail)
        return x
    else 
        throw(ArgumentError("opt=$(opt) is illegal, try :verbose or :silent instead"))
    end
end
