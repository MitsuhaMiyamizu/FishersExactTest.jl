using SpecialFunctions: logabsgamma
const KF_GAMMA_EPS = 1e-14
const KF_TINY = 1e-290

# complementary error function
# \frac{2}{\sqrt{\pi}} \int_x^{\infty} e^{-t^2} dt
# AS66, 2nd algorithm, http://lib.stat.cmu.edu/apstat/66

function kf_erfc(x::Float64)
    p0 = 220.2068679123761
	p1 = 221.2135961699311
	p2 = 112.0792914978709
	p3 = 33.912866078383
	p4 = 6.37396220353165
	p5 = 0.7003830644436881
	p6 = 0.03526249659989109
	q0 = 440.4137358247522
	q1 = 793.8265125199484
	q2 = 637.3336333788311
	q3 = 296.5642487796737
	q4 = 86.78073220294608
	q5 = 16.06417757920695
	q6 = 1.755667163182642
	q7 = 0.08838834764831844
    M_SQRT2 = sqrt(2)
    expntl = float(0)
    z = float(0)
    p = float(0)
    z = abs(x) * M_SQRT2
    if z > 37.0 
        return x > 0. ? x = 0. : x = 2.
    end
    expntl = exp(z * z * - 0.5)
    if z < (10.0 / M_SQRT2)
        p = expntl * ((((((p6 * z + p5) * z + p4) * z + p3) * z + p2) * z + p1) * z + p0) / (((((((q7 * z + q6) * z + q5) * z + q4) * z + q3) * z + q2) * z + q1) * z + q0)
	else 
        p = expntl / 2.506628274631001 / (z + 1. / (z + 2. / (z + 3. / (z + 4. / (z + .65)))))
    end
    return x > 0. ? x = 2. * p : x = 2. * (1. - p)
end

"""
# The following computes regularized incomplete gamma functions.
# Formulas are taken from Wiki, with additional input from Numerical
# Recipes in C (for modified Lentz's algorithm) and AS245
# (http://lib.stat.cmu.edu/apstat/245).
#
# A good online calculator is available at:
#
#   http://www.danielsoper.com/statcalc/calc23.aspx
#
# It calculates upper incomplete gamma function, which equals
# kf_gammaq(s,z)*tgamma(s).
"""

function _kf_gammap(s::Float64, z::Float64)
    sum = float(0)
    x = float(0)
    k = 1
    sum = x = 1.
    while k < 100
        sum += (x *= z / (s + k))
        if ((x / sum) < KF_GAMMA_EPS) break end
	k += 1
    end
    return exp(s * log(z) - z - logabsgamma(s + 1.)[1] + log(sum))
end

function _kf_gammaq(s::Float64, z::Float64)
    j = Int64(1)
    D = float(0)
    #f = float(0)
    f = 1. + z - s
    C = f
    for j in 1:99
        a = j * (s - j)
        b = (j << 1) + 1 + z - s
        D = b + a * D
        D < KF_TINY ? D = KF_TINY : D 
        C = b + a / C
        C < KF_TINY ? C = KF_TINY : C
        D = 1. / D
        d = C * D
        f *= d
        if (abs(d - 1.)) < KF_GAMMA_EPS 
            break
        end
    end
    return exp(s * log(z) - z - logabsgamma(s)[1] - log(f))
end

function kf_gammap(s::Float64, z::Float64)
    if z <= 1. || z < s 
        return _kf_gammap(s, z)
    else 
        return 1. - _kf_gammaq(s, z)
    end
end

function kf_gammaq(s::Float64, z::Float64)
    if z <= 1. || z < s
        return 1. - _kf_gammap(s, z) 
    else 
        return _kf_gammaq(s, z)
    end
end

"""
Regularized incomplete beta function. The method is taken from
Numerical Recipe in C, 2nd edition, section 6.4. The following web
page calculates the incomplete beta function, which equals
kf_betai(a,b,x) * gamma(a) * gamma(b) / gamma(a+b):

http://www.danielsoper.com/statcalc/calc36.aspx
"""

function kf_betai_aux(a::Float64, b::Float64, x::Float64)
    D = float(0)
    if x == 0. return 0. end
    if x == 1. return 1. end
    f = 1.
    C = f
    # Modified Lentz's algorithm for computing continued fraction
    for j in 1:199
        m = j >> 1
        (j & 1 == 1) ? aa = - (a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1)) : aa = m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m))
        D = 1. + aa * D
        D < KF_TINY ? D = KF_TINY : D
        C = 1. + aa / C
        C < KF_TINY ? C = KF_TINY : C
        D = 1. / D
        d = C * D
        f *= d
        if abs(d - 1.) < KF_GAMMA_EPS break end
    end
    return exp(logabsgamma(a + b)[1] - logabsgamma(a)[1] - logabsgamma(b)[1] + a * log(x) + b * log(1. - x)) / a / f
end

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
