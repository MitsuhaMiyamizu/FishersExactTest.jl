# FishersExactTest.jl
An accurate and Julia-native implementation of Fisher's Exact Test for Julia. (No need to ```ccall``` external C library)

# Motivation
We and [someone else](https://blog.goo.ne.jp/r-de-r/e/f206d5a1bfb6ed1f60d9bd47c0865299) found an accuracy-related [bug](https://github.com/JuliaStats/HypothesisTests.jl/issues/148) in the implementation of Fisher's Exact Test provided by [JuliaStats](https://github.com/JuliaStats/HypothesisTests.jl), it seems that due to technical limitations behind their algorithm, the bug cannot be fixed in a short period of time.

To ensure the accuracy of the implementation in our package, we compare our result to known facts in [Wikipedia](https://en.wikipedia.org/wiki/Fisher%27s_exact_test)
and use ```fisher.test``` in R to validate our result, which yields a spot-on match.

# How to use
Simply by querying
```
import Pkg; Pkg.add("FishersExactTest")
using FishersExactTest
FisherExact2x2Test(a, b, c, d);

k = FisherExact2x2Test(a, b, c, d, opt = :silent);
k.p          # p_value
k.left_tail  # left_tail
k.right_tail # right_tail
k.two_tail   # two_sides

# Options for printing
FisherExact2x2Test(a, b, c, d, opt = :verbose);
# Options for suppress print
k = FisherExact2x2Test(a, b, c, d, opt = :silent);
```
where a, b, c, d are distributed in
| a  |  b |
|---|---|
|  c |  d |

# Reference
Fisher's Exact Test was extracted from [bedtools](https://github.com/arq5x/bedtools2/tree/89b94dce487097e60bbd6d77c2515085c6e80431/src/fisher)

# Courtesy
[@Ionizing](https://github.com/Ionizing)

[@atlasmir](https://github.com/atlasmir)
