# df

Automatic differentiation is being built into programming languages such
as Swift and Julia and added as libraries to languages to Python. This
is changing the way many optimization problems are expressed and the
sophistication with which deep learning models are being made.

This repo is an attempt to explain and understand automatic differentiation
at the conceptual level. I may choose to elaborate it later on with
algorithms for efficient implementation, but, as of this writing,
the main intention is illustration, understanding and exploration.

## Blog posts on the topic

1. [AD - differentiable functions](http://sriku.org/blog/2019/03/08/automatic-differentiation/) - by treating functions of a single variable as numbers.
2. [AD - Dual numbers and Taylor numbers](http://sriku.org/blog/2019/03/12/automatic-differentiation-dual-numbers--taylor-numbers/)
3. [AD - Higher ranked beings](http://sriku.org/blog/2019/03/13/automatic-differentiation-higher-ranked-beings/) - generalized vectors and tensors.

## Source code index

1. [dmath.hs](https://github.com/srikumarks/df/blob/master/dmath.hs) - Differentiable functions
2. [dual.hs](https://github.com/srikumarks/df/blob/master/dual.hs) - Dual numbers
3. [taylor.hs](https://github.com/srikumarks/df/blob/master/taylor.hs) - Taylor numbers. Not exactly a known term, but this is definition using which all derivatives at a point for a given function can be calculated .. meaning we can produce arbitrary Taylor series approximations using this approach and hence I'm calling these "Taylor numbers".
4. [taylorz.hs](https://github.com/srikumarks/df/blob/master/taylorz.hs) - Same as `taylor.hs` above, but with special support for zero to prune the expression tree.
5. [dvmath.hs](https://github.com/srikumarks/df/blob/master/dvmath.hs) - Applying the approach of `dmath.hs` but for functions of vectors, tensors and such higher ranked objects. This introduces more operations such as inner/outer products, convolutions, slices, etc. .. compared to the usual number protocol.


