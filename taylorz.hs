module TaylorZ where

-- This is pretty much a copy of taylor.hs, except that the Taylor type
-- definition has been augmented with a `Zero`.  The reason for doing this is
-- that in the original Taylor definition, which is nevertheless correct and
-- works, all the compound expressions have "expanding" constructs in the
-- derivative part - meaning the resultant is always more complex than the
-- parts. "a + b" will always be more complx than "a" and "b". This is quite
-- unnecessary in many cases and having an explicit representation for `Zero`
-- helps identify many cases where we can express "contracting" constructs that
-- can simplify expression construction. Having a zero helps model constants,
-- for example, without forcing a calculation tree for the derivative part that
-- always needs to be evaluated even if only to get 0 as the final answer in all
-- the cases.
--
-- To feel the difference between the complexities of Taylor versus TaylorZ,
-- try evaluating the 25th derivative of `square x = x * x` at some `x` using
-- both. If you didn't feel the difference yet, try higher derivatives.

-- "Taylor" is a "Dual Number" type, only recursive - so we may call them
-- "Taylor numbers". The first part is a normal function value and the second
-- part is the value of the derivative of the function. The advantage of this
-- over the dmath.hs is that you should be able to pass Taylor numbers wherever
-- you write normal number functions and get derivatives also calculated. This
-- includes if-then-else conditionals as well.
data Taylor x = Zero | Taylor x (Taylor x)

-- The "Calculus in coinductive form" paper at -
-- http://lya.fciencias.unam.mx/favio/publ/cocalculus.pdf
-- talks about essentially the idea in this module.
-- The only remaining bit is to use the taylor series
-- calculated at a particular point of a function
-- to evaluate it around that value. Towards this end
-- I define a `teval` function that takes such a "taylor 
-- number" and uses it as the taylor series to calculate
-- nearby values for a given deviation x from that point,
-- up to k terms of the series.
--
-- For example, you can define e^x using the Taylor
-- expansion around 0 as `e = Taylor 1 e`, i.e. by
-- treating e^x as the solution to the differential
-- equation dy/dx = y .Now, you can e^x using
-- `teval 10 e x` which will calculate an approximation
-- using 10 terms of the series. Similarly, you can calculate
-- sin(x) also by treating sin(x) as the solution to
-- d^y/dx^2 = -y, with the initial conditions sin(0) = 0
-- and sin'(0) = 1 ... which can be expressed as --
-- `sin0 = Taylor 0.0 (Taylor 1.0 (-sin0))`
-- and evaluate sin(x) using `teval 10 sin0 x`.
teval 0 t x = f t
teval k Zero x = 0
teval k t x = (teval (k-1) t x) + (dfn k t) * (x ^ k) / fromIntegral (product [1..k])

-- f adn = the function part
-- df adn = the derivative part
f Zero = 0
f (Taylor f0 _) = f0
df Zero = Zero
df (Taylor _ df0) = df0

-- Recursive calculation of the nth derivative.
dfn 0 = f 
dfn n = dfn (n-1) . df

-- You can provide values for the x variable using this function.
-- If you have a normal mathematical function to calculate some
-- value from a given number x, you can evaluate the same function
-- by passing a "taylor number" as a value and have both the value
-- and the derivative calculated simultaneously.
x :: (Num a) => a -> Taylor a
x a = Taylor a 1

-- Basic arithmetic operations
instance (Num n) => Num (Taylor n) where
    Zero + b = b
    a + Zero = a
    a + b = Taylor (f a + f b) (df a + df b)
    a * Zero = 0
    Zero * b = 0
    a * b = Taylor (f a * f b) (df a * b + a * df b)
    negate Zero = 0
    negate a = Taylor (- f a) (- df a)
    abs Zero = 0
    abs a = Taylor (abs (f a)) (df a * signum a)
    signum Zero = 0
    signum a = Taylor (signum (f a)) 0
    fromInteger 0 = Zero
    fromInteger x = Taylor (fromInteger x) 0

instance (Eq n, Num n) => Eq (Taylor n) where
    a == b = f a == f b

instance (Eq n, Num n, Ord n) => Ord (Taylor n) where
    a <= b = f a <= f b

-- Reciprocal
instance (Fractional n) => Fractional (Taylor n) where
    fromRational 0 = Zero
    fromRational x = Taylor (fromRational x) 0
    recip a = Taylor (1 / f a) (- df a / (a * a))

-- Scientific functions
instance (Floating n) => Floating (Taylor n) where
    pi = Taylor pi 0
    exp a = Taylor (exp (f a)) (df a * exp a)
    log a = Taylor (log (f a)) (df a / a)
    sin a = Taylor (sin (f a)) (df a * cos a)
    cos a = Taylor (cos (f a)) (- df a * sin a)
    asin a = Taylor (asin (f a)) (df a / sqrt (1 - a * a))
    acos a = Taylor (acos (f a)) (- df a / sqrt (1 - a * a))
    atan a = Taylor (atan (f a)) (df a / (1 + a * a))
    sinh a = Taylor (sinh (f a)) (df a * cosh a)
    cosh a = Taylor (cosh (f a)) (df a * sinh a)
    asinh a = Taylor (asinh (f a)) (df a / sqrt (1 + a * a))
    acosh a = Taylor (acosh (f a)) (df a / sqrt (a * a - 1))
    atanh a = Taylor (atanh (f a)) (df a / (1 + a * a))
