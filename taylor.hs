module Taylor where

-- "Taylor" is a "Dual Number" type, only recursive - so we may call them
-- "Taylor numbers". The first part is a normal function value and the second
-- part is the value of the derivative of the function. The advantage of this
-- over the dmath.hs is that you should be able to pass Taylor numbers wherever
-- you write normal number functions and get derivatives also calculated. This
-- includes if-then-else conditionals as well.
data Taylor x = Taylor x (Taylor x)

-- Simple Show instance to display one term of the series.
instance (Show x) => Show (Taylor x) where
    show (Taylor v dv) = ":> " ++ (show v)

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
teval k t x = (teval (k-1) t x) + (dfn k t) * (x ^ k) / fromIntegral (product [1..k])

-- f adn = the function part
-- df adn = the derivative part
f (Taylor f0 _) = f0
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
    a + b = Taylor (f a + f b) (df a + df b)
    a * b = Taylor (f a * f b) (df a * b + a * df b)
    negate a = Taylor (- f a) (- df a)
    abs a = Taylor (abs (f a)) (df a * signum a)
    signum a = Taylor (signum (f a)) 0
    fromInteger x = Taylor (fromInteger x) 0

instance (Eq n) => Eq (Taylor n) where
    a == b = f a == f b

instance (Eq n, Ord n) => Ord (Taylor n) where
    a <= b = f a <= f b

-- Reciprocal
instance (Fractional n) => Fractional (Taylor n) where
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
