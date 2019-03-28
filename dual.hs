module Dual where

-- "Dual" is the "Dual Number" type. The first part is a normal value,
-- and the second part holds the derivative's value. By defining
-- operations on dual numbers carefully, we're able to simultaneously
-- calculate a function's value as well as its derivate's value at
-- a given point. The advantage of this over dmath.hs is that you can
-- use dual numbers with normal functions, though you cannot use them
-- to calculate higher order derivatives.
data Dual x = Dual x x deriving (Show)

-- f adn = the function part
-- df adn = the derivative part
f (Dual f0 _) = f0
df (Dual _ df0) = df0

-- You can provide values for the x variable using this function.
-- If you have a normal mathematical function to calculate some
-- value from a given number x, you can evaluate the same function
-- by passing a "dual number" as a value and have both the value
-- and the derivative calculated simultaneously.
x :: (Num a) => a -> Dual a
x a = Dual a 1

-- Basic arithmetic operations
instance (Num n) => Num (Dual n) where
    a + b = Dual (f a + f b) (df a + df b)
    a * b = Dual (f a * f b) (df a * f b + f a * df b)
    negate a = Dual (- f a) (- df a)
    abs a = Dual (abs (f a)) (df a * signum (f a))
    signum a = Dual (signum (f a)) 0
    fromInteger x = Dual (fromInteger x) 0

instance (Eq n) => Eq (Dual n) where
    a == b = f a == f b

instance (Eq n, Ord n) => Ord (Dual n) where
    a <= b = f a <= f b

-- Reciprocal
instance (Fractional n) => Fractional (Dual n) where
    fromRational x = Dual (fromRational x) 0
    recip a = Dual (1 / f a) (- df a / (f a * f a))

-- Scientific functions
instance (Floating n) => Floating (Dual n) where
    pi = Dual pi 0
    exp a = Dual (exp (f a)) (df a * exp (f a))
    log a = Dual (log (f a)) (df a / f a)
    sin a = Dual (sin (f a)) (df a * cos (f a))
    cos a = Dual (cos (f a)) (- df a * sin (f a))
    asin a = Dual (asin (f a)) (df a / sqrt (1 - f a * f a))
    acos a = Dual (acos (f a)) (- df a / sqrt (1 - f a * f a))
    atan a = Dual (atan (f a)) (df a / (1 + f a * f a))
    sinh a = Dual (sinh (f a)) (df a * cosh (f a))
    cosh a = Dual (cosh (f a)) (df a * sinh (f a))
    asinh a = Dual (asinh (f a)) (df a / sqrt (1 + f a * f a))
    acosh a = Dual (acosh (f a)) (df a / sqrt (f a * f a - 1))
    atanh a = Dual (atanh (f a)) (df a / (1 + f a * f a))
