module Taylor where

-- "Taylor" is the "Dual Number" type, only recursive - which we may call
-- "Taylor numbers". The first part is a normal function and the second part is
-- the derivative of the function. The advantage of this over the dmath.hs is
-- that you should be able to pass Taylor numbers wherever you write normal
-- number functions and get derivatives also calculated. This includes
-- if-then-else conditionals as well.
data Taylor x = Taylor x (Taylor x)

-- f adn = the function part
-- df adn = the derivative part
f (Taylor f0 _) = f0
df (Taylor _ df0) = df0

dfn 0 = f 
dfn n = dfn (n-1) . df

-- You can provide values for the x variable using this function.
-- If you have a normal mathematical function to calculate some
-- value from a given number x, you can evaluate the same function
-- by passing a "dual number" as a value and have both the value
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
