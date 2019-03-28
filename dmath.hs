module DMath where

-- "DF" is "Differentiable function". The first part is a normal function and
-- the second part is the derivative of the function, which is itself a
-- differentiable function.  This enables us to calculate as many derivatives
-- as we want.
data DF x = DF (x -> x) (DF x)

-- f adf = the function part
-- df adf = the derivative part
f (DF f0 _) = f0
df (DF _ df0) = df0

-- Recursive definition of the nth derivative
dfn 0 = id
dfn n = dfn (n-1) . df

-- The identity function
x :: (Real a) => DF a
x = DF id 1

-- Basic arithmetic operations
instance (Real n) => Num (DF n) where
    a + b = DF (\x -> f a x + f b x) (df a + df b)
    a * b = DF (\x -> f a x * f b x) (a * df b + df a * b)
    negate a = DF (negate . f a) (negate (df a))
    abs a = DF (abs . f a) (df a * signum a)
    signum a = DF (signum . f a) 0 
    fromInteger x = DF (\_ -> fromInteger x) 0

-- Reciprocal
instance (Real n, Fractional n) => Fractional (DF n) where
    fromRational x = DF (\_ -> fromRational x) 0
    recip a = DF (\x -> 1 / f a x) (- df a / (a * a))

-- Scientific functions
instance (Real n, Floating n) => Floating (DF n) where
    pi = DF (\_ -> pi) 0
    exp a = DF (exp . f a) (df a * exp a)
    log a = DF (log . f a) (df a / a)
    sin a = DF (sin . f a) (df a * cos a)
    cos a = DF (cos . f a) (- df a * sin a)
    asin a = DF (asin . f a) (df a / sqrt (1 - a * a))
    acos a = DF (acos . f a) (- df a / sqrt (1 - a * a))
    atan a = DF (atan . f a) (df a / (1 + a * a))
    sinh a = DF (sinh . f a) (df a * cosh a)
    cosh a = DF (cosh . f a) (df a * sinh a)
    asinh a = DF (asinh . f a) (df a / sqrt (1 + a * a))
    acosh a = DF (acosh . f a) (df a / sqrt (a * a - 1))
    atanh a = DF (atanh . f a) (df a / (1 + a * a))
