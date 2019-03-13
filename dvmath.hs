module DVMath where

-- "DF" is the same old "Differentiable Function", but now it is a function
-- from one vector space to another. We can model vector spaces loosely as a
-- mapping from an "index type" to a "field type" ... rather like a collection
-- type itself. For now, we'll stick to the function representation to
-- understand how it might work, but we'll pull it out as a typedef to indicate
-- that the representation need not be frozen. The "field", for the moment, is
-- also frozen to real numbers.
--
-- Since the derivative of such a vector function needs to be taken relative to
-- a dimension identified by an index, the second part is itself a function
-- from index to a `DF`.
--
-- Note: This whole module doesn't intend to talk about efficient algorithms
-- and only attempts to model the problem at hand conceptually. We'll need to
-- come up with efficient representations.
--
-- Disclaimer: This is still in exploratory mode.
type Vec i x = i -> x
data DF i j x = DF (Vec i x -> Vec j x) (i -> DF i j x)

-- f adf = the function part
-- df i adf = the derivative part
f (DF f0 _) = f0
df i (DF _ df0) = df0 i

-- It is useful to also be able to look at the derivative as though it
-- increases the "rank" of the function we're dealing with.
dft :: (Eq i, Eq j) => DF i j x -> DF i (j,i) x
dft a = DF f' df'
    where 
        f' x (j,i) = f (df i a) x j
        df' i = dft (df i a)

-- Recursive definition of nth derivative. The "nth derivative" is no longer
-- just specifiable with an "n", but needs a list of indices w.r.t. to take
-- the derivative in sequence. For the moment, we don't worry about modeling
-- the fact that such a derivative is independent of the specific sequence
-- of indices used. In other words, dfn ixs = dfn (sort ixs).
dfn [] = id
dfn (i:ixs) = dfn ixs . df i

dirac :: (Eq i, Num x) => i -> i -> x
dirac i j = if i == j then 1 else 0

-- By picking out a specific value of the input argument, we're in effect
-- turning a particular component of a vector into a scalar.
v :: (Eq i, Eq j, Real x) => i -> DF i j x
v i = DF (\x _ -> x i) (dirac i)

-- The whole input beast can itself be referred to, as usual, using a single
-- "x" definition.
x :: (Eq i, Real x) => DF i i x
x = DF id (\j -> DF (\_ -> dirac j) (\_ -> 0))

-- Utility to help take a list of values and make a constant DF from it.
fromList :: (Eq i, Real x) => [x] -> DF i Int x
fromList l = DF (\_ -> (l !!)) (\_ -> 0)

{--
 - Mapping functions over the underlying field onto vectors is straightforward,
 - and is essentially an expression of the `fmap` of a functor. The two-argument
 - operations need some attention though.
 --}

-- Basic arithmetic operations
instance (Eq i, Eq j, Real n) => Num (DF i j n) where
    a + b = DF (\x i -> f a x i + f b x i) (\i -> df i a + df i b)
    a * b = DF (\x i -> f a x i * f b x i) (\i -> a * df i b + df i a * b)
    negate a = DF (\x i -> negate (f a x i)) (\i -> negate (df i a))
    abs a = DF (\x i -> abs (f a x i)) (\i -> df i a * signum a)
    signum a = DF (\x i -> signum (f a x i)) (\i -> 0)
    fromInteger x = DF (\_ _ -> fromInteger x) (\_ -> 0)

-- Reciprocal
instance (Real n, Fractional n, Eq i, Eq j) => Fractional (DF i j n) where
    fromRational x = DF (\_ _ -> fromRational x) (\_ -> 0)
    recip a = DF (\x i -> 1 / f a x i) (\i -> - df i a / (a * a))

-- Scientific functions
instance (Real n, Floating n, Eq i, Eq j) => Floating (DF i j n) where
    pi = DF (\_ _ -> pi) (\_ -> 0)
    exp a = DF (\x i -> exp (f a x i)) (\i -> df i a * exp a)
    log a = DF (\x i -> log (f a x i)) (\i -> df i a / a)
    sin a = DF (\x i -> sin (f a x i)) (\i -> df i a * cos a)
    cos a = DF (\x i -> cos (f a x i)) (\i -> - df i a * sin a)
    asin a = DF (\x i -> asin (f a x i)) (\i -> df i a / sqrt (1 - a * a))
    acos a = DF (\x i -> acos (f a x i)) (\i -> - df i a / sqrt (1 - a * a))
    atan a = DF (\x i -> atan (f a x i)) (\i -> df i a / (1 + a * a))
    sinh a = DF (\x i -> sinh (f a x i)) (\i -> df i a * cosh a)
    cosh a = DF (\x i -> cosh (f a x i)) (\i -> df i a * sinh a)
    asinh a = DF (\x i -> asinh (f a x i)) (\i -> df i a / sqrt (1 + a * a))
    acosh a = DF (\x i -> acosh (f a x i)) (\i -> df i a / sqrt (a * a - 1))
    atanh a = DF (\x i -> atanh (f a x i)) (\i -> df i a / (1 + a * a))

-- While using if-then-else is not so straightforward and needs vectors
-- to be processed into scalars before doing that, we can parameterize
-- conditionals using a "region" function. While the region function takes
-- two arguments, it is free to ignore one if it so chooses.
cond :: (Vec i x -> j -> Bool) -> DF i j x -> DF i j x -> DF i j x
cond region a b = DF f' df'
    where
        f' x j = if region x j then f a x j else f b x j
        df' i = cond region (df i a) (df i b)

-- The famous relu function (rectified linear unit) can be expressed as a such
-- a conditional.
relu a = cond (\x j -> f a x j < 0) 0 a
                
-- There are many kinds of products we can create with vectors.  The outer
-- product increases the rank of the vectors and is a useful operation before
-- many kinds of reductions.
outer :: (Eq i, Eq j, Eq k, Real x) => DF i j x -> DF i k x -> DF i (j,k) x
outer a b = DF f' df'
    where
        f' x (i,j) = f a x i * f b x j
        df' i = outer (df i a) b + outer a (df i b)

-- The inner product usually ends up reducing the rank of the input vectors by
-- summing over a part of the index space. To generalize this idea, we just
-- parameterize the index range of the summation into an enumeration function
-- named `dot` in the argument.
inner :: (Eq i, Eq j, Eq k, Eq l, Real x) => (l -> Int -> Maybe (j,k)) -> DF i j x -> DF i k x -> DF i l x
inner dot a b = DF f' df'
    where
        f' x l = sum' x (dot l) 0 0
        sum' x dot ix result = case dot ix of
            Nothing -> result
            Just (j,k) -> sum' x dot (ix+1) (result + f a x j * f b x k)
        df' i = inner dot (df i a) b + inner dot a (df i b)

-- An inner product between two vectors can be computed as a summation
-- reduction of the outer product of the vectors. Here we code up "collapse",
-- which does such a summation reduction. Like `inner`, `collapse` also reduces
-- the rank of the input.
collapse :: (Eq i, Eq j, Eq k, Real x) => (k -> Int -> Maybe j) -> DF i j x -> DF i k x
collapse dot a = DF f' df'
    where
        f' x k = sum' x (dot k) 0 0
        sum' x dot ix result = case dot ix of
            Nothing -> result
            Just j -> sum' x dot (ix+1) (result + f a x j)
        df' i = collapse dot (df i a)

-- Vector fuction composition. The derivative is expressed using the chain rule.
chain :: (Eq i, Eq j, Eq k, Real x) => (Int -> Maybe j) -> DF j k x -> DF i j x -> DF i k x
chain js a b = DF f' df'
    where
        f' x = f a (f b x)
        df' i = inner (dot i) (chain js (dft a) b) (dft b)
        dot i k ix = case js ix of
            Nothing -> Nothing
            Just j -> Just ((k,j),(j,i))

type Slice j i = i -> Maybe j

-- A utility to take a slice of a vector. For simplicity, we model the slice
-- operation as a boolean selector over the index space.
slice :: (Eq i, Eq j, Eq k, Real x) => Slice j k -> DF i j x -> DF i k x
slice s a = DF f' df'
    where
        f' x k = case s k of
            Nothing -> 0
            Just j -> f a x j
        df' i = slice s (df i a)

stride start step i = Just (start + i * step)
range min max i = if i >= min && i <= max then Just i else Nothing

-- Sometimes, it is also useful to be able to change the shape of
-- a vector .. which basically means we change the way its dimensions
-- are addressed.
reshape :: (Eq i, Eq j, Eq k, Real x) => (k -> j) -> DF i j x -> DF i k x
reshape shaper a = DF f' df'
    where
        f' x k = f a x (shaper k)
        df' i = reshape shaper (df i a)

-- A dead simple notion of convolution as a reduction operation that does not
-- result in reduction of rank unlike the inner product.
conv :: (Eq i, Eq j, Eq k, Real x) => (k -> [(j,k)]) -> DF i j x -> DF i k x -> DF i k x
conv stride kernel a = DF f' df'
    where
        f' x k = sum (map (\(j,k) -> f kernel x j * f a x k) (stride k))
        df' i = conv stride (df i kernel) a + conv stride kernel (df i a)
        

