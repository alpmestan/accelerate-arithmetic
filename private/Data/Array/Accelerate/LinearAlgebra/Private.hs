{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE FlexibleContexts #-}
module Data.Array.Accelerate.LinearAlgebra.Private where

import qualified Data.Array.Accelerate.Utility.Loop as Loop
import qualified Data.Array.Accelerate.Utility.Lift.Exp as Exp
import qualified Data.Array.Accelerate.Utility.Arrange as Arrange
import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate
          (Acc, Array, Exp, Any(Any), All(All), Z(Z), (:.)((:.)))



type Scalar ix a = Acc (Array ix a)
type Vector ix a = Acc (Array (ix :. Int) a)
type Matrix ix a = Acc (Array (ix :. Int :. Int) a)

transpose ::
   (A.Shape ix, A.Slice ix, A.Elt a) =>
   Matrix ix a -> Matrix ix a
transpose m =
   A.backpermute
      (A.lift $ swapIndex $ matrixShape m)
      (A.lift . swapIndex . A.unlift)
      m

swapIndex ::
   Exp ix :. Exp Int :. Exp Int ->
   Exp ix :. Exp Int :. Exp Int
swapIndex (ix :. r :. c) = (ix :. c :. r)


numElems :: (A.Shape ix, A.Slice ix, A.Elt a) => Vector ix a -> Exp Int
numElems m = case vectorShape m of _ix :. n -> n

numRows :: (A.Shape ix, A.Slice ix, A.Elt a) => Matrix ix a -> Exp Int
numRows m = case matrixShape m of _ix :. rows :. _cols -> rows

numCols :: (A.Shape ix, A.Slice ix, A.Elt a) => Matrix ix a -> Exp Int
numCols m = case matrixShape m of _ix :. _rows :. cols -> cols

vectorShape ::
   (A.Shape ix, A.Slice ix, A.Elt a) =>
   Vector ix a -> Exp ix :. Exp Int
vectorShape m = A.unlift $ A.shape m

matrixShape ::
   (A.Shape ix, A.Slice ix, A.Elt a) =>
   Matrix ix a -> Exp ix :. Exp Int :. Exp Int
matrixShape m = A.unlift $ A.shape m

withVectorIndex ::
   (A.Shape ix, A.Slice ix, A.Lift Exp a) =>
   (Exp ix :. Exp Int -> a) ->
   (Exp (ix :. Int) -> Exp (A.Plain a))
withVectorIndex f = A.lift . f . A.unlift

withMatrixIndex ::
   (A.Shape ix, A.Slice ix, A.Lift Exp a) =>
   (Exp ix :. Exp Int :. Exp Int -> a) ->
   (Exp (ix :. Int :. Int) -> Exp (A.Plain a))
withMatrixIndex f = A.lift . f . A.unlift


outer ::
   (A.Shape ix, A.Slice ix, A.Num a) =>
   Vector ix a -> Vector ix a -> Matrix ix a
outer x y =
   A.zipWith (*)
      (A.replicate (A.lift $ Any :. All :. numElems y) x)
      (A.replicate (A.lift $ Any :. numElems x :. All) y)

multiplyMatrixVector ::
   (A.Shape ix, A.Slice ix, A.Num a) =>
   Matrix ix a ->
   Vector ix a ->
   Vector ix a
multiplyMatrixVector m v =
   case matrixShape m of
      (_ix :. rows :. _cols) ->
         A.fold1 (+) $
         A.zipWith (*) m
            (A.replicate (A.lift $ Any :. rows :. All) v)

multiplyMatrixMatrix ::
   (A.Shape ix, A.Slice ix, A.Num a) =>
   Matrix ix a ->
   Matrix ix a ->
   Matrix ix a
multiplyMatrixMatrix x y =
   case (matrixShape x, matrixShape y) of
      (_ :. rows :. _cols, _ :. _rows :. cols) ->
         A.fold1 (+) $ transpose $
         A.zipWith (*)
            (A.replicate (A.lift $ Any :. All :. All :. cols) x)
            (A.replicate (A.lift $ Any :. rows :. All :. All) y)

newtonInverseStep ::
   (A.Shape ix, A.Slice ix, A.Num a) =>
   Matrix ix a ->
   Matrix ix a ->
   Matrix ix a
newtonInverseStep a x =
   A.zipWith (-) (A.map (2*) x) $
   multiplyMatrixMatrix x $ multiplyMatrixMatrix a x

identity ::
   (A.Shape ix, A.Slice ix, A.Elt a, A.FromIntegral Int a) =>
   Exp (ix :. Int :. Int) -> Matrix ix a
identity sh =
   A.generate sh
      (withMatrixIndex $
       \(_ :. r :. c) -> A.fromIntegral $ A.boolToInt (r A.== c))

newtonInverse ::
   (A.Shape ix, A.Slice ix, A.Num a) =>
   Exp Int ->
   Matrix ix a ->
   Matrix ix a ->
   Matrix ix a
newtonInverse n seed a =
   Loop.nest n (newtonInverseStep a) seed



scaleRows ::
   (A.Slice ix, A.Shape ix, A.Num a) =>
   Vector ix a -> Matrix ix a -> Matrix ix a
scaleRows s x =
   zipScalarVectorWith (*) s x



zipScalarVectorWith ::
   (A.Slice ix, A.Shape ix, A.Elt a, A.Elt b, A.Elt c) =>
   (Exp a -> Exp b -> Exp c) ->
   Scalar ix a -> Vector ix b -> Vector ix c
zipScalarVectorWith f x ys =
   case vectorShape ys of
      _ix :. dim ->
         A.zipWith f (A.replicate (A.lift (Any :. dim)) x) ys

zipScalarMatrixWith ::
   (A.Slice ix, A.Shape ix, A.Elt a, A.Elt b, A.Elt c) =>
   (Exp a -> Exp b -> Exp c) ->
   Scalar ix a -> Matrix ix b -> Matrix ix c
zipScalarMatrixWith f x ys =
   case matrixShape ys of
      _ix :. rows :. cols ->
         A.zipWith f
            (A.replicate (A.lift (Any :. rows :. cols)) x) ys



columnFromVector ::
   (A.Shape ix, A.Slice ix, A.Elt a) =>
   Vector ix a -> Matrix ix a
columnFromVector a = A.reshape (Exp.indexCons (A.shape a) 1) a

{- |
input must be a matrix with exactly one column
-}
vectorFromColumn ::
   (A.Shape ix, A.Slice ix, A.Elt a) =>
   Matrix ix a -> Vector ix a
vectorFromColumn a = A.reshape (A.indexTail $ A.shape a) a



flattenMatrix, flattenMatrixReshape, flattenMatrixBackPermute ::
   (A.Slice ix, A.Shape ix, A.Elt a) =>
   Matrix ix a -> Vector ix a
flattenMatrix = flattenMatrixBackPermute

flattenMatrixReshape m =
   case matrixShape m of
      ix :. rows :. cols ->
         A.reshape (A.lift $ ix :. rows*cols) m

accDivMod :: Integral a => a -> a -> (a, a)
accDivMod x y = (div x y, mod x y)

flattenMatrixBackPermute m =
   case matrixShape m of
      ix :. rows :. cols ->
         A.backpermute
            (A.lift $ ix :. rows*cols)
            (withVectorIndex $
             \(vix :. n) -> case accDivMod n cols of (r,c) -> vix :. r :. c)
            m


restoreMatrix, restoreMatrixReshape, restoreMatrixBackPermute ::
   (A.Slice ix, A.Shape ix, A.Elt a) =>
   Exp Int -> Vector ix a -> Matrix ix a
restoreMatrix = restoreMatrixBackPermute

restoreMatrixReshape cols v =
   case vectorShape v of
      ix :. n ->
         A.reshape (A.lift $ ix :. div n cols :. cols) v

restoreMatrixBackPermute cols v =
   case vectorShape v of
      ix :. n ->
         A.backpermute
            (A.lift $ ix :. div n cols :. cols)
            (withMatrixIndex $ \(vix :. k :. j) -> vix :. k*cols+j)
            v



extrudeVector ::
   (A.Shape ix, A.Slice ix, A.Elt a) =>
   Exp ix -> Vector Z a -> Vector ix a
extrudeVector shape y =
   -- A.replicate (A.lift $ shape :. All) y
   A.backpermute
      (A.lift $ shape :. numElems y)
      (A.index1 . A.indexHead)
      y

extrudeMatrix ::
   (A.Shape ix, A.Slice ix, A.Elt a) =>
   Exp ix -> Matrix Z a -> Matrix ix a
extrudeMatrix shape y =
   A.backpermute
      (A.lift $ shape :. numRows y :. numCols y)
      (withMatrixIndex $ \(_:.r:.c) -> Z:.r:.c)
      y

zipExtrudedVectorWith ::
   (A.Slice ix, A.Shape ix, A.Elt a, A.Elt b, A.Elt c) =>
   (Exp a -> Exp b -> Exp c) ->
   Vector Z a ->
   Vector ix b ->
   Vector ix c
zipExtrudedVectorWith f x y =
   A.zipWith f (extrudeVector (A.indexTail $ A.shape y) x) y

zipExtrudedMatrixWith ::
   (A.Slice ix, A.Shape ix, A.Elt a, A.Elt b, A.Elt c) =>
   (Exp a -> Exp b -> Exp c) ->
   Matrix Z a ->
   Matrix ix b ->
   Matrix ix c
zipExtrudedMatrixWith f x y =
   A.zipWith f (extrudeMatrix (A.indexTail $ A.indexTail $ A.shape y) x) y

gatherFromVector ::
   (A.Shape ix, A.Elt a) =>
   Scalar ix Int -> Vector Z a -> Scalar ix a
gatherFromVector indices =
   Arrange.gather (A.map A.index1 indices)
