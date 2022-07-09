module Data.Array.Accelerate.Arithmetic.Interpolation (
   bisect,
   lookupInterval,
   Interpolator13, sampleBasisFunctions13,
   ) where

import qualified Data.Array.Accelerate.LinearAlgebra.Matrix.Sparse as Sparse
import qualified Data.Array.Accelerate.LinearAlgebra as LinAlg
import qualified Data.Array.Accelerate.Utility.Arrange as Arrange
import qualified Data.Array.Accelerate.Utility.Lift.Exp as Exp
import qualified Data.Array.Accelerate.Utility.Loop as Loop
import Data.Array.Accelerate.LinearAlgebra
          (Scalar, Vector, numElems, extrudeVector, )

import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate (Exp, Any(Any), Z(Z), (:.)((:.)), )

import Data.Ord.HT (limit, )


bisect ::
   (A.Slice ix, A.Shape ix, A.Ord a, A.Elt a) =>
   Vector ix a ->
   Scalar ix a ->
   Scalar ix (Int, Int) ->
   Scalar ix (Int, Int)
bisect nodes xs bounds =
   let centers =
          A.map
             (A.uncurry $ \lower upper -> div (lower+upper) 2)
             bounds
   in  A.zipWith3
          (\center interval leftBranch ->
              A.cond leftBranch
                 (Exp.mapSnd (const center) interval)
                 (Exp.mapFst (const center) interval))
          centers bounds $
       A.zipWith (A.<) xs $
       Arrange.gather (Arrange.mapWithIndex Exp.indexCons centers) nodes

lookupInterval ::
   (A.Slice ix, A.Shape ix, A.Ord a, A.Elt a) =>
   Vector ix a ->
   Scalar ix a ->
   Scalar ix Int
lookupInterval nodes x =
   A.map A.fst $
   Loop.nestLog2 (numElems nodes) (bisect nodes x) $
   A.fill (A.shape x) $
   A.lift (0 :: Exp Int, numElems nodes)


outerVector ::
   (A.Shape ix, A.Slice ix, A.Elt a, A.Elt b, A.Elt c) =>
   (Exp a -> Exp b -> Exp c) ->
   Scalar ix a -> Vector Z b -> Vector ix c
outerVector f x y =
   A.zipWith f
      (A.replicate (A.lift $ Any :. numElems y) x)
      (extrudeVector (A.shape x) y)


{- |
One node before index 0 and three nodes starting from index 0.
-}
type Interpolator13 a = (a,a) -> (a,a) -> (a,a) -> (a,a) -> a -> a

sampleBasisFunctions13 ::
   (A.Slice ix, A.Shape ix, A.Ord a, Num a) =>
   Interpolator13 (Exp a) ->
   Vector Z a -> Vector ix a -> Sparse.Rows ix a
sampleBasisFunctions13 interpolate nodes zs =
   Sparse.Rows (numElems nodes) $
   let indices = lookupInterval (extrudeVector (A.shape zs) nodes) zs
       minIx = 1
       maxIx = numElems nodes - 3
       limitIndices = A.map (limit (minIx, maxIx)) indices
       gatherFromNodes d =
          LinAlg.gatherFromVector (A.map (d+) limitIndices) nodes
   in  outerVector
          (A.lift2 $
           \(n, ln, z, x) (k, y) ->
              case (Exp.unliftQuadruple x, Exp.unliftQuadruple y) of
                 ((xm1,x0,x1,x2), (ym1,y0,y1,y2)) ->
                    (ln+k :: Exp Int,
                     A.cond (n A.< minIx) y0 $
                     A.cond (n A.> maxIx) y1 $
                     interpolate (xm1,ym1) (x0,y0) (x1,y1) (x2,y2) z))
          (A.zip4 indices limitIndices zs
             (A.zip4
                (gatherFromNodes (-1))
                (gatherFromNodes 0)
                (gatherFromNodes 1)
                (gatherFromNodes 2)))
          (A.use $
           A.fromList (Z:.4)
              [(-1, (1,0,0,0)), (0, (0,1,0,0)), (1, (0,0,1,0)), (2, (0,0,0,1))])
