{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE FlexibleContexts #-}
module Data.Array.Accelerate.LinearAlgebra.Matrix.Sparse (
   Columns(..),
   multiplyColumnsVector,
   transposeColumns,
   Rows(..),
   multiplyRowsVector,
   transposeRows,
   multiplyColumnsRows,
   -- realBandedGramian,
   scaleRowRows,
   ) where

import qualified Data.Array.Accelerate.LinearAlgebra.Matrix.Banded as BandMatrix
import qualified Data.Array.Accelerate.LinearAlgebra as LinAlg
import qualified Data.Array.Accelerate.Utility.Lift.Exp as Exp
import qualified Data.Array.Accelerate.Utility.Arrange as Arrange
import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate.Utility.Lift.Exp (expr, )

import Data.Array.Accelerate.LinearAlgebra (Matrix, Vector, matrixShape, )
import Data.Array.Accelerate (Exp, Any(Any), All(All), (:.)((:.)), (?), )


{- |
Sparse matrix with a definite number of non-zero entries per column.
-}
data Columns ix a =
        Columns {numRows :: Exp Int, columnMatrix :: Matrix ix (Int, a)}

realIndex ::
   (A.Shape ix, A.Slice ix, A.Elt a) =>
   Matrix ix (Int, a) ->
   Matrix ix (ix :. Int)
realIndex m =
   A.zipWith Exp.indexCons
      (A.generate (A.shape m) (A.indexTail . A.indexTail))
      (A.map A.fst m)

multiplyColumnsVector ::
   (A.Shape ix, A.Slice ix, A.Num a) =>
   Columns ix a ->
   Vector ix a ->
   Vector ix a
multiplyColumnsVector (Columns rows m) v =
   Arrange.scatter (+)
      (realIndex m)
      (case matrixShape m of
          sh :. _rows :. _cols -> A.fill (A.lift $ sh :. rows) 0) $
   A.zipWith (*)
      (A.map A.snd m)
      (A.replicate (A.lift $ Any :. LinAlg.numRows m :. All) v)

transposeColumns ::
   (A.Shape ix, A.Slice ix, A.Num a) =>
   Columns ix a ->
   Rows ix a
transposeColumns (Columns n x) =
   Rows n $ LinAlg.transpose x


{- |
Sparse matrix with a definite number of non-zero entries per row.
-}
data Rows ix a =
        Rows {numCols :: Exp Int, rowMatrix :: Matrix ix (Int, a)}

multiplyRowsVector ::
   (A.Shape ix, A.Slice ix, A.Num a) =>
   Rows ix a ->
   Vector ix a ->
   Vector ix a
multiplyRowsVector (Rows _cols m) v =
   A.fold1 (+) $
   A.zipWith (*) (A.map A.snd m) $
   Arrange.gather (realIndex m) v

transposeRows ::
   (A.Shape ix, A.Slice ix, A.Num a) =>
   Rows ix a ->
   Columns ix a
transposeRows (Rows n x) =
   (Columns n $ LinAlg.transpose x)

multiplyColumnsRows ::
   (A.Shape ix, A.Slice ix, A.Num a) =>
   Columns ix a ->
   Rows ix a ->
   Matrix ix a
multiplyColumnsRows (Columns rows x) (Rows cols y) =
   let (ixs,prods) = A.unzip $ matchMatrices x y
       global = A.indexTail . A.indexTail . A.indexTail
   in  Arrange.scatter (+)
          (Arrange.mapWithIndex
             (Exp.modify2 expr (expr,expr) $ \mix (k,j) ->
                global mix :. k :. j) $
           ixs)
          (A.fill (A.lift $ global (A.shape prods) :. rows :. cols) 0)
          prods

{- |
Compute x^T*x, given that it has a band structure.
You must pass the band-width as parameter
and you must make sure that the Gramian stays within this band.
Otherwise you cause out-of-bounds array accesses.
So far, only correct for real matrices.
-}
-- realBandedGramian ::
--    (A.Shape ix, A.Slice ix, A.Num a) =>
--    Exp Int ->
--    Rows ix a ->
--    BandMatrix.Symmetric ix a
-- realBandedGramian width (Rows cols y) =
--    let (ixs,prods) = A.unzip $ matchMatrices (LinAlg.transpose y) y
--        global = A.indexTail . A.indexTail . A.indexTail
--    in  BandMatrix.Symmetric $
--        Arrange.scatter (+)
--           (Arrange.mapWithIndex
--              (Exp.modify2 expr (expr,expr) $ \mix (k,j) ->
--                 k A.> j ? (A.ignore, A.lift $ global mix :. k :. j-k)) $
--            ixs)
--           (A.fill (A.lift $ global (A.shape prods) :. cols :. width) 0)
--           prods

matchMatrices ::
   (A.Shape ix, A.Slice ix, A.Num a) =>
   Matrix ix (Int, a) ->
   Matrix ix (Int, a) ->
   Matrix (ix :. Int) ((Int, Int), a)
matchMatrices x y =
   case (matrixShape x, matrixShape y) of
      (_ :. xRows :. _xCols, _ :. _yRows :. yCols) ->
         -- it must be xCols == yRows
         A.zipWith
            (Exp.modify2 (expr,expr) (expr,expr) $
             \(n,xi) (m,yi) -> ((n, m), xi*yi))
            (A.replicate (A.lift $ Any :. All :. All :. yCols) x)
            (A.replicate (A.lift $ Any :. xRows :. All :. All) y)


scaleRowRows ::
   (A.Slice ix, A.Shape ix, A.Num a) =>
   Vector ix a -> Rows ix a -> Rows ix a
scaleRowRows s (Rows n x) =
   Rows n $
   LinAlg.zipScalarVectorWith (\si xi -> Exp.mapSnd (si*) xi) s x
