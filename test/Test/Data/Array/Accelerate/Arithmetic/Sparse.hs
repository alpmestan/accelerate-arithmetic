{-# LANGUAGE Rank2Types #-}
{-# LANGUAGE TypeOperators #-}
module Test.Data.Array.Accelerate.Arithmetic.Sparse where

import Test.Data.Array.Accelerate.Arithmetic.Utility (arbitraryArray, (=!=), )

import qualified Data.Array.Accelerate.LinearAlgebra.Matrix.Banded as BandMatrix
import qualified Data.Array.Accelerate.LinearAlgebra.Matrix.Sparse as Sparse
import qualified Data.Array.Accelerate.LinearAlgebra as LinAlg

import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate (Exp, Z(Z), (:.)((:.)))

import qualified Test.QuickCheck as QC

import Control.Monad (liftM2, )

import Data.Word (Word32, )


data
   CRVTriple a =
      CRVTriple
         (Sparse.Columns Z a)
         (Sparse.Rows Z a)
         (LinAlg.Vector Z a)

instance (QC.Arbitrary a, A.Elt a) => QC.Arbitrary (CRVTriple a) where
   arbitrary = do
      k <- QC.choose (1,200)
      nc <- QC.choose (1,100)
      nr <- QC.choose (1,100)
      cc <- QC.choose (1,10)
      cr <- QC.choose (1,10)
      mc <-
         arbitraryArray (Z :. cc :. k) $
         liftM2 (,) (QC.choose (0,nc-1)) QC.arbitrary
      mr <-
         arbitraryArray (Z :. k :. cr) $
         liftM2 (,) (QC.choose (0,nr-1)) QC.arbitrary
      v <- arbitraryArray (Z :. nr) QC.arbitrary
      return $
         CRVTriple
            (Sparse.Columns (A.lift nc) (A.use mc))
            (Sparse.Rows (A.lift nr) (A.use mr))
            (A.use v)


multiplication :: CRVTriple Word32 -> Bool
multiplication (CRVTriple mc mr v) =
   LinAlg.multiplyMatrixVector (Sparse.multiplyColumnsRows mc mr) v
   =!=
   Sparse.multiplyColumnsVector mc (Sparse.multiplyRowsVector mr v)



data BandGramian a = BandGramian (Exp Int) (Sparse.Rows Z a)

instance (QC.Arbitrary a, A.Elt a) => QC.Arbitrary (BandGramian a) where
   arbitrary = do
      width <- QC.choose (1,10)
      rows <- QC.choose (1,100)
      cols <- QC.choose (width,100)

      m <-
         fmap (A.fromList (Z :. rows :. width) . concat) $
         QC.vectorOf rows $
         liftM2
            (\start row -> zip [start..] row)
            (QC.choose (0,cols-width))
            (QC.vectorOf width QC.arbitrary)

      return $
         BandGramian (A.lift width)
            (Sparse.Rows (A.lift cols) (A.use m))


bandedGramian :: BandGramian Word32 -> Bool
bandedGramian (BandGramian width m) =
   Sparse.multiplyColumnsRows (Sparse.transposeRows m) m
   =!=
   BandMatrix.flattenSymmetric (Sparse.realBandedGramian width m)
