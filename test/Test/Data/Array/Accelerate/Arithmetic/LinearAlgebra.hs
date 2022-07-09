module Test.Data.Array.Accelerate.Arithmetic.LinearAlgebra where

import qualified Data.Array.Accelerate.LinearAlgebra.Private as LinAlg
import qualified Data.Array.Accelerate as A

import Data.Array.Accelerate.LinearAlgebra.Private (Matrix, numCols, )
import Data.Array.Accelerate (Z(Z), (:.)((:.)),)

import Test.Data.Array.Accelerate.Arithmetic.Utility (arbitraryArray, (=!=), )

import qualified Test.QuickCheck as QC

import Data.Word (Word32, )


newtype ArbMatrix a = ArbMatrix (Matrix Z a)

instance (QC.Arbitrary a, A.Elt a) => QC.Arbitrary (ArbMatrix a) where
   arbitrary = do
      nc <- QC.choose (1,100)
      nr <- QC.choose (1,100)
      fmap (ArbMatrix . A.use) $
         arbitraryArray (Z :. nr :. nc) QC.arbitrary


flattenMatrix :: ArbMatrix Word32 -> Bool
flattenMatrix (ArbMatrix m) =
   LinAlg.flattenMatrixReshape m
   =!=
   LinAlg.flattenMatrixBackPermute m

restoreMatrix :: ArbMatrix Word32 -> Bool
restoreMatrix (ArbMatrix m) =
   let v = LinAlg.flattenMatrix m
   in  LinAlg.restoreMatrixReshape (numCols m) v
       =!=
       LinAlg.restoreMatrixBackPermute (numCols m) v

flattenRestoreMatrix :: ArbMatrix Word32 -> Bool
flattenRestoreMatrix (ArbMatrix m) =
   m =!= LinAlg.restoreMatrix (numCols m) (LinAlg.flattenMatrix m)
