{-# LANGUAGE FlexibleContexts #-}
module Data.Array.Accelerate.LinearAlgebra.Matrix.Banded (
   Symmetric(..),
   flattenSymmetric,
   ) where

import Data.Array.Accelerate.LinearAlgebra (Matrix, matrixShape)

import qualified Data.Array.Accelerate.Utility.Lift.Exp as Exp
import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate.Utility.Lift.Exp (expr)
import Data.Array.Accelerate ((:.)((:.)), (!), (?))


newtype Symmetric ix a = Symmetric (Matrix ix a)

flattenSymmetric ::
   (A.Slice ix, A.Shape ix, A.Num a) =>
   Symmetric ix a -> Matrix ix a
flattenSymmetric (Symmetric m) =
   case matrixShape m of
      (sh :. rows :. width) ->
         A.generate (A.lift $ sh :. rows :. rows) $
         Exp.modify (expr:.expr:.expr) $ \(ix:.k0:.j0) ->
            let k = min k0 j0
                j = max k0 j0 - k
            in  width A.> j ? (m ! A.lift(ix:.k:.j), 0)
