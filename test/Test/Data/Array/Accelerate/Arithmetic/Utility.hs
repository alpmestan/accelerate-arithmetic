{-# LANGUAGE Rank2Types #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE FlexibleInstances #-}
module Test.Data.Array.Accelerate.Arithmetic.Utility where

import qualified Data.Array.Accelerate.Interpreter as AI
import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate (Acc, Array, Z, (:.)((:.)), )

import qualified Test.QuickCheck as QC


infix 4 =!=

(=!=) ::
   (Eq sh, Eq e, A.Shape sh, A.Elt e) =>
   Acc (Array sh e) -> Acc (Array sh e) -> Bool
x =!= y  =
   let xi = AI.run x
       yi = AI.run y
   in  A.arrayShape xi == A.arrayShape yi
       &&
       A.toList xi == A.toList yi


class A.Shape sh => Shape sh where
   switchShape ::
      f Z ->
      (forall sh1. Shape sh1 => f (sh1 :. Int)) ->
      f sh

instance Shape Z where switchShape f _ = f
instance Shape sh => Shape (sh :. Int) where switchShape _ f = f

newtype ShapeSize sh = ShapeSize {getShapeSize :: sh -> Int}

shapeSize :: Shape sh => sh -> Int
shapeSize =
   getShapeSize $
   switchShape
      (ShapeSize $ const 1)
      (ShapeSize $ \(ix:.n) -> shapeSize ix * n)

arbitraryArray ::
   (Shape sh, A.Elt a) => sh -> QC.Gen a -> QC.Gen (Array sh a)
arbitraryArray sh gen =
   fmap (A.fromList sh) $
   QC.vectorOf (shapeSize sh) gen
