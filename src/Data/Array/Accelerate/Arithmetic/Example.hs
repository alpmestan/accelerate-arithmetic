module Data.Array.Accelerate.Arithmetic.Example where

import qualified Data.Array.Accelerate.Arithmetic.Interpolation as Ip
import qualified Data.Array.Accelerate.LinearAlgebra.Matrix.Sparse as Sparse
import Data.Array.Accelerate.LinearAlgebra (Vector, )

import qualified Data.Array.Accelerate.Interpreter as AI
import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate (Array, Z(Z), (:.)((:.)), )


exampleSparseColumnMatrix :: IO ()
exampleSparseColumnMatrix = do
   let m :: Sparse.Columns Z Double
       m =
          Sparse.Columns (A.lift (3::Int)) $
          A.use $ A.fromList (Z :. 2 :. 5) $
          (0,1) : (2,2) : (1,3) : (0,4) : (2,5) :
          (1,6) : (2,7) : (0,8) : (2,9) : (1,10) :
          []

       v :: Vector Z Double
       v = A.use $ A.fromList (Z :. 5) [1,10,100,1000,10000]

   print $ AI.run $ Sparse.multiplyColumnsVector m v

exampleSparseRowMatrix :: IO ()
exampleSparseRowMatrix = do
   let m :: Sparse.Rows Z Double
       m =
          Sparse.Rows (A.lift (5::Int)) $
          A.use $ A.fromList (Z :. 3 :. 2) $
          (0,1) : (0,2) :
          (3,3) : (1,4) :
          (3,5) : (4,6) :
          []

       v :: Vector Z Double
       v = A.use $ A.fromList (Z :. 5) [1,10,100,1000,10000]

   print $ AI.run $ Sparse.multiplyRowsVector m v

exampleLookup :: IO ()
exampleLookup = do
   let nodes :: Array A.DIM2 Double
       nodes = A.fromList (Z :. 3 :. 5) [0 ..]
       x :: Array A.DIM1 Double
       x = A.fromList (Z :. 3) [0.2, 6.7, 13.1]
   print $ AI.run1 (Ip.lookupInterval (A.use nodes)) x
