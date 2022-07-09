module Main where

import qualified Test.Data.Array.Accelerate.Arithmetic.Sparse as Sparse
import qualified Test.Data.Array.Accelerate.Arithmetic.LinearAlgebra as LinAlg

import qualified Test.QuickCheck.Modifiers as Mod
import Test.QuickCheck (quickCheck)


test :: IO ()
test = mapM_ (\(msg,act) -> putStr (msg++": ") >> act) $
   ("sparseMatrix", quickCheck (\(Mod.Blind x) -> Sparse.multiplication x)) :
   ("bandedGramian", quickCheck (\(Mod.Blind x) -> Sparse.bandedGramian x)) :
   ("flattenMatrix", quickCheck (\(Mod.Blind x) -> LinAlg.flattenMatrix x)) :
   ("restoreMatrix", quickCheck (\(Mod.Blind x) -> LinAlg.restoreMatrix x)) :
   ("flattenRestoreMatrix", quickCheck (\(Mod.Blind x) -> LinAlg.flattenRestoreMatrix x)) :
   []


main :: IO ()
main = test
