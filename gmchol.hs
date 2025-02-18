import Data.List (foldl')
import Data.Array
import Data.Maybe (fromMaybe)
import Control.Monad (forM_)

gmchol :: Array (Int, Int) Double -> Int -> Array (Int, Int) Double
gmchol a n = runSTUArray $ do
    r <- newArray ((0, 0), (n-1, n-1)) 0.0
    e <- newArray ((0, 0), (n-1, n-1)) 0.0

    let normA = maximum [sum [abs (a ! (i, j)) | i <- [0..n-1]] | j <- [0..n-1]]
        gamm = maximum [abs (a ! (i, i)) | i <- [0..n-1]]
        delta = max (epsilon * normA) epsilon

    forM_ [0..n-1] $ \j -> do
        let thetaJ = foldl' (\acc i -> let sumVal = sum [r ! (k, i) * r ! (k, j) | k <- [0..i-1]] in
                                          acc `max` (a ! (i, j) - sumVal)) 0.0 [0..n-1]
        forM_ [0..n-1] $ \i -> do
            let sumVal = sum [r ! (k, i) * r ! (k, j) | k <- [0..i-1]]
            writeArray r (i, j) ((a ! (i, j) - sumVal) / (r ! (i, i)))
            when (i > j) $ writeArray r (i, j) 0.0

        let sumVal = sum [r ! (k, j) * r ! (k, j) | k <- [0..j-1]]
            phiJ = a ! (j, j) - sumVal
            xiJ = if (j + 1) < n then maximum [abs (a ! (i, j)) | i <- [j+1..n-1]] else abs (a ! (n-1, j))
            betaJ = sqrt (max gamm (max (xiJ / fromIntegral n) epsilon))

        eVal <- readArray e (j, j)
        let eJ = if delta >= max (abs phiJ) ((thetaJ * thetaJ) / (betaJ * betaJ))
                 then delta - phiJ
                 else if abs phiJ >= max ((delta * delta) / (betaJ * betaJ)) delta
                      then abs phiJ - phiJ
                      else if ((thetaJ * thetaJ) / (betaJ * betaJ)) >= max delta (abs phiJ)
                           then ((thetaJ * thetaJ) / (betaJ * betaJ)) - phiJ
                           else 0.0
        writeArray e (j, j) eJ
        writeArray r (j, j) (sqrt (a ! (j, j) - sumVal + eJ))

    return r

epsilon :: Double
epsilon = 1.1102230246251565e-16

main :: IO ()
main = do
    let n = 3
        a = array ((0, 0), (n-1, n-1)) [((0, 0), 4), ((0, 1), 12), ((0, 2), -16),
                                          ((1, 0), 12), ((1, 1), 37), ((1, 2), -43),
                                          ((2, 0), -16), ((2, 1), -43), ((2, 2), 98)]
        r = gmchol a n

    putStrLn "R matrix:"
    forM_ [0..n-1] $ \i -> do
        forM_ [0..n-1] $ \j -> do
            putStr (show (r ! (i, j)) ++ " ")
        putStrLn ""

