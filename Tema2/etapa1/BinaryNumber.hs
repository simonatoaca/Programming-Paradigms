module BinaryNumber where

import Data.List
import Data.Tuple (swap)

import Debug.Trace

{-
    Reprezentarea unui număr binar ca listă finită sau infinită de biți.
    Primul bit este cel mai puțin semnificativ.

    De exemplu, 6 poate fi reprezentat ca [0, 1, 1, 0, 0, ...].

    type introduce un sinonim de tip, similar cu typedef din C.
-}
type BinaryNumber = [Int]

{-
    *** TODO ***

    Transformă un număr din reprezentarea binară finită în cea zecimală.
    
    Constrângeri: funcția trebuie definită
    * în stil point-free (fără explicitarea parametrului formal)
    * utilizând o funcțională (fără recursivitate explicită).

    Exemple:
    
    > toDecimal [0,1,1]
    6
-}
toDecimal :: BinaryNumber -> Int
-- toDecimal = foldr (\x acc -> 2 * acc + x) 0 - de aici am pornit (regula lui Horner)
toDecimal = foldr (flip ((+) . (2 *))) 0

{-
    *** TODO ***

    Transformă un număr din reprezentarea zecimală în cea binară infinită.

    Constrângeri: pt bonus 10p, funcția trebuie definită
    * în stil point-free (fără explicitarea parametrului formal)
    * utilizând funcționala unfoldr (fără recursivitate explicită).

    Spre deosebire de foldr, care împăturește o listă la o singură valoare,
    unfoldr desfășoară procesul invers, construind o listă pe baza unei valori
    inițiale.
    
    Hint: divMod.

    Exemple:
    
    > take 10 $ toBinary 6
    [0,1,1,0,0,0,0,0,0,0]
-}
toBinary :: Int -> BinaryNumber
-- toBinary = unfoldr (\x -> Just (snd $ divMod x 2, fst $ divMod x 2))
toBinary = unfoldr (Just . swap . (flip divMod 2))

{-
    *** TODO ***

    Incrementează cu 1 reprezentarea finită sau infinită a unui număr binar.
    Atenție la transport!

    Constrângeri: utilizați
    * recursivitate explicită
    * pattern matching.

    Exemple:

    > inc [0,1,1] 
    [1,1,1]

    > inc [1,1,1]
    [0,0,0,1]
-}
inc :: BinaryNumber -> BinaryNumber
inc (0 : xs) = 1 : xs
inc (1 : xs) = 0 : inc xs

{-
   *** TODO ***

    Decrementează cu 1 reprezentarea finită sau infinită a unui număr binar.
    Atenție la împrumut!

    Constrângeri: utilizați
    * recursivitate explicită
    * pattern matching.

    Exemple:

    > dec [1,1,1]
    [0,1,1]

    > dec [0,0,0,1]
    [1,1,1]
-}
dec :: BinaryNumber -> BinaryNumber
dec (1 : xs) = 0 : xs
dec (0 : xs) = 1 : dec xs

{-
    *** TODO ***

    Adună două numere binare, asumând reprezentări infinite, pentru ușurința
    aplicării funcționalelor.

    Constrângeri: utilizați
    * where sau let
    * pt bonus 10p, funcționala mapAccumL (fără recursivitate explicită).

    mapAccumL are tipul (a -> b -> (a, c)) -> a -> [b] -> (a, [c]).
    Așa cum sugerează numele, combină comportamentele funcționalelor:
    * map, transformând element cu element [b] în [c]
    * foldl, utilizând în același timp un acumulator de tipul a.

    Exemple:

    > take 10 $ add (toBinary 74) (toBinary 123)
    [1,0,1,0,0,0,1,1,0,0]
    
    > toDecimal $ take 10 $ add (toBinary 74) (toBinary 123)
    197
-}
add :: BinaryNumber -> BinaryNumber -> BinaryNumber
add bits1 bits2 = snd result ++ [fst result] -- se adauga bitul de carry la final
    where 
        operands = zipWith (+) bits1 bits2
        result = mapAccumL (\acc x -> case x + acc of
                                        0 -> (0, 0)
                                        1 -> (0, 1)
                                        2 -> (1, 0)
                                        3 -> (1, 1)
                                        otherwise -> (0, 0)) 0 operands

{-
    *** TODO ***

    În pregătirea operației de înmulțire a două numere binare, cu reprezentări
    infinite, stivuiește copii deplasate la dreapta ale lui bits1, înmulțite
    cu bit-ul curent din bits2. Deplasarea se face adăugând la stânga lui bits1
    un număr de 0-uri dat de numărul liniei curente. Întoarce o listă infinită
    de liste infinite.

    Vizual:

    0 1 1 0 ... *   <- bits1
    1 0 1 0 ...     <- bits2
    -----------
   |0 1 1 0 ...        înmulțire bits1 cu 1 și deplasare 0 poziții la dreapta
    0|0 0 0 0 ...      înmulțire bits1 cu 0 și deplasare 1 poziție la dreapta
    0 0|0 1 1 0 ...    înmulțire bits1 cu 1 și deplasare 2 poziții la dreapta

    Constrângeri:
    * Corpul funcției trebuie să fie un list comprehension.
    * Nu utilizați recursivitate explicită.

    Hint: iterate pt generarea secvențelor de deplasare spre dreapta cu biți 0.

    Exemple:

    (exemplul vizual)
    > take 3 $ map (take 6) $ stack (toBinary 6) (toBinary 5)
    [[0,1,1,0,0,0],[0,0,0,0,0,0],[0,0,0,1,1,0]]
-}
stack :: BinaryNumber -> BinaryNumber -> [BinaryNumber]
-- ia fiecare lista shiftata si o inmulteste cu bitul corespunzator din bits2
stack bits1 bits2 = [x | x <- zipWith (map . (*)) bits2 (iterate (0:) bits1)]

{-
    *** TODO ***

    Întoarce o listă infinită de numere binare, care pe poziția i >= 0 conține
    rezultatul înmulțirii lui bits1 cu numărul format din primii i biți
    ai lui bits2, i.e. suma primelor i linii întoarse de stack.

    Constrângeri:
    * Utilizați funcționala scanl (fără recursivitate explicită).

    Spre deosebire de foldl, care întoarce acumulatorul final, scanl întoarce
    o listă cu toate acumulatoarele intermediare.

    Exemple:
    
    > take 5 $ map (toDecimal . take 10) $ multiply (toBinary 6) (toBinary 5) 
    [0,6,6,30,30]
-}
multiply :: BinaryNumber -> BinaryNumber -> [BinaryNumber]
multiply bits1 bits2 = scanl (add) [0, 0 ..] bitStack
    where bitStack = stack bits1 bits2

{-
    *** TODO ***

    Întrebare de reflecție, la care veți răspunde la prezentarea temei.

    Având în vedere că liniile întoarse de stack care conțin exclusiv biți 0
    nu contribuie la rezultatul unei înmulțiri, să presupunem că modificăm
    definiția funcției stack astfel încât să construiască doar liniile utile;
    de exemplu, folosind filter sau pattern matching în list comprehension
    pt a păstra doar biții 1 din bits2. Ce probleme ar putea crea această
    abordare?
-}
