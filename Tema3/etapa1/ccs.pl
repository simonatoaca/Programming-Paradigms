:- dynamic detailed_mode_disabled/0.
:- ensure_loaded('files.pl').



% tile/2
% tile(Index, Tile)
%
% Fiecare soluție a predicatului tile este o corespondență între index
% (numărul piesei în lista din enunț) și reprezentarea internă a piesei
% respective.
%
% Puteți alege orice reprezentare doriți, în așa fel încât să puteți
% reprezenta toate piesele din enunț.
%
% Orice muchie a unei piese este o cetate, un drum, sau o pajiște.
% Pot exista cel mult 1 drum și cel mult 2 castele pe aceeași piesă.
%
% Reprezentarea trebuie să poată fi rotită (vezi predicatul ccw/3 mai
% jos) pentru a produce reprezentarea piesei rotite cu 90 de grade.
%
% Trebuie să definiți reprezentări pentru fiecare dintre cele 16 piese
% din enunțul temei.
%
% Exemplu: apelul tile(1, T1). trebuie să lege T1 la reprezentarea pe
% care o faceți pentru piesa 1. Această reprezentare poate fi transmisă
% celorlalte predicate din temă, pentru a întreba, de exemplu, ce se
% află pe muchia de nord a piesei 1, sau dacă piesa 1 se potrivește cu o
% altă piesă.

% Tile is a (center_type | list | number of citadels).
% Each position of the first list (N, E, S, W) holds the type of that area.
% The second list is (NE, SE, SW, NW)
tile(1, Tile) :- Tile = [c, [c, c, p, c], 1].
tile(2, Tile) :- Tile = [c, [c, c, d, c], 1].
tile(3, Tile) :- Tile = [p, [c, c, p, p], 1].
tile(4, Tile) :- Tile = [p, [c, c, p, p], 2].
tile(5, Tile) :- Tile = [p, [c, p, c, p], 2].
tile(6, Tile) :- Tile = [c, [c, p, c, p], 1].
tile(7, Tile) :- Tile = [p, [c, p, p, p], 1].
tile(8, Tile) :- Tile = [p, [c, c, d, d], 1].
tile(9, Tile) :- Tile = [p, [c, p, d, d], 1].
tile(10, Tile) :- Tile = [p, [c, d, d, p], 1].
tile(11, Tile) :- Tile = [d, [c, d, p, d], 1].
tile(12, Tile) :- Tile = [d, [c, d, d, d], 1].
tile(13, Tile) :- Tile = [d, [p, p, d, d], 1].
tile(14, Tile) :- Tile = [d, [p, d, p, d], 1].
tile(15, Tile) :- Tile = [d, [p, d, d, d], 1].
tile(16, Tile) :- Tile = [p, [d, d, d, d], 1].


% at/3
% at(+Tile, +Direction, ?What)
%
% Predicatul este adevărat dacă pe piesa Tile are pe muchia de pe
% direcția Direction o entitate de tipul What.
%
% Directions este una dintre n, e, s, w (vezi predicatul directions/1
% din utils.pl).
%
% Entitatea (What) este una dintre c, d, sau p. reprezentând cetate,
% drum, sau pajiște.
%
% De exemplu, piesa 4 are cetate în nord și în este, și pajiște în sud
% și vest. Iar piesa 10 are cetate în nord, drum în este și sud, și
% pajiște în vest.
%
% Dacă What nu este legat, trebuie legat la entitatea care se află pe
% muchia din direcția Dir.
at([_, [N, _, _, _], _], n, What) :- What = N.
at([_, [_, E, _, _], _], e, What) :- What = E.
at([_, [_, _, S, _], _], s, What) :- What = S.
at([_, [_, _, _, W], _], w, What) :- What = W.


% atL/3
% atL(+Tile, +Directions, +What)
%
% Predicatul este adevărat dacă piesa Tile are entitatea what pe toate
% direcțiile din lista Directions, cu aceleași valori pentru entități și
% direcții ca și la predicatul at/3.
%
% De exemplu, predicatul este adevărat pentru reprezentarea piesei 1,
% pentru lista [w,n,e], și pentru entitatea c. Este adevărat de asemenea
% pentru reprezentarea piesei 14, pentru lista [e,w], și pentru
% entitatea d.
%
% Atenție! Pentru ca predicatul să fie adevărat, nu este nevoie ca în
% Directions să fie *toate* direcțiile pe care se află entitatea
% respectivă, pot fi doar o submulțime a acestora.
% De exemplu, la piesa 14, predicatul este adevărat pentru entitatea d
% și pentru oricare dintre listele [w], [e] sau [e,w].

% E adevarat pt lista goala (in caz ca Rest de mai incolo e [])
atL(_, [], _) :- true.

% "Itereaza" prin Directions si verifica cu 'at'
atL(Tile, [X | Rest], What) :- at(Tile, X, What), atL(Tile, Rest, What).


% hasTwoCitadels/1
% hasTwoCitadels(+Tile)
%
% Predicatul întoarce adevărat dacă pe piesă există două cetăți diferite
% (ca în piesele 4 și 5).
hasTwoCitadels([_, _, 2]) :- true.


% ccw/3
% ccw(+Tile, +Rotation, -RotatedTile)
% Predicatul este adevărat dacă RotatedTile este reprezentarea piesei cu
% reprezentarea Tile, dar rotită de Rotation ori, în sens trigonometric.
%
% De exemplu, dacă T4 este reprezentarea piesei 4, atunci ccw(4, 1, R)
% va lega R la reprezentarea unei piese care are pajiște la nord și
% vest, și cetate la est și sud.
%
% Pentru piesele cu simetrie, reprezentarea unora dintre rotații este
% identică cu piesa inițială.
% De exemplu, la piesele 5, 6 și 14, rotirea cu Rotation=2 va duce la o
% reprezentare identică cu piesa inițială, și la fel rezultatele pentru
% Rotation=1 și Rotation=3 vor fi identice.
% La piesa 16, orice rotire trebuie să aibă aceeași reprezentare cu
% reprezentarea inițială.

% Helper function rotateOnce(+Tile, -RotatedTile)
rotateOnce([Center, [N, E, S, W], NumCitadels], RotatedTile) :- RotatedTile = [Center, [E, S, W, N], NumCitadels].

ccw(Tile, 0, RotatedTile) :- RotatedTile = Tile.
ccw(Tile, 1, RotatedTile) :- rotateOnce(Tile, RotatedTile).
ccw(Tile, 2, RotatedTile) :- ccw(Tile, 1, R1), rotateOnce(R1, RotatedTile).
ccw(Tile, 3, RotatedTile) :- ccw(Tile, 2, R2), rotateOnce(R2, RotatedTile).


% rotations/2
% rotations(+Tile, -RotationPairs)
%
% Predicatul leagă RotationPairs la o listă de perechi
% (Rotation, RotatedTile)
% în care Rotation este un număr de rotații între 0 și 3 inclusiv și
% RotatedTile este reprezentarea piesei Tile rotită cu numărul respectiv
% de rotații.
%
% Rezultatul trebuie întotdeauna să conțină perechea (0, Tile).
%
% IMPORTANT:
% Rezultatul nu trebuie să conțină rotații duplicate. De exemplu, pentru
% piesele 5,6 și 14 rezultatul va conține doar 2 perechi, iar pentru
% piesa 16 rezultatul va conține o singură pereche.
%
% Folosiți recursivitate (nu meta-predicate).

% addRotation(+Tile, +Acc, -RotationPairs, +NumRot)
addRotation(Tile, L, Res, NumRot) :- ccw(Tile, NumRot, R1), Tile \== R1, NextRot is NumRot + 1,
                        addRotation(R1, [(NumRot, R1) | L], Res, NextRot); Res = L. % Cand ajunge la RotationTile == Tile -> pune rezultatul in Res

rotations(Tile, RotationPairs) :- addRotation(Tile, [(0, Tile)], RotationPairs, 1).


% match/3
% match(+Tile, +NeighborTile, +NeighborDirection)
%
% Predicatul întoarce adevărat dacă NeighborTile poate fi pusă în
% direcția NeighborDirection față de Tile și se potrivește, adică muchia
% comună este de același fel.
%
% De exemplu, dacă T2 este reprezentarea piesei 2, iar T16 este
% reprezentarea piesei 16, atunci match(T2, T16, s) este adevărat.
%
% Similar, pentru piesele 8 și 10, este adevărat
% ccw(T8, 3, T8R), match(T8R, T10, w).
%
% Puteți folosi predicatul opposite/2 din utils.pl.

match(Tile, NeighborTile, NeighborDirection) :- opposite(NeighborDirection, TileDirection),
                                            at(Tile, NeighborDirection, What),
                                            at(NeighborTile, TileDirection, What).
                                            



% findRotation/3
% findRotation(+Tile, +Neighbors, -Rotation)
%
% Predicatul leagă Rotation la rotația (între 0 și 3 inclusiv) pentru
% piesa cu reprezentarea Tile, astfel încât piesa să se potrivească cu
% vecinii din Neighbors.
%
% Neighbors este o listă de perechi (NeighborTile, NeighborDirection) și
% specifică că pe direcția NeighborDirection se află piesa cu
% reprezentarea NeighborTile. Este posibil ca Neighbors să conțină mai
% puțin de 4 elemente.
%
% Se vor da toate soluțiile care duc la potrivire.
%
% De exemplu, pentru piesa 11, dacă la nord se află piesa 14 rotită o
% dată (drumul este vertical), iar la sud se află piesa 2 rotită de 2
% ori (drumul este spre nord), atunci posibilele rotații pentru piesa 11
% sunt 1 sau 3, deci findRotation trebuie să aibă 2 soluții, în care
% leagă R la 1, și la 3.
% În același exemplu, dacă am avea și piesa 1 ca vecin spre est, atunci
% soluția de mai sus s-ar reduce doar la rotația 3.
%
% Hint: Prolog face backtracking automat. Folosiți match/3.

% checkFit(+Tile, +Neighbors) e adevarata daca Tile se potriveste cu toti vecinii
checkFit(_, []) :- true.
checkFit(Tile, [(NeighborTile, NeighborDirection) | Rest]) :- 
                    match(Tile, NeighborTile, NeighborDirection),
                    checkFit(Tile, Rest).

findRotation(Tile, Neighbors, Rotation) :- checkFit(T1, Neighbors),
                                            ccw(Tile, Rotation, T1).





