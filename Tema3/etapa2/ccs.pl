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


%%%%%%%%%%%%%%%%%%%%%%%%%% Etapa 2


%% TODO
% emptyBoard/1
% emptyBoard(-Board)
%
% Leagă Board la reprezentarea unei table goale de joc (nu a fost
% plasată încă nicio piesă).
emptyBoard([]) :- true.



%% TODO
% boardSet/4
% boardSet(+BoardIn, +Pos, +Tile, -BoardOut)
%
% Predicatul întoarce false dacă se încearcă plasarea unei piese pe o
% poziție pe care este deja o piesă, pe o poziție fără muchie comună
% cu o piesă existentă, sau într-un loc unde piesa nu se potrivește cu
% vecinii săi.
%
% Pentru o tablă goală, predicatul reușește întotdeauna, și poziția Pos
% devine singura de pe tablă.
%
% Poziția este dată ca un tuplu (X, Y).

% isIsolated(+BoardIn, +(X, Y)) - intoarce adevarat daca piesa
% n ar avea muchii comune cu nicio alta piesa de pe tabla
isIsolated(BoardIn, (X, Y)) :- X1 is (X - 1), X2 is (X + 1), Y1 is (Y - 1), Y2 is (Y + 1),
    \+ member(X1/Y/_, BoardIn), \+ member(X/Y1/_, BoardIn),
    \+ member(X2/Y/_, BoardIn), \+ member(X/Y2/_, BoardIn).

boardSet([], (X, Y), Tile, BoardOut) :- BoardOut = [X/Y/Tile].
boardSet(BoardIn, (X, Y), Tile, BoardOut) :- 
                    \+ isIsolated(BoardIn, (X, Y)),
                    canPlaceTile(BoardIn, (X, Y), Tile),
                    BoardOut = [X/Y/Tile | BoardIn].


%% TODO
% boardGet/3
% boardGet(+Board, +Pos, -Tile)
%
% Predicatul leagă Tile la reprezentarea piesei de pe tabla Board, de la
% poziția Pos. Poziția este dată ca un tuplu (X, Y).
%
% Dacă la poziția Pos nu este nicio piesă, predicatul eșuează.

boardGet(Board, (X, Y), Tile) :- member(X/Y/Tile, Board).



%% TODO
% boardGetLimits/5
% boardGetLimits(+Board, -XMin, -Ymin, -XMax, -YMax)
%
% Predicatul leagă cele 4 argumente la coordonatele x-y extreme la
% care există piese pe tablă.
%
% Pentru o tablă goală, predicatul eșuează.
%
% Hint: max_list/2 și min_list/2
boardGetLimits([], _, _, _, _) :- false.
boardGetLimits(Board, Xmin, Ymin, Xmax, Ymax) :- findall(X, member(X/_/_, Board), XList),
                findall(Y, member(_/Y/_, Board), YList), max_list(XList, Xmax), max_list(YList, Ymax),
                min_list(XList, Xmin), min_list(YList, Ymin).

%% TODO
% canPlaceTile/3
% canPlaceTile(+Board, +Pos, +Tile)
%
% Întoarce adevărat dacă este o mișcare validă plasarea piese Tile la
% poziția Pos pe tabla Board. Poziția este dată ca un tuplu (X, Y).
%
% O mișcare este validă dacă tabla este goală sau dacă:
% - poziția este liberă;
% - poziția este adiacentă (are o muchie comună) cu o piesă deja
% existentă pe tablă;
% - piesa se potrivește cu toți vecinii deja existenți pe tablă.
%
% Hint: neighbor/3 și directions/1 , ambele din utils.pl

fitsOnDir(Board, (X, Y), Tile, Dir) :- neighbor((X, Y), Dir, (XN, YN)), \+ member(XN/YN/_, Board);
            neighbor((X, Y), Dir, (XN, YN)), member(XN/YN/NeighborTile, Board), match(Tile, NeighborTile, Dir).

fits(Board, (X, Y), Tile) :- fitsOnDir(Board, (X, Y), Tile, n), fitsOnDir(Board, (X, Y), Tile, e),
        fitsOnDir(Board, (X, Y), Tile, s), fitsOnDir(Board, (X, Y), Tile, w).

canPlaceTile([], _, _) :- true.
canPlaceTile(Board, (X, Y), Tile) :- \+ member(X/Y/_, Board), % pozitia e libera pe tabla
                                    \+ isIsolated(Board, (X, Y)), % piesa ar avea macar o muchie adiacenta altei piese
                                    fits(Board, (X, Y), Tile). % verifica daca piesa se potriveste cu piesele din jur



%% TODO
% getAvailablePositions/2
% getAvailablePositions(+Board, -Positions)
%
% Predicatul leagă Positions la o listă de perechi (X, Y)
% a tuturor pozițiilor de pe tabla Board unde se pot pune piese (poziții
% libere vecine pe o muchie cu piese existente pe tablă).
%
% Pentru o tablă goală, predicatul eșuează.
%
% Hint: between/3 (predefinit) și neighbor/3 din utils.pl
%
% Atenție! Și în afara limitelor curente există poziții disponibile.

% Verifica daca pozitia (X, Y) are cel putin un vecin
hasNeigh(Board, (X, Y)) :-  neighbor((X, Y), n, (XN, YN)), member(XN/YN/_, Board);
                            neighbor((X, Y), e, (XN, YN)), member(XN/YN/_, Board);
                            neighbor((X, Y), s, (XN, YN)), member(XN/YN/_, Board);
                            neighbor((X, Y), w, (XN, YN)), member(XN/YN/_, Board).


% Verifica daca pozitia e intre limite, e libera si are macar un vecin
availablePosition(Board, (X, Y)) :- boardGetLimits(Board, Xmin, Ymin, Xmax, Ymax),
                                    X1 is Xmin - 1, X2 is Xmax + 1, Y1 is Ymin - 1, Y2 is Ymax + 1, 
                                    between(X1, X2, X), between(Y1, Y2, Y),
                                    \+ member(X/Y/_, Board), hasNeigh(Board, (X, Y)).


getAvailablePositions(Board, Positions) :- \+ emptyBoard(Board),
                    findall((X, Y), availablePosition(Board, (X, Y)), PosWithDuplicates),
                    sort(PosWithDuplicates, Positions). % scoate duplicatele



%% TODO
% findPositionForTile/4
% findPositionForTile(+Board, +Tile, -Position, -Rotation)
%
% Predicatul are ca soluții toate posibilele plasări pe tabla Board ale
% piesei Tile, legând Position la o pereche (X, Y) care reprezintă
% poziția și Rotation la un număr între 0 și 3 inclusiv, reprezentând de
% câte ori trebuie rotită piesa ca să se potrivească.
%
% Unele piese se pot potrivi cu mai multe rotații pe aceeași poziție și
% acestea reprezintă soluții diferite ale predicatului, dar numai dacă
% rotațiile duc la rezultate diferite.
%
% Dacă tabla este goală, predicatul leagă Position la (0, 0) și Rotation
% la 0.
%
% De exemplu, dacă pe tablă se află doar piesa 11, la vest de ea piesa 9
% se potrivește cu rotația 1 sau 2 - două soluții diferite. Pentru
% plasarea la vest de piesa 11 a piesei 16 însă există o singură soluție
% - rotație 0.
%
% În ieșirea de la teste, rezultatele vor fi asamblate ca
% (X,Y):Rotation.

% getNeighbor(+Board, -Position, -Tile, - Dir)/4
% Genereaza lista de vecini de tip (NeighborTile, NeighborDir)
getNeighbor(Board, (X, Y), (Tile, Dir)) :- neighbor((X, Y), Dir, (XN, YN)),
                                           member(XN/YN/Tile, Board).

% Genereaza rezultatele
validMove(Board, Tile, Pos, Rot, Positions) :- member(Pos, Positions),
                                    findall((NeighTile, NeighDir), getNeighbor(Board, Pos, (NeighTile, NeighDir)), Neighbors),
                                    rotations(Tile, RotationPairs), % ca sa scot rotatiile duplicate
                                    findRotation(Tile, Neighbors, Rot),
                                    member((Rot, _), RotationPairs). % scot rotatiile duplicate

findPositionForTile([], _, Position, Rotation) :- Position = (0, 0), Rotation = 0.
findPositionForTile(Board, Tile, Position, Rotation) :-
                getAvailablePositions(Board, Positions),
                validMove(Board, Tile, Position, Rotation, Positions).






