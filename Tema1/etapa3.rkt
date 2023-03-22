#lang racket

(require "etapa2.rkt")

(provide (all-defined-out))

; TODO 1
; După modelul funcției stable-match?, implementați funcția
; get-unstable-couples care primește o listă de logodne
; engagements, o listă de preferințe masculine mpref și o 
; listă de preferințe feminine wpref, și întoarce lista
; tuturor cuplurilor instabile din engagements.
; Precizări (aspecte care se garantează, nu trebuie verificate):
; - fiecare cuplu din lista engagements are pe prima poziție
;   o femeie
; Nu este permisă recursivitatea pe stivă.
; Nu sunt permise alte funcții ajutătoare decât
; better-match-exists? și funcțiile de manipulare a listelor de
; preferințe definite în etapele anterioare.
; Nu este permisă apelarea multiplă a aceleiași funcții pe
; aceleași argumente.
; Folosiți una sau mai multe dintre expresiile let, let*, letrec,
; named let pentru a vă putea conforma acestor restricții.
(define (get-unstable-couples engagements mpref wpref)
  (let rev ((rev-eng '()) (eng-copy engagements)) ; inverseaza ordinea din fiecare logodna
    (if (null? eng-copy) ; dupa inversare se continua cautarea cuplurilor unstable cu un for pe coada
        (let for ((eng engagements) (result '()))
          (if (null? eng)
              result
              (let* ((x (car eng)) (woman (car x)) (man (cdr x)) (eng-rest (cdr eng)))
                (if (or (better-match-exists? man woman (get-pref-list mpref man) wpref engagements) ; se cauta un better match pt barbat
                        (better-match-exists? woman man (get-pref-list wpref woman) mpref rev-eng))  ; se cauta un better match pt femeie
                    (for eng-rest (cons x result))
                    (for eng-rest result)))))
        (rev (cons (cons (cdar eng-copy) (caar eng-copy)) rev-eng) (cdr eng-copy)))))

; TODO 2
; Implementați funcția engage care primește o listă free-men
; de bărbați încă nelogodiți, o listă de logodne parțiale 
; engagements (unde fiecare cuplu are pe prima poziție o femeie),
; o listă de preferințe masculine mpref și o listă de preferințe 
; feminine wpref, și întoarce o listă completă de logodne stabile,
; obținută conform algoritmului Gale-Shapley:
; - cât timp există un bărbat m încă nelogodit
;   - w = prima femeie din preferințele lui m pe care m nu a cerut-o încă
;   - dacă w este nelogodită, adaugă perechea (w, m) la engagements
;   - dacă w este logodită cu m'
;     - dacă w îl preferă pe m lui m'
;       - m' devine liber
;       - actualizează partenerul lui w la m în engagements
;     - altfel, repetă procesul cu următoarea femeie din lista lui m
; Folosiți named let pentru orice proces recursiv ajutător (deci nu
; veți defini funcții ajutătoare recursive).
; Folosiți let și/sau let* pentru a evita calcule duplicate.
(define (engage free-men engagements mpref wpref)
  (let iter-men ((men-list free-men) (engs engagements)) ; itereaza prin lista de barbati
    (if (null? men-list)
        engs
        (let ((man (car men-list)))
          (let iter-women-pref ((women-pref (get-pref-list mpref man))) ; itereaza prin lista de preferinte a fiecarui barbat
            (if (null? women-pref)
                (iter-men (cdr men-list) engs)
                (let ((woman (car women-pref)))
                  (let ((partner (get-partner engs woman))) ; cauta partenerul femeii
                    (if partner
                        (if (preferable? (get-pref-list wpref woman) man partner)
                            (iter-men (cons partner (cdr men-list)) (update-engagements engs woman man)) ; w il prefera pe m'
                            (iter-women-pref (cdr women-pref)))
                        (iter-men (cdr men-list) (append engs (list (cons woman man)))))))))))))


; TODO 3
; Implementați funcția gale-shapley care este un wrapper pentru
; algoritmul implementat de funcția engage. Funcția gale-shapley
; primește doar o listă de preferințe masculine mpref și o listă
; de preferințe feminine wpref și calculează o listă completă de
; logodne stabile conform acestor preferințe.
(define (gale-shapley mpref wpref)
  (engage (get-men mpref) '() mpref wpref))


; TODO 4
; Implementați funcția get-couple-members care primește o listă
; de perechi cu punct și întoarce o listă simplă cu toate elementele 
; care apar în perechi.
; Folosiți funcționale, fără recursivitate explicită.
(define (get-couple-members pair-list)
  (foldl (λ (x acc) (cons (car x) (cons (cdr x) acc))) null pair-list))
