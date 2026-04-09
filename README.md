# Testowanie generatorów liczb psuedolosowych

## Opis projektu

Projekt zawiera implementację i porównanie kilku generatorów liczb pseudolosowych z wykorzystaniem testów statystycznych pierwszego i drugiego poziomu. Łączy on teorię testowania generatorów (z pliku PDF) z praktyczną implementacją testów i wizualizacji statystycznych (skrypt Python).

Głównym celem projektu jest:

- wygenerowanie dużych prób liczb pseudolosowych,
- przeprowadzenie testu chi-kwadrat i testów zgodności,
- wizualizacja rozkładów i p-value,
- ocena jakości generatorów,
- porównanie wyników na poziomie testów 1. i 2. rzędu.
  
## Zawartość repozytorium
### 1. Testowanie_generatorów.pdf

Dokument omawia:
- podstawy teoretyczne generatorów liczb pseudolosowych (PRNG),
- wymagania na dobry generator,
- strukturę testów statystycznych:
- testy pierwszego stopnia – weryfikują zgodność rozkładu z oczekiwanym,
- testy drugiego stopnia – testują rozkład wartości p-value,
- przykłady testów takich jak:
  - test chi-kwadrat,
  - test Kołmogorowa-Smirnowa,
  - testy grupowe,
- interpretację p-value i typowe błędy w ich analizie.

Treść dokumentu stanowi kontekst teoretyczny do kodu, wyjaśniając dlaczego testy są wykonywane w opisany sposób.

### 2. Monte_Carlo_1_code.py

Skrypt implementuje i porównuje kilka generatorów:
- Xorshift128
- LCG (Linear Congruential Generator)
- random.random() (generator Pythona)
oraz generuje próbki liczb i wykonuje testy:
- testy chi-kwadrat,
- porównanie histogramów,
- budowa testów drugiego poziomu (testowanie rozkładu p-value),

Tworzy również wykresy:
- histogramy wartości generatorów,
- histogramy p-value,
- wykresy rozkładów Poissona i normalnych,
- wyświetla podsumowanie wyników testów.

Najważniejszą funkcją jest:

compare_generators_plots()

której wywołanie uruchamia wszystkie testy i generuje komplet wykresów.

## Testy statystyczne
Test pierwszego poziomu (test chi-kwadrat)

Każdy generator produkuje próbkę N liczb, która:

jest binowana,
porównywana z rozkładem oczekiwanym (jednostajnym),
przeliczana na statystykę χ² i p-value.
Test drugiego poziomu

Analizowany jest rozkład p-value z wielu testów pierwszego poziomu.
Jeśli p-value są równomiernie rozłożone, generator uznaje się za statystycznie poprawny.

 ## Wyniki

Skrypt wypisuje zbiorcze wyniki:

=== PODSUMOWANIE (p-value 2nd level) ===
Xorshift : p(2nd) = ...
LCG      : p(2nd) = ...
Python   : p(2nd) = ...

Interpretacja:

wysokie p-value (2nd level) → generator zgodny z rozkładem jednostajnym,
niskie p-value → generator ma defekty strukturalne.

## Uruchamianie
Wymagane biblioteki:
pip install numpy scipy matplotlib pandas
Uruchomienie testów:
python Monte_Carlo_1_code.py

Wygenerowane zostaną wszystkie wykresy i podsumowanie wyników.
