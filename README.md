# TAP - Projekt 1: Reaktor CSTR

## Jak wrzucic do Overleaf

1. Wrzucic `sprawozdanie1.tex` oraz `tabele_1c.tex` do Overleaf
2. Wrzucic caly folder `wykresy/` do Overleaf (tylko pliki PDF uzywane w sprawozdaniu, lista nizej)
3. Skompilowac - powinno dzialac od razu

### Wykresy uzywane w sprawozdaniu

```
wykresy/model_punkt_pracy.pdf
wykresy/model_odpowiedzi_skokowe.pdf
wykresy/1a_przestrzen_stanow.pdf
wykresy/1a_transmitancje.pdf
wykresy/1b_CAin.pdf
wykresy/1b_FC.pdf
wykresy/1b_Tin.pdf
wykresy/1b_TCin.pdf
wykresy/1d_jakosc_sterowania.pdf
wykresy/1d_jakosc_zaklocenia.pdf
wykresy/1d_model_dyskretny.pdf
wykresy/1d_weryfikacja_CAin.pdf
wykresy/1d_weryfikacja_FC.pdf
wykresy/1d_ss_vs_tf.pdf
```

## Struktura plikow

### Funkcja wspolna

| Plik | Opis |
|------|------|
| `punkt_pracy.m` | Oblicza dokladny punkt rownowagi (fsolve) i zwraca `[x0, u0, p]`. Uzywana przez wszystkie skrypty. |

### Skrypty MATLAB (kolejnosc uruchamiania)

| Plik | Zadanie | Co generuje |
|------|---------|-------------|
| `wyznacz_punkt_pracy.m` | Weryfikacja punktu pracy | Wypisuje dokladny punkt rownowagi w konsoli |
| `model.m` | Model nieliniowy | `model_punkt_pracy.pdf`, `model_odpowiedzi_skokowe.pdf` |
| `model_przestrzen_stanu.m` | Zad 1a: linearyzacja, transmitancje, dyskretyzacja | `1a_przestrzen_stanow.pdf`, `1a_transmitancje.pdf`, `1d_dyskretyzacja_*.pdf`, `1d_model_dyskretny.pdf` |
| `zad1b_porownanie.m` | Zad 1b: porownanie NL vs LIN (fan-ploty) | `1b_CAin.pdf`, `1b_FC.pdf`, `1b_Tin.pdf`, `1b_TCin.pdf` |
| `zad1c_jakosc_aproksymacji.m` | Zad 1c: RMSE aproksymacji liniowej | `tabele_1c.tex`, `wyniki_1c.csv` |
| `zad1d_dyskretyzacja.m` | Zad 1d: dyskretyzacja, dobor Tp, weryfikacja, SS vs TF | `1d_jakosc_sterowania.pdf`, `1d_jakosc_zaklocenia.pdf`, `1d_weryfikacja_CAin.pdf`, `1d_weryfikacja_FC.pdf`, `1d_ss_vs_tf.pdf` |

### Sprawozdanie

| Plik | Opis |
|------|------|
| `sprawozdanie1.tex` | Glowny plik LaTeX sprawozdania |
| `tabele_1c.tex` | Tabelki RMSE generowane przez `zad1c_jakosc_aproksymacji.m` (includowane w sprawozdaniu) |

### Inne pliki

| Plik | Opis |
|------|------|
| `model_zlinearyzowany.m` | Starszy skrypt z modelem zlinearyzowanym (nie jest uzywany w sprawozdaniu) |
| `TAP_proj_wprs.pdf` | Tresc projektu |
| `TAPpolecenie26L.pdf` | Polecenie - wariant 26L |
| `TAPzad1.pdf` | Tresc zadania 1 |

## Jak wygenerowac wykresy od nowa

Otworzyc MATLAB, ustawic working directory na folder projektu i uruchomic po kolei:

```matlab
model                       % wykresy modelu nieliniowego
model_przestrzen_stanu      % linearyzacja + dyskretyzacja (1a)
zad1b_porownanie            % fan-ploty NL vs LIN (1b)
zad1c_jakosc_aproksymacji   % tabelki RMSE (1c)
zad1d_dyskretyzacja         % dyskretyzacja + weryfikacja (1d)
```

Wykresy zapisuja sie automatycznie do folderu `wykresy/`.
