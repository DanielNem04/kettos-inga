# Kettős inga szimuláció – Technikai dokumentáció

A projekt egy **kettősinga inga** mozgásának szimulációját valósítja meg **4. rendű Runge–Kutta (RK4)** módszerrel. A cél a fizikai rendszer numerikus modellezése, a szögek időbeli alakulásának számolása és mentése. Illetve a rendszer vizsgálatára az első átfordulási idő meghatározását is meg lehet határozni.

A program háromféle futási módban használható:
- Teljes szimuláció és annak kimentése CSV fájlba (`csv`)
- Egy adott paraméterezésű kezdőállapot esetén az egyik szög első átfordulási idejének kiszámítása (`flip_single`)
- Kezdeti szögrácson történő átfordulási idők számítása, hőtérképes készítéséhez, ez is egy CSV fájlt generál  (`flip_grid`)

---

## Fájlszerkezet

```
.
│   CMakeLists.txt        # CMake konfiguráció
│   readme.md             # Ez a fájl
│   jegyzokonyv.pdf       # kiertekeles abrakkal
│   abrazolas.ipynb        # Jupyter Notebook az abra keszítéselhez
│
├───app/
│   └───main.cpp          # A fő program
│
├───figures/              # (A jegyzőkönyv ábrái, és az ahhoz használt adatfileok)
│
├───include/
│   └───double_pendulum/
│       ├───analysis.h    # Az első átfordulási idő meghatározása
│       ├───core.h        # Alapvető fizikai számítások (pl. deriváltak)
│       └───simulation.h  # Szimulációs függvények deklarációja
│
├───references/
│   ├───chen_2008_report.pdf  # Hivatkozott szakirodalom
│   └───Double.pdf            
│
└───src/
    ├───analysis.cpp     # Az átfordulási idő kiszámítása
    ├───core.cpp         # Fizikai deriváltak és koordináta-konvertálás
    └───simulation.cpp   # A teljes RK4-alapú szimuláció megvalósítása
```

---

## Programhasználat

A program `main.cpp` fájlja három különböző módon hívható:

```
./my_program <mode> [paraméterek...]
```

### `csv` mód

```
./my_program csv l1 l2 m1 m2 phi1_0 phi2_0 omega1_0 omega2_0 T h
```

Teljes mozgás szimulációja és kiírása `simulation_output.csv` fájlba.

**Kimeneti adatok:**
- idő
- pozíciók: x1, y1, x2, y2
- sebességek: v1, v2
- szögek: phi1, phi2
- szögsebességek: omega1, omega2

---

### `flip_single` mód

```
./my_program flip_single l1 l2 m1 m2 phi1_0 phi2_0 omega1_0 omega2_0 T h use_phi1
```

Egyetlen szög első átfordulási idejét számítja ki (amikor eléri az első (2k+1)π értéket).

- `use_phi1` = 1 → phi1 szög vizsgálata
- `use_phi1` = 0 → phi2 szög vizsgálata

---

### `flip_grid` mód

```
./my_program flip_grid l1 l2 m1 m2 omega1_0 omega2_0 T h min_phi max_phi resolution
```

Rácsszerűen vizsgálja a kezdőszögek hatását az első átfordulási időre, kétdimenziós CSV-fájlokat generál.

**Kimenetek:**
- `phi1_flip_...csv` és `phi2_flip_...csv` a `./app` mappába

---

##  áttekintése

### `core` (fizikai modell)

- `converter()` – szögekből és szögsebességekből pozíciókat és sebességeket számol
- `do1dt()`, `do2dt()` – a nemlineáris mozgásegyenletek jobb oldalát adják vissza

---

### `simulation`

- `run_simulation_full()` – RK4 lépésekkel végigszimulálja a teljes rendszert a megadott időhatárig
- `run_simulation_phi()` – csak a kiválasztott szög (`phi1` vagy `phi2`) időbeli változását számítja

---

### `analysis`

- `find_first_flip_time()` – meghatározza azt az időpontot, amikor a kiválasztott szög először eléri az (2k+1)π értéket
- Algoritmus lineáris interpolációval pontosítja az átfordulási időt

---

## Megjegyzések

A program CMake rendszerrel fordítható:

```bash
mkdir build
cd build
cmake ..
cmake --build .
```

- A szimulációk minden paramétert parancssorból várnak
- A `flip_grid` mód számításigényes, de nagy felbontású hőtérképek készítésére alkalmas
- A `simulation_output.csv` és a `flip_*.csv` fájlok Pythonból (pl. matplotlib) könnyen ábrázolhatók

---

## Függőségek

- C++17 ajánlott


## Kimenetek összefoglalása

| Mód           | Futtatási parancs  (linux)                                                                | Kimenet típusa            | Kimenet neve / helye                     | Tartalom                                                                 |
|----------------|------------------------------------------------------------------------------------------|----------------------------|------------------------------------------|--------------------------------------------------------------------------|
| `csv`          | `./my_program csv l1 l2 m1 m2 phi1_0 phi2_0 omega1_0 omega2_0 T h`                       | CSV fájl                   | `app/simulation_output.csv`             | idő, pozíciók (x1,y1,x2,y2), sebességek, szögek, szögsebességek          |
| `flip_single`  | `./my_program flip_single l1 l2 m1 m2 phi1_0 phi2_0 omega1_0 omega2_0 T h use_phi1`     | Konzolos szöveg (stdout)   | –                                        | első átfordulási idő (float, pl. `12.3652`)                             |
| `flip_grid`    | `./my_program flip_grid l1 l2 m1 m2 omega1_0 omega2_0 T h min_phi max_phi resolution`   | Két CSV fájl               | `app/phi1_flip_*.csv`, `phi2_flip_*.csv` | rácsos adatok: (phi1, phi2, flip_time), hőtérképekhez         |


---

## Kapcsolódó szakirodalom

- `references/chen_2008_report.pdf`
- `references/Double.pdf`
