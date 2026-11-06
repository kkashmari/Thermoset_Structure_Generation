# Reaction-template topology search

- Decision: **new template required**
- Reason: The complete force-field-typed pre/post fix bond/react topology has not been registered before.
- Typed template key: `642320ccba87c6d6d1ac1c322e34ab32ce8fd75bbc07c0f0f88d6051c1110ae4`

## Library topology changes

- remove 2:10–2:23
- remove 1:3–1:1
- add 2:10–1:1 (1)
- add 1:3–2:23 (1)

## Atom type and charge changes

- Atom 1: `c3h` → `c2`; charge +0.02700000 → +0.10480000
- Atom 2: `c3h` → `c1`; charge +0.08000000 → +0.08000000
- Atom 3: `o3e` → `oh`; charge -0.26600000 → -0.55710000
- Atom 13: `cp` → `cp`; charge -0.05000000 → +0.00000000
- Atom 17: `npc` → `n`; charge -0.73500000 → -0.65030000
- Atom 19: `h*` → `ho`; charge +0.39250000 → +0.42410000
- Atom 20: `h*` → `h*`; charge +0.39250000 → +0.43950000

## Bonds changes

- Remove [1, 2] type `c3h        c3h  via equivalence: c-c`
- Remove [1, 3] type `c3h        o3e  via equivalence: c-o`
- Remove [1, 6] type `c3h       hpan  via equivalence: c-h`
- Remove [1, 7] type `c3h       hpan  via equivalence: c-h`
- Remove [2, 3] type `c3h        o3e  via equivalence: c-o`
- Remove [2, 4] type `c2        c3h  via equivalence: c-c`
- Remove [2, 8] type `c3h       hpan  via equivalence: c-h`
- Remove [13, 17] type `cp        npc  via equivalence: cp-nh`
- Remove [17, 19] type `h*        npc  via equivalence: h*-nh`
- Remove [17, 20] type `h*        npc  via equivalence: h*-nh`
- Add [1, 2] type `c1         c2  via equivalence: c-c`
- Add [1, 6] type `c2       hpan  via equivalence: c-h`
- Add [1, 7] type `c2       hpan  via equivalence: c-h`
- Add [1, 17] type `c2          n  via equivalence: c-n`
- Add [2, 3] type `c1         oh  via equivalence: c-o`
- Add [2, 4] type `c1         c2  via equivalence: c-c`
- Add [2, 8] type `c1       hpan  via equivalence: c-h`
- Add [3, 19] type `ho         oh  via equivalence: h*-o`
- Add [13, 17] type `cp          n  via equivalence: cp_-n_`
- Add [17, 20] type `h*          n`

## Angles changes

- Remove [1, 2, 3] type `c3h        c3h        o3e  via equivalence: c-c-o`
- Remove [1, 2, 4] type `c2        c3h        c3h  via equivalence: c-c-c`
- Remove [1, 2, 8] type `c3h        c3h       hpan  via equivalence: c-c-h`
- Remove [1, 3, 2] type `c3h        o3e        c3h  via equivalence: c-o-c`
- Remove [2, 1, 3] type `c3h        c3h        o3e  via equivalence: c-c-o`
- Remove [2, 1, 6] type `c3h        c3h       hpan  via equivalence: c-c-h`
- Remove [2, 1, 7] type `c3h        c3h       hpan  via equivalence: c-c-h`
- Remove [2, 4, 5] type `c3h         c2         oc  via equivalence: c-c-o`
- Remove [2, 4, 9] type `c3h         c2       hpan  via equivalence: c-c-h`
- Remove [2, 4, 10] type `c3h         c2       hpan  via equivalence: c-c-h`
- Remove [3, 1, 6] type `hpan        c3h        o3e  via equivalence: h-c-o`
- Remove [3, 1, 7] type `hpan        c3h        o3e  via equivalence: h-c-o`
- Remove [3, 2, 4] type `c2        c3h        o3e  via equivalence: c-c-o`
- Remove [3, 2, 8] type `hpan        c3h        o3e  via equivalence: h-c-o`
- Remove [4, 2, 8] type `c2        c3h       hpan  via equivalence: c-c-h`
- Remove [6, 1, 7] type `hpan        c3h       hpan  via equivalence: h-c-h`
- Remove [12, 13, 17] type `cp         cp        npc  via equivalence: cp-cp-nh`
- Remove [13, 17, 19] type `cp        npc         h*  via equivalence: cp-nh-h*`
- Remove [13, 17, 20] type `cp        npc         h*  via equivalence: cp-nh-h*`
- Remove [14, 13, 17] type `cp         cp        npc  via equivalence: cp-cp-nh`
- Remove [19, 17, 20] type `h*        npc         h*  via wildcard: *3-np_-h_`
- Add [1, 2, 3] type `c2         c1         oh  via equivalence: c-c-o`
- Add [1, 2, 4] type `c2         c1         c2  via equivalence: c-c-c`
- Add [1, 2, 8] type `c2         c1       hpan  via equivalence: c-c-h`
- Add [1, 17, 13] type `c2          n         cp  via wildcard: *9-n_-c_`
- Add [1, 17, 20] type `c2          n         h*  via equivalence: c-n-h*`
- Add [2, 1, 6] type `c1         c2       hpan  via equivalence: c-c-h`
- Add [2, 1, 7] type `c1         c2       hpan  via equivalence: c-c-h`
- Add [2, 1, 17] type `c1         c2          n  via equivalence: c-c-n`
- Add [2, 3, 19] type `c1         oh         ho  via equivalence: c-o-h*`
- Add [2, 4, 5] type `c1         c2         oc  via equivalence: c-c-o`
- Add [2, 4, 9] type `c1         c2       hpan  via equivalence: c-c-h`
- Add [2, 4, 10] type `c1         c2       hpan  via equivalence: c-c-h`
- Add [3, 2, 4] type `c2         c1         oh  via equivalence: c-c-o`
- Add [3, 2, 8] type `hpan         c1         oh  via equivalence: h-c-o`
- Add [4, 2, 8] type `c2         c1       hpan  via equivalence: c-c-h`
- Add [6, 1, 7] type `hpan         c2       hpan  via equivalence: h-c-h`
- Add [6, 1, 17] type `hpan         c2          n  via equivalence: h-c-n`
- Add [7, 1, 17] type `hpan         c2          n  via equivalence: h-c-n`
- Add [12, 13, 17] type `cp         cp          n  via wildcard: *6-cp_-n_`
- Add [13, 17, 20] type `cp          n         h*  via wildcard: *3-n_-h_`
- Add [14, 13, 17] type `cp         cp          n  via wildcard: *6-cp_-n_`

## Dihedrals changes

- Remove [1, 2, 4, 5] type `c3h        c3h         c2         oc  via equivalence: c-c-c-o`
- Remove [1, 2, 4, 9] type `c3h        c3h         c2       hpan  via equivalence: c-c-c-h`
- Remove [1, 2, 4, 10] type `c3h        c3h         c2       hpan  via equivalence: c-c-c-h`
- Remove [1, 3, 2, 4] type `c2        c3h        o3e        c3h  via equivalence: c-c-o-c`
- Remove [1, 3, 2, 8] type `c3h        o3e        c3h       hpan  via equivalence: h-c-o-c`
- Remove [2, 3, 1, 6] type `c3h        o3e        c3h       hpan  via equivalence: h-c-o-c`
- Remove [2, 3, 1, 7] type `c3h        o3e        c3h       hpan  via equivalence: h-c-o-c`
- Remove [3, 1, 2, 4] type `c2        c3h        c3h        o3e  via equivalence: c-c-c-o`
- Remove [3, 1, 2, 8] type `hpan        c3h        c3h        o3e  via equivalence: h-c-c-o`
- Remove [3, 2, 1, 6] type `hpan        c3h        c3h        o3e  via equivalence: h-c-c-o`
- Remove [3, 2, 1, 7] type `hpan        c3h        c3h        o3e  via equivalence: h-c-c-o`
- Remove [3, 2, 4, 5] type `o3e        c3h         c2         oc  via equivalence: o-c-c-o`
- Remove [3, 2, 4, 9] type `hpan         c2        c3h        o3e  via equivalence: h-c-c-o`
- Remove [3, 2, 4, 10] type `hpan         c2        c3h        o3e  via equivalence: h-c-c-o`
- Remove [4, 2, 1, 6] type `c2        c3h        c3h       hpan  via equivalence: c-c-c-h`
- Remove [4, 2, 1, 7] type `c2        c3h        c3h       hpan  via equivalence: c-c-c-h`
- Remove [5, 4, 2, 8] type `hpan        c3h         c2         oc  via equivalence: h-c-c-o`
- Remove [6, 1, 2, 8] type `hpan        c3h        c3h       hpan  via equivalence: h-c-c-h`
- Remove [7, 1, 2, 8] type `hpan        c3h        c3h       hpan  via equivalence: h-c-c-h`
- Remove [8, 2, 4, 9] type `hpan         c2        c3h       hpan  via equivalence: h-c-c-h`
- Remove [8, 2, 4, 10] type `hpan         c2        c3h       hpan  via equivalence: h-c-c-h`
- Remove [11, 12, 13, 17] type `cp         cp         cp        npc  via equivalence: cp-cp-cp-nh`
- Remove [12, 13, 17, 19] type `cp         cp        npc         h*  via equivalence: cp-cp-nh-h*`
- Remove [12, 13, 17, 20] type `cp         cp        npc         h*  via equivalence: cp-cp-nh-h*`
- Remove [14, 13, 17, 19] type `cp         cp        npc         h*  via equivalence: cp-cp-nh-h*`
- Remove [14, 13, 17, 20] type `cp         cp        npc         h*  via equivalence: cp-cp-nh-h*`
- Remove [15, 14, 13, 17] type `c1         cp         cp        npc  via equivalence: c-cp-cp-nh`
- Remove [16, 14, 13, 17] type `c3         cp         cp        npc  via equivalence: c-cp-cp-nh`
- Remove [17, 13, 12, 18] type `c2         cp         cp        npc  via equivalence: c-cp-cp-nh`
- Add [1, 2, 3, 19] type `c2         c1         oh         ho  via equivalence: c-c-o-h*`
- Add [1, 2, 4, 5] type `c2         c1         c2         oc  via equivalence: c-c-c-o`
- Add [1, 2, 4, 9] type `c2         c1         c2       hpan  via equivalence: c-c-c-h`
- Add [1, 2, 4, 10] type `c2         c1         c2       hpan  via equivalence: c-c-c-h`
- Add [1, 17, 13, 12] type `c2          n         cp         cp  via wildcard: *-cp_-n_-*`
- Add [1, 17, 13, 14] type `c2          n         cp         cp  via wildcard: *-cp_-n_-*`
- Add [2, 1, 17, 13] type `c1         c2          n         cp  via wildcard: *-c_-n_-*`
- Add [2, 1, 17, 20] type `c1         c2          n         h*  via equivalence: c-c-n-h*`
- Add [3, 2, 1, 6] type `hpan         c2         c1         oh  via equivalence: h-c-c-o`
- Add [3, 2, 1, 7] type `hpan         c2         c1         oh  via equivalence: h-c-c-o`
- Add [3, 2, 1, 17] type `n         c2         c1         oh  via equivalence: n-c-c-o`
- Add [3, 2, 4, 5] type `oc         c2         c1         oh  via equivalence: o-c-c-o`
- Add [3, 2, 4, 9] type `hpan         c2         c1         oh  via equivalence: h-c-c-o`
- Add [3, 2, 4, 10] type `hpan         c2         c1         oh  via equivalence: h-c-c-o`
- Add [4, 2, 1, 6] type `c2         c1         c2       hpan  via equivalence: c-c-c-h`
- Add [4, 2, 1, 7] type `c2         c1         c2       hpan  via equivalence: c-c-c-h`
- Add [4, 2, 1, 17] type `c2         c1         c2          n  via equivalence: c-c-c-n`
- Add [4, 2, 3, 19] type `c2         c1         oh         ho  via equivalence: c-c-o-h*`
- Add [5, 4, 2, 8] type `hpan         c1         c2         oc  via equivalence: h-c-c-o`
- Add [6, 1, 2, 8] type `hpan         c1         c2       hpan  via equivalence: h-c-c-h`
- Add [6, 1, 17, 13] type `cp          n         c2       hpan  via wildcard: *-c_-n_-*`
- Add [6, 1, 17, 20] type `h*          n         c2       hpan  via equivalence: h-c-n-h*`
- Add [7, 1, 2, 8] type `hpan         c1         c2       hpan  via equivalence: h-c-c-h`
- Add [7, 1, 17, 13] type `cp          n         c2       hpan  via wildcard: *-c_-n_-*`
- Add [7, 1, 17, 20] type `h*          n         c2       hpan  via equivalence: h-c-n-h*`
- Add [8, 2, 1, 17] type `hpan         c1         c2          n  via equivalence: h-c-c-n`
- Add [8, 2, 3, 19] type `ho         oh         c1       hpan  via equivalence: h-c-o-h*`
- Add [8, 2, 4, 9] type `hpan         c1         c2       hpan  via equivalence: h-c-c-h`
- Add [8, 2, 4, 10] type `hpan         c1         c2       hpan  via equivalence: h-c-c-h`
- Add [11, 12, 13, 17] type `cp         cp         cp          n  via wildcard: *-cp_-cp_-*`
- Add [12, 13, 17, 20] type `cp         cp          n         h*  via wildcard: *-cp_-n_-*`
- Add [14, 13, 17, 20] type `cp         cp          n         h*  via wildcard: *-cp_-n_-*`
- Add [15, 14, 13, 17] type `c1         cp         cp          n  via wildcard: *-cp_-cp_-*`
- Add [16, 14, 13, 17] type `c3         cp         cp          n  via wildcard: *-cp_-cp_-*`
- Add [17, 13, 12, 18] type `c2         cp         cp          n  via wildcard: *-cp_-cp_-*`

## Impropers changes

- Remove [1, 2, 3, 4] type `c2        c3h        c3h        o3e  WARNING: zero fallback`
- Remove [1, 2, 3, 6] type `c3h        c3h       hpan        o3e  WARNING: zero fallback`
- Remove [1, 2, 3, 7] type `c3h        c3h       hpan        o3e  WARNING: zero fallback`
- Remove [1, 2, 3, 8] type `c3h        c3h       hpan        o3e  WARNING: zero fallback`
- Remove [1, 2, 4, 8] type `c2        c3h        c3h       hpan  WARNING: zero fallback`
- Remove [1, 2, 6, 7] type `c3h        c3h       hpan       hpan  WARNING: zero fallback`
- Remove [1, 3, 6, 7] type `hpan        c3h       hpan        o3e  WARNING: zero fallback`
- Remove [2, 3, 4, 8] type `c2        c3h       hpan        o3e  WARNING: zero fallback`
- Remove [2, 4, 5, 9] type `c3h         c2       hpan         oc  WARNING: zero fallback`
- Remove [2, 4, 5, 10] type `c3h         c2       hpan         oc  WARNING: zero fallback`
- Remove [2, 4, 9, 10] type `c3h         c2       hpan       hpan  WARNING: zero fallback`
- Remove [12, 13, 14, 17] type `cp         cp         cp        npc  via equivalence: cp-cp-cp-nh`
- Remove [13, 17, 19, 20] type `cp        npc         h*         h*  via wildcard: *-np_-*-*`
- Add [1, 2, 3, 4] type `c2         c1         c2         oh  WARNING: zero fallback`
- Add [1, 2, 3, 8] type `c2         c1       hpan         oh  WARNING: zero fallback`
- Add [1, 2, 4, 8] type `c2         c1         c2       hpan  WARNING: zero fallback`
- Add [1, 2, 6, 7] type `c1         c2       hpan       hpan  WARNING: zero fallback`
- Add [1, 2, 6, 17] type `c1         c2       hpan          n  WARNING: zero fallback`
- Add [1, 2, 7, 17] type `c1         c2       hpan          n  WARNING: zero fallback`
- Add [1, 6, 7, 17] type `hpan         c2       hpan          n  WARNING: zero fallback`
- Add [1, 13, 17, 20] type `c2          n         cp         h*  via wildcard: *-n_-*-*`
- Add [2, 3, 4, 8] type `c2         c1       hpan         oh  WARNING: zero fallback`
- Add [2, 4, 5, 9] type `c1         c2       hpan         oc  WARNING: zero fallback`
- Add [2, 4, 5, 10] type `c1         c2       hpan         oc  WARNING: zero fallback`
- Add [2, 4, 9, 10] type `c1         c2       hpan       hpan  WARNING: zero fallback`
- Add [12, 13, 14, 17] type `cp         cp         cp          n  via wildcard: *-cp_-*-*`
