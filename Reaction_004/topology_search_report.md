# Reaction-template topology search

- Decision: **new template required**
- Reason: The complete force-field-typed pre/post fix bond/react topology has not been registered before.
- Typed template key: `d08a2e44080466a01a85fea5c058b5cb34b587da82261a963062b88e9047928a`

## Library topology changes

- remove 1:53–1:67
- remove 2:28–2:37
- add 1:53–2:37 (1)
- add 2:28–1:67 (1)

## Atom type and charge changes

- Atom 14: `n` → `n`; charge -0.65030000 → -0.42160000
- Atom 16: `h*` → `ho`; charge +0.43950000 → +0.42410000
- Atom 19: `o3e` → `oh`; charge -0.26600000 → -0.55710000
- Atom 21: `c3h` → `c1`; charge +0.08000000 → +0.08000000
- Atom 24: `c3h` → `c2`; charge +0.02700000 → +0.10480000

## Bonds changes

- Remove [14, 16] type `h*          n`
- Remove [18, 21] type `c2        c3h  via equivalence: c-c`
- Remove [19, 21] type `c3h        o3e  via equivalence: c-o`
- Remove [19, 24] type `c3h        o3e  via equivalence: c-o`
- Remove [20, 24] type `c3h       hpan  via equivalence: c-h`
- Remove [21, 24] type `c3h        c3h  via equivalence: c-c`
- Remove [21, 26] type `c3h       hpan  via equivalence: c-h`
- Remove [24, 25] type `c3h       hpan  via equivalence: c-h`
- Add [14, 24] type `c2          n  via equivalence: c-n`
- Add [16, 19] type `ho         oh  via equivalence: h*-o`
- Add [18, 21] type `c1         c2  via equivalence: c-c`
- Add [19, 21] type `c1         oh  via equivalence: c-o`
- Add [20, 24] type `c2       hpan  via equivalence: c-h`
- Add [21, 24] type `c1         c2  via equivalence: c-c`
- Add [21, 26] type `c1       hpan  via equivalence: c-h`
- Add [24, 25] type `c2       hpan  via equivalence: c-h`

## Angles changes

- Remove [1, 14, 16] type `c2          n         h*  via equivalence: c-n-h*`
- Remove [10, 14, 16] type `cp          n         h*  via wildcard: *3-n_-h_`
- Remove [17, 18, 21] type `c3h         c2         oc  via equivalence: c-c-o`
- Remove [18, 21, 19] type `c2        c3h        o3e  via equivalence: c-c-o`
- Remove [18, 21, 24] type `c2        c3h        c3h  via equivalence: c-c-c`
- Remove [18, 21, 26] type `c2        c3h       hpan  via equivalence: c-c-h`
- Remove [19, 21, 24] type `c3h        c3h        o3e  via equivalence: c-c-o`
- Remove [19, 21, 26] type `hpan        c3h        o3e  via equivalence: h-c-o`
- Remove [19, 24, 20] type `hpan        c3h        o3e  via equivalence: h-c-o`
- Remove [19, 24, 21] type `c3h        c3h        o3e  via equivalence: c-c-o`
- Remove [19, 24, 25] type `hpan        c3h        o3e  via equivalence: h-c-o`
- Remove [20, 24, 21] type `c3h        c3h       hpan  via equivalence: c-c-h`
- Remove [20, 24, 25] type `hpan        c3h       hpan  via equivalence: h-c-h`
- Remove [21, 18, 22] type `c3h         c2       hpan  via equivalence: c-c-h`
- Remove [21, 18, 23] type `c3h         c2       hpan  via equivalence: c-c-h`
- Remove [21, 19, 24] type `c3h        o3e        c3h  via equivalence: c-o-c`
- Remove [21, 24, 25] type `c3h        c3h       hpan  via equivalence: c-c-h`
- Remove [24, 21, 26] type `c3h        c3h       hpan  via equivalence: c-c-h`
- Add [1, 14, 24] type `c2          n         c2  via equivalence: c-n-c`
- Add [10, 14, 24] type `c2          n         cp  via wildcard: *9-n_-c_`
- Add [14, 24, 20] type `hpan         c2          n  via equivalence: h-c-n`
- Add [14, 24, 21] type `c1         c2          n  via equivalence: c-c-n`
- Add [14, 24, 25] type `hpan         c2          n  via equivalence: h-c-n`
- Add [16, 19, 21] type `c1         oh         ho  via equivalence: c-o-h*`
- Add [17, 18, 21] type `c1         c2         oc  via equivalence: c-c-o`
- Add [18, 21, 19] type `c2         c1         oh  via equivalence: c-c-o`
- Add [18, 21, 24] type `c2         c1         c2  via equivalence: c-c-c`
- Add [18, 21, 26] type `c2         c1       hpan  via equivalence: c-c-h`
- Add [19, 21, 24] type `c2         c1         oh  via equivalence: c-c-o`
- Add [19, 21, 26] type `hpan         c1         oh  via equivalence: h-c-o`
- Add [20, 24, 21] type `c1         c2       hpan  via equivalence: c-c-h`
- Add [20, 24, 25] type `hpan         c2       hpan  via equivalence: h-c-h`
- Add [21, 18, 22] type `c1         c2       hpan  via equivalence: c-c-h`
- Add [21, 18, 23] type `c1         c2       hpan  via equivalence: c-c-h`
- Add [21, 24, 25] type `c1         c2       hpan  via equivalence: c-c-h`
- Add [24, 21, 26] type `c2         c1       hpan  via equivalence: c-c-h`

## Dihedrals changes

- Remove [2, 1, 14, 16] type `c1         c2          n         h*  via equivalence: c-c-n-h*`
- Remove [5, 1, 14, 16] type `h*          n         c2       hpan  via equivalence: h-c-n-h*`
- Remove [6, 1, 14, 16] type `h*          n         c2       hpan  via equivalence: h-c-n-h*`
- Remove [9, 10, 14, 16] type `cp         cp          n         h*  via wildcard: *-cp_-n_-*`
- Remove [11, 10, 14, 16] type `cp         cp          n         h*  via wildcard: *-cp_-n_-*`
- Remove [17, 18, 21, 19] type `o3e        c3h         c2         oc  via equivalence: o-c-c-o`
- Remove [17, 18, 21, 24] type `c3h        c3h         c2         oc  via equivalence: c-c-c-o`
- Remove [17, 18, 21, 26] type `hpan        c3h         c2         oc  via equivalence: h-c-c-o`
- Remove [18, 21, 19, 24] type `c2        c3h        o3e        c3h  via equivalence: c-c-o-c`
- Remove [18, 21, 24, 19] type `c2        c3h        c3h        o3e  via equivalence: c-c-c-o`
- Remove [18, 21, 24, 20] type `c2        c3h        c3h       hpan  via equivalence: c-c-c-h`
- Remove [18, 21, 24, 25] type `c2        c3h        c3h       hpan  via equivalence: c-c-c-h`
- Remove [19, 21, 18, 22] type `hpan         c2        c3h        o3e  via equivalence: h-c-c-o`
- Remove [19, 21, 18, 23] type `hpan         c2        c3h        o3e  via equivalence: h-c-c-o`
- Remove [19, 21, 24, 20] type `hpan        c3h        c3h        o3e  via equivalence: h-c-c-o`
- Remove [19, 21, 24, 25] type `hpan        c3h        c3h        o3e  via equivalence: h-c-c-o`
- Remove [19, 24, 21, 26] type `hpan        c3h        c3h        o3e  via equivalence: h-c-c-o`
- Remove [20, 24, 19, 21] type `c3h        o3e        c3h       hpan  via equivalence: h-c-o-c`
- Remove [20, 24, 21, 26] type `hpan        c3h        c3h       hpan  via equivalence: h-c-c-h`
- Remove [21, 19, 24, 25] type `c3h        o3e        c3h       hpan  via equivalence: h-c-o-c`
- Remove [22, 18, 21, 24] type `c3h        c3h         c2       hpan  via equivalence: c-c-c-h`
- Remove [22, 18, 21, 26] type `hpan         c2        c3h       hpan  via equivalence: h-c-c-h`
- Remove [23, 18, 21, 24] type `c3h        c3h         c2       hpan  via equivalence: c-c-c-h`
- Remove [23, 18, 21, 26] type `hpan         c2        c3h       hpan  via equivalence: h-c-c-h`
- Remove [24, 19, 21, 26] type `c3h        o3e        c3h       hpan  via equivalence: h-c-o-c`
- Remove [25, 24, 21, 26] type `hpan        c3h        c3h       hpan  via equivalence: h-c-c-h`
- Add [1, 14, 24, 20] type `c2          n         c2       hpan  via equivalence: h-c-n-c`
- Add [1, 14, 24, 21] type `c1         c2          n         c2  via equivalence: c-c-n-c`
- Add [1, 14, 24, 25] type `c2          n         c2       hpan  via equivalence: h-c-n-c`
- Add [2, 1, 14, 24] type `c1         c2          n         c2  via equivalence: c-c-n-c`
- Add [5, 1, 14, 24] type `c2          n         c2       hpan  via equivalence: h-c-n-c`
- Add [6, 1, 14, 24] type `c2          n         c2       hpan  via equivalence: h-c-n-c`
- Add [9, 10, 14, 24] type `c2          n         cp         cp  via wildcard: *-cp_-n_-*`
- Add [10, 14, 24, 20] type `cp          n         c2       hpan  via wildcard: *-c_-n_-*`
- Add [10, 14, 24, 21] type `c1         c2          n         cp  via wildcard: *-c_-n_-*`
- Add [10, 14, 24, 25] type `cp          n         c2       hpan  via wildcard: *-c_-n_-*`
- Add [11, 10, 14, 24] type `c2          n         cp         cp  via wildcard: *-cp_-n_-*`
- Add [14, 24, 21, 18] type `c2         c1         c2          n  via equivalence: c-c-c-n`
- Add [14, 24, 21, 19] type `n         c2         c1         oh  via equivalence: n-c-c-o`
- Add [14, 24, 21, 26] type `hpan         c1         c2          n  via equivalence: h-c-c-n`
- Add [16, 19, 21, 18] type `c2         c1         oh         ho  via equivalence: c-c-o-h*`
- Add [16, 19, 21, 24] type `c2         c1         oh         ho  via equivalence: c-c-o-h*`
- Add [16, 19, 21, 26] type `ho         oh         c1       hpan  via equivalence: h-c-o-h*`
- Add [17, 18, 21, 19] type `oc         c2         c1         oh  via equivalence: o-c-c-o`
- Add [17, 18, 21, 24] type `c2         c1         c2         oc  via equivalence: c-c-c-o`
- Add [17, 18, 21, 26] type `hpan         c1         c2         oc  via equivalence: h-c-c-o`
- Add [18, 21, 24, 20] type `c2         c1         c2       hpan  via equivalence: c-c-c-h`
- Add [18, 21, 24, 25] type `c2         c1         c2       hpan  via equivalence: c-c-c-h`
- Add [19, 21, 18, 22] type `hpan         c2         c1         oh  via equivalence: h-c-c-o`
- Add [19, 21, 18, 23] type `hpan         c2         c1         oh  via equivalence: h-c-c-o`
- Add [19, 21, 24, 20] type `hpan         c2         c1         oh  via equivalence: h-c-c-o`
- Add [19, 21, 24, 25] type `hpan         c2         c1         oh  via equivalence: h-c-c-o`
- Add [20, 24, 21, 26] type `hpan         c1         c2       hpan  via equivalence: h-c-c-h`
- Add [22, 18, 21, 24] type `c2         c1         c2       hpan  via equivalence: c-c-c-h`
- Add [22, 18, 21, 26] type `hpan         c1         c2       hpan  via equivalence: h-c-c-h`
- Add [23, 18, 21, 24] type `c2         c1         c2       hpan  via equivalence: c-c-c-h`
- Add [23, 18, 21, 26] type `hpan         c1         c2       hpan  via equivalence: h-c-c-h`
- Add [25, 24, 21, 26] type `hpan         c1         c2       hpan  via equivalence: h-c-c-h`

## Impropers changes

- Remove [1, 10, 14, 16] type `c2          n         cp         h*  via wildcard: *-n_-*-*`
- Remove [17, 18, 21, 22] type `c3h         c2       hpan         oc  WARNING: zero fallback`
- Remove [17, 18, 21, 23] type `c3h         c2       hpan         oc  WARNING: zero fallback`
- Remove [18, 19, 21, 24] type `c2        c3h        c3h        o3e  WARNING: zero fallback`
- Remove [18, 19, 21, 26] type `c2        c3h       hpan        o3e  WARNING: zero fallback`
- Remove [18, 21, 22, 23] type `c3h         c2       hpan       hpan  WARNING: zero fallback`
- Remove [18, 21, 24, 26] type `c2        c3h        c3h       hpan  WARNING: zero fallback`
- Remove [19, 20, 21, 24] type `c3h        c3h       hpan        o3e  WARNING: zero fallback`
- Remove [19, 20, 24, 25] type `hpan        c3h       hpan        o3e  WARNING: zero fallback`
- Remove [19, 21, 24, 25] type `c3h        c3h       hpan        o3e  WARNING: zero fallback`
- Remove [19, 21, 24, 26] type `c3h        c3h       hpan        o3e  WARNING: zero fallback`
- Remove [20, 21, 24, 25] type `c3h        c3h       hpan       hpan  WARNING: zero fallback`
- Add [1, 10, 14, 24] type `c2          n         c2         cp  via wildcard: *-n_-*-*`
- Add [14, 20, 21, 24] type `c1         c2       hpan          n  WARNING: zero fallback`
- Add [14, 20, 24, 25] type `hpan         c2       hpan          n  WARNING: zero fallback`
- Add [14, 21, 24, 25] type `c1         c2       hpan          n  WARNING: zero fallback`
- Add [17, 18, 21, 22] type `c1         c2       hpan         oc  WARNING: zero fallback`
- Add [17, 18, 21, 23] type `c1         c2       hpan         oc  WARNING: zero fallback`
- Add [18, 19, 21, 24] type `c2         c1         c2         oh  WARNING: zero fallback`
- Add [18, 19, 21, 26] type `c2         c1       hpan         oh  WARNING: zero fallback`
- Add [18, 21, 22, 23] type `c1         c2       hpan       hpan  WARNING: zero fallback`
- Add [18, 21, 24, 26] type `c2         c1         c2       hpan  WARNING: zero fallback`
- Add [19, 21, 24, 26] type `c2         c1       hpan         oh  WARNING: zero fallback`
- Add [20, 21, 24, 25] type `c1         c2       hpan       hpan  WARNING: zero fallback`
