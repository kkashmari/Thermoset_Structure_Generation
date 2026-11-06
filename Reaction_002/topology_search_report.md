# Reaction-template topology search

- Decision: **new template required**
- Reason: The complete force-field-typed pre/post fix bond/react topology has not been registered before.
- Typed template key: `eaa20a1d28c6fb418c538c9a5ad34f7db3658b2dccc81b9ed6794953b11c8a6d`

## Library topology changes

- remove 1:8–1:29
- remove 2:35–2:41
- add 1:8–2:41 (1)
- add 2:35–1:29 (1)

## Atom type and charge changes

- Atom 5: `n` → `n`; charge -0.86110000 → -0.63240000
- Atom 9: `h*` → `ho`; charge +0.43950000 → +0.42410000
- Atom 17: `n` → `n`; charge -0.42160000 → -0.63240000
- Atom 20: `c3h` → `c1`; charge +0.08000000 → +0.08000000
- Atom 23: `o3e` → `oh`; charge -0.26600000 → -0.55710000
- Atom 24: `c3h` → `c2`; charge +0.02700000 → +0.10480000

## Bonds changes

- Remove [5, 9] type `h*          n`
- Remove [18, 20] type `c2        c3h  via equivalence: c-c`
- Remove [19, 24] type `c3h       hpan  via equivalence: c-h`
- Remove [20, 23] type `c3h        o3e  via equivalence: c-o`
- Remove [20, 24] type `c3h        c3h  via equivalence: c-c`
- Remove [20, 25] type `c3h       hpan  via equivalence: c-h`
- Remove [23, 24] type `c3h        o3e  via equivalence: c-o`
- Remove [24, 26] type `c3h       hpan  via equivalence: c-h`
- Add [5, 24] type `c2          n  via equivalence: c-n`
- Add [9, 23] type `ho         oh  via equivalence: h*-o`
- Add [18, 20] type `c1         c2  via equivalence: c-c`
- Add [19, 24] type `c2       hpan  via equivalence: c-h`
- Add [20, 23] type `c1         oh  via equivalence: c-o`
- Add [20, 24] type `c1         c2  via equivalence: c-c`
- Add [20, 25] type `c1       hpan  via equivalence: c-h`
- Add [24, 26] type `c2       hpan  via equivalence: c-h`

## Angles changes

- Remove [2, 5, 9] type `c3          n         h*  via equivalence: c-n-h*`
- Remove [9, 5, 14] type `c2          n         h*  via equivalence: c-n-h*`
- Remove [17, 18, 20] type `c3h         c2          n  via equivalence: c-c-n`
- Remove [18, 20, 23] type `c2        c3h        o3e  via equivalence: c-c-o`
- Remove [18, 20, 24] type `c2        c3h        c3h  via equivalence: c-c-c`
- Remove [18, 20, 25] type `c2        c3h       hpan  via equivalence: c-c-h`
- Remove [19, 24, 20] type `c3h        c3h       hpan  via equivalence: c-c-h`
- Remove [19, 24, 23] type `hpan        c3h        o3e  via equivalence: h-c-o`
- Remove [19, 24, 26] type `hpan        c3h       hpan  via equivalence: h-c-h`
- Remove [20, 18, 21] type `c3h         c2       hpan  via equivalence: c-c-h`
- Remove [20, 18, 22] type `c3h         c2       hpan  via equivalence: c-c-h`
- Remove [20, 23, 24] type `c3h        o3e        c3h  via equivalence: c-o-c`
- Remove [20, 24, 23] type `c3h        c3h        o3e  via equivalence: c-c-o`
- Remove [20, 24, 26] type `c3h        c3h       hpan  via equivalence: c-c-h`
- Remove [23, 20, 24] type `c3h        c3h        o3e  via equivalence: c-c-o`
- Remove [23, 20, 25] type `hpan        c3h        o3e  via equivalence: h-c-o`
- Remove [23, 24, 26] type `hpan        c3h        o3e  via equivalence: h-c-o`
- Remove [24, 20, 25] type `c3h        c3h       hpan  via equivalence: c-c-h`
- Add [2, 5, 24] type `c2          n         c3  via equivalence: c-n-c`
- Add [5, 24, 19] type `hpan         c2          n  via equivalence: h-c-n`
- Add [5, 24, 20] type `c1         c2          n  via equivalence: c-c-n`
- Add [5, 24, 26] type `hpan         c2          n  via equivalence: h-c-n`
- Add [9, 23, 20] type `c1         oh         ho  via equivalence: c-o-h*`
- Add [14, 5, 24] type `c2          n         c2  via equivalence: c-n-c`
- Add [17, 18, 20] type `c1         c2          n  via equivalence: c-c-n`
- Add [18, 20, 23] type `c2         c1         oh  via equivalence: c-c-o`
- Add [18, 20, 24] type `c2         c1         c2  via equivalence: c-c-c`
- Add [18, 20, 25] type `c2         c1       hpan  via equivalence: c-c-h`
- Add [19, 24, 20] type `c1         c2       hpan  via equivalence: c-c-h`
- Add [19, 24, 26] type `hpan         c2       hpan  via equivalence: h-c-h`
- Add [20, 18, 21] type `c1         c2       hpan  via equivalence: c-c-h`
- Add [20, 18, 22] type `c1         c2       hpan  via equivalence: c-c-h`
- Add [20, 24, 26] type `c1         c2       hpan  via equivalence: c-c-h`
- Add [23, 20, 24] type `c2         c1         oh  via equivalence: c-c-o`
- Add [23, 20, 25] type `hpan         c1         oh  via equivalence: h-c-o`
- Add [24, 20, 25] type `c2         c1       hpan  via equivalence: c-c-h`

## Dihedrals changes

- Remove [3, 2, 5, 9] type `c1         c3          n         h*  via equivalence: c-c-n-h*`
- Remove [7, 2, 5, 9] type `c1         c3          n         h*  via equivalence: c-c-n-h*`
- Remove [9, 5, 14, 11] type `h*          n         c2       hpan  via equivalence: h-c-n-h*`
- Remove [9, 5, 14, 12] type `c1         c2          n         h*  via equivalence: c-c-n-h*`
- Remove [9, 5, 14, 16] type `h*          n         c2       hpan  via equivalence: h-c-n-h*`
- Remove [17, 18, 20, 23] type `n         c2        c3h        o3e  via equivalence: n-c-c-o`
- Remove [17, 18, 20, 24] type `c3h        c3h         c2          n  via equivalence: c-c-c-n`
- Remove [17, 18, 20, 25] type `hpan        c3h         c2          n  via equivalence: h-c-c-n`
- Remove [18, 20, 23, 24] type `c2        c3h        o3e        c3h  via equivalence: c-c-o-c`
- Remove [18, 20, 24, 19] type `c2        c3h        c3h       hpan  via equivalence: c-c-c-h`
- Remove [18, 20, 24, 23] type `c2        c3h        c3h        o3e  via equivalence: c-c-c-o`
- Remove [18, 20, 24, 26] type `c2        c3h        c3h       hpan  via equivalence: c-c-c-h`
- Remove [19, 24, 20, 23] type `hpan        c3h        c3h        o3e  via equivalence: h-c-c-o`
- Remove [19, 24, 20, 25] type `hpan        c3h        c3h       hpan  via equivalence: h-c-c-h`
- Remove [19, 24, 23, 20] type `c3h        o3e        c3h       hpan  via equivalence: h-c-o-c`
- Remove [20, 23, 24, 26] type `c3h        o3e        c3h       hpan  via equivalence: h-c-o-c`
- Remove [21, 18, 20, 23] type `hpan         c2        c3h        o3e  via equivalence: h-c-c-o`
- Remove [21, 18, 20, 24] type `c3h        c3h         c2       hpan  via equivalence: c-c-c-h`
- Remove [21, 18, 20, 25] type `hpan         c2        c3h       hpan  via equivalence: h-c-c-h`
- Remove [22, 18, 20, 23] type `hpan         c2        c3h        o3e  via equivalence: h-c-c-o`
- Remove [22, 18, 20, 24] type `c3h        c3h         c2       hpan  via equivalence: c-c-c-h`
- Remove [22, 18, 20, 25] type `hpan         c2        c3h       hpan  via equivalence: h-c-c-h`
- Remove [23, 20, 24, 26] type `hpan        c3h        c3h        o3e  via equivalence: h-c-c-o`
- Remove [23, 24, 20, 25] type `hpan        c3h        c3h        o3e  via equivalence: h-c-c-o`
- Remove [24, 23, 20, 25] type `c3h        o3e        c3h       hpan  via equivalence: h-c-o-c`
- Remove [25, 20, 24, 26] type `hpan        c3h        c3h       hpan  via equivalence: h-c-c-h`
- Add [2, 5, 24, 19] type `c3          n         c2       hpan  via equivalence: h-c-n-c`
- Add [2, 5, 24, 20] type `c1         c2          n         c3  via equivalence: c-c-n-c`
- Add [2, 5, 24, 26] type `c3          n         c2       hpan  via equivalence: h-c-n-c`
- Add [3, 2, 5, 24] type `c1         c3          n         c2  via equivalence: c-c-n-c`
- Add [5, 24, 20, 18] type `c2         c1         c2          n  via equivalence: c-c-c-n`
- Add [5, 24, 20, 23] type `n         c2         c1         oh  via equivalence: n-c-c-o`
- Add [5, 24, 20, 25] type `hpan         c1         c2          n  via equivalence: h-c-c-n`
- Add [7, 2, 5, 24] type `c1         c3          n         c2  via equivalence: c-c-n-c`
- Add [9, 23, 20, 18] type `c2         c1         oh         ho  via equivalence: c-c-o-h*`
- Add [9, 23, 20, 24] type `c2         c1         oh         ho  via equivalence: c-c-o-h*`
- Add [9, 23, 20, 25] type `ho         oh         c1       hpan  via equivalence: h-c-o-h*`
- Add [11, 14, 5, 24] type `c2          n         c2       hpan  via equivalence: h-c-n-c`
- Add [12, 14, 5, 24] type `c1         c2          n         c2  via equivalence: c-c-n-c`
- Add [14, 5, 24, 19] type `c2          n         c2       hpan  via equivalence: h-c-n-c`
- Add [14, 5, 24, 20] type `c1         c2          n         c2  via equivalence: c-c-n-c`
- Add [14, 5, 24, 26] type `c2          n         c2       hpan  via equivalence: h-c-n-c`
- Add [16, 14, 5, 24] type `c2          n         c2       hpan  via equivalence: h-c-n-c`
- Add [17, 18, 20, 23] type `n         c2         c1         oh  via equivalence: n-c-c-o`
- Add [17, 18, 20, 24] type `c2         c1         c2          n  via equivalence: c-c-c-n`
- Add [17, 18, 20, 25] type `hpan         c1         c2          n  via equivalence: h-c-c-n`
- Add [18, 20, 24, 19] type `c2         c1         c2       hpan  via equivalence: c-c-c-h`
- Add [18, 20, 24, 26] type `c2         c1         c2       hpan  via equivalence: c-c-c-h`
- Add [19, 24, 20, 23] type `hpan         c2         c1         oh  via equivalence: h-c-c-o`
- Add [19, 24, 20, 25] type `hpan         c1         c2       hpan  via equivalence: h-c-c-h`
- Add [21, 18, 20, 23] type `hpan         c2         c1         oh  via equivalence: h-c-c-o`
- Add [21, 18, 20, 24] type `c2         c1         c2       hpan  via equivalence: c-c-c-h`
- Add [21, 18, 20, 25] type `hpan         c1         c2       hpan  via equivalence: h-c-c-h`
- Add [22, 18, 20, 23] type `hpan         c2         c1         oh  via equivalence: h-c-c-o`
- Add [22, 18, 20, 24] type `c2         c1         c2       hpan  via equivalence: c-c-c-h`
- Add [22, 18, 20, 25] type `hpan         c1         c2       hpan  via equivalence: h-c-c-h`
- Add [23, 20, 24, 26] type `hpan         c2         c1         oh  via equivalence: h-c-c-o`
- Add [25, 20, 24, 26] type `hpan         c1         c2       hpan  via equivalence: h-c-c-h`

## Impropers changes

- Remove [2, 5, 9, 14] type `c2          n         c3         h*  via wildcard: *-n_-*-*`
- Remove [17, 18, 20, 21] type `c3h         c2       hpan          n  WARNING: zero fallback`
- Remove [17, 18, 20, 22] type `c3h         c2       hpan          n  WARNING: zero fallback`
- Remove [18, 20, 21, 22] type `c3h         c2       hpan       hpan  WARNING: zero fallback`
- Remove [18, 20, 23, 24] type `c2        c3h        c3h        o3e  WARNING: zero fallback`
- Remove [18, 20, 23, 25] type `c2        c3h       hpan        o3e  WARNING: zero fallback`
- Remove [18, 20, 24, 25] type `c2        c3h        c3h       hpan  WARNING: zero fallback`
- Remove [19, 20, 23, 24] type `c3h        c3h       hpan        o3e  WARNING: zero fallback`
- Remove [19, 20, 24, 26] type `c3h        c3h       hpan       hpan  WARNING: zero fallback`
- Remove [19, 23, 24, 26] type `hpan        c3h       hpan        o3e  WARNING: zero fallback`
- Remove [20, 23, 24, 25] type `c3h        c3h       hpan        o3e  WARNING: zero fallback`
- Remove [20, 23, 24, 26] type `c3h        c3h       hpan        o3e  WARNING: zero fallback`
- Add [2, 5, 14, 24] type `c2          n         c2         c3  via wildcard: *-n_-*-*`
- Add [5, 19, 20, 24] type `c1         c2       hpan          n  WARNING: zero fallback`
- Add [5, 19, 24, 26] type `hpan         c2       hpan          n  WARNING: zero fallback`
- Add [5, 20, 24, 26] type `c1         c2       hpan          n  WARNING: zero fallback`
- Add [17, 18, 20, 21] type `c1         c2       hpan          n  WARNING: zero fallback`
- Add [17, 18, 20, 22] type `c1         c2       hpan          n  WARNING: zero fallback`
- Add [18, 20, 21, 22] type `c1         c2       hpan       hpan  WARNING: zero fallback`
- Add [18, 20, 23, 24] type `c2         c1         c2         oh  WARNING: zero fallback`
- Add [18, 20, 23, 25] type `c2         c1       hpan         oh  WARNING: zero fallback`
- Add [18, 20, 24, 25] type `c2         c1         c2       hpan  WARNING: zero fallback`
- Add [19, 20, 24, 26] type `c1         c2       hpan       hpan  WARNING: zero fallback`
- Add [20, 23, 24, 25] type `c2         c1       hpan         oh  WARNING: zero fallback`
