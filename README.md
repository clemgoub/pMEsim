# pMEsim
polymorphic mobile element simulator

This repository hosts the base code to simulate polymporphic mobile elements insertions and deletions on the human chromosome 22 for validation of [`GraffiTE`](https://github.com/cgroza/GraffiTE).
It fist generate VCF files including simulated variants and background (random) insertion/deletion. The simulated VCF are next fed into [SimuG](https://github.com/yjx1217/simuG) to simulate alternative chromosome 22.

If time allows, I'm hoping to extend the script as a general-purpose pME simulator.
