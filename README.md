# pMEsim
polymorphic mobile element simulator

This repository hosts the base code to simulate polymporphic mobile elements insertions and deletions on the human chromosome 22 for validation of [`GraffiTE`](https://github.com/cgroza/GraffiTE).

It first generates VCF files including simulated variants and background (random) insertion/deletion. Next, the simulated VCFs are fed into [SimuG](https://github.com/yjx1217/simuG) to simulate alternative chromosome 22.

If time allows, I'm hoping to extend the script as a general-purpose pME simulator.
