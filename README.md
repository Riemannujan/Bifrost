This project is a toy experiment of the transciphering HE/FE framework Bifröst. The HE scheme chosen is BGV, and the experiments are run for a lattice-based IPFE from Mera et al. (2023). A preliminary version aimed at evaluating a lattice-based version of SPADE from Nuoskala et al. (2024), but this project is currently on hold.

The folder is organized as follows:
- utils.py contains the side functions and the polynomial class definition. This is not relevant to the experiments
- BGV.py contains the algorithms (KeyGen, Enc, Dec, Add, Mul) of the BGV encryption scheme. This is a short version of BGV as we only implement the algorithms useful to Bifröst
- IPFE.py contains the core algorithms (Setup, Enc, KeyDer, Dec) for the IPFE scheme.
- bifrost.py contains the implementation of the BGV evaluation algorithm on IPFE.
