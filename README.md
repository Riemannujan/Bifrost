This project is a toy experiment of the transciphering HE/FE framework Bifröst. The HE scheme chosen is BGV, and the experiments are run for two FE schemes: lattice-based IPFE from Mera et al. (2023), and a lattice-based version of SPADE from Nuoskala et al. (2024).

The folder is organized as follows:
- utils.py contains the side functions and the polynomial class definition. This is not relevant to the experiments
- BGV.py contains the algorithms (KeyGen, Enc, Dec, Add, Mul) of the BGV encryption scheme. This is a short version of BGV as we only implement the algorithms useful to Bifröst
- IPFE.py contains the core algorithms (Setup, Enc, KeyDer, Dec) for the IPFE scheme. No experiment is run in this file.
- SPADE.py contains the core algorithms (Setup, Enc, KeyDer, Dec) for the lwe-based SPADE scheme.
- bifrost.py contains the implementation of the BGV evaluation algorithm on SPADE and the BGV evaluation algorithm on IPFE (2 versions: vanilla and function-hiding)

