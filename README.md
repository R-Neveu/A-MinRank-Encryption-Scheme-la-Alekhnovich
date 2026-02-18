# MinRankPKE Proof of Concept Implementation

## Overview

This repository provides the proof of concept implementation from 'A Minrank-based Encryption Scheme à la Alekhnovich-Regev'. Timings are expressed as milliseconds (ms) to show that the computations are reasonable.

The codebase allows you to run a certain number of tests (NUM_OF_TEST = 300 by default), and then prints the average time it took.

Results from the paper were computed on a laptop with a CPU: 13th Gen Intel(R) Core(TM) i7-1365U, 32GB of RAM, 12 threads. This aims to show that the timings are reasonable on a laptop, i.e., the operations do not take several minutes to complete.

## Instructions

### Dependencies

This project depends on OpenSSL. It should be installed before running the implementation.

### Build

To build the executables, one should run the commands (each takes a bit of time to complete):

- gcc -O3 -march=native -I include -Wno-deprecated-declarations MRPKE_I.c include/*.c -lssl -lcrypto -o MRPKE_I

- gcc -O3 -march=native -I include -Wno-deprecated-declarations MRPKE_III.c include/*.c -lssl -lcrypto -o MRPKE_III

- gcc -O3 -march=native -I include -Wno-deprecated-declarations MRPKE_V.c include/*.c -lssl -lcrypto -o MRPKE_V

These commands are also available in the comments of the files MRPKE_I.c, MRPKE_III.c, MRPKE_V.c

### Comments

The explanation of the implementation can be found in the commented parts of MRPKE_I.c. The implementation uses several improvements to the scheme, namely, by using T = EH^\top instead of T = SG+E. Furthermore, we consider H=[I | H'], allowing a better efficiency. To check that the decryption worked, one can uncomment the 4 lines at the end of the computations to print the messages to see what was encrypted and what was recovered. Of course, the actual code is far from being clean. However, this is only a proof of concept aiming at showing that the scheme works and is reasonable.

### Expected results

The results obviously depend on the computer. However, the results should be consistent. We give here the timings when the computer is set in "performance mode", and when it is set in "powersave mode". This corresponds to a frequency of around 4 to 5 GHz in performance and approximately 2.2 GHz in powersave.

One should expect timings ranging (for level I) from 16 to 34 ms for KeyGen (depending on the computer), from 2 to 4 for Encryption (without the computation of the public key and its transpose), and decryption should be less than a millisecond. For level III, this should be 45 to 80 ms for KeyGen, 7 to 13 ms for Encryption, and still less than a millisecond to decrypt. Level V is roughly the same increase: 75 to 130 ms for KeyGen, 17 to 26 for Encryption, and less than a ms for Decryption (perhaps around a millisecond for slower computers).

For the computation of transpose(H), it ranges from 10 to 21 ms for level I, 27 to 50 for level III, and 40 to 75 for level V.

### Considered elements

In a previous version of the paper, the timings of KeyGen did not take into account the cost of transposing the matrix H, as we originally considered this could be done offline. This timing is now added in the KeyGen timings. Still, this keeps the scheme very practical as the KeyGen needs to be done only once. Furthermore, we also disregarded the re-generation of the public key in the encryption and considered encryption takes H (or rather, transpose(H) ) as an input.

Overall, the point of this implementation was to show that the scheme is practical: it is the case here as it takes only milliseconds to execute (with room for improvement. For instance, by using SHAKEx4 instead of AES in the generation of randomness, or by parallelising matrix multiplication or transposition).

Note that for all security levels, due to the construction of the library rbc-lib, we had to take different values of ELL1 than what we gave in the paper. The values are close, but induce a slight overhead in the timings.
