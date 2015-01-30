# pcg32
This is a tiny self-contained C++ immplementation of the PCG32 random number generator based on code available at http://www.pcg-random.org.

I decide to make my own version since the official small version lacks a C++ interface and various features (e.g. rewind support and floating point sample generation), and the big C++ version is extremely large and uses very recent C++ features that are not yet supported by all compilers.
