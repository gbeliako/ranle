/* =======================================================================
* This package is to generate randomly linear extensions of Boolean lattices B_n (typically up to n=17) by one of the 4 methods and also with
* baseline KK-walk.
* 
* The package is C++ header-only. It uses namespace ranle.
* 
* To use it: 
*  add 
* #include "ranle.h"
* Then either use the individual classes or the wrapper class LEGEN in the following way (see also the included main file for examples)
* ranle::LEGEN LEG(n);
* 
* provide space to stors LEs
* vector<ranle::int_64> Linext( 1ULL<<n );// for one LE
* 	for (int num = 0; num < Total; num++)
	{
		LEG.generate(Linext.data(), mode, 1, markov);
		// print Linext
	}
*  or:
* vector<ranle::int_64> LinextLarge( (1ULL<<n)*Total );// for many LEs
* LEG.generate(LinextLarge.data(), mode, Total, markov);
* 
* Alternatively one can use procedural interface
* 	ranle::GenerateLE(LinextLarge.data(), n, mode, Total, markov);
* 
* 
* The inputs are n, mode in {0,1,2,3,4,5}, and the number of additional markov steps to improve uniformity
* 
* mode=0 standard KK walk
* mode=1 merging incremental method (better for imbalanced LEs)
* mode=2 accelerated KK-walk from a balanced LE (better for nearly balanced LEs)
* mode=3 based on marginal distributions of k-tuple ranks (faster)
* mode=4 based on fitting the distance to balanced by method 2, but starts from a favourable balanced LE
* mode=5 is a mixture of methods 3 and 4
* 
* declarations:
* 	void LEGEN::generate(int_64 * L, int mode, int Total, int markov)
* or 
* int GenerateLE(int_64* L, int n, int mode, int Total, int markov)
* (returns 0 if success, 1 if n is deemed too large)
* 
* The class LEstats
* keeps track of the distance to balanced and ranks of k-tuples statistics, see example main for its usage
* typically call for every LE Stat.UpdateAllStats(LE)
* and at the end Stat.AverageCompute(); and print out its members
* it is used to assess the uniformity of the distribution of the sample
* reset (reinitialise) this class for each sample
* 
* 
* 
* Authot: Gleb Beliakov gleb@deakin.edu.au , 2025
======================================================================= */
