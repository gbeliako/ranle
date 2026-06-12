#pragma once

#include <cstdlib>
#include <queue>
#include <vector>
#include <set>
#include <algorithm>
#include <map>
#include<random>
#include<numeric>
#include <string>
#include <unordered_map>
#include <ctime>
#include <cmath>

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
* provide space to store LEs
* vector<ranle::int_64> Linext( 1ULL<<n );// for one LE
* 	for (int num = 0; num < Total; num++)
	{
		LEG.generate(Linext.data(), mode, 1, markov, rep);
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
* The inputs are n, mode in {0,1,2,3,4,5,6,7,8,9}, and the number of additional markov steps to improve uniformity. rep is relevant to mode=6-9
* 
* mode=0 standard KK walk
* mode=1 merging incremental method (better for imbalanced LEs)
* mode=2 accelerated KK-walk from a balanced LE (better for nearly balanced LEs)
* mode=3 based on marginal distributions of k-tuple ranks (faster)
* mode=4 based on fitting the distance to balanced by method 2, but starts from a favourable balanced LE
* mode=5 is a mixture of methods 3 and 4
* mode=6 uses cardinality arrangements, starts with additive LE+symmetric
* mode=7 uses cardinality arrangements, starts with mixture of 2 2-additive LE  + symmetric(small fraction). rep>0 overwrites the size of initial mixture
* mode=8 uses cardinality arrangements, starts with mixture of 2 WOWA LE . rep>0 overwrites the size of initial mixture
* mode=9 uses cardinality arrangements, starts with mixture of 2  additive LEs . rep>0 overwrites the size of initial mixture
* mode in {6..9} is the fastest method
*
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
* Authot: Gleb Beliakov gleb@deakin.edu.au , 2025-2026
======================================================================= */


using namespace std;


namespace ranle {

	using namespace std;

	double sqr(double a) { return a * a; };


 clock_t clockS, clockF;
 double TotalTime;
 double   ElapsedTime() {
	 clockF = clock();
	 double duration = (double)(clockF - clockS) / CLOCKS_PER_SEC;
	 TotalTime += duration;
	 clockS = clockF;
	 return TotalTime;
 } ;
 void   ResetTime() { TotalTime = 0; clockS = clock(); }


	typedef unsigned long long int_64;
	typedef unsigned long  uint;

::std::random_device r;

::std::mt19937 generator(r());
::std::uniform_int_distribution<int_64> distribution2(0, 1);
::std::normal_distribution<double> distributionN(0, 1);

::std::uniform_real_distribution<double> distribution1(0.0, 1.0);




#if defined(__clang__) 

inline double beta_stable(double x, double y) {
	return std::exp(std::lgamma(x) + std::lgamma(y) - std::lgamma(x + y));
}

double Bnorm(double a, double b)
{
	return beta_stable(a, b);
}
#else

double Bnorm(double a, double b) {	return ::std::beta(a, b); }
#endif

double Beta(double x, double a, double b)
{
	double r = pow(x, a - 1) * pow(1. - x, b - 1) / Bnorm(a, b);
	return r;
}
// beta binomial distribution (simple via gamma)
double BetaBinom(int n, int k, double a, double b)
{
	double den = ::std::lgamma(a) + ::std::lgamma(b) - ::std::lgamma(a + b);
	double num = ::std::lgamma(a + k) + ::std::lgamma(n - k + b) - ::std::lgamma(a + k + n - k + b);
	double ch = ::std::lgamma(n + 1) - ::std::lgamma(k + 1) - ::std::lgamma(n - k + 1);

	return exp(ch - den + num);
}
// beta binomial distribution
double BetaBinomCancel(int n, int k, double a, double b)
{
	double logC = 0.0;
	double log_alpha_k = 0.0;
	double log_beta_nk = 0.0;
	double log_ab_n = 0.0;
	int i;
	//	a -= 1;
	//	b -= 1;
	for (i = 1; i <= k; i++) {
		log_alpha_k += log1p((a - 1) / i);
	}
	for (i = 1; i <= n - k; i++) {
		log_beta_nk += log1p((b - 1) / i);
	}
	for (i = 1; i <= n; i++) {
		log_ab_n += log1p((a + b - 1) / i);
	}
	double r = exp(log_alpha_k + log_beta_nk - log_ab_n);

	return r;

	log_alpha_k = 0.0;
	log_beta_nk = 0.0;
	log_ab_n = 0.0;
	for (i = 1; i <= k; i++) {
		logC += log(n - i + 1.) - log(i);
	}
	//	logC= std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1);
	for (i = 0; i < k; i++) {
		log_alpha_k += log(a + i);
	}
	for (i = 0; i < n - k; i++) {
		log_beta_nk += log(b + i);
	}
	for (i = 0; i < n; i++) {
		log_ab_n += log(a + b + i);
	}
	cout << exp(logC + log_alpha_k + log_beta_nk - log_ab_n) << " " << r << endl;
	return exp(logC + log_alpha_k + log_beta_nk - log_ab_n);
}




/* Some predefined constants */
int_64 maxbound = 1ULL << 50;
double stdcoef = 0.1;

int meanlength0[] = { 2,2,2,2,2, 8,  28, 80,  360, 1580, 7090, 32000, 32000 * 4, 32000 * 4 * 4,
   32000 * 4 * 4 * 4, 32000 * 4 * 4 * 4, 32000 * 4 * 4 * 4, 32000 * 4 * 4 * 4, 32000 * 4 * 4 * 4, 32000 * 4 * 4 * 4 };

double AvalMore[100];
double BvalMore[100];

double Meanv4[] = { 0,2.7126, 7.53 ,15 - 2.7126 };
double Stdev4[] = { 0,1.46,  2.12, 1.46 };

double Meanv5[] = { 0,3.375, 10.7393 , 31 - 10.7393, 31 - 3.375 };
double Stdev5[] = { 0,2.0,  3.76, 3.76 , 2.0 };


double Meanv6[] = { 0,4.0456, 14.61751 , 31.5, 48.3793, 58.961 };
double Stdev6[] = { 0,2.564,  5.8673,  7.58269 , 5.86738,  2.56401 };

double Meanv7[] = { 0,4.71041, 19.1256 , 46.5694, 80.4314, 107.87, 122.283 };
double Stdev7[] = { 0,3.11938, 8.42794,  13.3162 ,13.3161,  8.42794, 3.11938 };

double Meanv8[] = { 0,5.4014 , 24.2570 , 66.0116 , 127.502, 255 - 66.0116, 255 - 24.2570, 255 - 5.4014 };
double Stdev8[] = { 0, 3.7047 , 11.4689 , 21.3819 , 26.1591 , 21.3819 , 11.4689 , 3.7047 };

double Meanv9[] = { 0,6.0864 , 29.9765 , 90.4233 , 193.778, 511 - 193.778, 511 - 90.4233 ,511 - 29.9765 , 511 - 6.0864 };
double Stdev9[] = { 0,4.2834 , 14.9345 , 32.1514 , 46.4923 , 46.4923, 32.1514 , 14.9345 , 4.2834 };

double Meanv10[] = { 0,6.7549 , 36.3676 , 120.3010 , 284.196 , 511.500 , 1023 - 284.196 ,1023 - 120.3010 ,1023 - 36.3676 ,1023 - 6.7549 };
double Stdev10[] = { 0,4.8242 , 18.9592 , 45.9986 , 76.9454 , 91.0463 ,76.9454 , 45.9986 , 18.9592 , 4.8242 };
int MaxN = 24;
int MaxN2 = MaxN * MaxN;
double CoefMean = 0.62;
double CoefStd = 0.53;
int shift2 = 0;  
double C1 = 0.7;

inline int coin() {	return 	(int)distribution2(generator); }
inline int biasedcoin(double prob)
{
	double t = distribution1(generator);
	if (t <= prob) return 1; else return 0;
}
double facfun(int i, int n)
{
	return 1.07 + (n * 0.002 * sqrt(n)) + 0.35 * sqrt(2 * fabs(n / 2. - i) / n);
	//return 1.00 + (C1 + 0.02 * n) * sqrt(2 * fabs(n / 2. - i) / n);
}

double shift1(int i, int n)
{	return  0.15 * (n / 2. - i) / n; }





#ifdef __GNUC__
unsigned int bitweight(int_64 i) {
	return __builtin_popcountl(i);
}

#elif _MSC_VER
#  include <intrin.h>

#ifdef  _WIN64
#  define __builtin_popcountl  __popcnt64  //_mm_popcnt_u64
unsigned int bitweight(int_64 i) {
	return (uint)__builtin_popcountl(i);
}
#else
#  define __builtin_popcountl  __popcnt  //_mm_popcnt_u64
inline unsigned int bitweight(int_64 i) {
	return __builtin_popcountl((uint32_t)(i >> 32)) + __builtin_popcountl((uint32_t)i);
}
#endif

#else 
uint bitweight(int_64 v) {
	v = v - ((v >> 1) & (int_64)~(int_64)0 / 3);                           // temp
	v = (v & (int_64)~(int_64)0 / 15 * 3) + ((v >> 2) & (int_64)~(int_64)0 / 15 * 3);      // temp
	v = (v + (v >> 4)) & (int_64)~(int_64)0 / 255 * 15;                      // temp
	unsigned int c = (int_64)(v * ((int_64)~(int_64)0 / 255)) >> (sizeof(int_64) - 1) * CHAR_BIT; // count
	return (unsigned int)c;
}
//#endif
/*
unsigned int bitweight(unsigned int i)
{
	 i = i - ((i >> 1) & 0x55555555);
	 i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
	 return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}*/
#endif
int_64 UnivSetTable[20] = { 0, 1, 3, 7, 15, 31, 63, 127, 255,511,1023,2047,4095,8191,16383,32767,65535,131071,262143 };
void RemoveFromSet(int_64* A, int i) { *A &= (~(int_64(1) << i)); }
void AddToSet(int_64* A, int i) { *A |= (int_64(1) << i); }
int  IsInSet(int_64 A, int i) { return int((A >> i) & 0x1); }
int  IsSubset(int_64 A, int_64 B) { return ((A & B) == B); } //
int_64 Setunion(int_64 A, int_64 B) { return (A | B); }
int_64 Setintersection(int_64 A, int_64 B) { return (A & B); }
int_64 Setdiff(int_64 A, int_64 B) { return (A & ~(A & B)); }
int IsOdd(int i) { return ((i & 0x1) ? 1 : 0); }
int Removei_th_bitFromSet(int_64* A, int i) {
	int_64 B = (*A) & (~(int_64(1) << i));
	if (B == *A) return 1; else *A = B; return 0;
}

int_64 UniversalSet(int n) {
	int_64 A = UnivSetTable[min(n, 19)];
	while (n > 19) { AddToSet(&A, --n); }
	return A;
};

inline int Cardinality(int_64 A) { return bitweight(A); }



int_64 swapbits1(int_64 a, vector<int>& perm, int n)
{
	int_64 b = 0;
	for (int i = 0; i < n; ++i)
	{
		if (IsInSet(a, i))
			AddToSet(&b, perm[i]);
	}

	return b;
}


int_64 BinomialCoefficient(const int n, const int k) {
//	assert(n >= k);
//	assert(k >= 0);

	if (k == 0) {
		return 1;
	}

	// Recursion ->  (n, k) = (n - 1, k - 1) * n / k
	int_64 step1 = n - k + 1;
	int_64 step0;
	for (int i = 1; i < k; ++i) {
		step1 = (step0 = step1) * (n - k + 1 + i) / (i + 1);
	}

	return step1;
}
int_64 choose(int i, int n) // choose i from n (binomial(n,i))
{
	if (i == 0) return 1;
	if (i == 1) return n;
	if (i == 2) return (int_64)(n * (n - 1)) / 2;
	if (i == 3) return (int_64)(n * (n - 1) * (n - 2)) / 6;
	if (i == 4) return (int_64)(n * (n - 1) * (n - 2) * (n - 3)) / 24;
	if (i == 5) return (int_64)(n * (n - 1) * (n - 2) * (n - 3) * (n - 4)) / 120;
	if (i == 6) return (int_64)(n * (n - 1) * (n - 2) * (n - 3) * (n - 4) * (n - 5)) / 120 / 6;
	if (i == 7) return (int_64)(n * (n - 1) * (n - 2) * (n - 3) * (n - 4) * (n - 5) * (n - 6)) / 120 / 42;
	if (i == 8) return (int_64)(n * (n - 1) * (n - 2) * (n - 3) * (n - 4) * (n - 5) * (n - 6) * (n - 7)) / 120 / 42 / 8;
//	if (i == 9) return (int_64)(n * (n - 1) * (n - 2) * (n - 3) * (n - 4) * (n - 5) * (n - 6) * (n - 7) * (n - 8)) / 120 / 42 / 72;

	return BinomialCoefficient(n, i);
	//(int_64)(m_factorials[n] / m_factorials[i] / m_factorials[n - i]);
}

// make a random balanced LE
void MakeBalanced(int_64* v, int_64 len, int n, int_64* posv)
{
	int_64 m = 1ULL << n;
	vector<int_64> starts(n + 1);

	int_64 start = 1;
	int_64 sz = 0;
	vector<int_64> temp(m);

	starts[0] = 0;
	starts[1] = 1;
	for (int i = 2; i < n; i++) {
		// i is cardinality
		int_64 sz = choose(i - 1, n);
		starts[i] = starts[i - 1] + sz;
	}

	for (int_64 j = 0; j < m; j++) temp[j] = j;
	::std::shuffle(temp.begin(), temp.end(), generator); // random

	for (int_64 j = 0; j < m; j++)
	{
		int c = Cardinality(temp[j]);
		v[starts[c]] = temp[j];
		if (posv != NULL) posv[temp[j]] = starts[c];
		starts[c]++;
	}
	v[0] = 0; 
	if (posv != NULL) posv[0] = 0;
	v[m - 1] = m - 1;
	if (posv != NULL) posv[m - 1] = m - 1;
}




class DistribKtuples {
public:
	int_64 Low, High, m;
	double mean;
	int Card;
	double height;
	vector<double> MeanvMore, StdevMore;
	int usebeta = 0;
	DistribKtuples() {};

	vector<double>  alpha, beta;
	vector<int>  betabin;

	int GetParam(int n, int_64 m, double mu, double sigma, double& al, double& bet)
	{ // return 1 is betabinomial else beta scaled
		if (usebeta) return 0;

		// use beta binomial
		m -= 1;

		double pi = mu / m;
		double s = pi * (1 - pi);

		double s2 = sqr(sigma) - m * s;
		s = m * m * s - sqr(sigma);
		s = s / s2;

		if (s > 1e-12) {
			// ensure s>0
			al = s * pi;
			bet = s * (1 - pi);
			return 1;
		}
		else {
			return 0;
		}
	}
	int_64 LowBorder(int card, int n)
	{
		if (card == 1) return 1;
		if (card == 2) return 3;
		// 0001111
		return (1ULL << (card)) - 1;
	}

	int_64 HighBorder(int card, int n)
	{
		if (card == 1) return (1ULL << (n - 1));

		int_64 N = 0ULL;
		for (int i = 1; i <= card; i++) N += (1ULL << (n - i));
		return N;
	}

	void Init(int n, int card)
	{
		m = 1ULL << n;
		Card = card;
		Low = LowBorder(card, n);
		High = HighBorder(card, n);

		double Fac = 1.0;
		double shift = 0;
		if (n > 10) {
			MeanvMore.resize(MaxN2, 0); StdevMore.resize(MaxN2, 0);
		};

		vector<double> mean(n + 1);
		vector<double> stdev(n + 1);
		alpha.resize(n + 1, 0);
		beta.resize(n + 1, 0);
		betabin.resize(n + 1);

		switch (n) {
		case 4:
			for (int i = 1; i < n; i++) {
				mean[i] = Meanv4[i];
				stdev[i] = Stdev4[i];
			}
			break;
		case 5:
			for (int i = 1; i < n; i++) {
				mean[i] = Meanv5[i];
				stdev[i] = Stdev5[i];
			}
			break;
		case 6:
			for (int i = 1; i < n; i++) {
				mean[i] = Meanv6[i];
				stdev[i] = Stdev6[i];
			}

			break;
		case 7:
			for (int i = 1; i < n; i++) {
				mean[i] = Meanv7[i];
				stdev[i] = Stdev7[i];
			}
			break;
		case 8:
			for (int i = 1; i < n; i++) {
				mean[i] = Meanv8[i];
				stdev[i] = Stdev8[i];
			}
			break;
		case 9:
			for (int i = 1; i < n; i++) {
				mean[i] = Meanv9[i];
				stdev[i] = Stdev9[i];
			}
			break;
		case 10:
			for (int i = 1; i < n; i++) {
				mean[i] = Meanv10[i];
				stdev[i] = Stdev10[i];
			}
			break;

		default: // predict stdev and mean for larger n

			// 0 is just a copy, 11 starts at 1

			for (int i = 1; i < 10; i++) { //just copy
				MeanvMore[MaxN * 0 + i] = Meanv10[i];
				StdevMore[MaxN * 0 + i] = Stdev10[i];
			}

			for (int j = 1; j <= n - 10; j++) {
				int locn = j + 10;
				MeanvMore[MaxN * j + 1] = MeanvMore[MaxN * (j - 1) + 1] + CoefMean;
				StdevMore[MaxN * j + 1] = StdevMore[MaxN * (j - 1) + 1] + CoefStd;
				for (int i = 2; i <= locn / 2; i++) {
					MeanvMore[MaxN * j + i] = MeanvMore[MaxN * (j - 1) + i] + MeanvMore[MaxN * (j - 1) + i - 1] + 0.5;
					StdevMore[MaxN * j + i] = StdevMore[MaxN * (j - 1) + i] + StdevMore[MaxN * (j - 1) + i - 1];
				}
				if (IsOdd(j)) { //repeat
					MeanvMore[MaxN * j + locn / 2 + 1] = (1ULL << (locn)) - 1 - MeanvMore[MaxN * j + locn / 2];
					StdevMore[MaxN * j + locn / 2 + 1] = StdevMore[MaxN * j + locn / 2];
				}
				else MeanvMore[MaxN * j + locn / 2] = ((1ULL << (locn)) - 1.) * 0.5;

			}
			int j = n - 10;
			for (int i = n / 2 + 1; i < n; i++) {
				MeanvMore[MaxN * j + i] = m - 1 - MeanvMore[MaxN * j + n - i];
				StdevMore[MaxN * j + i] = StdevMore[MaxN * (j)+n - i];
			}

			//			for (int i = 1; i < n; i++)
			//				cout << MeanvMore[i + MaxN * (n - 10)] << " " << StdevMore[i + MaxN * (n - 10)] << " ";
			//			cout << endl;

			for (int i = 1; i < n; i++) {
				mean[i] = MeanvMore[i + MaxN * (n - 10)];
				stdev[i] = StdevMore[i + MaxN * (n - 10)];
				//				cout << i << " " << AvalMore[i] << " " << BvalMore[i] << endl;
			}

		} // switch


		for (int i = 1; i < n; i++) {

			Fac = facfun(i, n);
			shift = shift1(i, n) / (m - shift2 - 1);

			int r = GetParam(n, m, mean[i] + shift, stdev[i] * (sqrt(Fac)), alpha[i], beta[i]);
			//			cout << i << " " << mean[i] << " " << stdev[i] << " " << alpha[i] << " " << beta[i] << endl;
			if (r) betabin[i] = 1; else
			{ // use standard beta
				double var = sqr(stdev[i]) / (m - shift2 - 1) / (m - shift2 - 1) * Fac;
				double me = mean[i] / (m - shift2 - 1) + shift;   // can be 

				alpha[i] = me * (me * (1 - me) / var - 1);
				beta[i] = (1 - me) * (me * (1 - me) / var - 1);
				betabin[i] = 0;
			}
		}
	};

	// returns probability of card-tuple at position pos1
	double GetProbPos(int_64 pos1, int card, int n)
	{
		//truncate
		if (pos1<Low || pos1> High) return 0;

		if (betabin[card])
		{
			//return BetaBinom(m-1, pos1, alpha[card], beta[card]);
			return BetaBinomCancel(int(m - 1), (int)pos1, alpha[card], beta[card]);
		}
		double pos = (double)pos1;

		double a = 1, b = 1, x;

		x = pos / ((m)-1);
		a = alpha[card]; b = beta[card];
		return Beta(x, a, b);
	};


};

// variants of KK-walk
void DoMarkovChainLong(int_64* v, int_64 len, int k)
{
	uniform_int_distribution<int_64> uni(1, len - 2); // guaranteed unbiased

	for (int j = 0; j < k; j++)
	{
		//if(coin(rng))
		{
			int_64 pos = uni(generator);
		L1:
			if (!IsSubset(v[pos + 1], v[pos])) {
				::std::swap(v[pos], v[pos + 1]);

				//if (0 && coin()) { pos++; goto L1; } does not help
			}

		}
	}
}
void DoMarkovChainLongAux(int_64* v, int_64 len, int k, int_64* posv)
{
	uniform_int_distribution<int_64> uni(1, len - 2); // guaranteed unbiased

	for (int j = 0; j < k; j++)
	{
		//if(coin(rng))
		{
			int_64 pos = uni(generator);
			//cout << pos << endl;
		L1:		if (!IsSubset(v[pos + 1], v[pos])) {
			//	cout << "swap  " << ShowValue(v[pos+1]) << " " << ShowValue(v[pos]) << " " << pos+1 << " " << pos << endl;

			::std::swap(posv[v[pos]], posv[v[pos + 1]]);
			::std::swap(v[pos], v[pos + 1]);
			//if ( coin()) { pos++; goto L1; }does not help
		}
		}
	}
}

struct Less_than_count {
	int_64 count;
	int operator()(const int_64& a, const int_64& b) {
		if (a > b) count++;
		return a < b;
	}
};
template <typename RandomIt>
void insertion_sort(RandomIt begin, RandomIt end) {
	for (auto it = begin + 1; it < end; ++it) {
		auto key = *it;
		auto j = it;
		while (j > begin && *(j - 1) > key) {
			*j = *(j - 1);
			--j;
		}
		*j = key;
	}
}
long long countInversionsInsertionSort(::std::vector<int_64>& arr) {
	long long inversion_count = 0;
	for (int i = 1; i < arr.size(); ++i) {
		int key = arr[i];
		int j = i - 1;

		// Count how many elements are greater than key to the left
		while (j >= 0 && arr[j] > key) {
			arr[j + 1] = arr[j];
			--j;
			++inversion_count;
		}
		arr[j + 1] = key;
	}
	return inversion_count;
}


typedef struct {
	int_64 A;
	float key;
	int rank;
} queue_el;

struct customLess0
{
	bool operator()(const queue_el l, const queue_el r) const { return l.key < r.key; }
};




// random walk from a balanced LE
class LESymmStart {

public:
	int n, k;
	int_64 m;
	vector<int_64> meas;
	vector<int_64> pos;
	vector<int_64> pool;
	vector<int_64> poolend, poolstart;

	int InitShuffleSteps, MainSteps;

	vector<bool> poolavail;

	::std::priority_queue<queue_el, vector<queue_el>, customLess0> available;

	LESymmStart(int N) {
		n = N;
		m = (int_64)1 << n;

		poolavail.resize(m, 1);
		meas.resize(m);
		pos.resize(m);
		pool.resize(m); //keeps measures in order given and for each cardinality
		poolend.resize(n + 1);
		poolstart.resize(n + 1);
		InitShuffleSteps = 0;
		MainSteps = 0;
	};

	void UpdateAvail(queue_el& A, int succ)
	{
		int c;
		long long int rank;
		switch (succ)
		{
		case 0: // cannot swap, remove from queue
			break;
		case 1: // success, add itself and subsequent el
		case 2:
			rank = A.rank;
			rank--;
			// can be <0
			c = Cardinality(A.A);

			queue_el El;
			if (rank >= 0) {
				int_64 Ap = pool[poolstart[c] + rank]; // measure 

				if (poolavail[Ap]) {
					El.A = Ap;
					El.rank = rank;
					El.key = distribution1(generator);
					poolavail[El.A] = 0;
					available.push(El);
				}
			}
			if (succ == 1) {
				A.key = distribution1(generator);
				poolavail[A.A] = 0;
				available.push(A);
			}
			else poolavail[A.A] = 1;
			break;
		}

	}


	int MakeStep()
	{
		if (available.empty()) return 0;

		queue_el El = available.top();
		available.pop();

		int_64 p = pos[El.A];

		if (!IsSubset(meas[p + 1], meas[p])) {
			::std::swap(pos[meas[p]], pos[meas[p + 1]]);
			::std::swap(meas[p], meas[p + 1]);

			if (1 && p < m - 4 && IsSubset(meas[p + 2], meas[p + 1]))// attempt to exchange the subsequent elements to avoid dead end
			{
				if (!IsSubset(meas[p + 3], meas[p + 2])) {
					::std::swap(pos[meas[p + 3]], pos[meas[p + 2]]);
					::std::swap(meas[p + 3], meas[p + 2]);
				}
			}


			if (IsSubset(meas[p + 2], meas[p + 1]) || Cardinality(meas[p + 2]) == Cardinality(meas[p + 1])) { // cannot swap and break order
				UpdateAvail(El, 2); // cannot continue swaps
			}
			else {
				UpdateAvail(El, 1); // can continue
			}
		}
		else
			UpdateAvail(El, 0); // cannot go further, no changes made
		return 1;
	};

	// this is standard walk from balanced
	void DoWalk(int K)
	{
		restart(1);
		for (k = 0; k < K; k++)
			if (!MakeStep()) break;

		if (k < K) // do generic walk
			DoMarkovChainLongAux(meas.data(), m, K - k, pos.data());
	};

	// restart for new generation from balanced random or from structured balanced
	void restart(int d)
	{
		available = ::std::priority_queue<queue_el, vector<queue_el>, customLess0>();

		fill(poolavail.begin(), poolavail.end(), 1);

		poolstart[0] = 0;
		poolstart[1] = 1;
		for (int i = 2; i < n; i++) {
			// i is cardinality
			int_64 sz = choose(i - 1, n);
			poolstart[i] = poolstart[i - 1] + sz;
		}
		poolstart[n] = m - 1;

		if (d) MakeBalanced(meas.data(), m, n, pos.data());

		vector<int_64> ranks(n + 1, 0);
		// fill pool with ranks
		for (int_64 j = 0; j < m; j++)
		{
			int c = Cardinality(meas[j]);
			int_64 pos1 = ranks[c] + poolstart[c];// just ranks , highest rank is poolsize
			pool[pos1] = meas[j];
			ranks[c]++;
		}

				// create first queue
		for (int i = 1; i < n - 1; i++) {
			queue_el El;
			int_64 sz = choose(i, n);
			El.rank = sz - 1;// start with 0
			El.key = distribution1(generator);
			El.A = meas[poolstart[i + 1] - 1];
			poolavail[El.A] = 0;
			available.push(El);
		}

	};


	void MakeBalancedStructured()
	{
		// pos will have staring positions of tuples for now
		vector<int_64> starts(n + 1);

		int_64 start = 1;
		int_64 sz = 0;
		vector<int_64> temp(m);

		starts[0] = 0;
		starts[1] = 1;
		for (int i = 2; i < n; i++) {
			// i is cardinality
			int_64 sz = choose(i - 1, n);
			starts[i] = starts[i - 1] + sz;
		}

		for (int_64 j = 0; j < m; j++) temp[j] = j;
		for (int_64 j = 1; j < m - 1; j++)
		{
			int c = Cardinality(temp[j]);
			meas[starts[c]] = temp[j];
			starts[c]++;
		}


		// recalculate
		starts[1] = 1;
		for (int i = 2; i < n + 1; i++) {
			// i is cardinality
			int_64 sz = choose(i - 1, n);
			starts[i] = starts[i - 1] + sz;
		}

		// sort
		for (int k = 1; k < n; k++) {
			int_64 curpos = 0;

			sort(&(meas[starts[k]]), &(meas[starts[k + 1]]));
			//starts[k];
		}
		meas[0] = 0;
		meas[m - 1] = m - 1;

		// reshuffle steps

		for (int i = 0; i < InitShuffleSteps; i++)
		{
			// select cardinality
			int k = (n)*distribution1(generator);
			
			if (k <= 0) k = 1;
			if (k >= n) k = n - 1;

			double r = distribution1(generator);
			int_64 pos1 = r * choose(k, n);
			r = distribution1(generator);
			int_64 pos2 = r * choose(k, n);

			swap(meas[pos1 + starts[k]], meas[pos2 + starts[k]]);
		}

		for (int_64 j = 0; j < m; j++)	pos[meas[j]] = j;
	}

	void Symmeterise()
	{
		// as usual permutations like in fmtools

		vector<int> perm;
		for (int ii = 0; ii < n; ++ii) perm.push_back(ii);
		vector<int_64> bitswapped;
		for (int_64 ii = 0; ii < m; ++ii) bitswapped.push_back(ii);
		int_64 N = UniversalSet(n);

		::std::shuffle(perm.begin(), perm.end(), generator);
		for (int_64 ii = 0; ii < m; ++ii) bitswapped[ii] = swapbits1(ii, perm, n); // initialise

		for (int_64 ii = 0; ii < m; ++ii) meas[ii] = bitswapped[meas[ii]];
	}

	// walk from specified structured balanced LE
	void DoWalkStructured(long int K)
	{

		InitShuffleSteps = (n - 4) * n * 2;
		//	MainSteps = chi2(generator) * meanlength[n] /6;

		::std::normal_distribution<double> distributionNn(0, stdcoef * meanlength0[n]);
		double r = distributionNn(generator);
		MainSteps = meanlength0[n] + r;

		MakeBalancedStructured();
		restart(0);

		for (k = 0; k < MainSteps; k++)
			if (!MakeStep()) break;

		//			cout << k << "steps" << endl;

		DoMarkovChainLongAux(meas.data(), m, (MainSteps - k) * n * n, pos.data());

		K = m * n * n;

		DoMarkovChainLongAux(meas.data(), m, K, pos.data());

		if (coin())  TakeDual(meas.data(), 0, m - 1, n);

		Symmeterise();
	};

	void TakeDual(int_64* measdest, int_64 start, int_64 end, int curv)
	{
		int_64 N = UniversalSet(curv);
		for (int_64 i = start; i <= end; i++) pos[i] = measdest[i];// copy

		for (int_64 i = start; i <= end; i++) measdest[end - i] = ~(pos[i]) & N;
	}
};


// generate based on k-tuples positions
class LECardStart {

public:
	int n, k;
	int_64 m;
	int_64 curr, curr1;
	vector<int_64> meas;
	vector<int_64> pos;

	vector<int> order;

	vector< set<int_64>> IndicesAvail; // indices for each cardinality

	vector<int_64> poolavail;
	vector<bool> notinserted;
	int predecessors_notinserted, successors_notinserted;
	int_64 argminpred, argminsucc;

	vector<DistribKtuples> Probabil;

	vector<::std::discrete_distribution<int>> distributions;

	LECardStart(int N) {
		n = N;
		m = 1ULL << n;
		IndicesAvail.resize(n); // 0 means singletons, up to n-1 tuples

		poolavail.resize(m, 1);
		meas.resize(m);
		pos.resize(m);
		notinserted.resize(m, true);
		//predecessors_notinserted.resize(m, 0); // cn it be shorter
		order.resize(n);

		vector<double> T(m);

		Probabil.resize(n);
		for (int i = 0; i < n - 1; i++) Probabil[i].Init(n, i + 1);

		for (int i = 0; i < m - 1; i++) {
			// include the first one which is just 1
			vector<double> init;
			if (i == 0) init.push_back(1); // placeholder for emptysepace
			else for (int c = 1; c < n; c++) {
				double r = Probabil[c - 1].GetProbPos(i, c, n);
				if (i <= 1 && c == 1) r = 1;
				init.push_back(r);
			}
						//
						//			cout << i << ": ";
						//			for (auto a : init) { cout << a << " "; }
						//				cout << endl;
			::std::discrete_distribution<int> Dist(init.begin(), init.end());
			distributions.push_back(Dist);
		}

	};


	int CheckPredecessors(int_64 A, int_64 it)
	{
		int minpred = n + 1;
		predecessors_notinserted = 0;
		for (int k = 0; k < n; k++) // check predecessors
		{
			int_64 B = A;
			RemoveFromSet(&B, k);
			if (B != A) { // k was in A, check B
				if (notinserted[B]) predecessors_notinserted++;
			}
		}

		if (predecessors_notinserted == 0) // can be inserted
		{		return 1;	}

		if (minpred > predecessors_notinserted) {
			minpred = predecessors_notinserted;
			argminpred = it; // the set to use next
		}

		return 0;
	};

	int CheckSuccessors(int_64 A, int_64 it)
	{
		int minpred = n + 1;
		successors_notinserted = 0;
		for (int k = 0; k < n; k++) // check successors
		{
			int_64 B = A;
			AddToSet(&B, k);
			if (B != A) { // k was in A, check B  can break here
				if (notinserted[B]) successors_notinserted++;
			}
		}

		if (successors_notinserted == 0) // can be inserted
		{		return 1;	}

		if (minpred > successors_notinserted) {
			minpred = successors_notinserted;
			argminsucc = it; // the set to use next
		}

		return 0;
	};


	int CheckPredecessorsRecursive(int_64 A, int_64 it, int card)
	{
		if (card == -1) return 1;

		for (int k = 0; k < n; k++) // check predecessors
		{
			int_64 B = A;
			RemoveFromSet(&B, order[k]); // randomise order
			if (B != A) { // k was in A, check B
				if (notinserted[B]) {
					CheckPredecessorsRecursive(B, pos[B], card - 1);
				}
			}
		}
		// insert me
		meas[curr] = A;
		notinserted[A] = 0;
		IndicesAvail[card].erase(it);
		curr++;

		return 1;
	};

	// 0-based, returns the cardinality selected randomly at position pos. card is 0-based (singletons==0)
	int RandomCard(int pos)
	{		return distributions[pos](generator);  // +1  = card
	}

	int InsertElement()
	{
		// into position curr
		if (curr1 < curr) {			return 0;		}
		if (curr == m - 1) {
			meas[curr] = m - 1;
			return 0;
		}

		// this is 1- card
		int c;

		c = RandomCard(curr);

				// if all ktuples inserted then all predecessors inserted, try larger sets
		while (IndicesAvail[c].empty() && c < n - 1) c++;
		if (c == n - 1) {
			// assert curr==m-1
			if(curr==m-1)
				meas[curr] = m - 1;
			return 0;
		}
	L1:
		for (set<int_64>::iterator it = IndicesAvail[c].begin(); it != IndicesAvail[c].end(); ++it)
		{

			int_64 A = poolavail[*it];
			int res = CheckPredecessors(A, *it);
			if (res == 1)
			{
				// insert A
				meas[curr] = A;
				notinserted[A] = 0;
				IndicesAvail[c].erase(*it);
				curr++;
				return 1;
			}


		}

		// now I have the best one to insert
		int_64 A = poolavail[argminpred];
		shuffle(order.begin(), order.end(), generator);
		CheckPredecessorsRecursive(A, argminpred, c);

		return 1;
	};

	int InsertElementBack()
	{
		// into position curr
		if (curr1 <= curr) {
			return 0;
		}

		// this is 1- card
		int c;
		c = RandomCard(curr1);
		while (c >= n-1)c--;
				// if all ktuples inserted then all predecessors inserted, try larger sets
		while (IndicesAvail[c].empty() && c > 0) c--;

		if (c == 0) {
			// assert curr==m-1
			return 0;
		}
	L1:
		for (set<int_64>::iterator it = IndicesAvail[c].begin(); it != IndicesAvail[c].end(); ++it)
		{
			int_64 A = poolavail[*it];
			int res = CheckSuccessors(A, *it);
			if (res == 1)
			{
				// insert A
				meas[curr1] = A;
				notinserted[A] = 0;
				IndicesAvail[c].erase(*it);
				curr1--;
				return 1;
			}


		}
		// none of the c tuples can be inserted, try c+1
		{
			if (c < n - 2)
			{
				c++;
				goto L1;
			}
		}
		// now I have the best one to insert
//		int_64 A = poolavail[argminpred];
//		shuffle(order.begin(), order.end(), generator);
//		CheckPredecessorsRecursive(A, argminpred, c);

		return 1;
	};

	// use extra K steps of KK-walk
	void GenerateMerge(int K)
	{
		restart();
		int ret = 1;
		while (ret) {
			ret = InsertElement();
			ret |= InsertElementBack();
		}
		// will exit after reaching m
		if ( coin())
			TakeDual(meas.data(), 0, m - 1, n);
		// do generic walk
		DoMarkovChainLong(meas.data(), m, K);
	};

	void TakeDual(int_64* measdest, int_64 start, int_64 end, int curv)
	{
		int_64 N = UniversalSet(curv);
		for (int_64 i = start; i <= end; i++) pos[i] = measdest[i];// copy

		for (int_64 i = start; i <= end; i++) measdest[end - i] = ~(pos[i]) & N;
	}

	void restart()
	{
		// preparation

		MakeBalanced(poolavail.data(), m, n, pos.data());
		for (int c = 0; c < n; c++) IndicesAvail[c].clear();
		fill(notinserted.begin(), notinserted.end(), 1);


		// just reshuffle i, but usually the tree is self-balancing, so no need
		for (int_64 i = 1; i < m - 1; i++)
		{
			int c = Cardinality(poolavail[i]);
			IndicesAvail[c - 1].insert(i);  // but can randomise to make more balanced tree
		}

		::std::iota(order.begin(), order.end(), 0);
		meas[0] = 0;
		notinserted[0] = 0;
		curr = 1;
		curr1 = m - 2;
		meas[curr1 + 1] = m - 1;
		notinserted[m - 1] = 0;
		// we have distributions of s,p,t,...
	};
};



// generate based on merge operation, adding one singleton at a time, produces mostly imbalanced LE
class LEgenerator {
public:
	int n;
	int_64 m;
	vector<int_64> meas;
	vector<int_64> measdest;
	vector<int_64> pos, posdest;
	int markov;

	int levels;
	int_64 leaves;

	int curvar, curlevel, curpos, gap, curelement;
	vector<double> Twork, Tpwork;
	LEgenerator(int N)
	{
		n = N;
		m = (int_64)1 << n;

		levels = n - 2;
		leaves = m / 4;
		meas.resize(m);
		measdest.resize(m);
		Twork.resize(m); Tpwork.resize(m);
		pos.resize(m); posdest.resize(m);
	}




	void Symmeterise()
	{
		vector<int> perm;
		for (int ii = 0; ii < n; ++ii) perm.push_back(ii);
		vector<int_64> bitswapped;
		for (int_64 ii = 0; ii < m; ++ii) bitswapped.push_back(ii);

		int_64 N = UniversalSet(n);

		::std::shuffle(perm.begin(), perm.end(), generator);
		for (int_64 ii = 0; ii < m; ++ii) bitswapped[ii] = swapbits1(ii, perm, n); // initialise

		for (int_64 ii = 0; ii < m; ++ii) meas[ii] = bitswapped[meas[ii]];
	}
	//same as symmeterise but for a fixed current n <n
	void SymmeteriseN(int_64 start, int_64 end, int_64* meas, int curn)
	{
		// as usual permutations like in fmtools

		vector<int> perm;
		for (int ii = 0; ii < curn; ++ii) perm.push_back(ii);
		int_64 mm = 1ULL << curn;
		vector<int_64> bitswapped;
		for (int_64 ii = 0; ii < mm; ++ii) bitswapped.push_back(ii);


		int_64 N = UniversalSet(curn);

		::std::shuffle(perm.begin(), perm.end(), generator);
		for (int_64 ii = 0; ii < mm; ++ii) bitswapped[ii] = swapbits1(ii, perm, curn); // initialise

		for (int_64 ii = 0; ii < mm; ++ii) meas[start + ii] = bitswapped[meas[start + ii]];
	}
	void TakeDual(int_64* measdest, int_64 start, int_64 end, int curv)
	{
		int_64 N = UniversalSet(curv);
		for (int_64 i = start; i <= end; i++) meas[i] = measdest[i];// copy

		for (int_64 i = start; i <= end; i++) measdest[end - i] = ~(meas[i]) & N;
	}

	// main algorithm
	int MergeLE2(int_64* measure, int_64* measuredest, int_64 start, int_64 end, int levelk, int_64 setsize)
	{

		if (biasedcoin(double(curvar - 1) / n))
			MakeBalanced(&(measure[start + setsize]), setsize, curvar - 1, &(pos[start + setsize]));

		for (int_64 i = start + setsize; i <= end; i++) AddToSet(&(measure[i]), curvar - 1);

		// here make curposset the expected position of the k-tuple), read from a table
		measuredest[start + 0] = measure[start];// emptyset

		posdest[0 + start] = 0;
		curpos = 1; curelement = setsize; // 0 based

		int_64 curposset1 = curpos;
		int_64 staypos;
//		vector<double> probs(m);
		while (curpos + start < end) {
			staypos = setsize;

			int card = Cardinality(measure[start + curelement]);

		L1:
			int_64 orset = measure[start + curelement];
			RemoveFromSet(&orset, curvar - 1);

			int_64 pos1 = pos[start + orset];
			if (pos1 < curposset1) pos1 = curposset1;

			double prob =  0;

			// fill 
			while (curposset1 <= pos1 && curposset1 < staypos) {   //
				measuredest[start + curpos] = measure[start + curposset1];
				posdest[start + measure[start + curposset1]] = curpos;
				curpos++; curposset1++;
			}

			// which position to go
			int_64 diff = staypos - pos1;
	//		CalculateSwapProbabilitiesRange(diff, 0.75, probs);

			::std::exponential_distribution<> d(1);
			double ranposd = d(generator);

			int_64 ranpos;//	(distribution3(generator) * diff) / maxbound; // can be 0,.. diff-1

			ranpos = ranposd * diff / 4;

			if (ranpos > diff) ranpos = diff;

			for (int_64 i = 0; i <= ranpos; i++) {  // was <
				prob = 0.5;

				if (curposset1 < staypos && biasedcoin(prob)) {
					measuredest[start + curpos] = measure[start + curposset1];
					posdest[start + measure[start + curposset1]] = curpos;
					curpos++; curposset1++;
				}
				else break;
			}

			measuredest[start + curpos] = measure[start + curelement];
			posdest[start + measure[start + curelement]] = curpos;
			curpos++; curelement++;

			if (curelement + start >= end) break;

		}

		while (curposset1 < staypos)
		{
			measuredest[start + curpos] = measure[start + curposset1];
			posdest[start + measure[start + curposset1]] = curpos;
			curpos++; curposset1++;
		}

		measuredest[end] = measure[end];// last one always
		posdest[measure[end] + start] = end - start + 1;

		return 0;
	}

	int MergeLEall(int K)
	{
		markov = K;
		// prepare the leaves

		for (int_64 i = 0; i < leaves; i++) {
			meas[i * 4] = 0; meas[i * 4 + 3] = 3; // 12
			if (coin()) { meas[i * 4 + 1] = 1; meas[i * 4 + 2] = 2; }
			else { meas[i * 4 + 1] = 2; meas[i * 4 + 2] = 1; }

			pos[0 + i * 4] = 0; pos[i * 4 + 3] = 3; pos[i * 4 + meas[i * 4 + 1]] = 1; pos[i * 4 + meas[i * 4 + 2]] = 2;
		}

		for (int level = n - 2; level > 0; level--)
		{
			//	cout << "level " << level << endl;
			curlevel = level;
			curvar = n - level + 1;
			int_64 setsize = (int_64)1 << (n - level);
			int_64 numsets = m / setsize;

			for (int_64 s = 0; s < numsets / 2; s++) {
				int_64 startpos = s * setsize * 2;
				int_64 endpos = startpos + setsize * 2 - 1;  // inclusive

				MergeLE2(meas.data(), measdest.data(), startpos, endpos, level, setsize);

				int marsteps = 1 << int(curvar + 4);
				if (marsteps > 1e5) marsteps = 1e5;
				//marsteps = 10000;
				if (curvar > 4 && curvar < n) DoMarkovChainLongAux(&(measdest[startpos]), setsize * 2, marsteps, &(posdest[startpos]));
				//if (curvar > 4)  DoMarkovChainLongJumps(&(measdest[startpos]), setsize * 2, curvar, 100, &(posdest[startpos]) );
			}
			meas = measdest;// copy
			pos = posdest;
		}

		DoMarkovChainLong(&(meas[0]), m, markov);
		Symmeterise();

		if ( coin())
			TakeDual(&(meas[0]), 0, m - 1, n);

		return 0;
	}
};




Less_than_count less_than1;


double* valuesptr1;

int fm_sortvals1(size_t i1, size_t i2)
{
	return (valuesptr1[i1] < valuesptr1[i2]);
}

#include "wowa.cpp"

// generator based on cardinality arrangements, very fast
class AdditiveFM {
public:
	int n;
	int_64 m, d;
	int K, Rep; // Rep defaults to 0, and if specified rep, then overwrites number in the mixture
	vector<int_64>  beg, LEsym;

	vector<short> Arrangement;
	vector<bool> notinserted;

	vector<double> FM, FM2, FM3;  // for random walk
	vector<int_64> LE;
	double a, b;
	vector<double> W, WT, W2;
	// A gamma distribution with alpha = 1, and beta = 2
// approximates an exponential distribution.

	int fm_arraysize_2add(int n)
	{
		// calculates the number of parameers needed in cardinal representation for 2-additive capacity
		// no 0 included !!!
		return (int)(choose(2, n)) + n;
	}

	AdditiveFM(int N)
	{
		n = N;
		m = 1ULL << n;
		d = m - 2;
		Rep = 0;
		K = 10;

		W.resize(n + 1);
		WT.resize(n + 1);
		FM.resize(m); FM2.resize(m); FM3.resize(m);
		LE.resize(m);

		W2.resize(n * n);

		LEsym.resize(m);
		Arrangement.resize(m);
		notinserted.resize(m);
		beg.resize(n + 1); //n

		beg[0] = 0;
		beg[1] = 1;
		for (int i = 2; i < n; i++) {
			// i is cardinality
			int_64 sz = choose(i - 1, n);
			beg[i] = beg[i - 1] + sz;
		}
		beg[n] = m;


	}

	void GenerateAdditive()
	{// FM will be additive
		GenerateOnSimplex();

		for (int_64 A = 1; A < m - 1; A++) {
			FM[A] = 0;
			for (int i = 0; i < n; i++) if (IsInSet(A, i)) FM[A] += W[i];
		}
		FM[m - 1] = 1; FM[0] = 0;
	};

	void GenerateOnSimplex1(vector<double>& w, int_64 dim)
	{
		for (int_64 i = 0; i < dim - 1; i++) w[i] = distribution1(generator);
		sort(w.begin(), w.begin() + dim - 1, less<double>());
		w[dim - 1] = 1;
		for (int_64 i = dim - 1; i > 0; i--) w[i] = w[i] - w[i - 1];
	}
	int generate_fm_2additive(int n, double* vv)
	{   // simply take convex combination of all vertices
		//option==1 means starting with 1 (include emptyset in the output), otherwise not , start with 0
		int length = fm_arraysize_2add(n);
		vector<double> w(length);
		vector<double> values(length);
		// take 0 or not? for now without 0

		{
			GenerateOnSimplex1(w, length);
			// decide on =ve and -ve by random flip
			for (int i = 0; i < n; i++)
				values[i] = w[i];
			int u = n;
			for (int i = 0; i < n - 1; i++)
				for (int j = i + 1; j < n; j++)
				{
					double r = distribution1(generator);
					if (r < 0.5) {
						values[u] = -w[u];
						values[i] += w[u];
						values[j] += w[u];
					}
					else {
						values[u] = w[u];
					}
					u++;
				}

			for (int_64 A = 0; A < length; A++) vv[A] = values[A];
			// no empty set

		}
		return (length + 1); //returns the length of each FM
	}


	void generate_wowa(int type)
	{
		vector<double> w(n);
		vector<double> p(n);
		vector<double> p1(n);
		vector<double> p2(n);
		vector<double> p3(n);
		vector<double> x(n);
		vector<double> temp(12 * (n + 1));

		GenerateOnSimplex1(w, n);
		GenerateOnSimplex1(p, n);
		GenerateOnSimplex1(p1, n);
		GenerateOnSimplex1(p2, n);
		for (int i = 0; i < n; i++) p[i] = 0.5 * (w[i] + p[i]);


		GenerateOnSimplex1(w, n);
		GenerateOnSimplex1(p1, n);
		for (int i = 0; i < n; i++) w[i] = 0.5 * (w[i] + p1[i]);

		int T;

		if (type == 0) {
			weightedOWAQuantifierBuild(p.data(), w.data(), n, temp.data(), T);
		}
		double coef = 0.6;

		int L = 5;
		FM[0] = 0;
		for (int_64 A = 1; A < m - 1; A++) {
			for (int i = 0; i < n; i++)
				if (IsInSet(A, i)) x[i] = 1; else x[i] = 0;

			switch (type) {
			case 0:
				FM[A] = weightedOWAQuantifier(x.data(), p.data(), w.data(), n, temp.data(), T);
				break;
			case 1:
				FM[A] = ImplicitWOWA(x.data(), p.data(), w.data(), n);
				break;
			case 2:
				FM[A] = weightedf(x.data(), p.data(), w.data(), n, &OWA, L);
				break;

			}
		}

		random_coefficients(n, WT); // sorted
		for (int_64 A = 1; A < m - 1; A++)
			FM[A] = coef * FM[A] + (1. - coef) * WT[Cardinality(A)];


		FM[m - 1] = 1;
	}


	void Generate2Additive()
	{// FM will be 2-additive
		generate_fm_2additive(n, W2.data());


		for (int_64 A = 1; A < m - 1; A++) {
			FM[A] = 0;
			for (int i = 0; i < n; i++) if (IsInSet(A, i)) FM[A] += W2[i];

			int u = n;
			for (int i = 0; i < n - 1; i++)
				for (int j = i + 1; j < n; j++)
				{
					if (IsInSet(A, i) && IsInSet(A, j)) FM[A] += W2[u];
					u++;
				}
		}
		FM[m - 1] = 1; FM[0] = 0;
	};

	void random_coefficients(int n, vector<double>& c)
		//Generates a vector of n random real numbers 0<=X1<=X1<=...<=Xn<=1
	{
		c[0] = 0; c[n] = 1;
		for (int i = 1; i < n; i++) {
			c[i] = (double)distribution1(generator);

		}
		sort(c.begin(), c.end());
	}
	void GenerateOnSimplex()
	{
		WT[n] = 1;
		random_coefficients(n, WT); //
		for (int i = 0; i < n; i++) {
			W[i] = WT[i + 1] - WT[i];
		}
	};


	void GenerateMixtureAdditive() {
		// just additive and add symmetric


		GenerateOnSimplex();
		double r = 0;
		double WW = 0.7;
		double w = WW * 1. / n;
		random_coefficients(n, WT); // sorted
		double coef = (double)distribution1(generator) * 2. / (n + 1);

		/*
		for (int i = 0; i < n; i++) {
			W[i] += w;
			r += W[i];
		}

		r = 1. / r;
		for (int i = 0; i < n; i++) {
			W[i] *= r;
		}
		*/

		for (int_64 A = 1; A < m - 1; A++) {
			FM[A] = 0;
			for (int i = 0; i < n; i++) if (IsInSet(A, i)) FM[A] += coef * W[i];

			FM[A] += (1 - coef) * WT[Cardinality(A)];
		}
		FM[m - 1] = 1; FM[0] = 0;

	};



	void GenerateMixtureAdditive2() {
		// just additive and add symmetric
		GenerateOnSimplex();

		int t = n / 2 - 1;
		if (Rep > 0) t = Rep;

		//if (n > 8) t += 1;
		//t = 3;
		double w = 1. / t;

		for (int_64 i = 0; i < n; i++) FM2[i] = w * W[i];
		for (int k = 0; k < t - 1; k++) {
			GenerateOnSimplex();
			for (int_64 i = 0; i < n; i++) FM2[i] += w * W[i];
		}
		for (int_64 i = 0; i < n; i++) W[i] = FM2[i];



		w = 1 / sqrt(n + 0.0);
		w = 1. / n;
		for (int i = 0; i < n; i++) {
			W[i] += w;
		}
		double r = 0;
		for (int i = 0; i < n; i++) r += W[i];

		r = 1. / r;
		for (int i = 0; i < n; i++) {
			W[i] *= r;
		}

		//		for (int_64 i = 0; i < n; i++)cout << W[i] << " ";
		//		cout << endl;

		for (int_64 A = 1; A < m - 1; A++) {
			FM[A] = 0;
			for (int i = 0; i < n; i++) if (IsInSet(A, i)) FM[A] += W[i];
		}
		FM[m - 1] = 1; FM[0] = 0;
	};


	void GenerateArrangement() {

		for (int_64 i = 1; i < m - 1; i++)
		{
			Arrangement[i] = Cardinality(LE[i] + 0);
			//		cout <<Arrangement[i] << " ";
		}
		//	cout << endl;
		Arrangement[0] = 0;
		Arrangement[m - 1] = n;
	};

	void GenerateFromArrangement() {

		fill(notinserted.begin(), notinserted.end(), 1);
		fill(LE.begin(), LE.end(), 0); //for debug
		MakeBalanced(LEsym.data(), m, n, NULL);
		notinserted[0] = 0;
		int r = 0;

		for (int_64 i = 1; i <= d; i++)
		{
			int res = 0;
			short c = Arrangement[i];
L1:
			for (int_64 j = beg[c]; j < beg[c + 1]; j++)
			{
				// check if can be inserted
				if (CanBeInserted(i, LEsym[j])) {
					LE[i] = LEsym[j];
					notinserted[LEsym[j]] = 0; // important
					res = 1;
					break;
				}
			}
			// if we are here, no tuple found, should be an error
			if (!res && i <= d) {
//				cout << i << " "<<c<<endl;
				r = RecursiveInsert(i, c);
				if (!r && c < n) {
					c++; goto L1;// all c-tuples are already in
					}
			}
		}

		LE[0] = 0;
		LE[m - 1] = m - 1;
	}

	void MarkovStepsArrangement() {
		uniform_int_distribution<int_64> uni(3, m - 2); // guaranteed unbiased
		for (int i = 0; i < K; i++)
		{
			int_64 pos = uni(generator);
			if (Arrangement[pos] > Arrangement[pos + 1] /* || Arrangement[pos] == 1 || Arrangement[pos + 1] == n - 1*/)
				::std::swap(Arrangement[pos], Arrangement[pos + 1]);
		}
		Arrangement[0] = 0;
		Arrangement[m - 1] = n;
	};

	int CanBeInserted(int_64 pos, int_64 A)
	{

		if (!notinserted[A]) return 0; // already used

		// check all predecessors (which are in notinserted)
		int predecessors_notinserted = 0;
		for (int k = 0; k < n; k++) // check predecessors
		{
			int_64 B = A;
			RemoveFromSet(&B, k);
			if (B != A) { // k was in A, check B
				if (notinserted[B]) {
					predecessors_notinserted++;
					break; // return 0
				}
			}
		}

		if (predecessors_notinserted == 0) // can be inserted
		{
			return 1;
		}

		return 0;
	}

	int RecursiveInsert(int_64& i, int c) {
		int r = 0;
		for (int_64 j = beg[c]; j < beg[c + 1]; j++)
		{
			if (notinserted[LEsym[j]]) {
				RecursiveInsert1(i, LEsym[j]);
				LE[i] = LEsym[j];
				notinserted[LEsym[j]] = 0;
	//			cout << "inserted " << i << " " << LE[i] << endl;
				r = 1;
				break;
			}
		}
		return r; // 0 means nothing was inserted, all are in
	};
	void RecursiveInsert1(int_64& i, int_64 A) {
		// insert all predecessors
		int_64 B;
		for (int k = 0; k < n; k++) {
			B = A;
			RemoveFromSet(&B, k);
			if (B != A) {
				if (notinserted[B]) {
					RecursiveInsert1(i, B);
					notinserted[B] = 0;
					LE[i] = B;
//					cout << i << " " << B << endl;
					i++;
				}
			}
		}
	}
	// from FM get LE
	void DistillLE()
	{
		// sort
		valuesptr1 = FM.data();

		::iota(LE.begin(), LE.end(), 0);
		// we use as key the values in FM
		stable_sort(LE.begin(), LE.end(), fm_sortvals1);
	}

	// use K steps, start anywhere
	void GenerateRandomLE(int K, int mode, int rep=0)
	{
		if (rep > 0) Rep = rep;

		// mode: 1 mixture add, 2 2add+add, 3 wowa mix, 4 2add, 5 2inter, 6 mix 2add
		double C = 0.5;
		int num = 2;
		if (Rep > 0) num = Rep;
		double w = 1.0 / num;

		//cout << Rep << " "<<mode<< endl;

		switch (mode) {
		case 0:
			GenerateMixtureAdditive();
			break;

		case 3:
			GenerateMixtureAdditive2();
			break;

		case 2:
			for (int i = 0;i < num;i++) {
				generate_wowa(0);
				if(i>0)
					for (int_64 A = 1; A < m - 1; A++) FM2[A] += FM[A] * w;
				else 
					for (int_64 A = 1; A < m - 1; A++) FM2[A] = FM[A] * w;
			}
			//generate_wowa(0);
			for (int_64 A = 1; A < m - 1; A++) FM[A] =  FM2[A];
			break;


		case 1:
			
		//	int len = ::generate_fm_2additive(num, n, 0, FM2.data());

			random_coefficients(n, WT); // sorted
			double coef = (double)distribution1(generator) * n / (n + 1);
			coef = 0.9;

			w = coef * w;
			//num = 1;
			for (int_64 A = 1; A < m - 1; A++) FM2[A] = 0;

			for (int i = 0;i < num;i++) {
				Generate2Additive();
				for (int_64 A = 1; A < m - 1; A++) FM2[A] += FM[A] * w;

			}

			for (int_64 A = 1; A < m - 1; A++) {
				FM[A] = FM2[A] +(1 - coef) * WT[Cardinality(A)];
			}

			FM[m - 1] = 1; FM[0] = 0;
			break;

		}


		FM[0] = 0; FM[m - 1] = 1;


		DistillLE(); // done, LE has the output

		GenerateArrangement();

		MarkovStepsArrangement();

		GenerateFromArrangement();

		vector<int_64> Pos(m);

		for (int_64 i = 0; i < m; i++)
			Pos[LE[i]] = i;


		if (1 && coin())
			TakeDual(LE.data(), 0, m - 1, n);
		//DoMarkovChainLongAuxEqualised(LE.data(), LE.size(), n, K, Pos.data());
		DoMarkovChainLong(LE.data(), LE.size(), K);


	}


	void TakeDual(int_64* measdest, int_64 start, int_64 end, int curv)
	{
		//
		vector<int_64> pos(m);
		int_64 N = UniversalSet(curv);
		for (int_64 i = start; i <= end; i++) pos[i] = measdest[i];// copy

		for (int_64 i = start; i <= end; i++) measdest[end - i] = ~(pos[i]) & N;

	}
};

// this class is for measuring some statistics of the generated LEs Calculates distance to balanced (Kendall's tau) and average tuple positions
class LEstats {
public:

	int n, Tot, current;
	int_64 m;
	vector<double> AvPos; // average pos if elements
	vector<double> PosProb;
	vector<int_64> Count;

	vector<double> AvPosStd; // Need to keep for all measures, hard

	vector<int_64> DistBalanced;

	double AvDistBalanced, TotInvertedPairs;
	double AvDistBalancedStd, TotInvertedPairsStd;

	LEstats(int N, int_64 M, int total) {
		n = N; m = M; Tot = total;
		AvPos.resize(n + 1, 0); 
		Count.resize(n + 1, 0); AvPosStd.resize(n + 1, 0);
		DistBalanced.resize(total, 0);
		PosProb.resize(m * (n), 0);
		current = 0;
	}

	void UpdateAllStats(int_64* measure)
	{
		current++;
		UpdateAvPos(measure);
		UpdatePosProb(measure);
		UpdateDistBalanced(measure);
	}

	void UpdateAvPos(int_64* measure)
	{
		for (int_64 i = 1; i < m - 1; i++) {
			int c = Cardinality(measure[i]);

			Count[c]++; // at the end we need AvPos/Count

			double delta = i - AvPos[c];
			AvPos[c] += delta / Count[c];
			AvPosStd[c] += delta * (i - AvPos[c]);
		}
	}
	void UpdatePosProb(int_64* measure)
	{
		for (int_64 i = 1; i < m - 1; i++) {
			int c = Cardinality(measure[i]);
			if(c<n)
			PosProb[c * m + i] += 1;
		}
	}


	void UpdateDistBalanced(int_64* measure)
	{
		int d = CalculateDistBalanced(measure);
		DistBalanced[current - 1] = d;
	}


	void AverageCompute()
	{
		for (int_64 i = 1; i < n; i++) {
			AvPosStd[i] /= ::std::max(int_64(1), Count[i] - 1);
			AvPosStd[i] = sqrt(AvPosStd[i]);
		}

		for (int_64 c = 1; c < n; c++)  for (int_64 i = 0; i < m; i++) PosProb[c * m + i] /= (m - 2) * Tot;// msut be choose(i,n)

		double avd = 0;
		AvDistBalancedStd = 0;
		for (int_64 i = 0; i < Tot; i++) {
			double delta = DistBalanced[i] - avd;

			avd += delta / (i + 1);
			AvDistBalancedStd += delta * (DistBalanced[i] - avd);
		}
		AvDistBalanced = avd;
		AvDistBalancedStd /= (Tot - 1);
		AvDistBalancedStd = sqrt(AvDistBalancedStd);
	}

	int_64 CalculateDistBalanced(int_64* measure)
	{
		// by sorting
		vector<int_64> ranks(m);
		for (int_64 i = 0; i < m; i++)  ranks[i] = Cardinality(measure[i]); // target

		int_64 c = countInversionsInsertionSort(ranks);
		return c;
	}
};


// this is a wrapper class for all the methods
class LEGEN {
public:
	int n;
	// method 1
	ranle::LEgenerator *LE;

	// method 2,4,5
	ranle::LESymmStart *LES;

	// method 3,5
	ranle::LECardStart *LEC;

	//methods 6-9
	ranle::AdditiveFM* LEA;

	LEGEN(int N) {
		n = N;

		if (n > 17) n = 17; // that would be a stretch for safety

		LE = new LEgenerator(n);
		LES = new LESymmStart(n);
		LEC = new LECardStart(n);
		LEA = new AdditiveFM(n);
	};
	~LEGEN() {
		delete LE; delete LES; delete LEC; delete LEA;
	};

	void generate(int_64 * L, int mode, int Total, int markov, int rep=0)
	{
		int_64 p = 0;

		// add logic for handling n
		if (n > 12 && (mode == 4 || mode == 5)) mode=3 ; //  otherwise will be too long

		int_64 m = 1ULL << n;

		switch (mode) {
		case 0:

			for (int num = 0; num < Total; num++) {
				MakeBalanced(LE->meas.data(), m, n, LE->pos.data());
				DoMarkovChainLong(LE->meas.data(), LE->m, markov);
				copy(LE->meas.begin(), LE->meas.end(), L + p);
				p += m;
			}

			break;
		case 1:

			for (int num = 0; num < Total; num++) {
				LE->MergeLEall(markov);
				copy(LE->meas.begin(), LE->meas.end(), L + p);
				p += m;
			}

			break;
		case 2:

			for (int num = 0; num < Total; num++) {
				LES->DoWalk(markov);
				copy(LES->meas.begin(), LES->meas.end(), L + p);
				p += m;
			}

			break;
		case 3:

			for (int num = 0; num < Total; num++) {
				LEC->GenerateMerge(markov);
				copy(LEC->meas.begin(), LEC->meas.end(), L + p);
				p += m;
			}

			break;
		case 4:

			for (int num = 0; num < Total; num++) {
				LES->DoWalkStructured(markov);
				copy(LES->meas.begin(), LES->meas.end(), L + p);
				p += m;
			}

			break;

		case 5: // mixture of 3 and 4


			for (int num = 0; num < Total; num++) {
				if (ranle::coin()) {
					LES->DoWalkStructured(markov);
					copy(LES->meas.begin(), LES->meas.end(), L + p);
					p += m;
				}
				else {
					LEC->GenerateMerge(markov);
					copy(LEC->meas.begin(), LEC->meas.end(), L + p);
					p += m;
				}
			}


			break;
		case 6:
		case 7:
		case 8:
		case 9:
			int rr = mode - 6; // will be 0 =add 1 =2add 2 =wowa
			for (int num = 0; num < Total; num++) {
				LEA->GenerateRandomLE(markov, rr, rep);
				copy(LEA->LE.begin(), LEA->LE.end(), L + p);
				p += m;
			}

			break;
		}

	};

	void generateFM(double* FM, int mode, int Total, int markov, int rep = 0) {
		int_64 m = 1ULL << n;
		int_64 p = 0;

		vector<int_64> L(m);
		vector<double> F(m);
		vector<double> F0(m);
		for (int num = 0; num < Total; num++) {
			generate(L.data(), mode, 1, markov, rep);
			LEA->random_coefficients(m, F0); F0[m - 1] = 1;
			for (auto A = 0;A < m;A++) F[L[A]] = F0[A]; // enforce the order to L
			copy(F.begin(), F.end(), FM + p);
			p += m;
		}
	}
};

// wrapper function for all methods
int GenerateLE(int_64* L, int n, int mode, int Total, int markov, int rep=0)
{
	int retval = 0;
	if (n > 12 && (mode == 4 || mode == 5)) {
		n = 12; retval = 1;
	} // 
	if (n > 17) {
		retval = 1;  n = 17;
	} // that would be a stretch

	int_64 p = 0;
	int_64 m = 1ULL << n;

	// method 1
	LEgenerator *LE;

	// method 2,4,5
	LESymmStart *LES;

	// method 3,5
	LECardStart *LEC;

	AdditiveFM* LEA;

	switch (mode) {
	case 0:
		LE = new LEgenerator(n);
		for (int num = 0; num < Total; num++) {
			MakeBalanced(LE->meas.data(), m, n, LE->pos.data());
			DoMarkovChainLong(LE->meas.data(), LE->m, markov);
			copy(LE->meas.begin(), LE->meas.end(), L + p);
			p += m;
		}
		delete LE;
		break;
	case 1:
		LE = new LEgenerator(n);
		for (int num = 0; num < Total; num++) {
			LE->MergeLEall(markov);
			copy(LE->meas.begin(), LE->meas.end(), L + p);
			p += m;
		}
		delete LE;
		break;
	case 2:
		LES = new LESymmStart(n);
		for (int num = 0; num < Total; num++) {
			LES->DoWalk(markov);
			copy(LES->meas.begin(), LES->meas.end(), L + p);
			p += m;
		}
		delete LES;
		break;
	case 3:
		LEC = new LECardStart(n);
		for (int num = 0; num < Total; num++) {
			LEC->GenerateMerge(markov);
			copy(LEC->meas.begin(), LEC->meas.end(), L + p);
			p += m;
		}
		delete LEC;
		break;
	case 4:
		LES = new LESymmStart(n);
		for (int num = 0; num < Total; num++) {
			LES->DoWalkStructured(markov);
			copy(LES->meas.begin(), LES->meas.end(), L + p);
			p += m;
		}
		delete LES;
		break;

	case 5: // mixture of 3 and 4
		LEC = new LECardStart(n);
		LES = new LESymmStart(n);

		for (int num = 0; num < Total; num++) {
			if (ranle::coin()) {
				LES->DoWalkStructured(markov);
				copy(LES->meas.begin(), LES->meas.end(), L + p);
				p += m;
			}
			else {
				LEC->GenerateMerge(markov);
				copy(LEC->meas.begin(), LEC->meas.end(), L + p);
				p += m;
			}
		}

		delete LEC;
		delete LES;
		break;
	case 6:
	case 7:
	case 8:
	case 9:
		LEA = new AdditiveFM(n);
		int rr = mode - 6; // will be 0 =add 1 =2add 2 =wowa
		for (int num = 0; num < Total; num++) {
			LEA->GenerateRandomLE(markov, rr, rep);
			copy(LEA->LE.begin(), LEA->LE.end(), L + p);
			p += m;
		}
		delete LEA;
	}

	return retval;

}
// wrapper function for all methods, this is not too effective as it initialises the internal classes every time
int GenerateFM(double* FM, int n, int mode, int Total, int markov, int rep = 0) {
	int_64 m = 1ULL << n;
	int_64 p = 0;
	int block = 128;
	int curr = 0;
	if (Total > block) block = Total;
	AdditiveFM* LEA;
	LEA = new AdditiveFM(n);

	vector<int_64> L(m*block);
	vector<double> F(m);
	vector<double> F0(m);
	for (auto num = 0; num < Total; num+=block) {
		curr += block;
		int needed = block;
		if (Total < curr + block) needed = Total - curr;
		if(needed>0)
			GenerateLE(L.data(), mode, needed, markov, rep);

		for (auto i = 0;i < needed; i++) {
			LEA->random_coefficients(m, F0); F0[m - 1] = 1;
			for (auto A = 0;A < m;A++) F[ L[A+i*m] ] = F0[A]; // enforce the order to L
			copy(F.begin(), F.end(), FM + p);
			p += m;
		}
	}
	delete LEA;
}

};



