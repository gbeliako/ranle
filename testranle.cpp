#include <cstdlib>
#include <iostream>
#include <fstream>

#include "ranle.h"

/* Example program for generating random LEs of Boolean lattices B_n  for n<=16  
* The test part exemplifies the usage of individual components, and the illustration part gives an example of practical generation

*/

// comment out to generate LEs and save to file. If active will generate LEs and test the distribution
#define TEST1

using namespace std;


/*
* call with arguments  testranle.exe  n Total KKsteps method output1 output2  
*/
int  main(int argc, char* argv[])
{
	int n = 8; int kadd = 4;
	ranle::int_64 num = 1000;
	ranle::int_64 m;
	int mode = 0, Total = 1;
	int markov = 100;
	n = 5; Total = 10;
	int rep = 2; // mixture

	string out, stats, stats1;

	// testranle.exe  n  total markov  mode repeats  outputfile statsfile statfile1
	// mode in {0,9}, 0 means KK walk, then methods 1,2,3,4,5,6-9 as per ranle.h
	// repeats=0 uses default values for mode in {6,9}, otherwise it is the size of the mixture of additive/2-additive/WOWA


	for (int i1 = 1; i1 < argc; i1++) {

		n = atoi(argv[i1]);
		i1++;
		Total = atoi(argv[i1]);
		i1++;
		markov = atoi(argv[i1]);
		i1++;
		mode = atoi(argv[i1]);
		i1++;
		rep = atoi(argv[i1]);
		i1++;

		out = (argv[i1]);
		i1++;
		stats = (argv[i1]);
		i1++;
		stats1 = (argv[i1]);
		i1++;
	}


	if (n > 12 && (mode==4 || mode==5)) n = 12;// 
	if (n > 17) n = 17; // that would be a stretch

	m = 1ULL << n;
	ranle::LEstats  Stat(n, m, Total);


/* =================
* test part: testing different methods
* 
=================*/
#ifdef TEST1


	// method 1
	ranle::LEgenerator LE(n);

	// method 2,4,5
	ranle::LESymmStart LES(n);

	// method 3,5
	ranle::LECardStart LEC(n);

	// methods 6-9
	ranle::AdditiveFM LEA(n);

	fstream fstats;
	fstats.open(stats.c_str(), ios::out);
	if (!fstats.is_open()) return 0;

	fstream fstats1;
	fstats1.open(stats1.c_str(), ios::out);
	if (!fstats1.is_open()) return 0;




	ranle::ResetTime();

	for (int num = 0; num < Total; num++) {
		switch (mode) {
		case 0:
			ranle::MakeBalanced(LE.meas.data(), m, n, LE.pos.data());
			ranle::DoMarkovChainLong(LE.meas.data(), LE.m, markov);
			break;
		case 1:
			LE.MergeLEall(markov);
			break;
		case 2:
			LES.DoWalk(markov);
			break;
		case 3:
			LEC.GenerateMerge(markov);
			break;
		case 4:
			LES.DoWalkStructured(markov);
			break;
		case 5: // mixture of 3 and 4
			if (ranle::coin()) {
				LES.DoWalkStructured(markov);
				LE.meas = LES.meas;
			}
			else {
				LEC.GenerateMerge(markov);
				LE.meas = LEC.meas;
			}
			break;

		case 6:
		case 7:
		case 8:
		case 9:

			int rr = mode - 6; // will be 0 =add 1 =2add 2 =wowa
			
			LEA.GenerateRandomLE(markov, rr, rep);
			LE.meas = LEA.LE;
			
			break;
		}
		

		// update statistics for the method chosen
		switch (mode) {
		case 3:
			Stat.UpdateAllStats(LEC.meas.data());
			break;
		case 4:
		case 2:
			Stat.UpdateAllStats(LES.meas.data());
			break;
		default:
			Stat.UpdateAllStats(LE.meas.data());
		}


	} // num

	// print out statistics

	double time = ranle::ElapsedTime();
	Stat.AverageCompute();
	fstats << "n,m,tot, mode markov time, avpos(i) (std),  avdist_bal (std) " << endl;
	fstats << n << " " << m << " " << Total << " " << mode << " " << markov << " " << time << endl;

	for (int i = 1; i < n; i++) fstats << Stat.AvPos[i] << " " << Stat.AvPosStd[i] << " ";
	fstats << endl;

	fstats << Stat.AvDistBalanced << " " << Stat.AvDistBalancedStd << " " << " " << endl;

	for (int c = 1; c < n; c++) {
		for (ranle::int_64 i = 1; i < m - 1; i++)
			fstats1 << Stat.PosProb[m * c + i] << " ";
		fstats1 << endl;
	}

#else

/*
 This is production run, prints the LEs (in binary represenation) to the specified file
 */
fstream fout;
fout.open(out.c_str(), ios::out);
if (!fout.is_open()) return 0;

	vector<ranle::int_64> Linext(m);

	// optional: if FM are needed : 
//	vector<double> FM(m);

	ranle::LEGEN LEG(n);  // just keeps all the classes with internal structures
	ranle::ResetTime();

	for (int num = 0; num < Total; num++)
	{
		LEG.generate(Linext.data(), mode, 1, markov, rep);
		// optional: 
		//LEG.generateFM(FM.data(), mode, 1, markov, rep);

		for (ranle::int_64 i = 0; i < m; i++) fout << (Linext[i]) << " ";

		// optional: 
		//for (ranle::int_64 i = 0; i < m; i++) fout << (FM[i]) << " ";

		fout << endl;

		// optional
		Stat.CalculateDistBalanced(Linext.data());
	}

	// another way: careful with large m and Total, procedural interface
	vector<ranle::int_64> LinextLarge(m*Total);
	ranle::GenerateLE(LinextLarge.data(), n, mode, Total, markov);

	// optional:
	//vector<double> FMLarge(m* Total);
	//ranle::GenerateFM(FMLarge.data(), n, mode, Total, markov);

	double time = ranle::ElapsedTime();
#endif



	return 0;

}





