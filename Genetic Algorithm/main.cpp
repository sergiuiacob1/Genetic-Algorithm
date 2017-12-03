#include <iostream>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#define DMAX 100
#define P_MUTATATION 0.01
#define P_CROSS 0.25
#define INF 2000000000
#define PI 3.1415
#define EPSILON 0.1
#define C 10000.0

using namespace std;

double prob[DMAX], probCumulated[DMAX];
int lgReprez[DMAX];
int discreteFactor;

double DeJong(double[], int);
double SixHump(double[], int);
double Schwefel7(double[], int);
double Rastrigin(double[], int);

void BuildLgReprez(int, double[], int);
double GeneticAlgorithm(double(*)(double[], int), int, double[][2]);
double RandomDouble(const double &, const double &);
int BinaryToInt(bool[], int);
void Copy(double[], double[], int);
void ValsToDouble(double[], bool[], int, double[]);
void BuildPartialProbabilities(double, double[], int);
double Fitness(double(*testFunction)(double[], int), double[], int);
void GenerateRandomSolution(bool[][DMAX * sizeof(int)], int, int, double[][2]);
void CreateRandomChromosome(bool[], int);
bool PopulationIsEvoluating(int, int);
double EvaluatePopulation(bool[][DMAX * sizeof(int)], int, int, double(*testFunction)(double[], int), double[][2]);
void SelectNextPopulation(bool[][DMAX * sizeof(int)], int &, int);
void MutateChromosomes();
void CrossChromosomes();

int main() {
	double fRez;
	double acceptedVals[DMAX][2];
	int nrDims, discreteFactor;

	srand((unsigned int)time(NULL));

	nrDims = 2; discreteFactor = 2;
	acceptedVals[0][0] = -5.12; acceptedVals[0][1] = 5.12;
	acceptedVals[1][0] = -5.12; acceptedVals[1][1] = 5.12;
	for (int i = 0; i < nrDims; ++i)
		BuildLgReprez(i, acceptedVals[i], discreteFactor);
	fRez = GeneticAlgorithm(Rastrigin, nrDims, acceptedVals);
	cout << "Rastrigin: " << fRez << '\n';

	cin >> fRez;//sa apara consola

	return 0;
}

void BuildLgReprez(int pos, double acceptedVals[2], int discreteFactor) {
	int lg;
	for (lg = 0; (acceptedVals[1] - acceptedVals[0]) * pow(10, discreteFactor) > pow(2, lg) - 1;)
		++lg;

	lgReprez[pos] = lg;
}

double GeneticAlgorithm(double(*testFunction)(double[], int), int nrDims, double acceptedVals[][2]) {
	double bestSol, generationResult;
	int nrIterations, popSize, currentGeneration, lastBestGeneration;
	bool chromosomes[DMAX][DMAX * sizeof(int)];

	nrIterations = 100;
	popSize = 100;
	bestSol = INF;
	for (int i = 0; i < nrIterations; ++i) {
		GenerateRandomSolution(chromosomes, popSize, nrDims, acceptedVals);
		currentGeneration = lastBestGeneration = 0;

		while (PopulationIsEvoluating(lastBestGeneration, currentGeneration)) {
			generationResult = EvaluatePopulation(chromosomes, popSize, nrDims, testFunction, acceptedVals);

			if (generationResult < bestSol) {
				bestSol = generationResult;
				lastBestGeneration = currentGeneration;
			}

			SelectNextPopulation(chromosomes, popSize, nrDims);
			MutateChromosomes();
			CrossChromosomes();

			++currentGeneration;
		}
	}

	return bestSol;
}

double EvaluatePopulation(bool chromosomes[DMAX][DMAX * sizeof(int)], int popSize, int nrDims, double(*testFunction)(double[], int), double acceptedVals[][2]) {
	double chromosomeFitness[DMAX], doubleVals[DMAX];
	double bestSol = INF, sumFitness = 0;

	for (int i = 0; i < popSize; ++i) {
		ValsToDouble(doubleVals, chromosomes[i], nrDims, acceptedVals[i]);

		chromosomeFitness[i] = Fitness(testFunction, doubleVals, nrDims);
		bestSol = min(bestSol, chromosomeFitness[i]);
		sumFitness += chromosomeFitness[i];
	}

	BuildPartialProbabilities(sumFitness, chromosomeFitness, popSize);

	return bestSol;
}

double Fitness(double(*testFunction)(double[], int), double doubleVals[], int nrDims) {
	if (testFunction == &Schwefel7)
		return abs(Schwefel7(doubleVals, nrDims) + C);

	return 1 / (testFunction(doubleVals, nrDims) + EPSILON);
}

void BuildPartialProbabilities(double sumFitness, double chromosomeFitness[DMAX], int popSize) {
	for (int i = 0; i < popSize; ++i)
		prob[i] = chromosomeFitness[i] / sumFitness;

	probCumulated[0] = prob[0];
	for (int i = 1; i < popSize; ++i)
		probCumulated[i] = probCumulated[i - 1] + prob[i];
}

void ValsToDouble(double doubleVals[DMAX], bool chromosome[DMAX * sizeof(int)], int nrDims, double acceptedVals[]) {
	int aux, sumLgReprez = 0;
	for (int i = 0; i < nrDims; ++i) {
		aux = BinaryToInt(chromosome + sumLgReprez, lgReprez[i]);
		doubleVals[i] = acceptedVals[0] + aux + (acceptedVals[1] - acceptedVals[0]) / (pow(2, lgReprez[i]) - 1);

		sumLgReprez += lgReprez[i];
	}
}

void GenerateRandomSolution(bool chromosomes[DMAX][DMAX * sizeof(int)], int popSize, int nrDims, double acceptedVals[DMAX][2]) {
	for (int i = 0; i < popSize; ++i) {
		CreateRandomChromosome(chromosomes[i], nrDims);
	}
}

void CreateRandomChromosome(bool chromosome[], int nrDims) {
	int i = 0, sumLgReprez = 0;
	for (int i = 0; i < nrDims; ++i) {
		for (int j = sumLgReprez; j < sumLgReprez + lgReprez[i]; ++j)
			chromosome[j] = rand() % 2;
		sumLgReprez += lgReprez[i];
	}
}

bool PopulationIsEvoluating(int lastBestGeneration, int currentGeneration) {
	if (currentGeneration - lastBestGeneration > 100)
		return false;
	return true;
}

void SelectNextPopulation(bool chromosomes[][DMAX * sizeof(int)], int &popSize, int nrDims) {
	double randomNr;
	bool chosen[DMAX];
	memset(chosen, 0, sizeof(chosen));

	randomNr = RandomDouble(0, 1);
	if (randomNr <= prob[0]) chosen[0] = true;

	for (int i = 1; i < popSize; ++i) {
		randomNr = RandomDouble(0, 1);
		if (probCumulated[i - 1] < randomNr && randomNr <= probCumulated[i])
			chosen[i] = true;
	}

	int i, aux = popSize;
	popSize = 0;
	/*for (i = 0; i < aux; ++i)
	if (chosen[i])
	Copy (chromosomes[)*/
}

void MutateChromosomes() {

}

void CrossChromosomes() {

}

double DeJong(double v[], int nrVals) {
	double functionValue = 0;

	for (int i = 0; i < nrVals; ++i) {
		functionValue += v[i] * v[i];
		if (v[i] < -5.12 || v[i] > 5.12)
		{
			cout << "Bad input for DeJong!\n";
			return 0;
		}
	}

	return functionValue;
}

double SixHump(double v[], int nrVals) {
	double functionValue = 0;

	if (v[0] < -3 || v[0] > 3 || v[1] < -2 || v[1] > 2)
	{
		cout << "Bad input for SixHump!\n";
		return 0;
	}

	functionValue = (4 - 2.1*v[0] * v[0] + v[0] * v[0] * v[0] * v[0] / 3)*v[0] * v[0] + v[0] * v[1] + (-4 + 4 * v[1] * v[1])*v[1] * v[1];

	return functionValue;
}

double Schwefel7(double v[], int nrVals) {
	double result = 0;
	for (int i = 0; i < nrVals; ++i)
		if (v[i] < -500 || v[i] > 500) {
			cout << "Bad input for Schwefel 7!\n";
			return 0;
		}

	for (int i = 0; i < nrVals; ++i)
		result += -v[i] * sin(sqrt(abs(v[i])));

	return result;
}

double Rastrigin(double v[], int nrVals) {
	double result = 0;

	for (int i = 0; i < nrVals; ++i)
		if (v[i] < -5.12 || v[i] > 5.12) {
			cout << "Bad input for Rastrigin!\n";
			return 0;
		}

	result = 10 * nrVals;
	for (int i = 0; i < nrVals; ++i)
		result += v[i] * v[i] - 10 * cos(2 * PI*v[i]);

	return result;
}

int BinaryToInt(bool binary[], int lg) {
	int rez, p2, i;

	for (rez = 0, p2 = 1, i = 0; i < lg; ++i, p2 *= 2)
		rez += p2*binary[i];

	return rez;
}

void Copy(double unde[], double deUnde[], int cate) {
	for (int i = 0; i < cate; ++i)
		unde[i] = deUnde[i];
}

void Copy(bool unde[], bool deUnde[], int cate) {
	for (int i = 0; i < cate; ++i)
		unde[i] = deUnde[i];
}

double RandomDouble(const double &minVal, const double &maxVal) {
	double f;
	f = (double)rand() / RAND_MAX;
	return minVal + f * (maxVal - minVal);
}