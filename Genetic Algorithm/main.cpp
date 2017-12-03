#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#define DMAX 5000
#define P_MUTATION 0.01
#define P_CROSS 0.25
#define INF 2000000000
#define PI 3.1415
#define EPSILON 0.1
#define C 10000.0

using namespace std;

struct pack {
	double fitness, functionValue;
};

double prob[DMAX], probCumulated[DMAX];
int lgReprez[DMAX];
int discreteFactor;

double DeJong(double[], int);
double SixHump(double[], int);
double Schwefel7(double[], int);
double Rastrigin(double[], int);

void BuildLgReprez(int, double[2], int);
double GeneticAlgorithm(double(*)(double[], int), int, double[][2]);
double RandomDouble(const double &, const double &);
int BinaryToInt(bool[], int);
void Copy(double[], double[], int);
void Copy(bool[], bool[], int);
void ValsToDouble(double[], bool[], int, double[][2]);
void BuildPartialProbabilities(double, double[], int);
struct pack Fitness(double(*testFunction)(double[], int), double[], int);
void GenerateRandomSolution(bool[][DMAX * sizeof(int)], int, int, double[][2]);
void CreateRandomChromosome(bool[], int);
bool PopulationIsEvoluating(int, int);
double EvaluatePopulation(bool[][DMAX * sizeof(int)], int, int, double(*testFunction)(double[], int), double[][2]);
void SelectNextPopulation(bool[][DMAX * sizeof(int)], int &, int);
void MutateChromosomes(bool[][DMAX * sizeof(int)], int, int);

void CrossChromosomes(bool[][DMAX * sizeof(int)], int, int);
void CrossIndividualChromosomes(bool[][DMAX * sizeof(int)], int, int, int);

int main() {
	double fRez;
	double acceptedVals[DMAX][2];
	int nrDims, discreteFactor;

	srand((unsigned int)time(NULL));

	nrDims = 2; discreteFactor = 2;
	for (int i = 0; i < nrDims; ++i) {
		acceptedVals[i][0] = -5.12;
		acceptedVals[i][1] = 5.12;
	}
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
	for (lg = 0; (acceptedVals[1] - acceptedVals[0]) * pow(10, discreteFactor) > pow(2, lg);)
		++lg;

	lgReprez[pos] = lg;
}

double GeneticAlgorithm(double(*testFunction)(double[], int), int nrDims, double acceptedVals[][2]) {
	double bestSol, generationResult;
	int nrIterations, popSize, currentGeneration, lastBestGeneration;
	bool chromosomes[DMAX][DMAX * sizeof(int)];

	nrIterations = 100;
	popSize = 1000;
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
			MutateChromosomes(chromosomes, popSize, nrDims);
			CrossChromosomes(chromosomes, popSize, nrDims);

			++currentGeneration;

			//cout << "Current Best solution: " << bestSol << '\n';
		}
	}

	return bestSol;
}

double EvaluatePopulation(bool chromosomes[][DMAX * sizeof(int)], int popSize, int nrDims, double(*testFunction)(double[], int), double acceptedVals[][2]) {
	double chromosomeFitness[DMAX], doubleVals[DMAX];
	double bestSol = INF, sumFitness = 0;
	struct pack result;

	for (int i = 0; i < popSize; ++i) {
		ValsToDouble(doubleVals, chromosomes[i], nrDims, acceptedVals);

		result = Fitness(testFunction, doubleVals, nrDims);
		chromosomeFitness[i] = result.fitness;
		bestSol = min(bestSol, result.functionValue);
		sumFitness += chromosomeFitness[i];
	}

	BuildPartialProbabilities(sumFitness, chromosomeFitness, popSize);

	return bestSol;
}

struct pack Fitness(double(*testFunction)(double[], int), double doubleVals[], int nrDims) {
	struct pack rez;

	rez.functionValue = testFunction(doubleVals, nrDims);

	if (testFunction == &Schwefel7)
		rez.fitness = abs(rez.functionValue + C);
	else
		rez.fitness = (1 / (rez.functionValue + EPSILON));
	return rez;
}

void BuildPartialProbabilities(double sumFitness, double chromosomeFitness[], int popSize) {
	for (int i = 0; i < popSize; ++i)
		prob[i] = chromosomeFitness[i] / sumFitness;

	probCumulated[0] = prob[0];
	for (int i = 1; i < popSize; ++i)
		probCumulated[i] = probCumulated[i - 1] + prob[i];
}

void ValsToDouble(double doubleVals[], bool chromosome[], int nrDims, double acceptedVals[][2]) {
	int aux, sumLgReprez = 0;
	for (int i = 0; i < nrDims; ++i) {
		aux = BinaryToInt(chromosome + sumLgReprez, lgReprez[i]);
		doubleVals[i] = acceptedVals[i][0] + aux * (acceptedVals[i][1] - acceptedVals[i][0]) / (pow(2, lgReprez[i]) - 1);

		sumLgReprez += lgReprez[i];
	}
}

void GenerateRandomSolution(bool chromosomes[][DMAX * sizeof(int)], int popSize, int nrDims, double acceptedVals[][2]) {
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
	for (i = 0; i < aux; ++i)
		if (chosen[i])
			Copy(chromosomes[popSize++], chromosomes[i], DMAX * sizeof(int));
}

void MutateChromosomes(bool chromosomes[][DMAX * sizeof(int)], int popSize, int nrDims) {
	int sumLgReprez = 0;
	double randomNr;

	for (int i = 0; i < nrDims; ++i)
		sumLgReprez += lgReprez[i];

	for (int i = 0; i < popSize; ++i)
		for (int j = 0; j < sumLgReprez; ++j) {
			randomNr = RandomDouble(0, 1);
			if (randomNr <= P_MUTATION)
				chromosomes[i][j] = 1 - chromosomes[i][j];
		}
}

void CrossChromosomes(bool chromosomes[][DMAX * sizeof(int)], int popSize, int nrDims) {
	vector <int> chosenForCrossing;
	double randomNr;

	for (int i = 0; i < popSize; ++i) {
		randomNr = RandomDouble(0, 1);
		if (randomNr <= P_CROSS)
			chosenForCrossing.push_back(i);
	}

	if (chosenForCrossing.size() % 2)
		chosenForCrossing.pop_back();

	for (unsigned int i = 0; i < chosenForCrossing.size(); i += 2) {
		CrossIndividualChromosomes(chromosomes, nrDims, chosenForCrossing[i], chosenForCrossing[i + 1]);
	}

}

void CrossIndividualChromosomes(bool chromosomes[][DMAX * sizeof(int)], int nrDims, int pos1, int pos2) {
	int randomPoz, sumLgReprez = 0;
	for (int i = 0; i < nrDims; ++i)
		sumLgReprez += lgReprez[i];

	randomPoz = rand() % sumLgReprez;
	for (int j = randomPoz; j < sumLgReprez; ++j)
		swap(chromosomes[pos1][j], chromosomes[pos2][j]);
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