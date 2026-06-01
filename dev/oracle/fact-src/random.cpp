#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include "random.h"
#include "stdlib.h"
#include <iostream>
#include <math.h>

using namespace std;

Random Random::RANDOM((unsigned) 0);

Random::Random(double seed) :
	seed_((unsigned) (RAND_MAX * seed)) {
}
Random::Random(unsigned seed) :
	seed_(seed) {
}

unsigned Random::seed() {
	return seed_;
}

unsigned Random::NextUnsigned() {
	return rand();
}

int Random::NextInt(int n) {
	// return (int) (n * ((double) rand() / RAND_MAX));
	return rand() % n;
}

int Random::NextInt(int min, int max) {
	return rand() % (max - min) + min;
}

double Random::NextDouble(double min, double max) {
	return (double) rand() / RAND_MAX * (max - min) + min;
}

double Random::NextDouble() {
	return (double) rand() / RAND_MAX;
}

double Random::NextGaussian() {
	double u = NextDouble(), v = NextDouble();
	return sqrt(-2 * log(u)) * cos(2 * M_PI * v);
}

int Random::NextCategory(const vector<double>& category_probs) {
	return GetCategory(category_probs, NextDouble());
}

int Random::GetCategory(const vector<double>& category_probs, double rand_num) {
	int c = 0;
	double sum = category_probs[0];
	while (sum < rand_num) {
		c++;
		sum += category_probs[c];
	}
	return c;
}
