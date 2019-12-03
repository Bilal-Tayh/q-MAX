#ifndef QMAX_H
#define QMAX_H
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <random>

class QMax
{
	int *_A;
	int _curIdx;
	int _q;
	int _actualsize;
	int _actualsizeMinusOne;
	int _qMinusOne;
	float _gamma;
	int _nminusq;
	int _phi;
    std::default_random_engine _generator;
	void maintenance();
	int PartitionAroundPivot(int left, int right, int pivot_idx, int* nums);
public:
	void reset();
	int findKthLargestAndPivot();
	QMax(int q, float gamma);
	void insert(int id);
    int GenerateRandom(int min,int max);
	int* largestQ();
	void print();
    int checkPivot(int value, double psi);
    int findValueIndex(int valuex);
};
#endif

