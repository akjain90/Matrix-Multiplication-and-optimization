#ifndef MATRIXMULT_H
#define MATRIXMULT_H
#include "Timer.h"
#include <iostream>
#include <fstream>
#include <vector>
extern "C" {
#include <cblas.h>
}
#ifdef USE_LIKWID
extern "C" {
#include <likwid.h>
}
#endif 

using namespace std;
class matrixmult
{
	public:
		matrixmult() {};
		void matmult_naive(int,char*[]);
		void matmult_padding(int ,char*[]);
		void matmult_block(int ,char*[]);
		void matmult_optimized_ic(int ,char*[]);
		void matmult_cblas(int ,char*[]);
};
#endif
