#include "matrixmult.h"
#include "matrixmult.cpp"
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
int main(int argc,char *argv[])
{
	
	siwir::Timer timer;
	cout<<"Please enter the desired operation to be performed: "<<endl;
	cout<<"a : For Naive Matrix-Matrix multiplication"<<endl;
	cout<<"b : For Matrix-Matrix multiplication with padding"<<endl;
	cout<<"c : For Matrix-Matrix multiplication with blocking"<<endl;
	cout<<"d : For Optimized Matrix-Matrix multiplication "<<endl;
	cout<<"e : For BLAS Matrix-Matrix multiplication"<<endl;
	char a;
	cin>>a;
	matrixmult mm;
	switch (a)
	{
		case 'a':
		cout<<" Naive Matrix matrix multiplication"<<endl;
		mm.matmult_naive (argc,argv);
		break;
		
		case 'b':
		cout<<"matrix matrix multiplication with padding"<<endl;
		mm.matmult_padding(argc,argv);
		break;
		
		case 'c':
		cout<<" matrix matrix multiplication with blocking"<<endl;
		mm.matmult_block(argc,argv);
		break;
		
		case 'd':
		cout<<" Optimized matrix matrix multiplication with:"<<endl<<"a)Padding"<<"	"<<"b)Loop Interchange"<<"	"<<"c) blocking"<<"	"<<"d) unrolling"<<endl;
		mm.matmult_optimized_ic(argc,argv);
		break;
		
		case 'e':
		cout<<" BLAS Matrix matrix multiplication"<<endl;
		mm.matmult_cblas(argc,argv);
		break;
	}	
	return 0;
}
