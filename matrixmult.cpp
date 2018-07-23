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
#include "matrixmult.h"

using namespace std;
						//////////*********///////
						//////////Naive///////
						/////////*********///////
void matrixmult::matmult_naive(int argc,char *argv[])
{
if(argc==4)
{
	vector<double> A;
	vector<double> B;
	int m;
	int n;
	int o;
	int p;
	double time_n=100;

	////file input for A/////

ifstream fino;
fino.open(argv[1]);
fino>>m>>n;
double zz=0;
for(int j=0;j<m*n;j++)
{
	fino>>zz;
	A.push_back(zz);
}
fino.close();

////file input for B/////

fino.open(argv[2]);
fino>>o>>p;
double y=0;
for(int j=0;j<o*p;j++)
{

		fino>>y;
		B.push_back(y);
}
fino.close();

vector<double> C(m*p);
siwir::Timer timer;
cout<<"Performance measure of Naive matmult"<< endl;
#ifdef USE_LIKWID
   likwid_markerInit();
   likwid_markerStartRegion( "naive" );
#endif

    timer.reset();
for(int i=0;i<m;i++)
	{
		for(int j=0;j<p;j++)
		{
			for(int k=0;k<n;k++)
			{
				C[j+(i*p)]+=A[k+(i*n)]*B[(p*k)+j];
			}
		}
	}
#ifdef USE_LIKWID
   likwid_markerStopRegion( "naive" );
   likwid_markerClose();
#endif

time_n = std::min(time_n, timer.elapsed());
   
cout << "Calculation for naive matmult took " << time_n << " seconds\n";
ofstream fout;
fout.open(argv[3]);
fout<<m<<" "<<p<<endl;
for(int i=0;i<m*p;i++)
	{
		fout<<C[i]<<endl;
	}
fout.close();
A.clear();
B.clear();
C.clear();
}
}
										/////////*********///////
										//////////padding///////
										/////////*********///////

void matrixmult::matmult_padding(int argc,char *argv[])
{
if (argc==4)
{
vector<double> A;
vector<double> B;
int m;
int n;
int o;
int p;
int padd_a = 0;
int padd_b = 0;
double time=100;
ifstream fin;
fin.open(argv[1]);
fin>>m>>n;
double z=0;
/////File input and Padding for Matrix A/////
if (n%512 ==0)
padd_a = 16;
else
	padd_a=0;
for(int j=0;j<m;j++)
{
	for(int i=0;i<n;i++)
	{
		fin>>z;
		A.push_back(z);
	}
	for(int k=0;k<padd_a;k++)
		A.push_back(0);
}
fin.close();

fin.open(argv[2]);
fin>>o>>p;
double y=0;
/////File input and Padding for Matrix B/////
if (p%512 ==0)
padd_b = 16;
else
	padd_b=0;
for(int j=0;j<o;j++)
{
	for(int i=0;i<p;i++)
	{
		fin>>y;
		B.push_back(y);
	}
	for(int k=0;k<padd_b;k++)
		B.push_back(0);
	}
fin.close();

vector<double> C(m*p);
cout<<"Performance of matmult with padding"<< endl;
#ifdef USE_LIKWID
    likwid_markerInit();
    likwid_markerStartRegion( "Padding" );
#endif
siwir::Timer timer;
timer.reset();

if(n%512 ==0)
{
	for(int i=0;i<m;i++)
		{
			for(int j=0;j<p;j++)
			{
				for(int k=0;k<n;k++)
				{
					C[j+i*p]+=A[k+i*(n+padd_a)]*B[(p+padd_b)*k+j];
				}
			
			}
			
		}
}
else
{	
	for(double i=0;i<m;i++)
	{
		for(double j=0;j<p;j++)
		{
			for(double k=0;k<n;k++)
			{
				C[j+i*p]+=A[k+i*n]*B[p*k+j];
			}
		}
	}
}
#ifdef USE_LIKWID
   likwid_markerStopRegion( "Padding" );
   likwid_markerClose();
#endif 

 time = std::min(time, timer.elapsed());
 
 cout << "Calculation for Matmult with Padding took " << time << " seconds\n"; 

ofstream fout;
fout.open(argv[3]);
fout<<m<<" "<<p<<endl;
for(int i=0;i<m*p;i++)
	{
		fout<<C[i]<<endl;

	}
fout.close();
}
}
					/////////*********///////
					//////////blocking///////
					/////////*********///////

void matrixmult::matmult_block(int argc,char *argv[])
{
if (argc==4)
{
vector<double> A;
vector<double> B;
int m;
int n;
int o;
int p;

double time=100;
    
////file input for A/////

ifstream fino;
fino.open(argv[1]);
fino>>m>>n;
double zz=0;

for(int j=0;j<m*n;j++)
{
    fino>>zz;
	A.push_back(zz);
}
fino.close();

////file input for B/////

fino.open(argv[2]);
fino>>o>>p;
double y=0;
for(int j=0;j<o*p;j++)
{
		fino>>y;
		B.push_back(y);
}

fino.close();

vector<double> C(m*p);

cout<<"Performance of blocked code"<< endl;
#ifdef USE_LIKWID
    likwid_markerInit();
    likwid_markerStartRegion( "Blocking" );
#endif
siwir::Timer timer;
timer.reset();

////Optimization with blocking/////

int cls=96;  	
	
        for(int ib=0;ib<m;ib+=cls)
        {   int is=ib;
            int ie=min(ib+cls,m);
    		
            for(int jb=0;jb<p;jb+=cls)
            {   int js=jb;
                int je=min(jb+cls,p);
                 for( int kb=0;kb<n;kb+=cls)
                {
                    int ks=kb;
                    int ke=min(kb+cls,n);
                    for(int i=is;i<ie;i++)
                    {   int i_p = i*p;
			           for(int j=js;j<je;j++)
                        {
							int j_0_i_p =   j+i_p;
                            for(int k=ks;k<ke;k++)
                            {
                                C[j_0_i_p]+=A[k+i*n]*B[p*k+j];
                            }
                       }
                    }
                }   
            }    
	    }

#ifdef USE_LIKWID
   likwid_markerStopRegion( "Blocking" );
   likwid_markerClose();
#endif

 time = std::min(time, timer.elapsed());
 
 cout << "Calculation of Matmult with blocking took " << time << " seconds\n"; 

ofstream fout;
fout.open(argv[3]);
fout<<m<<" "<<p<<endl;
for(int i=0;i<m*p;i++)
	{
		fout<<C[i]<<endl;
	}
fout.close();
A.clear();
B.clear();
C.clear();
}
}
			//////***********/////			
			////// Optimized /////
			//////***********/////

void matrixmult::matmult_optimized_ic(int argc,char *argv[])
{
if (argc==4)
{
	vector<double> A;
	vector<double> B;
	int m;
	int n;
	int o;
	int p;
	int padd_a = 0;
	int padd_b = 0;
	double time=100;
/////File input and Padding for Matrix A/////
ifstream fin;
fin.open(argv[1]);
fin>>m>>n;
double z=0;
if (n%512 ==0)
{padd_a = 16;
}
else
{	padd_a=0;
}
for(int j=0;j<m;j++)
{
	for(int i=0;i<n;i++)
	{
		fin>>z;
		A.push_back(z);
	}
	for(int k=0;k<padd_a;k++)
		A.push_back(0);
}
fin.close();
/////File input and Padding for Matrix B/////
fin.open(argv[2]);
fin>>o>>p;
double y=0;
if (p%512 ==0)
{
padd_b = 16;
}
else
{	padd_b=0;
}
for(int j=0;j<o;j++)
{
	for(int i=0;i<p;i++)
	{
		fin>>y;
		B.push_back(y);
	}
	for(int k=0;k<padd_b;k++)
		B.push_back(0);
	}
fin.close();
	int n_padd_a = n+padd_a;
	int p_padd_b = p+padd_b;
vector<double> C(m*p);
#ifdef USE_LIKWID
   likwid_markerInit();
   likwid_markerStartRegion( "Optimized" );
#endif

siwir::Timer timer;
timer.reset();

int cls=512; //cache line size// 
if(p%4 ==0)
	{
        for(int ib=0;ib<m;ib+=cls)
        {   int is=ib;
            int ie=min(ib+cls,m);
    		
           for(int jb=0;jb<p;jb+=cls)
            {   int js=jb;
                int je=min(jb+cls,p);
                for( int kb=0;kb<n;kb+=cls)
                {
                    int ks=kb;
                    int ke=min(kb+cls,n);
                    for(int i=is;i<ie;i++)	
                    {
			int i_p=i*p;
			int i_n=i*n_padd_a;
			for(int k=ks;k<ke;k+=8)	//Interchanging the loop so as to improve the temporal locality of cache
                         {	int i_k =k+i_n;
				int i_k1 =k+1+i_n;
				int i_k2 =k+2+i_n;
				int i_k3 =k+3+i_n;
				int i_k4 =k+4+i_n;
				int i_k5 =k+5+i_n;
				int i_k6 =k+6+i_n;
				int i_k7 =k+7+i_n;
			         for(int j=js;j<je;j++)
                                  {	int k_p=k*p_padd_b;
					int k1_p=(k+1)*p_padd_b;
					int k2_p=(k+2)*p_padd_b;
					int k3_p=(k+3)*p_padd_b;
					int k4_p=(k+4)*p_padd_b;
					int k5_p=(k+5)*p_padd_b;
					int k6_p=(k+6)*p_padd_b;
					int k7_p=(k+7)*p_padd_b;
					int k_j =k_p+j;
					int k1_j =k1_p+j;
					int k2_j =k2_p+j;
					int k3_j =k3_p+j;
					int k4_j =k4_p+j;
					int k5_j =k5_p+j;
					int k6_j =k6_p+j;
                                	int k7_j =k7_p+j;
					int i_j =j+i_p;
						C[i_j]+=A[i_k]*B[k_j];;
                               			C[i_j]+=A[i_k1]*B[k1_j];
                                		C[i_j]+=A[i_k2]*B[k2_j];
						C[i_j]+=A[i_k3]*B[k3_j];
						C[i_j]+=A[i_k4]*B[k4_j];
						C[i_j]+=A[i_k5]*B[k5_j];
						C[i_j]+=A[i_k6]*B[k6_j];
						C[i_j]+=A[i_k7]*B[k7_j];
                            }
                       
                        }
                    }
                   
                }
	        }
	    }
			
    }
else
	{
        for(int ib=0;ib<m;ib+=cls)
        {   int is=ib;
            int ie=min(ib+cls,m);
    
	        for(int jb=0;jb<p;jb+=cls)
            {   int js=jb;
                int je=min(jb+cls,p);
                for( int kb=0;kb<n;kb+=cls)
                {
                    int ks=kb;
                    int ke=min(kb+cls,n);
                    for(int i=is;i<ie;i++)
                    {
                        for(int k=ks;k<ke;k++)
                        {                                             
                            for(int j=js;j<je;j++)
                            {
                                C[j+(i*p_padd_b)]+=A[k+(i*n_padd_a)]*B[(p_padd_b*k)+j];
                            }
                       
                        }
                    }
                   
                }
	        }   
	    }	
	}
#ifdef USE_LIKWID
   likwid_markerStopRegion( "Optimized" );
 likwid_markerClose();
#endif
 time = std::min(time, timer.elapsed());
   
cout << "Calculation of optimized case took " << time << " seconds\n";

ofstream fout;
fout.open(argv[3]);
fout<<m<<" "<<p<<endl;
for(int i=0;i<m*p;i++)
	{
		fout<<C[i]<<endl;

	}
fout.close();
A.clear();
B.clear();
C.clear();
}	
}
							/////////******/////////
							//////////BLAS//////////
							/////////*********//////

void matrixmult::matmult_cblas(int argc,char *argv[])
{
if (argc==4)
{
vector<double> A;
vector<double> B;
int m;
int n;
int o;
int p;
double time=100;
ifstream fin;
fin.open(argv[1]);
fin>>m>>n;
double z=0;

for(int j=0;j<m*n;j++)
{
	
		fin>>z;
	A.push_back(z);
}
fin.close(); 

fin.open(argv[2]);
fin>>o>>p;
double y=0;
for(int j=0;j<o*p;j++)
{

		fin>>y;
		B.push_back(y);
}
fin.close();

const int lda = n;
const int ldb = p;
const int ldc = p;
vector<double> C(m*p);
siwir::Timer timer;
cout<<"Performance of Blas"<< endl;
#ifdef USE_LIKWID
   likwid_markerInit();
   likwid_markerStartRegion( "blas" );
#endif

      timer.reset();
      cblas_dgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, m, p, n, 1.0, A.data(), lda, B.data(), ldb, 0.0, C.data(), ldc );
      time = std::min(time, timer.elapsed());
      cout << "Calculation for blas took " << time << " seconds\n"; 

#ifdef USE_LIKWID
   likwid_markerStopRegion( "blas" );
   likwid_markerClose();
#endif
ofstream fout;
fout.open(argv[3]);
fout<<m<<" "<<p<<endl;
for(int i=0;i<m*p;i++)
        {
                fout<<C[i]<<endl;
        }
fout.close();
A.clear();
B.clear();
C.clear();
}
}
