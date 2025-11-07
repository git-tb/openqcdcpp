#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <fstream>
#include <cmath>

int getIntervalIdx(double& x, std::vector<double> xsamples)	{
	auto it		= std::find_if(xsamples.begin(), xsamples.end(), [x](double x_){return x<x_;});
	int idx	= std::distance(xsamples.begin(), it)-1;
	return idx;
}

int getIntervalIdxSmart(double& x, const std::vector<double> &xsamples, uint loweridx, uint upperidx)	{
	if(loweridx >= upperidx)	{
		std::cerr << "ERROR: loweridx >= upperidx, not sensible!" << std::endl;
		abort();
	} else if(loweridx + 1 == upperidx)	{
		if((x>=xsamples[loweridx]) && (x<xsamples[loweridx+1]))	{
			return loweridx;
		} else {
			std::cerr << "ERROR: requested x (" << x << ") is not in the intervall [" << xsamples[loweridx] << ", " << xsamples[loweridx+1] << ") of size " << upperidx-loweridx << "!" << std::endl;
			abort();
		}
	} else if(loweridx + 2 == upperidx)	{
		if((x>=xsamples[loweridx]) && (x<xsamples[loweridx+1]))	{
			return loweridx;
		} else if((x>=xsamples[loweridx+1]) && (x<xsamples[loweridx+2]))	{
			return loweridx+1;
		} else {
			std::cerr << "ERROR: requested x (" << x << ") is not in the intervall [" << xsamples[loweridx] << ", " << xsamples[loweridx+2] << ") of size " << upperidx-loweridx << "!" << std::endl;
			abort();
		}
	} else	{
		int mididx = (upperidx+loweridx)/2;
		if(x >= xsamples[mididx]) 	return getIntervalIdxSmart(x, xsamples, mididx, upperidx);
		else						return getIntervalIdxSmart(x, xsamples, loweridx, mididx);
	}
}

int main()	{

	double 	xmin	= 0.0,
			xmax 	= 10.0;

	std::ofstream fileout("test.dat");
	fileout << "listlen;time algo1;time algo2" << std::endl;

	uint Nsizes		= 50;
	uint SizeMin	= 10;
	uint SizeMax	= 1000;

	double logsizemin = std::log((double)SizeMin);
	double logsizemax = std::log((double)SizeMax);

	std::vector<double> sizes(Nsizes);
	for(uint i = 0; i < Nsizes; i++)	{
		double logsize	= logsizemin + (logsizemax-logsizemin)*(double)i/(double)(Nsizes-1);
		uint size 		= (uint)std::ceil(std::exp(logsize));
		sizes[i]		= size;
	}

	for(uint k = 0; k < sizes.size(); k++)	{
		{
			std::cout << "\33[2K\r" << std::flush;
			std::cout << "[\33[32m";
			for(int j = 0; j < 100; j++)	{
				if((int)(100*(double)(k)/(double)(sizes.size()-1)) >= j)	{
					std::cout << "\u2589";
				} else	{
					std::cout << "\u2591";
				}
			}
			std::cout << "\33[0m]" << std::flush;
		}


		uint Ntest = sizes[k];
		std::vector<double> xsamples(Ntest);
		for(uint i = 0; i < Ntest; i++)	{
			xsamples[i] = xmin + (xmax - xmin)*(double)i/(double)(Ntest-1);
		}

		double	total1	= 0.0,
				total2	= 0.0;
		for(uint i = 0; i < Ntest-1; i++)	{ 								///< run over the whole list of samples to get the average runtime			
			double x = xsamples[i] + 0.5 * (xsamples[i+1]-xsamples[i]); 	///< test the midpoints of each interval

			uint Nstat = 1000; ///< number of runs to average out fluctuations
			for(uint j = 0; j < Nstat; j++)	{
				uint idx1, idx2;
				{
					auto start		= std::chrono::high_resolution_clock::now();
					idx1 = getIntervalIdx(x, xsamples);
					auto end		= std::chrono::high_resolution_clock::now();
					auto duration	= std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
					total1			+= (double)duration.count()/(double)(Nstat*Ntest);
				}				
				{
					auto start		= std::chrono::high_resolution_clock::now();
					idx2 = getIntervalIdxSmart(x, xsamples,0,xsamples.size()-1);
					auto end		= std::chrono::high_resolution_clock::now();
					auto duration	= std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
					total2			+= (double)duration.count()/(double)(Nstat*Ntest);
				}	
				if(idx1 != idx2)	{
					std::cerr << "ERROR: algorithms get different results!" << std::endl;
					abort();
				}	
			}
		}

		fileout << Ntest << ";" << total1 << ";" << total2 << std::endl;
	}
	std::cout << std::endl;
	fileout.close();

	return 0;
}
