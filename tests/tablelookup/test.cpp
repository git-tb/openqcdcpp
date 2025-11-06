#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>

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
		if(x > xsamples[mididx]) return getIntervalIdxSmart(x, xsamples, mididxcd ..)
	}
}

class MyInterval	{
	public:
		MyInterval(double lowerbound_, double upperbound_)	:
			lowerbound(lowerbound_), upperbound(upperbound_)	{};
	
		bool contains(double x)	{
			return (x >= lowerbound) && (x < upperbound);
		}
	protected:
		double lowerbound, upperbound;
};

class MyIntervalTree : MyInterval	{
	public:
		MyIntervalTree(double lowerbound_, double upperbound_) :
			MyInterval(lowerbound_, upperbound_) {};
		
		void getAddress(double x, std::vector<int> &address)	{
			if(not contains(x))	{
				std::cerr << "ERROR: cannot find address of point outside my interval!" << std::endl;
				abort();
			} else if(children.empty())	{
				return;
			} else if(children.size() != 2)	{
				std::cerr << "ERROR: each MyIntervalTree must contain either 0 or 2 children!" << std::endl;
				abort();
			} else if(children[0].contains(x))	{
				address.push_back(0);
				children[0].getAddress(x, address);
			} else if (children[1].contains(x))	{
				address.push_back(1);
				children[1].getAddress(x,address);
			} else {
				std::cerr << "ERROR: If MyIntervalTree contains x, then exactly 1 of its children must also contain x!" << std::endl;
				abort();
			}
		}

		void addGridPoint(double x)	{
			if(not contains(x))	{
				std::cerr << "ERROR: cannot add point outside my interval!" << std::endl;
				abort();
			} else if ((x == lowerbound) or (x == upperbound)) {
				return;
			} else if(children.empty())	{
				children.push_back(MyIntervalTree(lowerbound, x));
				children.push_back(MyIntervalTree(x, upperbound));
			} else if(children.size() != 2)	{
				std::cerr << "ERROR: each MyIntervalTree must contain either 0 or 2 children!" << std::endl;
				abort();
			} else if(children[0].contains(x))	{
				children[0].addGridPoint(x);
			} else if (children[1].contains(x))	{
				children[1].addGridPoint(x);
			} else {
				std::cerr << "ERROR: If MyIntervalTree contains x, then exactly 1 of its children must also contain x!" << std::endl;
				abort();
			}
		}

		void print(uint indentationlevel, uint WIDTH = 10)	{
			for(uint i = 0; i < indentationlevel; i++) std::cout << "\t";
			std::cout 	<< std::setw(WIDTH) << lowerbound  
						<< " --- "
						<< std::setw(WIDTH) << upperbound
						<< std::endl;
			if(children.empty())	{
				return;
			} else if(children.size() != 2)	{
				std::cerr << "ERROR: each MyIntervalTree must contain either 0 or 2 children!" << std::endl;
				abort();
			} else {
				children[0].print(indentationlevel+1, WIDTH);
				children[1].print(indentationlevel+1, WIDTH);
			}
		}

		uint addressToIdx(std::vector<int> address, int level)	{
			if(level > address.size()-1)	{
				return 0;
			} else if(address[level] == 0)	{
				if(children.empty())	{
					return 0;
				} else {
					return children[0].addressToIdx(address, level+1);
				}
			} else {
				if (children.empty())	{
					return 1;
				} else {
					return children[0].totalSubintervals() + children[1].addressToIdx(address,level+1);
				}
			}			
		}

		uint totalSubintervals()	{
			if(children.empty())	{
				return 1;
			} else	{
				return children[0].totalSubintervals() + children[1].totalSubintervals();
			}
		}

	private:
		std::vector<MyIntervalTree> children;
};

int main()	{
	uint Nsize	= (uint)1e2;
	std::vector<double> xsamples(Nsize);
	double 	xmin	= 0.0,
			xmax	= 10.0;

	MyIntervalTree mytree(xmin, xmax);
	for(uint i = 0; i < Nsize; i++)	{
		xsamples[i]	= xmin + (xmax-xmin)*(double)i/(double)(Nsize-1);
		if((i != 0) and (i != Nsize-1)) mytree.addGridPoint(xsamples[i]);
	}

	mytree.print(0);

	uint Ntest = 20;
	for(uint i = 0; i < Ntest; i++)	{
		double r	= (double)std::rand()/(double)RAND_MAX;
		double x	= xmin + r*(xmax-xmin);

		std::vector<int> address;
		mytree.getAddress(x, address);

		std::cout 	<< std::setw(10) << x
					<< std::setw(4) << mytree.addressToIdx(address, 0)
					<< std::setw(4) << getIntervalIdx(x, xsamples)
					<< std::setw(10) << (bool)(mytree.addressToIdx(address, 0) == getIntervalIdx(x, xsamples))
					<< std::endl;
	}

	std::vector<int> address;
	mytree.getAddress(xmax,address);

	// uint Nrun	= (uint)1e5;
	// for(uint i = 0; i < Nrun; i++)	{
	// 	double r	= (double)std::rand()/(double)RAND_MAX;
	// 	double x	= xmin + r*(xmax-xmin);

	// 	auto it		= std::find_if(xsamples.begin(), xsamples.end(), [x](double x_){return x>x_;});
	// 	uint idx	= std::distance(xsamples.begin(), it);
	// }


	return 0;
}
