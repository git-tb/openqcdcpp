#include <iostream>
#include <functional>
#include <iomanip>

/// @brief emulate behaviour of a normal function f(x)
class MyFunction	{
	public:
		MyFunction(std::function<double(double)> func_) :
			func(func_)	{
		};

		MyFunction operator*(MyFunction other)	{
			return MyFunction([this,other](double x){return func(x)*other.func(x);});
		}

		double operator()(double x)	{
			return func(x);
		}

		std::function<double(double)> func;
};

/// @brief emulate behaviour of a plus distribution p+(x), such that (p+ * f)(x) = p+(x) * (f(1) - f(x))
class MyPseudoPlus	:	public MyFunction	{
	public:
		MyPseudoPlus(std::function<double(double)> func_)	: MyFunction(func_)	{};

		MyPseudoPlus operator*(MyFunction other)	{
			return MyPseudoPlus([this,other](double x){return func(x)*(other.func(1) - other.func(x));});
		}
};

int main()	{

	MyFunction		myf1([](double x){return x;});
	MyFunction		myf2([](double x){return 2*x;});
	MyPseudoPlus	myplus([](double x){return 1;});

	double x0 = 0.1;
	std::cout 	<< std::setw(25)	<< "myf1(x0)"					
				<< std::setw(10)	<<	myf1(x0)						<< std::endl;
	std::cout 	<< std::setw(25)	<< "myf2(x0)"					
				<< std::setw(10)	<<	myf2(x0)						<< std::endl;
	std::cout 	<< std::setw(25)	<< "(myf1*myf2)(x0)"			
				<< std::setw(10)	<<	(myf1*myf2)(x0)					<< std::endl;
	std::cout 	<< std::setw(25)	<< "myplus(x0)*myf1(x0)"		
				<< std::setw(10)	<<	myplus(x0)*myf1(x0)				<< std::endl;
	std::cout 	<< std::setw(25)	<< "(myplus*myf1)(x0)"			
				<< std::setw(10)	<<	(myplus*myf1)(x0)				<< std::endl;
	std::cout 	<< std::setw(25)	<< "(myplus*myf2)(x0)"			
				<< std::setw(10)	<<	(myplus*myf2)(x0)				<< std::endl;
	std::cout 	<< std::setw(25)	<< "((myplus*myf1)*myf2)(x0)"	
				<< std::setw(10)	<<	((myplus*myf1)*myf2)(x0)		<< std::endl;
	std::cout 	<< std::setw(25)	<< "(myplus*(myf1*myf2))(x0)"	
				<< std::setw(10)	<<	(myplus*(myf1*myf2))(x0)		<< std::endl;

	return 0;
}