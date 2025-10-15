#ifndef FLAGS_H
#define FLAGS_H

/**
 * @brief Define flags that can be set in the program. Each function can be passed a flag variable
 * (unsigned int) where each bit represents a certain option in the program that can be either on
 * or of. Flags can be combined via bitwise OR, i.e.
 * 
 * 	flags = FLAG1 | FLAG2 | ...
 * 
 * and be read with bitwise AND, i.e.
 * 
 * 	if(flags & FLAG1) {...}
 * 
 * This makes the code more readable and avoid long argument lists in the functional programming
 * approach, which is especially relevant in the case of default arguments
 * 
 * 	somefunction(double x, 
 * 		double Q2, 
 * 		bool option1 = false,
 * 		bool option2 = true,
 * 		bool option3 = true,
 * 		bool option4 = false) 
 * 	->
 * 	somefunction(double x, double Q2, FLAG2 | FLAG3)
 */

typedef unsigned int FLAG;

namespace FLAGS	{
	extern FLAG EXACT;
	extern FLAG FLAG2;
	extern FLAG FLAG3;
}

#endif