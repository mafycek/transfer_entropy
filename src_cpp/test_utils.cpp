
#include <string>
#include <gtest/gtest.h>
#include <iostream>

#include "utils.h"


TEST(Utils, ArrayExtraction)
{
	std::string test_string{"1.2 2 3 4 5 1e220"};
	std::vector<double> values;
	
	convert_array(test_string, values, std::stod);
	for (auto item: values)
	{
		std::cout << item << std::endl;
	}
}
