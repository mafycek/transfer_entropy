//
// Created by hynek on 5.3.25.
//
#pragma once

#include <string>
#include <vector>

template<typename T>
void convert_array(const std::string &array, std::vector<T>& values, T (* convertor) (const std::string &, std::size_t*) )
{
	try
	{
		std::string_view data{array};
		std::string::size_type position{0};
		std::string::size_type sz{0};
		for(; ;)
		{
			//double value = std::stod( std::string(data.substr(position)), & sz);
			const T value = convertor( std::string(data.substr(position)), & sz);
			position += sz;
			values.push_back(value);
		}
	}
	catch(std::invalid_argument & exp)
	{
		// conversion was unsucceful
	}
}

template<typename T>
void convert_array(const std::string &array, std::vector<T>& values, std::function<T (const std::string &, std::size_t*)> convertor )
{
	try
	{
		std::string_view data{array};
		std::string::size_type position{0};
		std::string::size_type sz{0};
		for(; ;)
		{
			//double value = std::stod( std::string(data.substr(position)), & sz);
			const T value = convertor( std::string(data.substr(position)), & sz);
			position += sz;
			values.push_back(value);
		}
	}
	catch(std::invalid_argument & exp)
	{
		// conversion was unsucceful
	}
}

