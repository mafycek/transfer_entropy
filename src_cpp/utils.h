//
// Created by hynek on 5.3.25.
//
#pragma once

#include <string>
#include <string_view>
#include <vector>

template<typename S, typename T>
void convert_array(const S &array, std::vector<T>& values, T (* convertor) (const std::string &, std::size_t*) )
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

template<typename S, typename T>
void convert_array(const S &array, std::vector<T>& values, std::function<T (const std::string &, std::size_t*)> convertor )
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

template<typename S, typename T>
void chop_string_arrays (const S &array, char chop_separator, std::vector<std::vector<T>>& values, T (* convertor) (const std::string &, std::size_t*) )
{
	size_t position{0};
	size_t location_separator{0};
	S data;
	for(; location_separator != array.npos ;)
	{
		data = array.substr( position );
		location_separator = data.find(chop_separator);
		size_t position_of_separator;
		if ( location_separator == data.npos )
		{
			position_of_separator = data.size();
		}
		else
		{
			position_of_separator = location_separator;
		}
		auto chop_substring = data.substr(0, position_of_separator);
		values.push_back(std::vector<T>());
		convert_array<std::string_view, T>(chop_substring, values.back(), convertor );
		position += position_of_separator + 1;
	}
}

template<typename S, typename T>
void chop_string_arrays (const S &array, char chop_separator, std::vector<std::vector<T>>& values, std::function<T (const std::string &, std::size_t*)> convertor )
{
	size_t position{0};
	size_t location_separator{0};
	S data;
	for(; location_separator != array.npos ;)
	{
		data = array.substr( position );
		location_separator = data.find(chop_separator);
		size_t position_of_separator;
		if ( location_separator == data.npos )
		{
			position_of_separator = data.size();
		}
		else
		{
			position_of_separator = location_separator;
		}
		auto chop_substring = data.substr(0, position_of_separator);
		values.push_back(std::vector<T>());
		convert_array<std::string_view, T>(chop_substring, values.back(), convertor );
		position += position_of_separator + 1;
	}
}

