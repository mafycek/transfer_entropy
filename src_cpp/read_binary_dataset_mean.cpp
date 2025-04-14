//
// Created by hynek on 9.4.25.
//

#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <string>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <gtest/gtest.h>

#include "renyi_entropy.h"

int main(int argc, char *argv[])
{
    std::map<std::tuple<unsigned int, double>, double> result;
    boost::filesystem::path myFile =
        boost::filesystem::current_path() / "myfile.dat";
    boost::filesystem::ifstream ifs(myFile);
    boost::archive::binary_iarchive iarch(ifs);
    iarch >> result;
    std::cout << result.size() << std::endl;

    for (auto item : result) {
        std::cout << std::get<0>(item.first) << " " << std::get<1>(item.first)
                  << " " << item.second << std::endl;
    }
}
