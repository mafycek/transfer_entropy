
#include "msgpack.hpp"

#include <gtest/gtest.h>


TEST ( MSGPACK, Test )
{
	std::string output {"CRE.bin"};
	boost::filesystem::path output_file =
	boost::filesystem::path(directory) / boost::filesystem::path(output);
	boost::filesystem::ofstream output_file_handler ( output_file );

	msgpack::object_handle oh =
        msgpack::unpack(str.data(), str.size());

    // deserialized object is valid during the msgpack::object_handle instance is alive.
    msgpack::object deserialized = oh.get();

}
