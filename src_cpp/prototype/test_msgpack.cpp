#include <vector>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/zstd.hpp>
#include <boost/iostreams/filter/lzma.hpp>
#include <boost/iostreams/filter/bzip2.hpp>

#include "msgpack.hpp"

#include <gtest/gtest.h>

namespace bio = boost::iostreams;

struct zstd_ostream : boost::iostreams::filtering_ostream
{
    zstd_ostream(std::ostream& os)
    {
        bio::zstd_params zstd_params ( 13 );
        bio::filtering_ostream::push( bio::zstd_compressor{zstd_params} );
        bio::filtering_ostream::push( os );
    }
};

TEST ( MSGPACK, Test )
{
	std::string directory{""};
	std::string output {"CRE.bin"};
	boost::filesystem::path output_file =
	boost::filesystem::path(directory) / boost::filesystem::path(output);
	boost::filesystem::ofstream output_file_handler ( output_file );

	//msgpack::object_handle oh = msgpack::unpack(str.data(), str.size());

    // deserialized object is valid during the msgpack::object_handle instance is alive.
	//msgpack::object deserialized = oh.get();
}

TEST ( MSGPACK, TestVector )
{
	std::string directory {"."};
	std::string file {"CRE.bin"};
	std::vector<int> result{1,2,3,4,5};
	std::stringstream ss;
	msgpack::pack(ss, result); //renyi_entropy::renyi_entropy<calculation_type>::storage_RTE(collection_result_RTE, processing_RTE)
	std::cout << ss.str().size() << std::endl;
	boost::filesystem::path output_file =
	boost::filesystem::path(directory) / boost::filesystem::path(file);
	boost::filesystem::ofstream output_file_handler ( output_file );
	zstd_ostream zstd_compression_stream{output_file_handler};

	zstd_compression_stream << ss.str();
}

TEST ( MSGPACK, TestList )
{
	std::string directory {"."};
	std::string file {"CRE.bin"};
	std::list<int> result{1,2,3,4,5};
	std::stringstream ss;
	msgpack::pack(ss, result); //renyi_entropy::renyi_entropy<calculation_type>::storage_RTE(collection_result_RTE, processing_RTE)
	std::cout << ss.str().size() << std::endl;
	boost::filesystem::path output_file =
	boost::filesystem::path(directory) / boost::filesystem::path(file);
	boost::filesystem::ofstream output_file_handler ( output_file );
	zstd_ostream zstd_compression_stream{output_file_handler};

	zstd_compression_stream << ss.str();
}
