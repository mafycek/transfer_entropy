
#include <chrono>

#include <bsoncxx/builder/basic/array.hpp>
#include <bsoncxx/builder/basic/document.hpp>
#include <bsoncxx/builder/basic/kvp.hpp>
#include <bsoncxx/types.hpp>

#include <mongocxx/client.hpp>
#include <mongocxx/instance.hpp>
#include <mongocxx/uri.hpp>

using bsoncxx::builder::basic::kvp;
using bsoncxx::builder::basic::make_array;
using bsoncxx::builder::basic::make_document;

int main() {
  try {
    // Create an instance.
    mongocxx::instance inst{};

    // Replace the connection string with your MongoDB deployment's connection
    // string.
    const auto uri =
        mongocxx::uri{"mongodb://127.0.0.1:27017/"
                      "localhost?directConnection=true&"
                      "serverSelectionTimeoutMS=2000&appName=mongosh+2.12.0"};

    // Set the version of the Stable API on the client.
    mongocxx::options::client client_options;
    const auto api = mongocxx::options::server_api{
        mongocxx::options::server_api::version::k_version_1};
    client_options.server_api_opts(api);

    // Setup the connection and get a handle on the "admin" database.
    mongocxx::client conn{uri, client_options};
    mongocxx::database db = conn["admin"];

    // Ping the database.
    const auto ping_cmd = bsoncxx::builder::basic::make_document(
        bsoncxx::builder::basic::kvp("ping", 1));
    db.run_command(ping_cmd.view());
    std::cout
        << "Pinged your deployment. You successfully connected to MongoDB!"
        << std::endl;

    // The mongocxx::instance constructor and destructor initialize and shut
    // down the driver, respectively. Therefore, a mongocxx::instance must be
    // created before using the driver and must remain alive for as long as the
    // driver is in use.
    auto db_test = conn["test"];

    // TODO: fix dates

    bsoncxx::document::value restaurant_doc = make_document(
        kvp("address",
            make_document(kvp("street", "2 Avenue"), kvp("zipcode", 10075),
                          kvp("building", "1480"),
                          kvp("coord", make_array(-73.9557413, 40.7720266)))),
        kvp("borough", "Manhattan"), kvp("cuisine", "Italian"),
        kvp("grades",
            make_array(
                make_document(kvp("date",
                                  bsoncxx::types::b_date{
                                      std::chrono::milliseconds{12323}}),
                              kvp("grade", "A"), kvp("score", 11)),
                make_document(kvp("date",
                                  bsoncxx::types::b_date{
                                      std::chrono::milliseconds{121212}}),
                              kvp("grade", "B"), kvp("score", 17)))),
        kvp("name", "Vella"), kvp("restaurant_id", "41704620"));

    // We choose to move in our document here, which transfers ownership to
    // insert_one()
    auto res = db_test["restaurants"].insert_one(std::move(restaurant_doc));

  } catch (const std::exception &e) {
    // Handle errors.
    std::cout << "Exception: " << e.what() << std::endl;
  }

  return 0;
}
