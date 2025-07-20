
#include <iostream>
#include <map>
#include <array>
#include <variant>

#include <msgpack.hpp>
//#include <matplotlibcpp.h>

#include <boost/stacktrace.hpp>
#include <boost/exception/all.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/zstd.hpp>
#include <boost/iostreams/filter/lzma.hpp>
#include <boost/iostreams/filter/bzip2.hpp>

#include <cpptrace/from_current.hpp>
#include <cpptrace/formatting.hpp>
#include <cpptrace/cpptrace.hpp>

#include "cpptrace_helper.h"
#include "renyi_entropy.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
#include "matplotlib_wrapper.h"
#pragma GCC diagnostic pop

namespace po = boost::program_options;
namespace bio = boost::iostreams;
//namespace plt = matplotlibcpp;
namespace pyw = python_wrappers;
using namespace pybind11::literals;

pyw::matplotlib_wrapper plt;

constexpr std::string left_latex_math_curly_bracket = "\\{";
constexpr std::string right_latex_math_curly_bracket = "\\}";

typedef double calculation_type;
typedef std::variant<double, long, unsigned int, std::string, std::vector<std::string>, std::array<std::string, 2>, std::array<unsigned int, 2>, std::vector<std::tuple<std::string, std::string> > > figure_parameter_t; //
typedef std::map<std::string, figure_parameter_t> parameter_t;

std::string create_label_from_index(const std::vector<unsigned int> & index, bool brackets = false)
{
	std::string result{(brackets ? left_latex_math_curly_bracket : "")};
	for (unsigned int i = 0; i < index.size(); ++i)
	{
		result += std::to_string( index[i] ) + ( i < index.size() - 1 ? "," : "" );
	}
	return (result + (brackets? right_latex_math_curly_bracket : ""));
}

template<typename TYPE>
void convert_figure_parameters (parameter_t &parameters, std::string selector, TYPE & destination)
{
	if(parameters.contains(selector))
	{
		auto &parameter = parameters[selector];
		if (std::holds_alternative<TYPE>(parameter) )
		{
			destination = std::get<TYPE>(parameter);
		}
	}
}

void simple_result_image(renyi_entropy::renyi_entropy<calculation_type>::type_average_result_conditional_information_transfer_type & processed_results, parameter_t & parameters)
{
	long number_of_columns_in_legent = 6L;
	convert_figure_parameters(parameters, "cols", number_of_columns_in_legent);
	long dpi = 600L;
	convert_figure_parameters(parameters, "dpi", dpi);
	std::array<unsigned int, 2> limits{0, 100};
	convert_figure_parameters(parameters, "limits", limits);
	std::string selector{"mean"};
	convert_figure_parameters(parameters, "selector", selector);
	std::string title{"{\\Large \\textrm{RTE\\ dependence\\ on\\ }$\\alpha$ }"};
	convert_figure_parameters(parameters, "title", title);
	std::string xlabel{"$\\alpha$"};
	convert_figure_parameters(parameters, "xlabel", xlabel);
	std::array<std::string, 2> RTE_source_names{"X", "Y"};
	convert_figure_parameters(parameters, "source_names", RTE_source_names);
	std::string output_base_name{"RTE.png"};
	convert_figure_parameters(parameters, "output_base_name", output_base_name);
	std::vector<std::string> suffices{"png"};
	convert_figure_parameters(parameters, "suffices", suffices);
	std::string color_map{"jet"};
	convert_figure_parameters(parameters, "color_map", color_map);
	std::string average_type_selector{renyi_entropy::average_runs};
	convert_figure_parameters(parameters, "average_type_selector", average_type_selector);

	const std::array<std::string, 4> RTE_type_extra_names{"effective", "nonlinear", "ballance", "reverse"};

	for (auto& suffix: suffices)
	{
		auto collection_result_RTE = processed_results[average_type_selector][selector];
		auto collection_key_views = std::views::keys(collection_result_RTE);
		std::vector<renyi_entropy::renyi_entropy<calculation_type>::collection_conditional_information_transfer_key_type> collection_keys{ collection_key_views.begin(), collection_key_views.end() };

		plt.matplotlib_rcParams["text.usetex"] = py::bool_(true);

		for ( auto &[history_key, collection_of_data_for_history] :collection_result_RTE )
		{
			//const auto history_key = (*(collection_result_RTE.begin())).first;
			const auto first_history = std::get<1>( history_key );
			const auto first_future = std::get<0>( history_key );
			const auto second_history = std::get<2>( history_key );
			const std::string first_history_label{create_label_from_index(first_history, true)}, future_first_label{create_label_from_index(first_future, true)}, second_history_label{create_label_from_index(second_history, true)};
			const auto filename_history_suffix{create_label_from_index(first_history, false) + "_" + create_label_from_index(first_future, false) + "_" + create_label_from_index(second_history, false)};

			//const auto selector_directon_type_sample = std::make_tuple(false, true, false, 0);
			//auto collection_of_data = (*(collection_result_RTE.begin())).second[ selector_directon_type_sample ];
			//auto collection_of_data = collection_of_data_for_history[ selector_directon_type_sample ];
			for (auto &[selector_directon_type_sample, collection_of_data]: collection_of_data_for_history)
			{
				auto [RTE_balance, RTE_reverse_orientation, RTE_shuffled, RTE_surrogate, RTE_sample] = selector_directon_type_sample;

				std::string &source_timeseries{RTE_source_names[0]}, &destination_timeseries{RTE_source_names[1]};
				if (RTE_reverse_orientation)
				{
					source_timeseries = RTE_source_names[1];
					destination_timeseries = RTE_source_names[0];
				}

				const auto prefix_filename = (RTE_reverse_orientation ? RTE_type_extra_names[3]:"") + (RTE_reverse_orientation && (RTE_balance||RTE_shuffled||RTE_surrogate) ? "_" : "") + (RTE_balance ? RTE_type_extra_names[2]:"") + (RTE_balance&&(RTE_shuffled||RTE_surrogate)  ? "_" : "") + (RTE_shuffled ? RTE_type_extra_names[0]:"") + (RTE_shuffled && RTE_surrogate ? "_" : "") + (RTE_surrogate ? RTE_type_extra_names[1]:"") + (RTE_reverse_orientation || RTE_balance || RTE_shuffled || RTE_surrogate ? "_" : "");
				const auto file_output_name{prefix_filename + output_base_name + "_" + filename_history_suffix + "." + suffix};
				std::cout << file_output_name << std::endl;
				const std::string RTE_type_extra = std::format("{}", std::string{"\\textrm{"} + (RTE_shuffled ? RTE_type_extra_names[0]:"") + (RTE_shuffled && RTE_surrogate ? "," : "") + (RTE_surrogate ? RTE_type_extra_names[1]:"")  + "}");
				const auto RTE_type = "{" +std::format("\\left( R, {} \\right)", RTE_type_extra) + "}";
				const auto RTE_direction = "{" + std::format("\\alpha: {}\\longrightarrow {}", source_timeseries, destination_timeseries) + "}";
				const auto ylabel = std::format("$T^{}_{}\\left({},{},{} \\right)$", RTE_type, RTE_direction, first_history_label, future_first_label, second_history_label);

				auto colormap = plt.matplotlib_colormaps [py::str(color_map)];
				plt.matplotlib_pyplot_xlabel(py::str(xlabel));
				plt.matplotlib_pyplot_ylabel(py::str(ylabel));
				plt.matplotlib_pyplot_title(py::str(title));

				std::vector<std::tuple<std::vector<double>, std::vector<double>>> raw_results(limits[1]);
				auto CRE_key_views = std::views::keys(collection_of_data);
				std::vector<renyi_entropy::renyi_entropy<calculation_type>::renyi_key_type> CRE_keys{ CRE_key_views.begin(), CRE_key_views.end() };
				for (auto & CRE_key: CRE_keys)
				{
					auto [neighbour, alpha] = CRE_key;
					const auto RTE_value = collection_of_data[CRE_key];
					std::get<0>(raw_results[neighbour]).push_back(alpha);
					std::get<1>(raw_results[neighbour]).push_back(RTE_value);
				}

				for (auto neighbour = limits[0]; neighbour < limits[1]; ++ neighbour)
				{
					const auto color_position_within_colormap = static_cast<double>(neighbour) / ( static_cast<double>((limits[1] <= 1 ? 2 : limits[1]) - 1 ));
					auto color = colormap( color_position_within_colormap );
					py::dict kwargs = py::dict("color"_a = color);
					if (average_type_selector == renyi_entropy::average_runs)
					{
						kwargs [py::str("label")] = py::str(std::format("$k={}$", neighbour));
					}
					plt.matplotlib_pyplot_plot(std::get<0>(raw_results[neighbour]), std::get<1>(raw_results[neighbour]), **kwargs);
				}
				py::dict kwargs = py::dict("loc"_a = py::str("best"), "ncols"_a = py::int_(number_of_columns_in_legent), "fontsize"_a = py::str("xx-small"));

				if (average_type_selector == renyi_entropy::average_runs)
				{
					plt.matplotlib_pyplot_legend(**kwargs);
				}
				//plt::show();
				py::dict kwargs_savefig = py::dict("dpi"_a = py::int_(dpi));
				plt.matplotlib_pyplot_savefig(file_output_name, **kwargs_savefig);
				plt.matplotlib_pyplot_close();
			}
		}
	}
}

void quantile_result_image(renyi_entropy::renyi_entropy<calculation_type>::type_average_result_conditional_information_transfer_type & processed_results, parameter_t & parameters)
{
	long number_of_columns_in_legent = 6L;
	convert_figure_parameters(parameters, "cols", number_of_columns_in_legent);
	long dpi = 600L;
	convert_figure_parameters(parameters, "dpi", dpi);
	std::array<unsigned int, 2> limits{0, 100};
	convert_figure_parameters(parameters, "limits", limits);
	std::vector<std::string> selectors{"median"};
	convert_figure_parameters(parameters, "selectors", selectors);
	std::string color_map_selector{"jet"};
	convert_figure_parameters(parameters, "color_map_selector", color_map_selector);
	std::vector<std::tuple<std::string, std::string>> selector_fills{{"quantile_06", "quantile_04"}};
	convert_figure_parameters(parameters, "selector_fills", selector_fills);
	std::string title{"{\\Large \\textrm{RTE\\ dependence\\ on\\ }$\\alpha$ }"};
	convert_figure_parameters(parameters, "title", title);
	std::string xlabel{"$\\alpha$"};
	convert_figure_parameters(parameters, "xlabel", xlabel);
	std::array<std::string, 2> RTE_source_names{"X", "Y"};
	convert_figure_parameters(parameters, "source_names", RTE_source_names);
	std::string output_base_name{"RTE.png"};
	convert_figure_parameters(parameters, "output_base_name", output_base_name);
	std::vector<std::string> suffices{"png"};
	convert_figure_parameters(parameters, "suffices", suffices);
	std::string color_map{"jet"};
	convert_figure_parameters(parameters, "color_map", color_map);
	std::string average_type_selector{renyi_entropy::average_runs};
	convert_figure_parameters(parameters, "average_type_selector", average_type_selector);
	double linewidth{1};
	convert_figure_parameters(parameters, "linewidth", linewidth);
	double alpha{1};
	convert_figure_parameters(parameters, "alpha", alpha);

	const std::array<std::string, 4> RTE_type_extra_names{"effective", "nonlinear", "ballance", "reverse"};

	for (const auto& suffix: suffices)
	{
		auto collection_result_RTE = processed_results[average_type_selector][std::get<0>(selector_fills[0])];
		auto collection_key_views = std::views::keys(collection_result_RTE);
		std::vector<renyi_entropy::renyi_entropy<calculation_type>::collection_conditional_information_transfer_key_type> collection_keys{ collection_key_views.begin(), collection_key_views.end() };

		plt.matplotlib_rcParams["text.usetex"] = py::bool_(true);

		for ( auto &[history_key, collection_of_data_for_history] :collection_result_RTE )
		{
			const auto first_history = std::get<1>( history_key );
			const auto first_future = std::get<0>( history_key );
			const auto second_history = std::get<2>( history_key );
			const std::string first_history_label{create_label_from_index(first_history, true)}, future_first_label{create_label_from_index(first_future, true)}, second_history_label{create_label_from_index(second_history, true)};
			const auto filename_history_suffix{create_label_from_index(first_history, false) + "_" + create_label_from_index(first_future, false) + "_" + create_label_from_index(second_history, false)};

			for (const auto &[selector_directon_type_sample, collection_of_data]: collection_of_data_for_history)
			{
				auto [RTE_balance, RTE_reverse_orientation, RTE_shuffled, RTE_surrogate, RTE_sample] = selector_directon_type_sample;

				std::string &source_timeseries{RTE_source_names[0]}, &destination_timeseries{RTE_source_names[1]};
				if (RTE_reverse_orientation)
				{
					source_timeseries = RTE_source_names[1];
					destination_timeseries = RTE_source_names[0];
				}

				const auto prefix_filename = (RTE_reverse_orientation ? RTE_type_extra_names[3]:"") + (RTE_reverse_orientation && (RTE_balance||RTE_shuffled||RTE_surrogate) ? "_" : "") + (RTE_balance ? RTE_type_extra_names[2]:"") + (RTE_balance&&(RTE_shuffled||RTE_surrogate)  ? "_" : "") + (RTE_shuffled ? RTE_type_extra_names[0]:"") + (RTE_shuffled && RTE_surrogate ? "_" : "") + (RTE_surrogate ? RTE_type_extra_names[1]:"") + (RTE_reverse_orientation || RTE_balance || RTE_shuffled || RTE_surrogate ? "_" : "");
				const auto file_output_name{prefix_filename + output_base_name + "_" + filename_history_suffix + "." + suffix};
				std::cout << file_output_name << std::endl;
				const std::string RTE_type_extra = std::format("{}", std::string{"\\textrm{"} + (RTE_shuffled ? RTE_type_extra_names[0]:"") + (RTE_shuffled && RTE_surrogate ? "," : "") + (RTE_surrogate ? RTE_type_extra_names[1]:"")  + "}");
				const auto RTE_type = "{" +std::format("\\left( R, {} \\right)", RTE_type_extra) + "}";
				const auto RTE_direction = "{" + std::format("\\alpha: {}\\longrightarrow {}", source_timeseries, destination_timeseries) + "}";
				const auto ylabel = std::format("$T^{}_{}\\left({},{},{} \\right)$", RTE_type, RTE_direction, first_history_label, future_first_label, second_history_label);

				auto fill_colormap = plt.matplotlib_colormaps [py::str(color_map)];
				auto selector_colormap = plt.matplotlib_colormaps [py::str(color_map_selector)];
				plt.matplotlib_pyplot_xlabel(py::str(xlabel));
				plt.matplotlib_pyplot_ylabel(py::str(ylabel));
				plt.matplotlib_pyplot_title(py::str(title));

				for (const auto &[raw_color_index, selector_fill] : std::views::enumerate(selector_fills))
				{
					std::vector<std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>> raw_results(1);
					auto upper_datasset = processed_results[average_type_selector][std::get<0>(selector_fill)][history_key][selector_directon_type_sample];
					auto lower_datasset = processed_results[average_type_selector][std::get<1>(selector_fill)][history_key][selector_directon_type_sample];
					auto CRE_key_views = std::views::keys(upper_datasset);
					std::vector<renyi_entropy::renyi_entropy<calculation_type>::renyi_key_type> CRE_keys{ CRE_key_views.begin(), CRE_key_views.end() };
					for (auto & CRE_key: CRE_keys)
					{
						auto [neighbour, alpha] = CRE_key;
						const auto upper_RTE_value = upper_datasset[CRE_key];
						const auto lower_RTE_value = lower_datasset[CRE_key];
						std::get<0>(raw_results[neighbour]).push_back(alpha);
						std::get<1>(raw_results[neighbour]).push_back(upper_RTE_value);
						std::get<2>(raw_results[neighbour]).push_back(lower_RTE_value);
					}

					for (auto neighbour = limits[0]; neighbour < limits[1]; ++ neighbour)
					{
						const auto color_position_within_colormap = static_cast<double>(raw_color_index) / ( static_cast<double>(selector_fills.size() - 1 ));
						auto color = fill_colormap( py::float_(color_position_within_colormap ));
						py::dict kwargs ("color"_a = color, "linewidth"_a = py::float_(linewidth), "alpha"_a = py::float_(alpha));
						if (average_type_selector == renyi_entropy::average_runs)
						{
							kwargs[py::str("label")] = py::str( std::format("$k={}$", neighbour) );
						}
						plt.matplotlib_pyplot_fill_between(std::get<0>(raw_results[neighbour]), std::get<1>(raw_results[neighbour]), std::get<2>(raw_results[neighbour]), **kwargs);
					}
				}

				for (const auto &[raw_color_index, selector] : std::views::enumerate(selectors))
				{
					std::vector<std::tuple<std::vector<double>, std::vector<double>>> raw_results(1);
					auto dataset = processed_results[average_type_selector][selector][history_key][selector_directon_type_sample];
					auto CRE_key_views = std::views::keys(dataset);
					std::vector<renyi_entropy::renyi_entropy<calculation_type>::renyi_key_type> CRE_keys{ CRE_key_views.begin(), CRE_key_views.end() };
					for (auto & CRE_key: CRE_keys)
					{
						auto [neighbour, alpha] = CRE_key;
						const auto data_RTE_value = dataset[CRE_key];
						std::get<0>(raw_results[neighbour]).push_back(alpha);
						std::get<1>(raw_results[neighbour]).push_back(data_RTE_value);
					}

					for (auto neighbour = limits[0]; neighbour < limits[1]; ++ neighbour)
					{
						const auto color_position_within_colormap = static_cast<double>(raw_color_index) / ( static_cast<double>(selectors.size() - 1 ));
						const auto color = selector_colormap( color_position_within_colormap );
						py::dict kwargs ("color"_a = color, "linewidth"_a = py::float_(linewidth), "alpha"_a = py::float_(alpha));
						if (average_type_selector == renyi_entropy::average_runs)
						{
							kwargs ["label"] = py::str (std::format("$k={}$", neighbour));
						}
						plt.matplotlib_pyplot_plot(std::get<0>(raw_results[neighbour]), std::get<1>(raw_results[neighbour]), **kwargs);
					}
				}
				py::dict kwargs_savefig = py::dict("dpi"_a = py::int_(dpi));
				plt.matplotlib_pyplot_savefig(file_output_name, **kwargs_savefig);
				plt.matplotlib_pyplot_close();
			}
		}
	}
}

void deviation_result_image(renyi_entropy::renyi_entropy<calculation_type>::type_average_result_conditional_information_transfer_type & processed_results, parameter_t & parameters)
{
	long number_of_columns_in_legent = 6L;
	convert_figure_parameters(parameters, "cols", number_of_columns_in_legent);
	long dpi = 600L;
	convert_figure_parameters(parameters, "dpi", dpi);
	std::array<unsigned int, 2> limits{0, 100};
	convert_figure_parameters(parameters, "limits", limits);
	std::vector<std::string> selectors{"mean"};
	convert_figure_parameters(parameters, "selectors", selectors);
	std::string color_map_selector{"jet"};
	convert_figure_parameters(parameters, "color_map_selector", color_map_selector);
	std::vector<std::tuple<std::string, std::string>> selector_fills{{"mean", "standard_deviation"}};
	convert_figure_parameters(parameters, "selector_fills", selector_fills);
	std::string title{"{\\Large \\textrm{RTE\\ dependence\\ on\\ }$\\alpha$ }"};
	convert_figure_parameters(parameters, "title", title);
	std::string xlabel{"$\\alpha$"};
	convert_figure_parameters(parameters, "xlabel", xlabel);
	std::array<std::string, 2> RTE_source_names{"X", "Y"};
	convert_figure_parameters(parameters, "source_names", RTE_source_names);
	std::string output_base_name{"RTE.png"};
	convert_figure_parameters(parameters, "output_base_name", output_base_name);
	std::vector<std::string> suffices{"png"};
	convert_figure_parameters(parameters, "suffices", suffices);
	std::string color_map{"jet"};
	convert_figure_parameters(parameters, "color_map", color_map);
	std::string average_type_selector{renyi_entropy::average_runs};
	convert_figure_parameters(parameters, "average_type_selector", average_type_selector);
	double linewidth{1};
	convert_figure_parameters(parameters, "linewidth", linewidth);
	double alpha{1};
	convert_figure_parameters(parameters, "alpha", alpha);

	const std::array<std::string, 4> RTE_type_extra_names{"effective", "nonlinear", "ballance", "reverse"};

	for (const auto& suffix: suffices)
	{
		auto collection_result_RTE = processed_results[average_type_selector][std::get<0>(selector_fills[0])];
		auto collection_key_views = std::views::keys(collection_result_RTE);
		std::vector<renyi_entropy::renyi_entropy<calculation_type>::collection_conditional_information_transfer_key_type> collection_keys{ collection_key_views.begin(), collection_key_views.end() };

		plt.matplotlib_rcParams["text.usetex"] = py::bool_(true);

		for ( auto &[history_key, collection_of_data_for_history] :collection_result_RTE )
		{
			const auto first_history = std::get<1>( history_key );
			const auto first_future = std::get<0>( history_key );
			const auto second_history = std::get<2>( history_key );
			const std::string first_history_label{create_label_from_index(first_history, true)}, future_first_label{create_label_from_index(first_future, true)}, second_history_label{create_label_from_index(second_history, true)};
			const auto filename_history_suffix{create_label_from_index(first_history, false) + "_" + create_label_from_index(first_future, false) + "_" + create_label_from_index(second_history, false)};

			for (const auto &[selector_directon_type_sample, collection_of_data]: collection_of_data_for_history)
			{
				auto [RTE_balance, RTE_reverse_orientation, RTE_shuffled, RTE_surrogate, RTE_sample] = selector_directon_type_sample;

				std::string &source_timeseries{RTE_source_names[0]}, &destination_timeseries{RTE_source_names[1]};
				if (RTE_reverse_orientation)
				{
					source_timeseries = RTE_source_names[1];
					destination_timeseries = RTE_source_names[0];
				}

				const auto prefix_filename = (RTE_reverse_orientation ? RTE_type_extra_names[3]:"") + (RTE_reverse_orientation && (RTE_balance||RTE_shuffled||RTE_surrogate) ? "_" : "") + (RTE_balance ? RTE_type_extra_names[2]:"") + (RTE_balance&&(RTE_shuffled||RTE_surrogate)  ? "_" : "") + (RTE_shuffled ? RTE_type_extra_names[0]:"") + (RTE_shuffled && RTE_surrogate ? "_" : "") + (RTE_surrogate ? RTE_type_extra_names[1]:"") + (RTE_reverse_orientation || RTE_balance || RTE_shuffled || RTE_surrogate ? "_" : "");
				const auto file_output_name{prefix_filename + output_base_name + "_" + filename_history_suffix + "." + suffix};
				std::cout << file_output_name << std::endl;
				const std::string RTE_type_extra = std::format("{}", std::string{"\\textrm{"} + (RTE_shuffled ? RTE_type_extra_names[0]:"") + (RTE_shuffled && RTE_surrogate ? "," : "") + (RTE_surrogate ? RTE_type_extra_names[1]:"")  + "}");
				const auto RTE_type = "{" +std::format("\\left( R, {} \\right)", RTE_type_extra) + "}";
				const auto RTE_direction = "{" + std::format("\\alpha: {}\\longrightarrow {}", source_timeseries, destination_timeseries) + "}";
				const auto ylabel = std::format("$T^{}_{}\\left({},{},{} \\right)$", RTE_type, RTE_direction, first_history_label, future_first_label, second_history_label);

				auto fill_colormap = plt.matplotlib_colormaps [py::str(color_map)];
				auto selector_colormap = plt.matplotlib_colormaps [py::str(color_map_selector)];
				plt.matplotlib_pyplot_xlabel(py::str(xlabel));
				plt.matplotlib_pyplot_ylabel(py::str(ylabel));
				plt.matplotlib_pyplot_title(py::str(title));

				for (const auto &[raw_color_index, selector_fill] : std::views::enumerate(selector_fills))
				{
					std::vector<std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>> raw_results(1);
					auto mean_dataset = processed_results[average_type_selector][std::get<0>(selector_fill)][history_key][selector_directon_type_sample];
					auto std_dataset = processed_results[average_type_selector][std::get<1>(selector_fill)][history_key][selector_directon_type_sample];
					auto CRE_key_views = std::views::keys(mean_dataset);
					std::vector<renyi_entropy::renyi_entropy<calculation_type>::renyi_key_type> CRE_keys{ CRE_key_views.begin(), CRE_key_views.end() };
					for (auto & CRE_key: CRE_keys)
					{
						auto [neighbour, alpha] = CRE_key;
						const auto mean_value = mean_dataset[CRE_key];
						const auto std_value = mean_dataset[CRE_key];
						std::get<0>(raw_results[neighbour]).push_back(alpha);
						std::get<1>(raw_results[neighbour]).push_back(mean_value - std_value);
						std::get<2>(raw_results[neighbour]).push_back(mean_value + std_value);
					}

					for (auto neighbour = limits[0]; neighbour < limits[1]; ++ neighbour)
					{
						const auto color_position_within_colormap = static_cast<double>(raw_color_index) / ( static_cast<double>(selector_fills.size() - 1 ));
						auto color = fill_colormap( py::float_(color_position_within_colormap ));
						py::dict kwargs ("color"_a = color, "linewidth"_a = py::float_(linewidth), "alpha"_a = py::float_(alpha));
						if (average_type_selector == renyi_entropy::average_runs)
						{
							kwargs [py::str("label")] = py::str(std::format("$k={}$", neighbour));
						}
						plt.matplotlib_pyplot_fill_between(std::get<0>(raw_results[neighbour]), std::get<1>(raw_results[neighbour]), std::get<2>(raw_results[neighbour]), **kwargs);
					}
				}

				for (const auto &[raw_color_index, selector] : std::views::enumerate(selectors))
				{
					std::vector<std::tuple<std::vector<double>, std::vector<double>>> raw_results(1);
					auto dataset = processed_results[average_type_selector][selector][history_key][selector_directon_type_sample];
					auto CRE_key_views = std::views::keys(dataset);
					std::vector<renyi_entropy::renyi_entropy<calculation_type>::renyi_key_type> CRE_keys{ CRE_key_views.begin(), CRE_key_views.end() };
					for (auto & CRE_key: CRE_keys)
					{
						auto [neighbour, alpha] = CRE_key;
						const auto data_RTE_value = dataset[CRE_key];
						std::get<0>(raw_results[neighbour]).push_back(alpha);
						std::get<1>(raw_results[neighbour]).push_back(data_RTE_value);
					}

					for (auto neighbour = limits[0]; neighbour < limits[1]; ++ neighbour)
					{
						const auto color_position_within_colormap = static_cast<double>(raw_color_index) / ( static_cast<double>(selectors.size() - 1 ));
						const auto color = selector_colormap( color_position_within_colormap );
						py::dict kwargs ("color"_a = color, "linewidth"_a = py::float_(linewidth), "alpha"_a = py::float_(alpha));
						if (average_type_selector == renyi_entropy::average_runs)
						{
							kwargs["label"] = std::format("$k={}$", neighbour);
						}
						plt.matplotlib_pyplot_plot(std::get<0>(raw_results[neighbour]), std::get<1>(raw_results[neighbour]), **kwargs);
					}
				}
				py::dict kwargs_savefig = py::dict("dpi"_a = py::int_(dpi));
				plt.matplotlib_pyplot_savefig(file_output_name, **kwargs_savefig);
				plt.matplotlib_pyplot_close();
			}
		}
	}
}

int main ( int argc, char *argv[] )
{
	cpptrace::absorb_trace_exceptions(false);
	cpptrace::use_default_stderr_logger();
	cpptrace::register_terminate_handler();
	warmup_cpptrace();
	segfault_handler_cpptrace ();

	po::options_description desc ( "Allowed options" );
	desc.add_options() ( "help,h", "produce help message" )
	("directory,d", po::value<std::string>()->default_value("."), "Folder to export results" )
	("file,f", po::value<std::string>()->default_value("CRE.bin.zstd"), "Output file" );

	po::variables_map vm;
	po::store ( po::parse_command_line ( argc, argv, desc ), vm );
	po::notify ( vm );

	CPPTRACE_TRY
	{
		std::string output, directory;
		std::string color_map{"jet"};

		if ( vm.count ( "help" ) )
		{
				std::cout << desc << std::endl;
				cpptrace::generate_trace().print();
				return 1;
		}
		if ( vm.count ( "directory" ) )
		{
				directory = vm["directory"].as<std::string>();
		}
		output = vm["file"].as<std::string>();

		std::string str;
		boost::filesystem::path output_file =
		boost::filesystem::path(directory) / boost::filesystem::path(output);
		boost::filesystem::ifstream input_file_handler ( output_file );

		if (input_file_handler.is_open())
		{
			bio::filtering_istream istream;
			istream.push ( bio::zstd_decompressor{} );
			istream.push ( input_file_handler );

			std::string buffer(std::istreambuf_iterator<char>{istream}, {});
			istream >> buffer;
			msgpack::object_handle object_handle_buffer = msgpack::unpack(buffer.data(), buffer.size());

			// deserialized object is valid during the msgpack::object_handle instance is alive.
			msgpack::object deserialized = object_handle_buffer.get();

			//renyi_entropy::renyi_entropy<calculation_type>::storage_RTE
			renyi_entropy::renyi_entropy<calculation_type>::type_average_result_conditional_information_transfer_type result_conditional_information_transfer;

			deserialized.convert(result_conditional_information_transfer);
			auto &processed_results = result_conditional_information_transfer; //std::get<1>(result_conditional_information_transfer);

			{
				parameter_t parameters_fill{{"output_base_name", "std_RTE"}, {"color_map", "plasma"}, {"average_type_selector", renyi_entropy::average_runs_neighbors}, {"selector_fills", std::vector<std::tuple<std::string, std::string>>{std::tuple<std::string, std::string>{"mean", "standard_deviation"}}}, {"selector", "mean"}, {"color_map_selector", "jet"}, {"alpha", 0.8}, {"limits", std::array<unsigned int, 2>{0, 1}}};
				deviation_result_image (processed_results, parameters_fill);
			}

			{
				parameter_t parameters { {"selector", "mean"}, {"color_map", "plasma"}, {"output_base_name", "mean_RTE_R"}, {"average_type_selector", renyi_entropy::average_runs}, {"alpha", 0.8}, {"limits", std::array<unsigned int, 2>{1, 100}} };
				simple_result_image (processed_results, parameters);
				parameters["selector"] = "standard_deviation";
				parameters["output_base_name"] = "variance_RTE_R";
				simple_result_image (processed_results, parameters);
			}

			{
				parameter_t parameters{ {"selector", "mean"}, {"color_map", "plasma"}, {"output_base_name", "mean_RTE"}, {"average_type_selector", renyi_entropy::average_runs_neighbors}, {"alpha", 0.8}, {"limits", std::array<unsigned int, 2>{0, 1}} };
				simple_result_image (processed_results, parameters);
				parameters["selector"] = "standard_deviation";
				parameters["output_base_name"] = "variance_RTE";
				simple_result_image (processed_results, parameters);
			}

			{
				parameter_t parameters_fill{{"output_base_name", "quantile_RTE"}, {"color_map", "plasma"}, {"average_type_selector", renyi_entropy::average_runs_neighbors}, {"selector_fills", std::vector<std::tuple<std::string, std::string>>{std::tuple<std::string, std::string>{"quantile_09", "quantile_01"}, std::tuple<std::string, std::string>{"quantile_08", "quantile_02"}, std::tuple<std::string, std::string>{"quantile_07", "quantile_03"}, std::tuple<std::string, std::string>{"quantile_06", "quantile_04"}}}, {"selector", "median"}, {"color_map_selector", "jet"}, {"alpha", 0.8}, {"limits", std::array<unsigned int, 2>{0, 1}}};
				quantile_result_image (processed_results, parameters_fill);
			}
		}
		else
		{
			std::cerr << "Cannot open input file" << std::endl;
		}
	}
	CPPTRACE_CATCH (std::exception& exc)
	{
		std::cerr << "Exception: "<< exc.what() << std::endl;
		cpptrace::from_current_exception().print();
	}

}
