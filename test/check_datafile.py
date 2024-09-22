import pprint
import datetime
import numpy as np
import pandas as pd
from pathlib import Path

if __name__ == "__main__":
    path = Path(
        "/home/hynek/work/skola/prog/python/transfer_entropy/data/test/refined_conditional_information_transfer-visa_visa.bin")  # refined_conditional_information_transfer-visa_visa.bin conditional_information_transfer-visa_visa.bin
    print(
        f"{datetime.datetime.now().isoformat()} Processing input file {path}"
    )
    table = pd.read_pickle(path)

    if isinstance(table, dict):
        refined_table_converted_multiindices = {}
        for key, value in table.items():
            refined_table_converted_multiindices[key] = {}
            for key_q, items in value.items():
                refined_key = key
                for index, value in enumerate(items):
                    refined_table_converted_multiindices[refined_key][
                        (index, key_q)
                    ] = value

        frame_refined = pd.DataFrame.from_dict(refined_table_converted_multiindices)
        frame_refined.index.set_names(["k", "q"], inplace=True)

    col = table.columns

    multi_sample_columns = {}
    for item in col[:-1]:
        (
            variable,
            history_first,
            future_first,
            history_second,
            statistical_value,
            shuffled_dataset,
            reversed_order,
            sample_number,
            remark,
            sample_dataset,
        ) = item
        if (
                not (
                            variable,
                            history_first,
                            future_first,
                            history_second,
                            shuffled_dataset,
                            reversed_order,
                            remark,
                    )
                    in multi_sample_columns
        ):
            multi_sample_columns[
                (
                    variable,
                    history_first,
                    future_first,
                    history_second,
                    shuffled_dataset,
                    reversed_order,
                    remark,
                )
            ] = 1
        else:
            multi_sample_columns[
                (
                    variable,
                    history_first,
                    future_first,
                    history_second,
                    shuffled_dataset,
                    reversed_order,
                    remark,
                )
            ] += 1

    pprint.pprint(multi_sample_columns)

    for key, value in multi_sample_columns.items():
        if value > 300:
            pprint.pprint(f"{key} {value}")

    time_set_name_variable = table.xs(
        "effective_transfer_entropy", level="Name of variable", axis=1
    )
    time_set_remark = time_set_name_variable.xs(
        "Avaraging over nearest neighbours starting from 5", level="Remark", axis=1
    )
    time_set_selection_hf = time_set_remark.xs(
        "0,1,2,3,4", level="History first", axis=1
    )
    time_set_selection_ff = time_set_selection_hf.xs(
        "1", level="Future first", axis=1
    )
    time_set_selection_hs = time_set_selection_ff.xs(
        "0,1,2,3,4", level="History second", axis=1
    )
    pprint.pprint(time_set_selection_hs)
