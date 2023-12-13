#!/usr/bin/env python3
# -*- coding: utf-8 -*-


def process_CLI_arguments(arguments, separator=(",", "'", "/", "|"), separator_of_groups=(";", ":", "~")):
    processed_arguments = []

    new_group = [[]]
    new_set = new_group[0]
    for item in arguments:

        if item in separator_of_groups:
            processed_arguments.append(new_group)
            new_group = [[]]
            new_set = new_group[0]
        elif item in separator:
            new_group.append([])
            new_set = new_group[-1]
        else:
            new_set.append(int(item))

    processed_arguments.append(new_group)
    return processed_arguments
