#!/usr/bin/python

# TBD: read about Factory Method, variable in docstring, sage documentation,
# print calc. to output file
# delete separation for twisted_part and untwisted_part
# decide about printing option

import os
import sys

import itertools as it
import re
import numpy as np


attach("cable_signature.sage")
attach("my_signature.sage")



class Config(object):
    def __init__(self):
        self.f_results = os.path.join(os.getcwd(), "results.out")

        # knot_formula is a schema for knots which signature function
        # will be calculated
        self.knot_formula = "[[k[0], k[1], k[3]], " + \
                             "[-k[1], -k[3]], " + \
                             "[k[2], k[3]], " + \
                             "[-k[0], -k[2], -k[3]]]"

        # self.knot_formula = "[[k[0], k[1], k[4]], [-k[1], -k[3]], \
        #                      [k[2], k[3]], [-k[0], -k[2], -k[4]]]"
        #
        # self.knot_formula = "[[k[3]], [-k[3]], \
        #                      [k[3]], [-k[3]] ]"
        #
        # self.knot_formula = "[[k[3], k[2], k[0]], [-k[2], -k[0]], \
        #                      [k[1], k[0]], [-k[3], -k[1], -k[0]]]"
        #
        # self.knot_formula = "[[k[0], k[1], k[2]], [k[3], k[4]], \
        #                      [-k[0], -k[3], -k[4]], [-k[1], -k[2]]]"
        # self.knot_formula = "[[k[0], k[1], k[2]], [k[3]],\
        #                          [-k[0], -k[1], -k[3]], [-k[2]]]"
        self.limit = 3

        self.verbose = True
        # self.verbose = False

def main(arg=None):
    try:
        limit = int(arg[1])
    except (IndexError, TypeError):
        limit = None

    global cable  # , cab_2, cab_1
    # self.knot_formula = "[[k[0], k[1], k[3]], " + \
    #                      "[-k[1], -k[3]], " + \
    #                      "[k[2], k[3]], " + \
    #                      "[-k[0], -k[2], -k[3]]]"

    # knot_formula = config.knot_formula
    # q_vector = (3, 5, 7, 13)
    # q_vector = (3, 5, 7, 11)

    # q_vector = (5, 13, 19, 41,\
    #             5, 17, 23, 43)

    formula_1 = "[[k[0], k[5], k[3]], " + \
                         "[-k[1], -k[3]], " + \
                         "[k[2], k[3]], " + \
                         "[-k[0], -k[2], -k[3]]]"
    formula_2 = "[[k[4], k[1], k[7]], " + \
                         "[-k[5], -k[7]], " + \
                         "[k[6], k[7]], " + \
                         "[-k[4], -k[6], -k[7]]]"
    q_vector = TorusCable.get_q_vector(formula_1[:-1] + ", " + formula_2[1:])
    cab_1 = TorusCable(knot_formula=formula_1, q_vector=q_vector)
    cab_2 = TorusCable(knot_formula=formula_2, q_vector=q_vector)
    cable = cab_1 + cab_2

if __name__ == '__main__':
    global config
    config = Config()
    if '__file__' in globals():
        # skiped in interactive mode as __file__ is not defined
        main(sys.argv)
    else:
        pass
        # main()


"""
This script calculates signature functions for knots (cable sums).

The script can be run as a sage script from the terminal
or used in interactive mode.

A knot (cable sum) is encoded as a list where each element (also a list)
corresponds to a cable knot, e.g. a list
[[1, 3], [2], [-1, -2], [-3]] encodes
T(2, 3; 2, 7) # T(2, 5) # -T(2, 3; 2, 5) # -T(2, 7).

To calculate the number of characters for which signature function vanish use
the function eval_cable_for_null_signature as shown below.

sage: eval_cable_for_null_signature([[1, 3], [2], [-1, -2], [-3]])

T(2, 3; 2, 7) # T(2, 5) # -T(2, 3; 2, 5) # -T(2, 7)
Zero cases: 1
All cases: 1225
Zero theta combinations:
(0, 0, 0, 0)

sage:

The numbers given to the function eval_cable_for_null_signature are k-values
for each component/cable in a direct sum.

To calculate signature function for a knot and a theta value, use function
get_signature_as_function_of_theta (see help/docstring for details).

About notation:
Cables that we work with follow a schema:
    T(2, q_1; 2, q_2; 2, q_4) # -T(2, q_2; 2, q_4) #
            # T(2, q_3; 2, q_4) # -T(2, q_1; 2, q_3; 2, q_4)
In knot_formula each k[i] is related with some q_i value, where
q_i = 2*k[i] + 1.
So we can work in the following steps:
1) choose a schema/formula by changing the value of knot_formula
2) set each q_i all or choose range in which q_i should varry
3) choose vector v / theata vector.
"""
