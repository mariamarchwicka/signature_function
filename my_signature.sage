#!/usr/bin/python

# TBD: read about Factory Method, variable in docstring, sage documentation
# move settings to sep file

import os
import sys

import collections
# import inspect
import itertools as it
import numpy as np
import re
# try:
#     from cable_signature import SignatureFunction, TorusCable
# except ModuleNotFoundError:
os.system('sage --preparse cable_signature.sage')
os.system('mv cable_signature.sage.py cable_signature.py')
from cable_signature import SignatureFunction, TorusCable

class Config(object):
    def __init__(self):
        self.f_results = os.path.join(os.getcwd(), "results.out")

        # knot_formula is a schema for knots which signature function
        # will be calculated
        self.knot_formula = "[[k[0], k[1], k[3]], [-k[1], -k[3]], \
                             [k[2], k[3]], [-k[0], -k[2], -k[3]]]"

        # self.knot_formula = "[[k[3], k[2], k[0]], [-k[2], -k[0]], \
        #                      [k[1], k[0]], [-k[3], -k[1], -k[0]]]"



        # self.knot_formula = "[[k[0], k[1], k[2]], [k[3], k[4]], \
        #                      [-k[0], -k[3], -k[4]], [-k[1], -k[2]]]"
        # self.knot_formula = "[[k[0], k[1], k[2]], [k[3]],\
        #                          [-k[0], -k[1], -k[3]], [-k[2]]]"
        self.limit = 3

        # in search for large sigma, for 1. checked knot q_1 = 3 + start_shift
        self.start_shift =  0

        self.verbose = True
        self.verbose = False

        self.print_results = True
        # self.print_results = False

        self.print_calculations_for_large_sigma = True
        self.print_calculations_for_large_sigma = False

        # is the ratio restriction for values in k_vector taken into account
        # False flag is usefull to make quick script tests
        self.only_slice_candidates = True
        self.only_slice_candidates = False


    # range for a_i, v = [a_1, a_2, a_3, a_4], for sigma calculations
    def get_list_of_ranges(self, q):
        list_of_ranges = [
            # all characters a_1, a_2, a_3, a_4 != 0
            it.product(range(1, q), range(1, q), range(1, q), range(1, 2)),

            # a_1 == 0, a_2, a_3, a_4 != 0
            it.product(range(1), range(1, q), range(1, q), range(1, 2)),
            # a_2 == 0, a_1, a_3, a_4 != 0
            it.product(range(1, q), range(1), range(1, q), range(1, 2)),
            # a_3 == 0, a_1, a_2, a_4 != 0
            it.product(range(1, q), range(1, q), range(1), range(1, 2)),
            # a_4 == 0, a_1, a_2, a_3 != 0
            it.product(range(1, q), range(1, q), range(1, 2), range(1)),

            # a_1 == 0, a_2 == 0, a_3, a_4 != 0
            it.product(range(1), range(1), range(1, q), range(1, 2)),
            # a_1 == 0, a_3 == 0, a_2, a_4 != 0
            it.product(range(1), range(1, q), range(1), range(1, 2)),
            # a_1 == 0, a_4 == 0, a_3, a_2 != 0
            it.product(range(1), range(1, q), range(1, 2), range(1)),
            # a_2 == 0, a_3 == 0, a_1, a_4 != 0
            it.product(range(1, q), range(1), range(1), range(1, 2)),
            # a_2 == 0, a_4 == 0, a_1, a_3 != 0
            it.product(range(1, q), range(1), range(1, 2), range(1)),
            # a_3 == 0, a_4 == 0, a_1, a_2 != 0
            it.product(range(1, q), range(1, 2), range(1), range(1)),
            ]
        # list_of_ranges = [
        #     # all characters a_1, a_2, a_3, a_4 != 0
        #     # 1, 1, 1, 1
        #     it.product(range(1, 2), range(1, 2), range(1, 2), range(1, 2)),
        #
        #     # -1, -1, -1, 1
        #     it.product(range(q - 1, q), range(q - 1, q), range(q - 1, q), range(1, 2)),
        #
        #     # 1, -1, -1, 1
        #     it.product(range(1, 2), range(q - 1, q), range(q - 1, q), range(1, 2)),
        #     # -1 , -1, 1, 1
        #     it.product(range(q - 1, q), range(q - 1, q), range(1, 2), range(1, 2)),
        #     # -1, 1, -1, 1
        #     it.product(range(q - 1, q), range(1, 2), range(q - 1, q), range(1, 2)),
        #
        #     # 1, 1, -1, 1
        #     it.product(range(1, 2), range(1, 2), range(q - 1, q), range(1, 2)),
        #     # 1, -1, 1, 1
        #     it.product(range(1, 2), range(q - 1, q), range(1, 2), range(1, 2)),
        #     # -1, 1, 1, 1
        #     it.product(range(q - 1, q), range(1, 2), range(1, 2), range(1, 2)),
        #
        # ]

        return list_of_ranges




def main(arg):
    if arg[1]:
        limit = int(arg[1])
    else:
        limit = None
    knots_with_large_sigma = search_for_large_signature_value(limit=limit)
    # search_for_null_signature_value(limit=limit)



# searching for sigma > 5 + #(v_i != 0) over given knot schema
def __search_for_large_signature_value(knot_formula, limit,
                                        verbose, print_results):
    # number of k_i (q_i) variables to substitute
    k_size = extract_max(knot_formula) + 1
    combinations = it.combinations_with_replacement(range(0, limit + 1), k_size)
    P = Primes()
    good_knots = []
    # iterate over q-vector
    for c in combinations:
        q = list(c)
        q[0] = P.unrank(q[0] + 1 + config.start_shift)
        q[1] = P.next(q[0] * 4 + q[1])
        q[2] = P.next(q[1] * 4 + q[2])
        q[3] = P.next(q[2] * 4 + q[3])
        cable = TorusCable(knot_formula=knot_formula, q_vector=q)
        list_of_ranges = config.get_list_of_ranges(cable.q_vector[-1])
        if cable.eval_cable_for_large_sigma(list_of_ranges, verbose=verbose,
                                            print_results=print_results):
            good_knots.append(cable)
    return good_knots

# searching for sigma > 5 + #(v_i != 0) over given knot schema
def search_for_large_signature_value(knot_formula=None, limit=None,
                                     verbose=None, print_results=None):
    if limit is None:
        limit = config.limit
    if knot_formula is None:
        knot_formula = config.knot_formula
    if verbose is None:
        vebose = config.verbose
    if print_results is None:
        print_results = config.print_results

    k_vector_size = extract_max(knot_formula) + 1
    limit = max(limit, k_vector_size)
    if config.only_slice_candidates:
        return __search_for_large_signature_value(knot_formula, limit, verbose,
                                                    print_results)

    # number of k_i (q_i) variables to substitute
    combinations = it.combinations(range(1, limit + 1), k_vector_size)
    P = Primes()
    good_knots = []

    # iterate over q-vector
    for c in combinations:
        k = [(P.unrank(i + config.start_shift) - 1)/2 for i in c]
        cable = TorusCable(knot_formula=knot_formula, k_vector=k)
        list_of_ranges = config.get_list_of_ranges(cable.q_vector[-1])
        if cable.eval_cable_for_large_sigma(list_of_ranges, verbose=verbose,
                                            print_results=print_results):
            good_knots.append(cable)
    return good_knots

# searching for signature == 0
def search_for_null_signature_value(knot_formula=None, limit=None):
    if limit is None:
        limit = config.limit
    if knot_formula is None:
        knot_formula = config.knot_formula
    print_results = config.print_results
    verbose = config.verbose

    k_vector_size = extract_max(knot_formula) + 1
    combinations = it.combinations_with_replacement(range(1, limit + 1),
                                                    k_vector_size)
    with open(config.f_results, 'w') as f_results:
        for k in combinations:
            if config.only_slice_candidates and k_vector_size == 5:
                k = get_shifted_combination(k)
            cable = TorusCable(knot_formula, k_vector=k)
            if is_trivial_combination(cable.knot_sum):
                print(cable.knot_sum)
                continue
            result = cable.eval_cable_for_null_signature(verbose=verbose,
                                                    print_results=print_results)
            if result is not None:
                null_comb, all_comb = result
                line = (str(k) + ", " + str(null_comb) + ", " +
                        str(all_comb) + "\n")
                f_results.write(line)

def extract_max(string):
    numbers = re.findall('\d+', string)
    numbers = map(int, numbers)
    return max(numbers)

def is_trivial_combination(knot_sum):
    # for now is applicable only for schema that are sums of 4 cables
    if len(knot_sum) == 4:
        oposit_to_first = [-k for k in knot_sum[0]]
        if oposit_to_first in knot_sum:
            return True
    return False


search_for_null_signature_value.__doc__ = \
    """
    This function calculates signature functions for knots constracted
    accordinga a schema for a cable sum. The schema is given as an argument
    or defined in the class Config.
    Results of calculations will be writen to a file and the stdout.
    limit is the upper bound for the first value in k_vector,
    i.e k[0] value in a cable sum, where q_0 = 2 * k[0] + 1.

    (the number of knots that will be constracted depends on limit value).
    For each knot/cable sum the function eval_cable_for_null_signature
    is called.
    eval_cable_for_null_signature calculetes the number of all possible thetas
    (characters) and the number of combinations for which signature function
    equeles zero. In case the first number is larger than squere of the second,
    eval_cable_for_null_signature returns None (i.e. the knot can not be slice).
    Data for knots that are candidates for slice knots are saved to a file.
    """

main.__doc__ = \
    """
    This function is run if the script was called from the terminal.
    It calls another function, search_for_null_signature_value,
    to calculate signature functions for a schema
    of a cable sum defined in the class Config.
    Optionaly a parameter (a limit for k_0 value) can be given.
    Thought to be run for time consuming calculations.
    """


extract_max.__doc__ = \
    """
    Return:
        maximum of absolute values of numbers from given string
    Examples:
        sage: extract_max("([1, 3], [2], [-1, -2], [-10])")
        10
        sage: extract_max("3, 55, ewewe, -42, 3300, 50")
        3300
    """


if __name__ == '__main__':
    global config
    config = Config()
    if '__file__' in globals():
        # skiped in interactive mode as __file__ is not defined
        main(sys.argv)

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
