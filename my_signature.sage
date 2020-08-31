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

def is_trivial_combination(knot_sum):
    # for now is applicable only for schema that are sums of 4 cables
    if len(knot_sum) == 4:
        oposit_to_first = [-k for k in knot_sum[0]]
        if oposit_to_first in knot_sum:
            return True
    return False

def get_shifted_combination(combination):
    # for now applicable only for schama
    # "[[k[0], k[1], k[2]], [k[3], k[4]],
    # [-k[0], -k[3], -k[4]], [-k[1], -k[2]]]"
    # shift the combination so that the knot can be a candidate for slice
    combination = [combination[0], 4 * combination[0] + combination[1],
         4 * (4 * combination[0] + combination[1]) + combination[2],
         4 * combination[0] + combination[3],
         4 * (4 * combination[0] + combination[3]) + combination[4]]
    return combination

def get_blanchfield_for_pattern(k_n, theta):
    if theta == 0:
        a = get_untwisted_signature_function(k_n)
        return a.square_root() + a.minus_square_root()

    results = []
    k = abs(k_n)
    ksi = 1/(2 * k + 1)

    # lambda_odd, i.e. (theta + e) % 2 != 0
    for e in range(1, k + 1):
        if (theta + e) % 2 != 0:
            results.append((e * ksi, 1 * sgn(k_n)))
            results.append((1 - e * ksi, -1 * sgn(k_n)))

    # lambda_even
    # print("normal")
    for e in range(1, theta):
        if (theta + e) % 2 == 0:
            results.append((e * ksi, 1 * sgn(k_n)))
            results.append((1 - e * ksi, -1 * sgn(k_n)))
    # print("reversed")
    for e in range(theta + 1, k + 1):
        if (theta + e) % 2 != 0:
            continue
        results.append((e * ksi, -1 * sgn(k_n)))
        results.append((1 - e * ksi, 1 * sgn(k_n)))
    return SignatureFunction(values=results)

def get_signature_summand_as_theta_function(*arg):
    def get_signture_function(theta):
        # TBD: another formula (for t^2) description

        k_n = abs(arg[-1])
        if theta > k_n:
            msg = "k for the pattern in the cable is " + str(arg[-1]) + \
                  ". Parameter theta should not be larger than abs(k)."
            raise ValueError(msg)

        # twisted part
        cable_signature = get_blanchfield_for_pattern(arg[-1], theta)

        # untwisted part
        for i, k in enumerate(arg[:-1][::-1]):
            ksi = 1/(2 * k_n + 1)
            power = 2^i
            a = get_untwisted_signature_function(k)
            shift = theta * ksi * power
            b = a >> shift
            c = a << shift
            for _ in range(i):
                b = b.double_cover()
                c = c.double_cover()
            cable_signature += b + c
            test = b - c
            test2 = -c + b
            assert test == test
        return cable_signature
    get_signture_function.__doc__ = get_signture_function_docsting
    return get_signture_function

def get_untwisted_signature_function(j):
    # return the signature function of the T_{2,2k+1} torus knot
    k = abs(j)
    w = ([((2 * a + 1)/(4 * k + 2), -1 * sgn(j)) for a in range(k)] +
         [((2 * a + 1)/(4 * k + 2), 1 * sgn(j))
         for a in range(k + 1, 2 * k + 1)])
    return SignatureFunction(values=w)

def extract_max(string):
    numbers = re.findall('\d+', string)
    numbers = map(int, numbers)
    return max(numbers)

def mod_one(n):
    return n - floor(n)


get_blanchfield_for_pattern.__doc__ = \
    """
    Arguments:
        k_n:    a number s.t. q_n = 2 * k_n + 1, where
                T(2, q_n) is a pattern knot for a single cable from a cable sum
        theta:  twist/character for the cable (value form v vector)
    Return:
        SignatureFunction created for twisted signature function
        for a given cable and theta/character
    Based on:
        Proposition 9.8. in Twisted Blanchfield Pairing
        (https://arxiv.org/pdf/1809.08791.pdf)
    """

TorusCable.get_number_of_combinations_of_theta.__doc__ = \
    """
    Arguments:
        arbitrary number of lists of numbers, each list encodes a single cable
    Return:
        number of possible theta values combinations that could be applied
        for a given cable sum,
        i.e. the product of q_j for j = {1,.. n},
        where n is a number of direct components in the cable sum,
        and q_j is the last q parameter for the component (a single cable)
    """

TorusCable.get_knot_descrption.__doc__ = \
    """
    Arguments:
        arbitrary number of lists of numbers, each list encodes a single cable.
    Examples:
        sage: get_knot_descrption([1, 3], [2], [-1, -2], [-3])
        'T(2, 3; 2, 7) # T(2, 5) # -T(2, 3; 2, 5) # -T(2, 7)'
    """

mod_one.__doc__ = \
    """
    Argument:
        a number
    Return:
        the fractional part of the argument
    Examples:
        sage: mod_one(9 + 3/4)
        3/4
        sage: mod_one(-9 + 3/4)
        3/4
        sage: mod_one(-3/4)
        1/4
    """

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

TorusCable.eval_cable_for_null_signature.__doc__ = \
    """
    This function calculates all possible twisted signature functions for
    a knot that is given as an argument. The knot should be encoded as a list
    of its direct component. Each component schould be presented as a list
    of integers. This integers correspond to the k - values in each component/
    cable. If a component is a mirror image of a cable the minus sign should
    be written before each number for this component. For example:
    eval_cable_for_null_signature([[1, 8], [2], [-2, -8], [-2]])
    eval_cable_for_null_signature([[1, 2], [-1, -2]])

    sage: eval_cable_for_null_signature([[1, 3], [2], [-1, -2], [-3]])

    T(2, 3; 2, 7) # T(2, 5) # -T(2, 3; 2, 5) # -T(2, 7)
    Zero cases: 1
    All cases: 1225
    Zero theta combinations:
    (0, 0, 0, 0)

    sage:
    The numbers given to the function eval_cable_for_null_signature
    are k-values for each component/cable in a direct sum.
    """

TorusCable.get_signature_as_function_of_theta.__doc__ = \
    """
    Function intended to construct signature function for a connected
    sum of multiple cables with varying theta parameter values.
    Accept arbitrary number of arguments (depending on number of cables in
    connected sum).
    Each argument should be given as list of integer representing
    k - parameters for a cable: parameters k_i (i=1,.., n-1) for satelit knots
    T(2, 2k_i + 1) and - the last one - k_n for a pattern knot T(2, 2k_n + 1).
    Returns a function that will take theta vector as an argument and return
    an object SignatureFunction.

    To calculate signature function for a cable sum and a theta values vector,
    use as below.

    sage: signature_function_generator = get_signature_as_function_of_theta(
                                             [1, 3], [2], [-1, -2], [-3])
    sage: sf = signature_function_generator(2, 1, 2, 2)
    sage: print(sf)
    0: 0
    5/42: 1
    1/7: 0
    1/5: -1
    7/30: -1
    2/5: 1
    3/7: 0
    13/30: -1
    19/42: -1
    23/42: 1
    17/30: 1
    4/7: 0
    3/5: -1
    23/30: 1
    4/5: 1
    6/7: 0
    37/42: -1

    Or like below.
    sage: print(get_signature_as_function_of_theta([1, 3], [2], [-1, -2], [-3]
                                                )(2, 1, 2, 2))
    0: 0
    1/7: 0
    1/6: 0
    1/5: -1
    2/5: 1
    3/7: 0
    1/2: 0
    4/7: 0
    3/5: -1
    4/5: 1
    5/6: 0
    6/7: 0
    """

get_signature_summand_as_theta_function.__doc__ = \
    """
    Argument:
        n integers that encode a single cable, i.e.
        values of q_i for T(2,q_0; 2,q_1; ... 2, q_n)
    Return:
        a function that returns SignatureFunction for this single cable
        and a theta given as an argument
    """
SignatureFunction.__doc__ = \
    """
    This simple class encodes twisted and untwisted signature functions
    of knots. Since the signature function is entirely encoded by its signature
    jump, the class stores only information about signature jumps
    in a dictionary self.cnt_signature_jumps.
    The dictionary stores data of the signature jump as a key/values pair,
    where the key is the argument at which the functions jumps
    and value encodes the value of the jump. Remember that we treat
    signature functions as defined on the interval [0,1).
    """
get_signture_function_docsting = \
    """
    This function returns SignatureFunction for previously defined single
    cable T_(2, q) and a theta given as an argument.
    The cable was defined by calling function
    get_signature_summand_as_theta_function(*arg)
    with the cable description as an argument.
    It is an implementaion of the formula:
        Bl_theta(K'_(2, d)) =
            Bl_theta(T_2, d) + Bl(K')(ksi_l^(-theta) * t)
            + Bl(K')(ksi_l^theta * t)
    """

signature_as_function_of_theta_docstring = \
    """
    Arguments:

    Returns object of SignatureFunction class for a previously defined
    connected sum of len(arg) cables.
    Accept len(arg) arguments: for each cable one theta parameter.
    If call with no arguments, all theta parameters are set to be 0.
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
