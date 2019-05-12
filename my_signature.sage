#!/usr/bin/env python

# TBD: remove part of the description to readme/example
"""
This script calculates signature functions for knots (cable sums).

The script can be run as a sage script from the terminal or used in inetactive
mode.


To calculate the number of characters for which signature function vanish use
the function eval_cable_for_thetas as shown below:

sage: eval_cable_for_thetas([[1, 3], [2], [-1, -2], [-3]])

T(2, 3; 2, 7) # T(2, 5) # -T(2, 3; 2, 5) # -T(2, 7)
Zero cases: 69
All cases: 1225
Zero theta combinations:
(0, 0, 0, 0)
(1, 1, 1, 1)
(1, 2, 2, 1)
(2, 1, 1, 2)
(2, 2, 2, 2)
(3, 0, 0, 3)
('T(2, 3; 2, 7) # T(2, 5) # -T(2, 3; 2, 5) # -T(2, 7)', 69, 1225)
sage:


The numbers given to the function eval_cable_for_thetas are k-values for each
component/cable in a direct sum.


To calculate signature function for a knot and a theta value, use function
get_function_of_theta_for_sum as follow:

sage: signature_function_generator = get_function_of_theta_for_sum([1, 3], [2], [-1, -2], [-3])
sage: sf = signature_function_generator(2, 1, 2, 2)
sage: print sf
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
sage:

or like below:

sage: print  get_function_of_theta_for_sum([1, 3], [2], [-1, -2], [-3])(2, 1, 2, 2)
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
sage:
"""

import os
import sys

import collections
import inspect
import itertools as it
import pandas as pd
import re


class MySettings(object):
    def __init__(self):
        self.f_results = os.path.join(os.getcwd(), "results.out")

        # is the ratio restriction for values in k_vector taken into account
        # False flag is usefull to make quick script tests
        self.only_slice_candidates = True
        # self.only_slice_candidates = False

        # knot_sum_formula is a schema for knots which signature function
        # will be calculated
        self.knot_sum_formula = "[[k[0], k[1], k[2]], [k[3], k[4]], \
                                 [-k[0], -k[3], -k[4]], [-k[1], -k[2]]]"
        # self.knot_sum_formula = "[[k[0], k[1], k[2]], [k[3]],\
        #                          [-k[0], -k[1], -k[3]], [-k[2]]]"
        self.default_limit = 3


class SignatureFunction(object):
    """
    This simple class encodes twisted and untwisted signature functions
    of knots. Since the signature function is entirely encoded by its signature
    jump, the class stores only information about signature jumps
    in a dictionary self.data.
    The dictionary stores data of the signature jump as a key/values pair,
    where the key is the argument at which the functions jumps
    and value encodes the value of the jump. Remember that we treat
    signature functions as defined on the interval [0,1).
    """
    def __init__(self, values=[]):
        # We will store data of signature jumps here.
        self.data = collections.defaultdict(int)
        # values contain initial data of singature jumps
        for jump_arg, jump in values:
            assert 0 <= jump_arg < 1, \
                "Signature function is defined on the interval [0, 1)."
            self.data[jump_arg] = jump


    def sum_of_absolute_values(self):
        return sum([abs(i) for i in self.data.values()])

    def double_cover(self):
        # to read values for t^2
        new_data = []
        for jump_arg, jump in self.data.items():
            new_data.append((mod_one(jump_arg/2), jump))
            new_data.append((mod_one(1/2 + jump_arg/2), jump))
        return SignatureFunction(new_data)

    def __lshift__(self, shift):
        # A shift of the signature functions corresponds to the rotation.
        return self.__rshift__(-shift)

    def __rshift__(self, shift):
        new_data = []
        for jump_arg, jump in self.data.items():
            new_data.append((mod_one(jump_arg + shift), jump))
        return SignatureFunction(new_data)

    def __neg__(self):
        # we can perform arithmetic operations on signature functions.
        new_data = []
        for jump_arg, jump in self.data.items():
            new_data.append(jump_arg, -jump)
        return SignatureFunction(new_data)

    def __add__(self, other):
        new_signature_function = SignatureFunction()
        new_data = collections.defaultdict(int)
        for jump_arg, jump in other.data.items():
            new_data[jump_arg] = jump + self.data.get(jump_arg, 0)
        for jump_arg, jump in self.data.items():
            if jump_arg not in new_data.keys():
                new_data[jump_arg] = self.data[jump_arg]
        new_signature_function.data = new_data
        return new_signature_function

    def __sub__(self, other):
        return self + other.__neg__()

    def __str__(self):
        return '\n'.join([str(jump_arg) + ": " + str(jump)
                          for jump_arg, jump in sorted(self.data.items())])


def main(arg):
    """
    This function is run if the script was called from the terminal.
    It calls another function to calculate signature functions for a schema
    of a cable sum defined in the class MySettings.
    Optionaly a parameter (a limit for k_0 value) can be given.
    Thought to be run for time consuming calculations.
    """

    try:
        new_limit = int(arg[1])
    except:
        new_limit = None
    perform_calculations(limit=new_limit)


def perform_calculations(knot_sum_formula=None, limit=None):
    """
    This function calculates signature functions for knots constracted
    accordinga a schema for a cable sum. The schema is given as an argument
    or defined in the class MySettings.
    Results of calculations will be writen to a file and to the stdout.
    limit is the upper bound for the first value in k_vector, i.e first k value
    in a cable sum (the number of knots that will be constracted depends
    on limit value).
    For each knot/cable sum the function eval_cable_for_thetas is called.
    eval_cable_for_thetas calculetes the number of all possible thetas
    (characters) and the number of combinations for which signature function
    equeles zero. In case the first number is larger than squere of the second,
    eval_cable_for_thetas returns None (i.e. the knot can not be slice).
    Data for knots that are candidates for slice knots are saved to a file.
    """

    settings = MySettings()
    if limit is None:
        limit = settings.default_limit
    if knot_sum_formula is None:
        knot_sum_formula = settings.knot_sum_formula

    k_vector_size = extract_max(knot_sum_formula) + 1
    combinations = it.combinations_with_replacement(range(1, limit + 1),
                                                    k_vector_size)

    with open(settings.f_results, 'w') as f_results:
        for k in combinations:
            # print
            print k
            # TBD:  maybe the following condition or the function
            # get_shifted_combination should be redefined to a dynamic version
            if settings.only_slice_candidates and k_vector_size == 5:
                k = get_shifted_combination(k)
            print k
            knot_sum = eval(knot_sum_formula)

            if is_trivial_combination(knot_sum):
                continue
            result = eval_cable_for_thetas(knot_sum, print_results=False)
            if result is not None:
                knot_description, null_comb, all_comb = result
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
    """
    This function calculates a twisted signature function for a given cable
    and theta/character. It returns object of class SignatureFunction.
    It is based on Proposition 9.8. in Twisted Blanchfield Pairing.
    """
    if theta == 0:
        return get_untwisted_signature_function(k_n)
    results = []
    k = abs(k_n)
    ksi = 1/(2 * k + 1)
    # lambda_odd (theta + e) % 2 == 0:
    for e in range(1, k + 1):
        if (theta + e) % 2 != 0:
            results.append((e * ksi, 1 * sgn(k_n)))
            results.append((1 - e * ksi, -1 * sgn(k_n)))
    # lambda_even
    # print "normal"
    for e in range(1, theta):
        if (theta + e) % 2 == 0:
            results.append((e * ksi, 1 * sgn(k_n)))
            results.append((1 - e * ksi, -1 * sgn(k_n)))
    # print "reversed"
    for e in range(theta + 1, k + 1):
        if (theta + e) % 2 != 0:
            continue
        results.append((e * ksi, -1 * sgn(k_n)))
        results.append((1 - e * ksi, 1 * sgn(k_n)))
    return SignatureFunction(results)


def get_cable_signature_as_theta_function(*arg):
    """
    This function takes as an argument a single cable and returns another
    function that alow to calculate signature function for previously defined
    cable and a theta given as an argument.
    """
    def get_signture_function(theta):
        """
        This function returns SignatureFunction for previously defined cable
        and a theta given as an argument.
        It is an implementaion of the formula:
        Bl_theta(K'_(2, d)) =
        Bl_theta(T_2, d) + Bl(K')(ksi_l^(-theta) * t)
        + Bl(K')(ksi_l^theta * t)
        """
        if theta > abs(arg[-1]):
            print "k for the pattern in the cable is " + str(arg[-1])
            print "theta shouldn't be larger than this"
            return None
        cable_signature = get_blanchfield_for_pattern(arg[-1], theta)

        for i, k in enumerate(arg[:-1][::-1]):
            ksi = 1/(2 * abs(k) + 1)
            power = 2^i
            a = get_untwisted_signature_function(k)
            shift = theta * ksi * power
            b = a >> shift
            c = a << shift
            for _ in range(i):
                b = b.double_cover()
                c = c.double_cover()
            cable_signature += b + c
        return cable_signature
    return get_signture_function


def get_untwisted_signature_function(j):
    """This function returns the signature function of the T_{2,2k+1}
    torus knot."""
    k = abs(j)
    w = ([((2 * a + 1)/(4 * k + 2), -1 * sgn(j)) for a in range(k)] +
         [((2 * a + 1)/(4 * k + 2), 1 * sgn(j))
         for a in range(k + 1, 2 * k + 1)])
    return SignatureFunction(w)


def get_function_of_theta_for_sum(*arg):
    """
    Function intended to calculate signature function for a connected
    sum of multiple cables with varying theta parameter values.
    Accept arbitrary number of arguments (depending on number of cables in
    connected sum).
    Each argument should be given as list of integer representing
    k - parameters for a cable: parameters k_i (i=1,.., n-1) for satelit knots
    T(2, 2k_i + 1) and - the last one - k_n for a pattern knot T(2, 2k_n + 1).
    Returns a function that will take theta vector as an argument and return
    an object SignatureFunction.
    """

    def signature_function_for_sum(*thetas):
        """
        Returns object of SignatureFunction class for a previously defined
        connercted sum of len(arg) cables.
        Accept len(arg) arguments: for each cable one theta parameter.
        If call with no arguments, all theta parameters are set to be 0.
        """

        la = len(arg)
        lt = len(thetas)
        if lt == 0:
            return signature_function_for_sum(*(la * [0]))
        if lt != la:
            msg = "This function takes exactly " + str(la) + \
                  " arguments or no argument at all (" + str(lt) + " given)."
            raise TypeError(msg)
        sf = SignatureFunction([(0, 0)])
        for i, knot in enumerate(arg):
            sf += (get_cable_signature_as_theta_function(*knot))(thetas[i])
        return sf
    return signature_function_for_sum


def eval_cable_for_thetas(knot_sum, print_results=True):
    """
    This function calculates all possible twisted signature functions for
    a knot that is given as an argument. The knot should be encoded as a list
    of its direct component. Each component schould be presented as a list
    of integers. This integers correspond to the k - values in each component/
    cable. If a component is a mirror image of a cable the minus sign should
    be written before each number for this component. For example:
    eval_cable_for_thetas([[1, 8], [2], [-2, -8], [-2]])
    eval_cable_for_thetas([[1, 2], [-1, -2]])
    """
    f = get_function_of_theta_for_sum(*knot_sum)
    knot_description = get_knot_descrption(*knot_sum)
    all_combinations = get_number_of_combinations(*knot_sum)

    null_combinations = 0
    zero_theta_combinations = []

    ranges_list = [range(abs(knot[-1]) + 1) for knot in knot_sum]
    for v_theta in it.product(*ranges_list):

        if f(*v_theta).sum_of_absolute_values() == 0:
            zero_theta_combinations.append(v_theta)
            m = len([theta for theta in v_theta if theta != 0])
            null_combinations += 2^m
        # else:
        #     assert sum(v_theta) != 0

    if print_results:
        print
        print knot_description
        print "Zero cases: " + str(null_combinations)
        print "All cases: " + str(all_combinations)
        if zero_theta_combinations:
            print "Zero theta combinations: "
            for el in zero_theta_combinations:
                print el
    if null_combinations^2 >= all_combinations:
        return knot_description, null_combinations, all_combinations
    return None


def get_knot_descrption(*arg):
    description = ""
    for knot in arg:
        if knot[0] < 0:
            description += "-"
        description += "T("
        for k in knot:
            description += "2, " + str(2 * abs(k) + 1) + "; "
        description = description[:-2] + ") # "
    return description[:-3]


def get_number_of_combinations(*arg):
    number_of_combinations = 1
    for knot in arg:
        number_of_combinations *= (2 * abs(knot[-1]) + 1)
    return number_of_combinations


def extract_max(string):
     """This function returns maximal number from given string."""
     numbers = re.findall('\d+', string)
     numbers = map(int, numbers)
     return max(numbers)


def mod_one(n):
    """This function returns the fractional part of some number."""
    return n - floor(n)


if __name__ == '__main__' and '__file__' in globals():
    main(sys.argv)
