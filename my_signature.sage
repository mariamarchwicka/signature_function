#!/usr/bin/python

def calculate_form(x, y, q4):
    x1, x2, x3, x4 = x
    y1, y2, y3, y4 = y
    form = (x1 * y1 - x2 * y2 + x3 * y3 - x4 * y4) % q_4
    # TBD change for ring modulo q_4
    return form

def check_condition(v, q4):
    form = calculate_form(v, v, q4)
    if form:
        return False
    return True

def find_v(q4):
    results = []
    for i in range(q4):
        for j in range(q4):
            for k in range(q4):
                for m in range(q4):
                    if check_condition([i, j, k, m], q_4):
                        results.add(v)
    return results

def check_inequality(q, v):
    a1, a2, a3, a4 = v
    q1, q2, q3, q4 = q
    pattern = [q1, q2, q4],[-q2, -q4],[q3, q4],[-q1, -q3, -q4]
    signature_function_generator = get_function_of_theta_for_sum(pattern)
    signature_function_for_sum = signature_function_generator(a1, a2, a3, a4)

    # sigma_v = sigma(q4, a1) - s(a2) + s(a3) - s(a4)


"""
This script calculates signature functions for knots (cable sums).

The script can be run as a sage script from the terminal
or used in interactive mode.

A knot (cable sum) is encoded as a list where each element (also a list)
corresponds to a cable knot, e.g. a list
[[1, 3], [2], [-1, -2], [-3]] encodes
T(2, 3; 2, 7) # T(2, 5) # -T(2, 3; 2, 5) # -T(2, 7).

To calculate the number of characters for which signature function vanish use
the function eval_cable_for_thetas as shown below.

sage: eval_cable_for_thetas([[1, 3], [2], [-1, -2], [-3]])

T(2, 3; 2, 7) # T(2, 5) # -T(2, 3; 2, 5) # -T(2, 7)
Zero cases: 1
All cases: 1225
Zero theta combinations:
(0, 0, 0, 0)

sage:


The numbers given to the function eval_cable_for_thetas are k-values for each
component/cable in a direct sum.

To calculate signature function for a knot and a theta value, use function
get_function_of_theta_for_sum (see help/docstring for details).
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
        """
        About notation:
        Cables that we work with follow a schema:
        T(2, q_0; 2, q_1; 2, q_2) # T(2, q_1; 2, q_2) #
                # -T(2, q_3; 2, q_2) # -T(2, q_0; 2, q_3; 2, q_2)
        In knot_sum_formula each k[i] is related with some q_i value, where
        q_i = 2*k[i] + 1.
        So we can work in the following steps:
        1) choose a schema/formula by changing the value of knot_sum_formula
        2) set each q_i all or choose range in which q_i should varry
        3) choose vector v / theata vector.

        """
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

    def square_root(self):
        # to read values for t^(1/2)
        new_data = []
        for jump_arg, jump in self.data.items():
            if jump_arg < 1/2:
                new_data.append((2 * jump_arg, jump))
        return SignatureFunction(new_data)

    def get_signture_jump(self, t):
        return self.data.get(t, 0)

    def minus_square_root(self):
        # to read values for t^(1/2)
        new_data = []
        for jump_arg, jump in self.data.items():
            if jump_arg >= 1/2:
                new_data.append((mod_one(2 * jump_arg), jump))
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
        return ''.join([str(jump_arg) + ": " + str(jump) + "\n"
                          for jump_arg, jump in sorted(self.data.items())])


def main(arg):
    """
    This function is run if the script was called from the terminal.
    It calls another function, perform_calculations,
    to calculate signature functions for a schema
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
    Results of calculations will be writen to a file and the stdout.
    limit is the upper bound for the first value in k_vector,
    i.e k[0] value in a cable sum, where q_0 = 2 * k[0] + 1.

    (the number of knots that will be constracted depends on limit value).
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
            # print k
            # TBD:  maybe the following condition or the function
            # get_shifted_combination should be redefined to a dynamic version
            if settings.only_slice_candidates and k_vector_size == 5:
                k = get_shifted_combination(k)
            # print k
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

    # TBD: k_n explanation

    if theta == 0:
        a = get_untwisted_signature_function(k_n)
        return a.square_root() + a.minus_square_root()
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
    This function takes as an argument a single cable T_(2, q), i.e.
    arbitrary number of integers that encode the cable,
    and returns another function that alow to calculate signature function
    for this single cable and a theta given as an argument.
    """
    def get_signture_function(theta):
        """
        This function returns SignatureFunction for previously defined single
        cable T_(2, q) and a theta given as an argument.
        The cable was defined by calling function
        get_cable_signature_as_theta_function(*arg)
        with the cable description as an argument.
        It is an implementaion of the formula:
        Bl_theta(K'_(2, d)) =
        Bl_theta(T_2, d) + Bl(K')(ksi_l^(-theta) * t)
        + Bl(K')(ksi_l^theta * t)
        """

        # TBD: another formula (for t^2) description

        k_n = abs(arg[-1])
        if theta > k_n:
            msg = "k for the pattern in the cable is " + str(arg[-1]) + \
                  ". Parameter theta should not be larger than abs(k)."
            raise ValueError(msg)
        cable_signature = get_blanchfield_for_pattern(arg[-1], theta)
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

    sage: signature_function_generator = get_function_of_theta_for_sum(
                                             [1, 3], [2], [-1, -2], [-3])
    sage: sf = signature_function_generator(2, 1, 2, 2)
    sage: print sf
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
    sage: print get_function_of_theta_for_sum([1, 3], [2], [-1, -2], [-3]
                                                )(2, 1, 2, 2)
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

    def signature_function_for_sum(*thetas, **kwargs):
        """
        Returns object of SignatureFunction class for a previously defined
        connected sum of len(arg) cables.
        Accept len(arg) arguments: for each cable one theta parameter.
        If call with no arguments, all theta parameters are set to be 0.
        """
        if 'verbose' in kwargs:
            verbose = kwargs['verbose']
        else:
            verbose = False

        la = len(arg)
        lt = len(thetas)

        # call with no arguments
        if lt == 0:
            return signature_function_for_sum(*(la * [0]))

        if lt != la:
            msg = "This function takes exactly " + str(la) + \
                  " arguments or no argument at all (" + str(lt) + " given)."
            raise TypeError(msg)

        sf = SignatureFunction([(0, 0)])

        # for each cable in cable sum apply theta
        for i, knot in enumerate(arg):
            try:
                sf += (get_cable_signature_as_theta_function(*knot))(thetas[i])
            # in case wrong theata value was given
            except ValueError as e:
                print "ValueError: " + str(e.args[0]) +\
                      " Please change " + str(i + 1) + ". parameter."
                return None
        if verbose:
            print
            print str(*thetas)
            print sf
        return sf
    return signature_function_for_sum


def eval_cable_for_thetas(knot_sum, print_results=True, verbose=False):
    """
    This function calculates all possible twisted signature functions for
    a knot that is given as an argument. The knot should be encoded as a list
    of its direct component. Each component schould be presented as a list
    of integers. This integers correspond to the k - values in each component/
    cable. If a component is a mirror image of a cable the minus sign should
    be written before each number for this component. For example:
    eval_cable_for_thetas([[1, 8], [2], [-2, -8], [-2]])
    eval_cable_for_thetas([[1, 2], [-1, -2]])

    sage: eval_cable_for_thetas([[1, 3], [2], [-1, -2], [-3]])

    T(2, 3; 2, 7) # T(2, 5) # -T(2, 3; 2, 5) # -T(2, 7)
    Zero cases: 1
    All cases: 1225
    Zero theta combinations:
    (0, 0, 0, 0)

    sage:
    The numbers given to the function eval_cable_for_thetas are k-values for each
    component/cable in a direct sum.


    """
    f = get_function_of_theta_for_sum(*knot_sum)
    knot_description = get_knot_descrption(*knot_sum)
    all_combinations = get_number_of_combinations(*knot_sum)

    null_combinations = 0
    zero_theta_combinations = []

    ranges_list = [range(abs(knot[-1]) + 1) for knot in knot_sum]
    if verbose:
        print
        print knot_description
    for v_theta in it.product(*ranges_list):
        if f(*v_theta, verbose=verbose).sum_of_absolute_values() == 0:
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


def check_squares(a, k):
    print
    p = 2 * k + 1
    k_0 = (p^2 - 1)/2
    knot_sum = [[a, k], [k_0], [-a, -k_0], [-k]]
    print get_knot_descrption(*knot_sum)
    if a * 4 >= p or is_trivial_combination(knot_sum):
        if a * 4 >= p:
            print str(knot_sum)
            print "a * 4 >= p"
        else:
            print "Trivial " + str(knot_sum)
        return None

    eval_cable_for_thetas(knot_sum)


def get_number_of_combinations(*arg):
    """
    Arguments:
        arbitrary number of lists of numbers, each list encodes a single cable.
    Return:
        number of possible theta values combinations that could be applied
        for a given cable sum,
        i.e. the product of q_j for j = {1,.. n},
        where n is a number of direct components in the cable sum,
        and q_j is the last q parameter for the component (a single cable).
    """
    number_of_combinations = 1
    for knot in arg:
        number_of_combinations *= (2 * abs(knot[-1]) + 1)
    return number_of_combinations


def extract_max(string):
    """
    Return:
        maximum of absolute values of numbers from given string
    Examples:
        sage: extract_max("([1, 3], [2], [-1, -2], [-10])")
        10
        sage: extract_max("3, 55, ewewe, -42, 3300, 50")
        3300
    """
    numbers = re.findall('\d+', string)
    numbers = map(int, numbers)
    return max(numbers)


def mod_one(n):
    """
    Argument:
        a number
    Return:
        the fractional part of a number
    Examples:
        sage: mod_one(9 + 3/4)
        3/4
        sage: mod_one(-9 + 3/4)
        3/4
        sage: mod_one(-3/4)
        1/4
    """
    return n - floor(n)


def get_knot_descrption(*arg):
    """
    Arguments:
        arbitrary number of lists of numbers, each list encodes a single cable.
    Examples:
        sage: get_knot_descrption([1, 3], [2], [-1, -2], [-3])
        'T(2, 3; 2, 7) # T(2, 5) # -T(2, 3; 2, 5) # -T(2, 7)'
    """
    description = ""
    for knot in arg:
        if knot[0] < 0:
            description += "-"
        description += "T("
        for k in knot:
            description += "2, " + str(2 * abs(k) + 1) + "; "
        description = description[:-2] + ") # "
    return description[:-3]


if __name__ == '__main__' and '__file__' in globals():
    # not called in interactive mode as __file__ is not defined
    main(sys.argv)
