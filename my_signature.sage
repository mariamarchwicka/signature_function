#!/usr/bin/env python

import collections
import sys
import os
import inspect
import pandas as pd
import itertools as it
import re


class MySettings(object):
    def __init__(self):
        self.f_results = os.path.join(os.getcwd(), "results.out")
        self.only_slice_candidates = True
        self.only_slice_candidates = False
        self.default_limit = 3


def main(arg):
    try:
        limit = int(arg[1])
    except:
        limit = None
    second_knot_sum_formula = "[[k[0], k[1], k[2]], [k[3], k[4]], \
                                [-k[0], -k[3], -k[4]], [-k[1], -k[2]]]"
    first_knot_sum_formula = "[[k[0], k[1], k[2]], [k[3]],\
                                [-k[0], -k[1], -k[3]], [-k[2]]]"
    tests(second_knot_sum_formula, limit)
    tests(first_knot_sum_formula, limit)


def is_trivial_combination(knot_sum):
    if len(knot_sum) == 4:
        oposit_to_first = [-k for k in knot_sum[0]]
        if oposit_to_first in knot_sum:
            return True
    return False


def tests(knot_sum_formula, limit=None):
    settings = MySettings()
    if limit is None:
        limit = settings.default_limit

    k_vector_size = extract_max(knot_sum_formula) + 1

    with open(settings.f_results, 'w') as f_results:
        combinations = it.combinations_with_replacement(range(1, limit + 1),
                                                        k_vector_size)
        for comb in combinations:
            if settings.only_slice_candidates and k_vector_size == 5:
                comb = [comb[0], 4 * comb[0] + comb[1],
                     4 * (4 * comb[0] + comb[1]) + comb[2],
                     4 * comb[0] + comb[3],
                     4 * (4 * comb[0] + comb[3]) + comb[4]]
            k = comb
            knot_sum = eval(knot_sum_formula)

            if is_trivial_combination(knot_sum):
                continue
            result = eval_cable_for_thetas(knot_sum)
            if result is not None:
                knot_description, null_comb, all_comb = result
                line = (str(k) + ", " + str(null_comb) + ", " +
                        str(all_comb) + "\n")
                f_results.write(line)


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

    def value(self, arg):
        # Compute the value of the signature function at the point arg.
        # This requires summing all signature jumps that occur before arg.
        assert 0 <= arg < 1, \
            "Signature function is defined on the interval [0, 1)."
        val = 0
        for jump_arg, jump in self.data.items():
            if jump_arg < arg:
                val += 2 * jump
            elif jump_arg == arg:
                val += jump
        return val

    def sum_of_absolute_values(self):
        return sum([abs(i) for i in self.data.values()])

    def double_cover(self):
        new_data = []
        for jump_arg, jump in self.data.items():
            new_data.append((mod_one(jump_arg/2), jump))
            new_data.append((mod_one(1/2 + jump_arg/2), jump))
        return SignatureFunction(new_data)

    def __lshift__(self, shift):
        # Shift of the signature functions correspond to the rotations.
        return self.__rshift__(-shift)

    def __rshift__(self, shift):
        new_data = []
        for jump_arg, jump in self.data.items():
            new_data.append((mod_one(jump_arg + shift), jump))
        return SignatureFunction(new_data)

    def __sub__(self, other):
        # we can perform arithmetic operations on signature functions.
        return self + other.__neg__()

    def __neg__(self):
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

    def __str__(self):
        return '\n'.join([str(jump_arg) + ": " + str(jump)
                          for jump_arg, jump in sorted(self.data.items())])

    # def __repr__(self):
    #     return self.__str__()


# Proposition 9.8.
def get_blanchfield_for_pattern(k_n, theta):
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


# Bl_theta(K'_(2, d) =
# Bl_theta(T_2, d) + Bl(K')(ksi_l^(-theta) * t)
# + Bl(K')(ksi_l^theta * t)
def get_cable_signature_as_theta_function(*arg):
    def signture_function(theta):
        if theta > abs(arg[-1]):
            print "k for pattern is " + str(arg[-1])
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
    return signture_function


def get_untwisted_signature_function(j):
    # Return the signature function of the T_{2,2k+1} torus knot.
    k = abs(j)
    w = ([((2 * a + 1)/(4 * k + 2), -1 * sgn(j)) for a in range(k)] +
         [((2 * a + 1)/(4 * k + 2), 1 * sgn(j))
         for a in range(k + 1, 2 * k + 1)])
    return SignatureFunction(w)


def get_function_of_theta_for_sum(*arg):
    """
    Function intended to calculate signature function for a connected
    sum of multiple cables with varying theta parameter values.
    Accept arbitrary number of arguments (number of cables in connected sum).
    Each argument should be given as list of integer representing
    k - parameters for a cable: parameters k_i (i=1,.., n-1) for satelit knots
    T(2, 2k_i + 1) and - the last one - k_n for a pattern knot T(2, 2k_n + 1).
    Returns a function described below.
    """
    def signature_function_for_sum(*thetas):
        # Returns object of SignatureFunction class for a previously defined
        # connercted sum of len(arg) cables.
        # Accept len(arg) arguments: for each cable one theta parameter.
        # If call with no arguments, all theta parameters are set to be 0.
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


def mod_one(n):
    """This function returns the fractional part of some number."""
    return n - floor(n)


def eval_cable_for_thetas(knot_sum):

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
        else:
            assert sum(v_theta) != 0

    if null_combinations^2 >= all_combinations:
        print
        print knot_description
        print "Zero cases: " + str(null_combinations)
        print "All cases: " + str(all_combinations)
        print "Zero theta combinations: "
        for el in zero_theta_combinations:
            print el
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
     numbers = re.findall('\d+',string)
     numbers = map(int,numbers)
     return max(numbers)


if __name__ == '__main__' and '__file__' in globals():
    main(sys.argv)

#
# def get_sigma(t, k):
#     p = 2
#     q = 2 * k + 1
#     sigma_set = get_sigma_set(p, q)
#     sigma = len(sigma_set) - 2 * len([z for z in sigma_set if t < z < 1 + t])
#     return sigma
#
#
# def get_sigma_set(p, q):
#     sigma_set = set()
#     for i in range(1, p):
#         for j in range(1, q):
#             sigma_set.add(j/q + i/p)
#     return sigma_set
