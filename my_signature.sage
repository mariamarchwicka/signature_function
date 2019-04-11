#!/usr/bin/env python

import collections
import sys
import inspect
import pandas as pd
import itertools as it


class MySettings(object):
    def __init__(self):
        k = 0

def main(arg):
    my_settings = MySettings()
    try:
        tests(int(arg[1]))
    except:
        tests()



def tests(limit=10):
    for comb in it.combinations_with_replacement(range(1, limit + 1), 5):
        knot_description, null_comb, all_comb = second_sum(*comb)
        if null_comb^2 >= all_comb:
            print "\n\nHURA!!"
            print comb
            print knot_description
            print "Zero cases: " + str(null_comb)
            print "All cases: " + str(all_comb)
    # for comb in it.combinations_with_replacement(range(1, limit + 1), 4):
    #     print comb
    #     print first_sum(*comb)


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
            b += c
            cable_signature += b
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


# ###################### TEMPORARY TESTS #########

# def first_sum(*arg):
#     k_0, k_1, k_2, k_3 = arg
#     F = get_function_of_theta_for_sum([k_3], [-k_2],
#                                       [-k_0, -k_1, -k_3],
#                                       [k_0, k_1, k_2])
#     all_combinations = (k_3 + 1) * (k_2 + 1) * (k_3 + 1) * (k_2 + 1)
#     null_combinations = 0
#     non_trivial_zeros = 0
#     for v_theta in it.product(range(k_3 + 1), range(k_2 + 1),
#                               range(k_3 + 1), range(k_2 + 1)):
#         f = F(*v_theta)
#         if f.sum_of_absolute_values() != 0 and sum(v_theta) == 0:
#                 print 4 * "\n" + "something wrong!!!!!!!!!!"
#                 print inspect.stack()[0][3]
#                 print arg
#                 print v_theta
#
#         if f.sum_of_absolute_values() == 0:
#             null_combinations += 1
#             if sum(v_theta) != 0:
#                 if len(arg) == len(set(arg)) and len(set(v_theta)) > 1:
#                     non_trivial_zeros += 1
#                     # print "\nNontrivial zero"
#                     # print inspect.stack()[0][3]
#                     print arg
#                     print v_theta
#                     print
#     return non_trivial_zeros, null_combinations, all_combinations


def get_knot_descrption(*arg):
    description = ""
    for knot in arg:
        if knot[0] < 0:
            description += "-"
        description += "T("
        for k in knot:
            description += "2, " + str(abs(k)) + "; "
        description = description[:-2]
        description += ") # "
    return description[:-3]

def get_number_of_combinations(*arg):
    number_of_combinations = 1
    for knot in arg:
        number_of_combinations *= (2 * knot[-1] + 1)
    return number_of_combinations

def second_sum(*arg):
    k_0, k_1, k_2, k_3, k_4 = arg
    knot_sum = [[k_0, k_1, k_2], [k_3, k_4], [-k_0, -k_3, -k_4], [-k_1, -k_2]]
    F = get_function_of_theta_for_sum(*knot_sum)
    knot_description = get_knot_descrption(*knot_sum)
    all_combinations = get_number_of_combinations(*knot_sum)
    null_combinations = 1
    # non_trivial_zeros = 0

    for v_theta in it.product(range(k_2 + 1), range(k_4 + 1),
                              range(k_4 + 1), range(k_2 + 1)):
        f = F(*v_theta)
        assert f.sum_of_absolute_values() == 0 or sum(v_theta) != 0
        if f.sum_of_absolute_values() == 0 and sum(v_theta) != 0:
            null_combinations += 2
            # if len(arg) == len(set(arg)) and len(set(v_theta)) > 1:
            #     non_trivial_zeros += 1
            #     print "\nNontrivial zero"
            #     print inspect.stack()[0][3]
            #     print arg
            #     print v_theta
            #     print
    return knot_description, null_combinations, all_combinations


if __name__ == '__main__' and '__file__' in globals():
    main(sys.argv)
