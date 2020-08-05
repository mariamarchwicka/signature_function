#!/usr/bin/python

# TBD: read about Factory Method, variable in docstring, sage documentation
# move settings to sep file

import os
import sys

import collections
import inspect
import itertools as it
import numpy as np
import re



class Config(object):
    def __init__(self):
        self.f_results = os.path.join(os.getcwd(), "results.out")

        # knot_formula is a schema for knots which signature function
        # will be calculated
        self.knot_formula = "[[k[0], k[1], k[3]], [-k[1], -k[3]], \
                             [k[2], k[3]], [-k[0], -k[2], -k[3]]]"

        # self.knot_formula = "[[k[0], k[1], k[2]], [k[3], k[4]], \
        #                      [-k[0], -k[3], -k[4]], [-k[1], -k[2]]]"
        # self.knot_formula = "[[k[0], k[1], k[2]], [k[3]],\
        #                          [-k[0], -k[1], -k[3]], [-k[2]]]"
        self.limit = 3

        self.verbose = True
        self.verbose = False

        self.print_calculations_for_small_signature = True
        self.print_calculations_for_small_signature = False


        self.print_calculations_for_large_signature = True
        self.print_calculations_for_large_signature = False


        # is the ratio restriction for values in k_vector taken into account
        # False flag is usefull to make quick script tests
        self.only_slice_candidates = True
        self.only_slice_candidates = False

        self.stop_after_firts_large_signature = True
        self.stop_after_firts_large_signature = False


class SignatureFunction(object):
    """
    This simple class encodes twisted and untwisted signature functions
    of knots. Since the signature function is entirely encoded by its signature
    jump, the class stores only information about signature jumps
    in a dictionary self.signature_jumps.
    The dictionary stores data of the signature jump as a key/values pair,
    where the key is the argument at which the functions jumps
    and value encodes the value of the jump. Remember that we treat
    signature functions as defined on the interval [0,1).
    """
    def __init__(self, values=None, counter=None):
        # set values of signature jumps
        if counter is None:
            counter = collections.Counter()
            if values is None:
                values = []
            for jump_arg, jump in values:
                assert 0 <= jump_arg < 1, \
                    "Signature function is defined on the interval [0, 1)."
                counter[jump_arg] = jump
        self.cnt_signature_jumps = counter
        self.signature_jumps = collections.defaultdict(int, counter)

    def sum_of_absolute_values(self):
        result = sum([abs(i) for i in self.signature_jumps.values()])
        test = sum([abs(i) for i in self.cnt_signature_jumps.values()])
        assert test == result
        return sum([abs(i) for i in self.cnt_signature_jumps.values()])

    def is_zero_everywhere(self):
        result = not any(self.signature_jumps.values())
        assert  result == (not any(self.cnt_signature_jumps.values()))
        if self.sum_of_absolute_values():
            assert result == False
        else:
            assert result == True
        return result

    def double_cover(self):
        # to read values for t^2
        new_data = []
        for jump_arg, jump in self.signature_jumps.items():
            new_data.append((jump_arg/2, jump))
            new_data.append((1/2 + jump_arg/2, jump))

        t_data = []
        for jump_arg, jump in self.cnt_signature_jumps.items():
            t_data.append((jump_arg/2, jump))
            t_data.append((1/2 + jump_arg/2, jump))

        sf = SignatureFunction(values=t_data)
        a = SignatureFunction(values=new_data)
        assert a == sf
        return sf

    def square_root(self):
        # to read values for t^(1/2)
        new_data = []
        for jump_arg, jump in self.signature_jumps.items():
            if jump_arg < 1/2:
                new_data.append((2 * jump_arg, jump))

        t_data = []
        for jump_arg, jump in self.cnt_signature_jumps.items():
            if jump_arg < 1/2:
                t_data.append((2 * jump_arg, jump))

        sf = SignatureFunction(values=t_data)
        a = SignatureFunction(values=new_data)
        assert a == sf
        return sf

    def minus_square_root(self):
        # to read values for t^(1/2)
        counter = collections.Counter()
        new_data = []
        for jump_arg, jump in self.cnt_signature_jumps.items():
            if jump_arg >= 1/2:
                counter[mod_one(2 * jump_arg)] = jump
                new_data.append((mod_one(2 * jump_arg), jump))
        t_data = []
        for jump_arg, jump in self.signature_jumps.items():
            if jump_arg >= 1/2:
                t_data.append((mod_one(2 * jump_arg), jump))
        print(t_data)
        a = SignatureFunction(values=t_data)
        sf = SignatureFunction(values=new_data)
        sf2 = SignatureFunction(counter=counter)
        print(new_data)
        print(counter.items())
        assert a == sf
        print("repr")
        print(repr(sf2))
        print(repr(a))
        assert a == sf2
        return sf

    def __lshift__(self, shift):
        # A shift of the signature functions corresponds to the rotation.
        return self.__rshift__(-shift)

    def __rshift__(self, shift):
        t_data = []
        for jump_arg, jump in self.signature_jumps.items():
            t_data.append((mod_one(jump_arg + shift), jump))
        new_data = []
        for jump_arg, jump in self.cnt_signature_jumps.items():
            new_data.append((mod_one(jump_arg + shift), jump))
        sf = SignatureFunction(values=new_data)
        a = SignatureFunction(values=t_data)
        assert a == sf
        return sf

    def __neg__(self):
        new_data = []
        print("neg")
        print("start values sign and cnt")
        print(self.signature_jumps.items())
        print(self.cnt_signature_jumps.items())
        for jump_arg, jump in self.signature_jumps.items():
            new_data.append((jump_arg, -jump))
        a = SignatureFunction(values=new_data)
        counter = collections.Counter()
        counter.subtract(self.cnt_signature_jumps)
        sf = SignatureFunction(counter=counter)
        assert a == sf
        return sf

    # TBD short
    def __add__(self, other):
        new_data = collections.defaultdict(int)
        for jump_arg, jump in other.signature_jumps.items():
            new_data[jump_arg] = jump + self.signature_jumps.get(jump_arg, 0)
        for jump_arg, jump in self.signature_jumps.items():
            if jump_arg not in new_data.keys():
                new_data[jump_arg] = self.signature_jumps[jump_arg]

        counter = copy(self.cnt_signature_jumps)
        counter.update(other.cnt_signature_jumps)
        assert collections.defaultdict(int, counter) == new_data
        return SignatureFunction(counter=counter)

    def __eq__(self, other):
        return self.cnt_signature_jumps == other.cnt_signature_jumps

    def __sub__(self, other):
        a = self + other.__neg__()
        counter = copy(self.cnt_signature_jumps)
        counter.subtract(other.cnt_signature_jumps)
        sf = SignatureFunction(counter=counter)
        assert a == sf
        return sf

    def __str__(self):
        result2 = ''.join([str(jump_arg) + ": " + str(jump) + "\n"
                for jump_arg, jump in sorted(self.signature_jumps.items())])
        result = ''.join([str(jump_arg) + ": " + str(jump) + "\n"
                for jump_arg, jump in sorted(self.cnt_signature_jumps.items())])
        assert result == result2

        return result

    def __repr__(self):
        result2 = ''.join([str(jump_arg) + ": " + str(jump) + ", "
                for jump_arg, jump in sorted(self.signature_jumps.items())])
        result = ''.join([str(jump_arg) + ": " + str(jump) + ", "
                for jump_arg, jump in sorted(self.cnt_signature_jumps.items())])


        assert result == result2

        return result[:-2] + "."

    def __call__(self, arg):
        # Compute the value of the signature function at the point arg.
        # This requires summing all signature jumps that occur before arg.
        arg = mod_one(arg)
        cnt = self.cnt_signature_jumps
        before_arg = [jump for jump_arg, jump in cnt.items() if jump_arg < arg]
        result = 2 * sum(before_arg) + cnt[arg]

        # TBD to delete
        val = 0
        for jump_arg, jump in self.signature_jumps.items():
          if jump_arg < arg:
              val += 2 * jump
          elif jump_arg == arg:
              val += jump
        assert result == val
        # end of to delete

        return result


def main(arg):
    try:
        new_limit = int(arg[1])
    except IndexError:
        new_limit = None
    search_for_large_signature_value(limit=new_limit)
    # search_for_null_signature_value(limit=new_limit)


# searching for signture > 5 + #(v_i != 0) over given knot schema
def search_for_large_signature_value(knot_formula=None,
                                     limit=None,
                                     verbose=None):
    if limit is None:
        limit = config.limit
    if knot_formula is None:
        knot_formula = config.knot_formula
    if verbose is None:
        vebose = config.verbose

    # number of k_i (q_i) variables to substitute
    k_vector_size = extract_max(knot_formula) + 1

    limit = max(limit, k_vector_size)
    combinations = it.combinations(range(1, limit + 1), k_vector_size)
    P = Primes()
    good_knots = []
    # with open(config.f_results, 'w') as f_results:
    for c in combinations:
        k = [(P.unrank(i) - 1)/2 for i in c]
        if config.only_slice_candidates:
            if not (k[3] > 4 * k[2] and
                    k[2] > 4 * k[1] and
                    k[1] > 4 * k[0]):
                if verbose:
                    print("Ratio-condition does not hold")
                continue
        result = eval_cable_for_large_signature(k_vector=k,
                                                knot_formula=knot_formula,
                                                print_results=False)
        good_knots.append(result)
    return good_knots



# searching for signture > 5 + #(v_i != 0)
def eval_cable_for_large_signature(k_vector=None,
                                   knot_formula=None,
                                   print_results=True,
                                   verbose=None,
                                   q_vector=None):
    if knot_formula is None:
        knot_formula = config.knot_formula
    if verbose is None:
        verbose = config.verbose
    if k_vector is None:
        if q_vector is None:
            # TBD docstring
            print("Please give a list of k (k_vector) or q values (q_vector).")
            return None
        else:
            k_vector = [(i - 1)/2 for i in q_vector]

    k = k_vector
    knot_sum = eval(knot_formula)
    knot_description = get_knot_descrption(*knot_sum)

    k_1, k_2, k_3, k_4 = [abs(i) for i in k]
    q_4 = 2 * k_4 + 1
    ksi = 1/q_4

    if verbose:
        print("\n\n")
        print(100 * "*")
        print("Searching for a large signature values for the cable sum: ")
        print(knot_description)

    if len(knot_sum) != 4:
        print("Wrong number of cable direct summands!")
        return None

    large_sigma_for_all_v_comninations = True
    good_knots = [("nic")]
    # iteration over all possible character combinations
    ranges_list = [range(abs(knot[-1]) + 1) for knot in knot_sum]
    for v_theta in it.product(*ranges_list):
        theta_squers = [i^2 for i in v_theta]
        condition = "(" + str(theta_squers[0]) \
                    + " - " + str(theta_squers[1]) \
                    + " + " + str(theta_squers[2]) \
                    + " - " + str(theta_squers[3]) \
                    + ") % " + str(q_4)
        # if verbose:
        #     print "\nChecking for characters: " + str(v_theta)
        if (theta_squers[0] - theta_squers[1] +
            theta_squers[2] - theta_squers[3]) % q_4:
            if verbose:
                print("The condition is not satisfied: " + \
                       str(condition) + " != 0.")
            continue
        if v_theta[0] == v_theta[1] == v_theta[2] == v_theta[3] == 0:
            print("\nSkip")
            continue
        if v_theta[0] == v_theta[1] == v_theta[2] == v_theta[3]:
            print("\nall v == a")


        # T(2, q_1; 2, q_2; 2, q_4) # -T(2, q_2; 2, q_4) #
        #         # T(2, q_3; 2, q_4) # -T(2, q_1; 2, q_3; 2, q_4)

        # "untwisted" part (Levine-Tristram signatures)
        sigma_q_1 = get_untwisted_signature_function(k_1)
        sigma_q_2 = get_untwisted_signature_function(k_2)
        sigma_q_3 = get_untwisted_signature_function(k_3)
        a_1, a_2, a_3, a_4 = v_theta
        untwisted_part = 2 * (sigma_q_2(ksi * a_1) +
                              sigma_q_1(ksi * a_1 * 2) -
                              sigma_q_2(ksi * a_2) +
                              sigma_q_3(ksi * a_3) -
                              sigma_q_3(ksi * a_4) -
                              sigma_q_1(ksi * a_4 * 2))

        # "twisted" part
        tp = [0, 0, 0, 0]
        for i, a in enumerate(v_theta):
            if a:
                tp[i] = -q_4 + 2 * a - 2 * (a^2/q_4)
        twisted_part = tp[0] - tp[1] + tp[2] - tp[3]
        assert twisted_part == int(twisted_part)

        sigma_v = untwisted_part + twisted_part
        if abs(sigma_v) > 5 + np.count_nonzero(v_theta):
            if config.print_calculations_for_large_signature:
                print("*" * 100)
                print("\n\nLarge signature value\n")
                print(knot_description)
                print("\nv_theta: ", end="")
                print(v_theta)
                print("k values: ", end="")
                print(str(k_1) + " " + str(k_2) + " " + \
                      str(k_3) + " " + str(k_4))
                print(condition)
                print("non zero value in v_theta: " + \
                      str(np.count_nonzero(v_theta)))
                print("sigma_v: " + str(sigma_v))
                print("\ntwisted_part: ", end="")
                print(twisted_part)
                print("untwisted_part: ", end="")
                print(untwisted_part)
                print("\n\nCALCULATIONS")
                print("*" * 100)
                print_results_LT(v_theta, knot_description,
                                 ksi, untwisted_part,
                                 k, sigma_q_1, sigma_q_2, sigma_q_3)
                print_results_sigma(v_theta, knot_description, tp, q_4)
                print("*" * 100 + "\n" * 5)
            else:
                print(knot_description + "\t" + str(v_theta) +\
                      "\t" + str(sigma_v) + "\t" + str(2 * sigma_q_1(2 * ksi * a_4)))
            if config.stop_after_firts_large_signature:
                break
        else:
            if config.print_calculations_for_small_signature:
                print("\n" * 5 + "*" * 100)
                print("\nSmall signature value\n")
                print(knot_description)
                print_results_LT(v_theta, knot_description, ksi, untwisted_part,
                                 k, sigma_q_1, sigma_q_2, sigma_q_3)
                print_results_sigma(v_theta, knot_description, tp, q_4)
                print("*" * 100 + "\n" * 5)
            else:
                print(knot_description + "\t" + str(v_theta) +\
                      "\t" + str(sigma_v) + "\t" + str(2 * sigma_q_1(2 * ksi * a_4)))


            large_sigma_for_all_v_comninations = False
            print("ojojojoj")
            break

    if large_sigma_for_all_v_comninations:
        print("\n\n\nHura hura")
        good_knots.append((knot_description, v_theta))

        # else:
        #     print "\n\tSmall signature value"
        #     print knot_description
        #     print "v_theta: " + str(v_theta)
        #     print condition
        #     print "non zero value in v_theta: " + str(np.count_nonzero(v_theta))
        #     print "signature at 1/2: " + str(y)
    return good_knots


def print_results_LT(v_theta, knot_description, ksi, untwisted_part,
                     k, sigma_q_1, sigma_q_2, sigma_q_3):
            a_1, a_2, a_3, a_4 = v_theta
            k_1, k_2, k_3, k_4 = [abs(i) for i in k]
            print("\n\nLevine-Tristram signatures for the cable sum:  ")
            print(knot_description)
            print("and characters:\n" + str(v_theta) + ",")
            print("ksi = " + str(ksi))
            print("\n\n2 * (sigma_q_2(ksi * a_1) + " + \
                    "sigma_q_1(ksi * a_1 * 2) - " +\
                    "sigma_q_2(ksi * a_2) + " +\
                    "sigma_q_3(ksi * a_3) - " +\
                    "sigma_q_3(ksi * a_4) - " +\
                    "sigma_q_1(ksi * a_4 * 2))" +\
                    \
                    " = \n\n2 * (sigma_q_2(" + \
                    str(ksi) + " * " + str(a_1) + \
                    ") + sigma_q_1(" + \
                    str(ksi) + " * " + str(a_1) + " * 2" + \
                    ") - sigma_q_2(" + \
                    str(ksi) + " * " + str(a_2) + \
                    ") + sigma_q_3(" + \
                    str(ksi) + " * " + str(a_3) + \
                    ") - sigma_q_3(" + \
                    str(ksi) + " * " + str(a_4) + \
                    ") - sigma_q_1(" + \
                    str(ksi) + " * " + str(a_4) + " * 2)) " + \
                    \
                    " = \n\n2 * (sigma_q_2(" + \
                    str(mod_one(ksi * a_1)) + \
                    ") + sigma_q_1(" + \
                    str(mod_one(ksi * a_1 * 2)) + \
                    ") - sigma_q_2(" + \
                    str(mod_one(ksi * a_2)) + \
                    ") + sigma_q_3(" + \
                    str(mod_one(ksi * a_3)) + \
                    ") - sigma_q_3(" + \
                    str(mod_one(ksi * a_4)) + \
                    ") - sigma_q_1(" + \
                    str(mod_one(ksi * a_4 * 2)) + \
                    \
                    ") = \n\n2 * ((" + \
                    str(sigma_q_2(ksi * a_1)) + \
                    ") + (" + \
                    str(sigma_q_1(ksi * a_1 * 2)) + \
                    ") - (" + \
                    str(sigma_q_2(ksi * a_2)) + \
                    ") + (" + \
                    str(sigma_q_3(ksi * a_3)) + \
                    ") - (" + \
                    str(sigma_q_3(ksi * a_4)) + \
                    ") - (" + \
                    str(sigma_q_1(ksi * a_4 * 2)) + ")) = " + \
                    "\n\n2 * (" + \
                    str(sigma_q_2(ksi * a_1) +
                    sigma_q_1(ksi * a_1 * 2) -
                    sigma_q_2(ksi * a_2) +
                    sigma_q_3(ksi * a_3) -
                    sigma_q_3(ksi * a_4) -
                    sigma_q_1(ksi * a_4 * 2)) + \
                    ") = " + str(untwisted_part))
            print("\nSignatures:")
            print("\nq_1 = " + str(2 * k_1 + 1) + ": " + repr(sigma_q_1))
            print("\nq_2 = " + str(2 * k_2 + 1) + ": " + repr(sigma_q_2))
            print("\nq_3 = " + str(2 * k_3 + 1) + ": " + repr(sigma_q_3))


def print_results_sigma(v_theta, knot_description, tp, q_4):
            a_1, a_2, a_3, a_4 = v_theta

            print("\n\nSigma values for the cable sum:  ")
            print(knot_description)
            print("and characters: " + str(v_theta))
            print("\nsigma(T_{2, q_4}, ksi_a) = " + \
                  "-q + (2 * a * (q_4 - a)/q_4) " +\
                  "= -q + 2 * a - 2 * a^2/q_4 if a != 0,\n\t\t\t" +\
                  " = 0 if a == 0.")
            print("\nsigma(T_{2, q_4}, chi_a_1) = ", end="")
            if a_1:
                print("- (" + str(q_4) + ") + 2 * " + str(a_1) + " + " +\
                      "- 2 * " + str(a_1^2) + "/" + str(q_4) + \
                      " = " + str(tp[0]))
            else:
                print("0")
            print("\nsigma(T_{2, q_4}, chi_a_2) = ", end ="")
            if a_2:
                print("- (" + str(q_4) + ") + 2 * " + str(a_2) + " + " +\
                      "- 2 * " + str(a_2^2) + "/" + str(q_4) + \
                      " = " + str(tp[1]))
            else:
                print("0", end="")
            print("\nsigma(T_{2, q_4}, chi_a_3) = ", end="")
            if a_3:
                print("- (" + str(q_4) + ") + 2 * " + str(a_3) + " + " +\
                      "- 2 * " + str(a_3^2) + "/" + str(q_4) + \
                      " = " + str(tp[2]))
            else:
                print("0", end="")
            print("\nsigma(T_{2, q_4}, chi_a_4) = ", end="")
            if a_4:
                print("- (" + str(q_4) + ") + 2 * " + str(a_4) + " + " +\
                      "- 2 * " + str(a_4^2) + "/" + str(q_4) + \
                      " = " + str(tp[3]))
            else:
                print("0")

            print("\n\nsigma(T_{2, q_4}, chi_a_1) " + \
                    "- sigma(T_{2, q_4}, chi_a_2) " + \
                    "+ sigma(T_{2, q_4}, chi_a_3) " + \
                    "- sigma(T_{2, q_4}, chi_a_4) =\n" + \
                    "sigma(T_{2, q_4}, " + str(a_1) + \
                    ") - sigma(T_{2, q_4}, " + str(a_2) + \
                    ") + sigma(T_{2, q_4}, " + str(a_3) + \
                    ") - sigma(T_{2, q_4}, " + str(a_4) + ") = " + \
                    str(tp[0] - tp[1] + tp[2] - tp[3]))

# searching for signature == 0
def search_for_null_signature_value(knot_formula=None, limit=None):
    if limit is None:
        limit = config.limit
    if knot_formula is None:
        knot_formula = config.knot_formula

    k_vector_size = extract_max(knot_formula) + 1
    combinations = it.combinations_with_replacement(range(1, limit + 1),
                                                    k_vector_size)

    with open(config.f_results, 'w') as f_results:

        for k in combinations:
            if config.only_slice_candidates and k_vector_size == 5:
                k = get_shifted_combination(k)
            knot_sum = eval(knot_formula)
            if is_trivial_combination(knot_sum):
                print(knot_sum)
                continue

            result = eval_cable_for_null_signature(knot_sum)
            if result is not None:
                knot_description, null_comb, all_comb = result
                line = (str(k) + ", " + str(null_comb) + ", " +
                        str(all_comb) + "\n")
                f_results.write(line)

# searching for signature == 0
def eval_cable_for_null_signature(knot_sum, print_results=False, verbose=None):
    # search for zero combinations
    if verbose is None:
        vebose = config.verbose
    f = get_signature_as_theta_function(*knot_sum, verbose=False)
    knot_description = get_knot_descrption(*knot_sum)
    all_combinations = get_number_of_combinations(*knot_sum)

    null_combinations = 0
    zero_theta_combinations = []

    ranges_list = [range(abs(knot[-1]) + 1) for knot in knot_sum]
    if verbose:
        print()
        print(knot_description)
    for v_theta in it.product(*ranges_list):
        if f(*v_theta, verbose=False).is_zero_everywhere():
            zero_theta_combinations.append(v_theta)
            m = len([theta for theta in v_theta if theta != 0])
            null_combinations += 2^m
        # else:
        #     assert sum(v_theta) != 0

    if print_results:
        print()
        print(knot_description)
        print("Zero cases: " + str(null_combinations))
        print("All cases: " + str(all_combinations))
        if zero_theta_combinations:
            print("Zero theta combinations: ")
            for el in zero_theta_combinations:
                print(el)
    if null_combinations^2 >= all_combinations:
        return knot_description, null_combinations, all_combinations
    return None


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


def get_signature_as_theta_function(*arg, **key_args):
    if 'verbose' in key_args:
        verbose_default = key_args['verbose']
    else:
        verbose_default = config.verbose
    def signature_as_theta_function(*thetas, **kwargs):
        verbose = verbose_default
        if 'verbose' in kwargs:
            verbose = kwargs['verbose']
        la = len(arg)
        lt = len(thetas)

        # call with no arguments
        if lt == 0:
            return signature_as_theta_function(*(la * [0]))

        if lt != la:
            msg = "This function takes exactly " + str(la) + \
                  " arguments or no argument at all (" + str(lt) + " given)."
            raise TypeError(msg)

        sf = SignatureFunction()

        # for each cable in cable sum apply theta
        for i, knot in enumerate(arg):
            try:
                sf += get_signature_summand_as_theta_function(*knot)(thetas[i])
            # in case wrong theata value was given
            except ValueError as e:
                print("ValueError: " + str(e.args[0]) +\
                      " Please change " + str(i + 1) + ". parameter.")
                return None
        if verbose:
            print()
            print(str(thetas))
            print(sf)
        return sf
    signature_as_theta_function.__doc__ = signature_as_theta_function_docstring
    return signature_as_theta_function


def get_number_of_combinations(*arg):
    number_of_combinations = 1
    for knot in arg:
        number_of_combinations *= (2 * abs(knot[-1]) + 1)
    return number_of_combinations


def extract_max(string):
    numbers = re.findall('\d+', string)
    numbers = map(int, numbers)
    return max(numbers)


def mod_one(n):
    return n - floor(n)


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

get_number_of_combinations.__doc__ = \
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

get_knot_descrption.__doc__ = \
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
    For each knot/cable sum the function eval_cable_for_null_signature is called.
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

eval_cable_for_null_signature.__doc__ = \
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
    The numbers given to the function eval_cable_for_null_signature are k-values for each
    component/cable in a direct sum.
    """

get_signature_as_theta_function.__doc__ = \
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

    sage: signature_function_generator = get_signature_as_theta_function(
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
    sage: print(get_signature_as_theta_function([1, 3], [2], [-1, -2], [-3]
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

signature_as_theta_function_docstring = \
    """
    Arguments:

    Returns object of SignatureFunction class for a previously defined
    connected sum of len(arg) cables.
    Acept len(arg) arguments: for each cable one theta parameter.
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

The numbers given to the function eval_cable_for_null_signature are k-values for each
component/cable in a direct sum.

To calculate signature function for a knot and a theta value, use function
get_signature_as_theta_function (see help/docstring for details).

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
