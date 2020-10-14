#!/usr/bin/python
import numpy as np
import itertools as it
from typing import Iterable
from collections import Counter
from sage.arith.functions import LCM_list
import warnings
import re

SIGNATURE = 0
SIGMA = 1

# 9.11 (9.8)
# 9.15 (9.9)


class SignatureFunction(object):

    def __init__(self, values=None, counter=None):
        # builed counter based on values of signature jumps
        if counter is None:
            counter = Counter()

            if values is None:
                values = []
            else:
                msg = "Signature function is defined on the interval [0, 1)."
                assert all(k < 1 for k, v in values), msg

            for k, v in values:
                counter[k] += v
        self.cnt_signature_jumps = counter
        # self.tikz_plot("bum.tex")


    def is_zero_everywhere(self):
        return not any(self.cnt_signature_jumps.values())

    def double_cover(self):
        # to read values for t^2
        items = self.cnt_signature_jumps.items()
        counter = Counter({(1 + k) / 2 : v for k, v in items})
        counter.update(Counter({k / 2 : v for k, v in items}))
        return SignatureFunction(counter=counter)

    def square_root(self):
        # to read values for t^(1/2)
        counter = Counter()
        for jump_arg, jump in self.cnt_signature_jumps.items():
            if jump_arg < 1/2:
                counter[2 * jump_arg] = jump
        return SignatureFunction(counter=counter)

    def minus_square_root(self):
        # to read values for t^(1/2)
        items = self.cnt_signature_jumps.items()
        counter = Counter({mod_one(2 * k) : v for k, v in items if k >= 1/2})
        return SignatureFunction(counter=counter)

    def extremum(self):
        max = 0
        current = 0
        items = sorted(self.cnt_signature_jumps.items())
        for arg, jump in items:
            current += 2 * jump
            assert current == self(arg) + jump
            if abs(current) > abs(max):
                max = current
                # if abs(max) > 9:
                #     return max
        return max

    def __rshift__(self, shift):
        # A shift of the signature functions corresponds to the rotation.
        counter = Counter({mod_one(k + shift) : v \
                          for k, v in self.cnt_signature_jumps.items()})
        return SignatureFunction(counter=counter)

    def __lshift__(self, shift):
        return self.__rshift__(-shift)

    def __neg__(self):
        counter = Counter()
        counter.subtract(self.cnt_signature_jumps)
        return SignatureFunction(counter=counter)

    def __add__(self, other):
        counter = copy(self.cnt_signature_jumps)
        counter.update(other.cnt_signature_jumps)
        return SignatureFunction(counter=counter)

    def __sub__(self, other):
        counter = copy(self.cnt_signature_jumps)
        counter.subtract(other.cnt_signature_jumps)
        return SignatureFunction(counter=counter)

    def __eq__(self, other):
        self_cnt = Counter({k : v for k, v in self.cnt_signature_jumps.items()
                           if v != 0})
        other_cnt = Counter({k : v for k, v in other.cnt_signature_jumps.items()
                            if v != 0})
        return self.cnt_signature_jumps == other.cnt_signature_jumps


    def __str__(self):
        result = ''.join([str(jump_arg) + ": " + str(jump) + "\n"
                for jump_arg, jump in sorted(self.cnt_signature_jumps.items())
                if jump != 0])
        return result

    def __repr__(self):
        result = ''.join([str(jump_arg) + ": " + str(jump) + ", "
                for jump_arg, jump in sorted(self.cnt_signature_jumps.items())])
        return result[:-2] + "."

    def __call__(self, arg):
        # return the value of the signature function at the point arg, i.e.
        # sum of all signature jumps that occur before arg
        items = self.cnt_signature_jumps.items()
        result = [jump for jump_arg, jump in items if jump_arg < mod_one(arg)]
        return 2 * sum(result) + self.cnt_signature_jumps[arg]

    def total_sign_jump(self):
        # Total signature jump is the sum of all jumps.
        return sum([j[1] for j in self.to_list()])

    def to_list(self):
        # Return signature jumps formated as a list
        return sorted(self.cnt_signature_jumps.items())

    def step_function_data(self):
        # Transform the signature jump data to a format understandable
        # by the plot function.
        l = self.to_list()
        assert l == sorted(self.cnt_signature_jumps.items())
        vals = ([(d[0], sum(2 * j[1] for j in l[:l.index(d)+1])) for d in l] +
              [(0,self.cnt_signature_jumps[0]), (1,self.total_sign_jump())])
        print("step_function_data")
        print(vals)
        counter = copy(self.cnt_signature_jumps)
        counter[0] = self.cnt_signature_jumps[0]
        counter[1] = self.total_sign_jump()
        print(sorted(counter.items()))
        return vals

    def plot(self):
        # plot the signture function
        plot_step_function(self.step_function_data())

    def tikz_plot(self, file_name):
        # Draw the graph of the signature and transform it into TiKz.
        # header of the LaTeX file

        with open(file_name, "w") as f:
            f.write("\\documentclass[tikz]{standalone}\n")
            f.write("\\usetikzlibrary{datavisualization, " +
                        "datavisualization.formats.functions}\n")
            f.write("\\begin{document}\n")
            f.write("\\begin{tikzpicture}\n")
            data = sorted(self.step_function_data())
            print("data")
            print(data)
            f.write("\\datavisualization[scientific axes, " +
                    "visualize as smooth line,\n")
            f.write("x axis={ticks={none,major={at={")
            f.write(", " + str(N(data[0][0],digits=4)) + " as \\(" + \
                                     str(data[0][0]) + "\\)")
            for jump_arg, jump in data:
                f.write(", " + str(N(jump_arg,digits=4)) +
                        " as \\(" + str(jump_arg) + "\\)")
                f.write("}}}}\n")
                f.write("  ]\n")
                f.write("data [format=function]{\n")
                f.write("var x : interval [0:1];\n")
                f.write("func y = \\value x;\n")
                f.write("};\n")
                # close LaTeX enviroments
                f.write("\\end{tikzpicture}\n")
                f.write("\\end{document}\n")


class TorusCable(object):
    def __init__(self, knot_formula, k_vector=None, q_vector=None):

        self._knot_formula = knot_formula
        # q_i = 2 * k_i + 1
        if k_vector is not None:
            self.k_vector = k_vector
        elif q_vector is not None:
            self.q_vector = q_vector
        else:
            msg = "Please give a list of k (k_vector) or q values (q_vector)."
            raise ValueError(msg)

        self._sigma_function = None
        self._signature_as_function_of_theta = None


    # SIGMA & SIGNATURE
    @property
    def sigma_function(self):
        if self._sigma_function is None:
            self._sigma_function = self.get_sigma_function()
        return self._sigma_function

    @property
    def signature_as_function_of_theta(self):
        if self._signature_as_function_of_theta is None:
            self._signature_as_function_of_theta = \
                    self.get_signature_as_function_of_theta()
        return self._signature_as_function_of_theta

    # KNOT ENCODING
    @property
    def knot_formula(self):
        return self._knot_formula
    # @knot_formula.setter
    # def knot_formula(self, knot_formula):
    #     self._knot_formula = knot_formula

    @property
    def knot_description(self):
        return self._knot_description

    @property
    def knot_sum(self):
        return self._knot_sum
    @knot_sum.setter
    def knot_sum(self, knot_sum):
        self._knot_sum = knot_sum
        self._knot_description = self.get_knot_descrption(knot_sum)
        self._last_k_list = [abs(i[-1]) for i in knot_sum]
        self._last_q_list = [2 * i + 1 for i in self._last_k_list]
        if any(n not in Primes() for n in self._last_q_list):
            msg = "Incorrect q-vector. This implementation assumes that" + \
                  " all last q values are prime numbers.\n" + \
                  str(self._last_q_list)
            raise ValueError(msg)

        self.q_order = LCM_list(self._last_q_list)

    @property
    def last_k_list(self):
        return self._last_k_list

    @property
    def last_q_list(self):
        return self._last_q_list

    @property
    def q_order(self):
        return self._q_order
    @q_order.setter
    def q_order(self, val):
        self._q_order = val

    @property
    def k_vector(self):
        return self._k_vector
    @k_vector.setter
    def k_vector(self, k):
        self._k_vector = k
        if self.extract_max(self.knot_formula) > len(k) - 1:
            msg = "The vector for knot_formula evaluation is to short!"
            msg += "\nk_vector " + str(k) + " \nknot_formula " \
                + str(self.knot_formula)
            raise IndexError(msg)
        self.knot_sum = eval(self.knot_formula)
        self._q_vector = [2 * k_val + 1 for k_val in k]

    @property
    def q_vector(self):
        return self._q_vector
    @q_vector.setter
    def q_vector(self, new_q_vector):
        self.k_vector = [(q - 1)/2 for q in new_q_vector]


    def add_with_shift(self, other):
        # print("*" * 100)
        # print("BEFORE")
        # print(self.knot_description)
        # print(self.knot_sum)
        # print("*" * 100)
        # print("BEFORE k_vectors self, other")
        # print(self.k_vector)
        # print(other.k_vector)

        shift = len(self.k_vector)
        formula = re.sub(r'\d+', lambda x: str(int(x.group()) + shift),
                  other.knot_formula)

        knot_formula = self.knot_formula[:-1] + ",\n" + formula[1:]
        k_vector = self.k_vector + other.k_vector
        cable = TorusCable(knot_formula, k_vector=k_vector)
        s_signature_as_function_of_theta = self.signature_as_function_of_theta
        o_signature_as_function_of_theta = other.signature_as_function_of_theta

        shift = len(self.knot_sum)
        shift = len(self.knot_sum)
        def signature_as_function_of_theta(*thetas, **kwargs):
            result = s_signature_as_function_of_theta(*thetas[shift:]) + \
                     o_signature_as_function_of_theta(*thetas[0:shift])
            return result
        cable._signature_as_function_of_theta = signature_as_function_of_theta
        # print("*" * 100)
        # print("AFTER")
        # print(self.knot_description)
        # print(self.knot_formula)
        # print(self.knot_sum)
        # print("*" * 100)
        # print("AFTER k_vector, q_vector")
        # print(self.k_vector)
        # print(self.q_vector)
        return cable

    def __add__(self, other):
        if self.k_vector != other.k_vector:
            msg = "k_vectors are different. k-vector preserving addition is " +\
                  "impossible. The function add_with_shift was called instead"
            warnings.warn(msg)
        # print("*" * 100)
        # print("BEFORE")
        # print(self.knot_description)
        # print(self.knot_sum)
        # print("*" * 100)
        # print("BEFORE k_vectors self, other")

        knot_formula = self.knot_formula[:-1] + ",\n" + other.knot_formula[1:]
        cable = TorusCable(knot_formula, k_vector=self.k_vector)
        s_signature_as_function_of_theta = self.signature_as_function_of_theta
        o_signature_as_function_of_theta = other.signature_as_function_of_theta
        # print("FUNCTIONS ")
        # print(s_signature_as_function_of_theta([1,1,1,2]))
        # print(o_signature_as_function_of_theta([1,1,1,2]))
        # print("FUNCTIONS 1111")
        # print(s_signature_as_function_of_theta([1,1,1,1]))
        # print(o_signature_as_function_of_theta([1,1,1,1]))

        shift = len(self.knot_sum)
        def signature_as_function_of_theta(*thetas, **kwargs):
            result = s_signature_as_function_of_theta(*thetas[shift:]) + \
                     o_signature_as_function_of_theta(*thetas[0:shift])
            return result
        cable._signature_as_function_of_theta = signature_as_function_of_theta
        # print("*" * 100)
        # print("AFTER")
        # print(self.knot_description)
        # print(self.knot_formula)
        # print(self.knot_sum)
        # print("*" * 100)
        # print("AFTER k_vector, q_vector")
        # print(self.k_vector)
        # print(self.q_vector)
        return cable



    @staticmethod
    def extract_max(string):
        numbers = re.findall(r'\d+', string)
        numbers = map(int, numbers)
        return max(numbers)

    @staticmethod
    def get_blanchfield_for_pattern(k_n, theta):
        if theta == 0:
            sf = TorusCable.get_untwisted_signature_function(k_n)
            return sf.square_root() + sf.minus_square_root()

        results = []
        k = abs(k_n)
        ksi = 1/(2 * k + 1)

        counter = Counter()
        # print("lambda_odd, i.e. (theta + e) % 2 != 0")
        for e in range(1, k + 1):
            if (theta + e) % 2 != 0:
                counter[e * ksi] = 1 * sgn(k_n)
                counter[1 - e * ksi] = -1 * sgn(k_n)

                results.append((e * ksi, 1 * sgn(k_n)))
                results.append((1 - e * ksi, -1 * sgn(k_n)))

        # for example for k = 9 (q = 19) from this part we get
        # for even theta
        # 2/19: 1
        # 4/19: 1
        # 6/19: 1
        # 8/19: 1
        # 11/19: -1
        # 13/19: -1
        # 15/19: -1
        # 17/19: -1
        #
        # for odd theta
        # 1/19: 1
        # 3/19: 1
        # 5/19: 1
        # 7/19: 1
        # 9/19: 1
        # 10/19: -1
        # 12/19: -1
        # 14/19: -1
        # 16/19: -1
        # 18/19: -1

        # print("lambda_even")
        # print("normal")
        for e in range(1, theta):
            if (theta + e) % 2 == 0:
                results.append((e * ksi, 1 * sgn(k_n)))
                results.append((1 - e * ksi, -1 * sgn(k_n)))
        # print("reversed")
        for e in range(theta + 1, k + 1):
            if (theta + e) % 2 == 0:
                results.append((e * ksi, -1 * sgn(k_n)))
                results.append((1 - e * ksi, 1 * sgn(k_n)))

        return SignatureFunction(values=results)

    @staticmethod
    def get_untwisted_signature_function(j):
        # return the signature function of the T_{2,2k+1} torus knot
        k = abs(j)
        q = 2 * k + 1
        values = ([((2 * a + 1)/(2 * q), -1 * sgn(j)) for a in range(k)] +
             [((2 * a + 1)/(2 * q), 1 * sgn(j))
             for a in range(k + 1, 2 * k + 1)])
        return SignatureFunction(values=values)

    @staticmethod
    def get_knot_descrption(knot_sum):
        description = ""
        for knot in knot_sum:
            if knot[0] < 0:
                description += "-"
            description += "T("
            for k in knot:
                description += "2, " + str(2 * abs(k) + 1) + "; "
            description = description[:-2] + ") # "
        return description[:-3]


    def get_signature_as_function_of_theta(self, **key_args):
        if 'verbose' in key_args:
            verbose_default = key_args['verbose']
        else:
            verbose_default = False
        def signature_as_function_of_theta(*thetas, **kwargs):
            verbose = verbose_default
            if 'verbose' in kwargs:
                verbose = kwargs['verbose']
            len_a = len(self.knot_sum)
            len_t = len(thetas)
            # call with no arguments
            if len_t == 0:
                return signature_as_function_of_theta(*(len_a * [0]))
            if len_t != len_a:
                if isinstance(thetas, Iterable):
                    if len(thetas[0]) == len_a:
                        thetas = thetas[0]
                else:
                    msg = "This function takes exactly " + str(len_a) + \
                          " arguments or no argument at all (" + str(len_t) + \
                          " given)."
                    raise TypeError(msg)
            sf = SignatureFunction()
            untwisted_part = SignatureFunction()
            # for each cable knot in cable sum apply theta

            # print(self.knot_sum)

            for i, knot in enumerate(self.knot_sum):
                try:
                    ssf = self.get_summand_signature_as_theta_function(*knot)
                    plus, _, up = ssf(thetas[i])
                    # sf += ssf(thetas[i])
                    sf += plus
                    untwisted_part += up
                # in case wrong theata value was given
                except ValueError as e:
                    print("ValueError: " + str(e.args[0]) +\
                          " Please change " + str(i + 1) + ". parameter.")
                    return None
            # a = thetas[0]
            # # last_q = abs (2 * self.knot_sum[-1][-1]) + 1
            # if all(i == thetas[0] for i in thetas):
            #     print()
            #     print("\n" + "*" * 100)
            #     print(self.knot_description)
            #     print("one vector " + str(thetas))
            #     print("max sf " + str(sf.extremum()))
            #     print()
            #     # assert untwisted_part.is_zero_everywhere()

            if verbose:
                print()
                print(str(thetas))
                print(sf)
            msg = "tota signature jump = " + str(sf.total_sign_jump())
            msg += "\nfunction\n" + str(sf)
            assert sf.total_sign_jump() == 0, msg

            return sf
        signature_as_function_of_theta.__doc__ =\
                                        signature_as_function_of_theta_docstring
        return signature_as_function_of_theta


    def get_summand_signature_as_theta_function(self, *knot_as_k_values):
        def get_summand_signture_function(theta):
            # TBD: another formula (for t^2) description
            # TBD if theata condition
            k_n = knot_as_k_values[-1]
            if theta > 2 * abs(k_n):
                msg = "k for the pattern in the cable is " + str(k_n) + \
                      ". Parameter theta should not be larger than abs(k)."
                raise ValueError(msg)

            # twisted part
            cable_signature = self.get_blanchfield_for_pattern(k_n, theta)
            twisted_part = self.get_blanchfield_for_pattern(k_n, theta)
            untwisted_part = SignatureFunction()
            # untwisted part
            # for each knot summand consider k values in reversed order
            # ommit last k = k_n value

            ksi = 1/(2 * abs(k_n) + 1)
            for i, k in enumerate(knot_as_k_values[:-1][::-1]):
                power = 2^i
                a = TorusCable.get_untwisted_signature_function(k)
                shift = theta * ksi * power
                b = a >> shift
                c = a << shift
                for _ in range(i):
                    b = b.double_cover()
                    c = c.double_cover()
                cable_signature += b + c
                untwisted_part += b + c
            return cable_signature, twisted_part, untwisted_part
        get_summand_signture_function.__doc__ = \
            get_summand_signture_function_docsting
        return get_summand_signture_function


    def get_number_of_combinations_of_theta(self):
        number_of_combinations = 1
        for knot in self.knot_sum:
            number_of_combinations *= (2 * abs(knot[-1]) + 1)
        return number_of_combinations

    # searching for signature == 0
    def check_for_null_theta_combinations(self, verbose=False):
        list_of_good_vectors= []
        number_of_null_comb = 0
        f = self.signature_as_function_of_theta
        range_list = [range(abs(knot[-1]) + 1) for knot in self.knot_sum]
        for theta_vector in it.product(*range_list):
            if f(*theta_vector).is_zero_everywhere():
                list_of_good_vectors.append(theta_vector)
                m = len([theta for theta in theta_vector if theta != 0])
                number_of_null_comb += 2^m
        return number_of_null_comb, list_of_good_vectors

    # searching for signature or sigma > 5 + #(v_i != 0)
    def check_all_combinations_in_ranges(self, list_of_ranges,
                                            sigma_or_sign,
                                            print_results=False):
        all_combinations_pass = True
        all_bad_vectors = []
        number_of_all_good_v = 0
        for i, range_product in enumerate(list_of_ranges):
            good_v, bad_v = self.check_combinations_in_range(range_product,
                                                            sigma_or_sign)
            number_of_all_good_v += len(good_v)
            all_bad_vectors = list(it.chain(all_bad_vectors, bad_v))
            if bad_v:
                all_combinations_pass = False
            # if print_results:
            #     print("good : bad:\t " + str(len(good_v)) +\
            #           " : " + str(len(bad_v)))
            #     if i in [0, 4,]:
            #         print()
            #     if bad_v:
            #         print(bad_v)

        if print_results:
            print("good : bad:\t " + str(number_of_all_good_v) +\
                  " : " + str(len(all_bad_vectors)))
        return all_combinations_pass

    # searching for signature or sigma > 5 + #(v_i != 0)
    def check_combinations_in_range(self, range_product, sigma_or_sign):
        bad_vectors = []
        good_vectors = []
        last_q = self.q_vector[-1]

        for vector in range_product:
            a_1, a_2, a_3, a_4 = vector
            if (a_1^2 - a_2^2 + a_3^2 - a_4^2) % last_q:
                continue
            # if all(a in [1, last_q - 1] for a in vector):
            #     pass
            # else:
            #     continue
            if self.is_value_for_vector_class_big(vector, sigma_or_sign):
                good_vectors.append(vector)
            else:
                # print(vector)
                bad_vectors.append(vector)
        return good_vectors, bad_vectors


    def is_metaboliser(self, theta):
        i = 1
        sum = 0
        for idx, el in enumerate(theta):
            to_add = i * el^2
            # print("i * el^2 " + str(i * el^2))
            to_add /= self.last_q_list[idx]
            sum += to_add
            # print("i * el^2 % q_4:  " + str(to_add))
            # print("sum ", sum)
            i *= -1
            # if sum is integer
            #     continue
            # if all(a in [1, last_q - 1] for a in vector):
            #     pass
            # else:
            #     continue
        # print(theta, end=" ")
        # print(sum)
        if sum.is_integer():
            print("#" * 100)
            print(theta)
            return True
        return False
        #     if self.is_value_for_vector_class_big(vector, sigma_or_sign):
        #         good_vectors.append(vector)
        #     else:
        #         # print(vector)
        #         bad_vectors.append(vector)
        # return good_vectors, bad_vectors




    # searching for signature == 0
    def eval_cable_for_null_signature(self, print_results=False, verbose=False):
        # search for zero combinations
        number_of_all_comb = self.get_number_of_combinations_of_theta()
        result = self.check_for_null_theta_combinations(verbose=verbose)
        number_of_null_comb, list_of_good_vectors = result

        if print_results:
            print()
            print(self.knot_description)
            print("Zero cases: " + str(number_of_null_comb))
            print("All cases: " + str(number_of_all_comb))
            if list_of_good_vectors:
                print("Zero theta combinations: ")
                for el in list_of_good_vectors:
                    print(el)
        if number_of_null_comb^2 >= number_of_all_comb:
            return number_of_null_comb, number_of_all_comb
        return None

    # check sigma or signature function value
    # for all v = s * [a_1, a_2, a_3, a_4] for s in [1, last_q - 1]
    def is_value_for_vector_class_big(self, theta_vector, sigma_or_sign):
        [a_1, a_2, a_3, a_4] = theta_vector
        k_4 = self.knot_sum[-1][-1]
        q_4 = 2 * k_4 + 1

        max_sigma = 0

        if sigma_or_sign == SIGNATURE:
            f = self.signature_as_function_of_theta
        else:
            f = self.sigma_function
        # print(theta_vector, end="\t")

        for shift in range(1, k_4 + 1):
            shifted_theta = [(shift * a) % q_4 for a in
                             [a_1, a_2, a_3, a_4]]
            if sigma_or_sign == SIGNATURE:
                sf = f(shifted_theta)
                sig_v = sf.extremum()
            else:
                sig_v = f(shifted_theta)
            print(sig_v, end=" ")
            if abs(sig_v) > abs(max_sigma):
                max_sig = sig_v
            if abs(sig_v) > 5 + np.count_nonzero(shifted_theta):
                # print("\tok " + str(sigma_v))
                return True
        # print("\tbad class " + str(max_sigma))
        return False

    # searching for sigma > 5 + #(v_i != 0)
    def eval_cable_for_large_values(self, list_of_ranges,
                                    sigma_or_sign,
                                    print_results=False,
                                    verbose=False):
        if print_results:
            print(self.knot_description) # , end="\t\t\t")
        if self.check_all_combinations_in_ranges(list_of_ranges,
                                                sigma_or_sign,
                                                print_results=print_results):
            return True
        return False


##############################################################################
    # sigma function

    def get_sigma_function(self):
        if len(self.k_vector) != 4:
            msg = "This function is not implemented for k_vectors " +\
                  "with len other than 4."
            raise IndexError(msg)
        k_1, k_2, k_3, k_4 = [abs(k) for k in self.k_vector]
        last_q = 2 * k_4 + 1
        ksi = 1/last_q
        sigma_q_1 = self.get_untwisted_signature_function(k_1)
        sigma_q_2 = self.get_untwisted_signature_function(k_2)
        sigma_q_3 = self.get_untwisted_signature_function(k_3)

        def sigma_function(theta_vector, print_results=False):
            # "untwisted" part (Levine-Tristram signatures)
            a_1, a_2, a_3, a_4 = theta_vector
            untwisted_part = 2 * (sigma_q_2(ksi * a_1) -
                                  sigma_q_2(ksi * a_2) +
                                  sigma_q_3(ksi * a_3) -
                                  sigma_q_3(ksi * a_4) +
                                  sigma_q_1(ksi * a_1 * 2) -
                                  sigma_q_1(ksi * a_4 * 2))
            # "twisted" part
            tp = [0, 0, 0, 0]
            for i, a in enumerate(theta_vector):
                if a:
                    tp[i] = -last_q + 2 * a - 2 * (a^2/last_q)
            twisted_part = tp[0] - tp[1] + tp[2] - tp[3]
            # if print_results:
            #     self.print_results_LT(theta_vector, untwisted_part)
            #     self.print_results_LT(theta_vector, twisted_part)

            sigma_v = untwisted_part + twisted_part
            return sigma_v
        return sigma_function



    def print_results_LT(self, theta_vector, untwisted_part):
        knot_description = self.knot_description
        k_1, k_2, k_3, k_4 = [abs(k) for k in self.k_vector]
        a_1, a_2, a_3, a_4 = theta_vector
        last_q = 2 * k_4 + 1
        ksi = 1/last_q
        sigma_q_1 = self.get_untwisted_signature_function(k_1)
        sigma_q_2 = self.get_untwisted_signature_function(k_2)
        sigma_q_3 = self.get_untwisted_signature_function(k_3)
        print("\n\nLevine-Tristram signatures for the cable sum:  ")
        print(knot_description)
        print("and characters:\n" + str(theta_vector) + ",")
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

    def print_results_sigma(self, theta_vector, twisted_part):
        a_1, a_2, a_3, a_4 = theta_vector
        knot_description = self.knot_description
        last_q = self.q_vector[-1]
        print("\n\nSigma values for the cable sum: ")
        print(knot_description)
        print("and characters: " + str(v_theta))
        print("\nsigma(T_{2, q_4}, ksi_a) = " + \
              "-q + (2 * a * (q_4 - a)/q_4) " +\
              "= -q + 2 * a - 2 * a^2/q_4 if a != 0,\n\t\t\t" +\
              " = 0 if a == 0.")
        print("\nsigma(T_{2, q_4}, chi_a_1) = ", end="")
        if a_1:
            print("- (" + str(last_q) + ") + 2 * " + str(a_1) + " + " +\
                  "- 2 * " + str(a_1^2) + "/" + str(last_q) + \
                  " = " + str(tp[0]))
        else:
            print("0")
        print("\nsigma(T_{2, q_4}, chi_a_2) = ", end ="")
        if a_2:
            print("- (" + str(last_q) + ") + 2 * " + str(a_2) + " + " +\
                  "- 2 * " + str(a_2^2) + "/" + str(last_q) + \
                  " = " + str(tp[1]))
        else:
            print("0", end="")
        print("\nsigma(T_{2, q_4}, chi_a_3) = ", end="")
        if a_3:
            print("- (" + str(last_q) + ") + 2 * " + str(a_3) + " + " +\
                  "- 2 * " + str(a_3^2) + "/" + str(last_q) + \
                  " = " + str(tp[2]))
        else:
            print("0", end="")
        print("\nsigma(T_{2, q_4}, chi_a_4) = ", end="")
        if a_4:
            print("- (" + str(last_q) + ") + 2 * " + str(a_4) + " + " +\
                  "- 2 * " + str(a_4^2) + "/" + str(last_q) + \
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

    # def __tmp_print_all_sigma_for_vector_class(self, theta_vector):
    #     print("\n")
    #     print(self.knot_description)
    #     print("vector = " + str(theta_vector))
    #     [a_1, a_2, a_3, a_4] = theta_vector
    #     last_q = self.q_vector[-1]
    #     for shift in range(1, last_q):
    #         shifted_theta = [(shift * a) % last_q for a in
    #                          [a_1, a_2, a_3, a_4]]
    #         print(str(shifted_theta) + "\t\t" + \
    #                 str(self.__sigma_function(shifted_theta)))
    #     print("\n")
    #
    # def __tmp_get_max_sigma_for_vector_class(self, theta_vector):
    #     max_sigma = (theta_vector, 0)
    #     [a_1, a_2, a_3, a_4] = theta_vector
    #     last_q = self.q_vector[-1]
    #     for shift in range(1, last_q):
    #         shifted_theta = [(shift * a) % last_q for a in
    #                          [a_1, a_2, a_3, a_4]]
    #         sigma = self.__sigma_function(shifted_theta)
    #         if abs(sigma) > abs(max_sigma[1]):
    #             max_sigma = (shifted_theta, sigma)
    #     return max_sigma[1]
    #


def mod_one(n):
    return n - floor(n)

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

get_summand_signture_function_docsting = \
    """
    This function returns SignatureFunction for previously defined single
    cable T_(2, q) and a theta given as an argument.
    The cable was defined by calling function
    get_summand_signature_as_theta_function(*arg)
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

TorusCable.get_blanchfield_for_pattern.__doc__ = \
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

TorusCable.get_summand_signature_as_theta_function.__doc__ = \
    """
    Argument:
        n integers that encode a single cable, i.e.
        values of q_i for T(2,q_0; 2,q_1; ... 2, q_n)
    Return:
        a function that returns SignatureFunction for this single cable
        and a theta given as an argument
    """
