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
            self.q_vector = self.get_q_vector(self.knot_formula)

        self._sigma_function = None
        self._signature_as_function_of_theta = None


    @property
    def signature_as_function_of_theta(self):
        if self._signature_as_function_of_theta is None:
            self._signature_as_function_of_theta = \
                    self.get_signature_as_function_of_theta()
        return self._signature_as_function_of_theta

    # knot encoding
    @property
    def knot_formula(self):
        return self._knot_formula
    # @knot_formula.setter
    # def knot_formula(self, knot_formula):
    #     self._knot_formula = knot_formula

    # knot encoding
    @property
    def knot_description(self):
        return self._knot_description

    # knot encoding
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


    def get_q_vector(knot_formula, slice=True):
        lowest_number = 2
        q_vector = [0] * (TorusCable.extract_max(knot_formula) + 1)
        P = Primes()
        for layer in TorusCable.get_layers_from_formula(knot_formula)[::-1]:
            for el in layer:
                q_vector[el] = P.next(lowest_number)
                lowest_number = q_vector[el]
            lowest_number *= 4
        return q_vector

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

    @staticmethod
    def get_layers_from_formula(knot_formula):
        layers = []
        k_indices = re.sub(r'k', '', knot_formula)
        k_indices = re.sub(r'-', '', k_indices)
        k_indices = re.sub(r'\n', '', k_indices)
        k_indices = re.sub(r'\[\d+\]', lambda x: x.group()[1:-1], k_indices)
        k_indices = eval(k_indices)
        number_of_layers = max(len(lst) for lst in k_indices)
        print(k_indices)
        layers = []
        for i in range(1, number_of_layers + 1):
            layer = set()
            for lst in k_indices:
                if len(lst) >= i:
                    layer.add(lst[-i])
            layers.append(layer)
        return layers


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

    def is_metabolizer(self, theta):
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
            # print("#" * 100)
            # print(theta)
            return True
        return False
        #     if self.is_value_for_vector_class_big(vector, sigma_or_sign):
        #         good_vectors.append(vector)
        #     else:
        #         # print(vector)
        #         bad_vectors.append(vector)
        # return good_vectors, bad_vectors


    def is_signature_big_in_ranges(self, ranges_list):
        for theta in it.product(*ranges_list):
            if not any(theta):
                continue
            we_have_a_problem = True
            if self.is_metabolizer(theta):
                for shift in range(1, self.q_order):
                    shifted_theta = [(shift * th) % self.last_q_list[i]
                                     for i, th in enumerate(theta)]
                    shifted_theta = [min(th, self.last_q_list[i] - th)
                                     for i, th in enumerate(shifted_theta)]
                    sf = self.signature_as_function_of_theta(*shifted_theta)
                    extremum = abs(sf.extremum())
                    if shift > 1:
                        print(shifted_theta, end=" ")
                        print(extremum)
                    if extremum > 5 + np.count_nonzero(shifted_theta):
                        # print("ok")
                        we_have_a_problem = False
                        break
                    elif shift == 1:
                        print("*" * 10)
                        print(shifted_theta, end=" ")
                        print(extremum)
                if we_have_a_problem:
                    print("\n" * 10 + "!" * 1000)
                    return False
        return True

    def is_signature_big_for_all_metabolizers(self):
        if len(self.knot_sum) == 8:
            for shift in range(0, 8, 4):
                ranges_list = 8 * [range(0, 1)]
                ranges_list[shift : shift + 3] = [range(0, i + 1) for i in \
                                             self.last_k_list[shift: shift + 3]]
                ranges_list[shift + 3] = range(0, 2)
                if not self.is_signature_big_in_ranges(ranges_list):
                    return False
                else:
                    print("\n\nok")
            return True

        elif len(self.knot_sum) == 4:
            print("\n\n\nhohohohohoho")
            upper_bounds = self.last_k_list[:3]
            ranges_list = [range(0, i + 1) for i in upper_bounds]
            ranges_list.append(range(0, 2))
            if not self.is_signature_big_in_ranges(ranges_list):
                return False
            return True

        msg = "Function implemented only for knots with 4 or 8 summands"
        raise ValueError(msg)


def mod_one(n):
    return n - floor(n)


TorusCable.get_knot_descrption.__doc__ = \
    """
    Arguments:
        arbitrary number of lists of numbers, each list encodes a single cable.
    Examples:
        sage: get_knot_descrption([1, 3], [2], [-1, -2], [-3])
        'T(2, 3; 2, 7) # T(2, 5) # -T(2, 3; 2, 5) # -T(2, 7)'
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
