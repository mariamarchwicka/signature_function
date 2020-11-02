#!/usr/bin/python
import numpy as np
import itertools as it
from typing import Iterable
from collections import Counter
from sage.arith.functions import LCM_list
import warnings
import re
import matplotlib.pyplot as plt
import inspect

# 9.11 (9.8)
# 9.15 (9.9)


class SignatureFunction():

    def __init__(self, values=None, counter=None):

        # counter of signature jumps
        if counter is None:
            counter = Counter()
            if values is None:
                values = []
            for k, v in values:
                counter[k] += v

        counter = Counter({k : v for k, v in counter.items() if v != 0})
        if any(k >= 1 for k in counter.keys()):
                msg = "Signature function is defined on the interval [0, 1)."
                raise ValueError(msg)

        counter[0] += 0
        counter[1] += 0
        self.jumps_counter = counter

    def __rshift__(self, shift):
        # A shift of the signature functions corresponds to the rotation.
        counter = Counter({mod_one(k + shift) : v \
                          for k, v in self.jumps_counter.items()})
        return SignatureFunction(counter=counter)

    def __lshift__(self, shift):
        return self.__rshift__(-shift)

    def __neg__(self):
        counter = Counter()
        counter.subtract(self.jumps_counter)
        return SignatureFunction(counter=counter)

    def __add__(self, other):
        counter = copy(self.jumps_counter)
        counter.update(other.jumps_counter)
        return SignatureFunction(counter=counter)

    def __sub__(self, other):
        counter = copy(self.jumps_counter)
        counter.subtract(other.jumps_counter)
        return SignatureFunction(counter=counter)

    def __eq__(self, other):
        return self.jumps_counter == other.jumps_counter

    def __str__(self):
        result = ''.join([str(jump_arg) + ": " + str(jump) + "\n"
                for jump_arg, jump in sorted(self.jumps_counter.items())])
        return result

    def __repr__(self):
        result = ''.join([str(jump_arg) + ": " + str(jump) + ", "
                for jump_arg, jump in sorted(self.jumps_counter.items())])
        return result[:-2] + "."

    def __call__(self, arg):
        # return the value of the signature function at the point arg, i.e.
        # sum of all signature jumps that occur before arg
        items = self.jumps_counter.items()
        result = [jump for jump_arg, jump in items if jump_arg < mod_one(arg)]
        return 2 * sum(result) + self.jumps_counter[arg]

    def is_zero_everywhere(self):
        return not any(self.jumps_counter.values())

    def double_cover(self):
        # to read values for t^2
        items = self.jumps_counter.items()
        counter = Counter({(1 + k) / 2 : v for k, v in items})
        counter.update(Counter({k / 2 : v for k, v in items}))
        return SignatureFunction(counter=counter)

    def square_root(self):
        # to read values for t^(1/2)
        counter = Counter()
        for jump_arg, jump in self.jumps_counter.items():
            if jump_arg < 1/2:
                counter[2 * jump_arg] = jump
        return SignatureFunction(counter=counter)

    def minus_square_root(self):
        # to read values for t^(1/2)
        items = self.jumps_counter.items()
        counter = Counter({mod_one(2 * k) : v for k, v in items if k >= 1/2})
        return SignatureFunction(counter=counter)

    def extremum(self, limit=None):
        max = 0
        current = 0
        items = sorted(self.jumps_counter.items())
        for arg, jump in items:
            current += 2 * jump
            assert current == self(arg) + jump
            if abs(current) > abs(max):
                max = current
                if limit is not None:
                    if abs(max) > limit:
                        break
        return max

    def total_sign_jump(self):
        # Total signature jump is the sum of all jumps.
        return sum([j[1] for j in sorted(self.jumps_counter.items())])

class SignatureWriter():

    def __init__(self, signature_function):
        self.sf = signature_function

    def plot(self, title=None, subplot=False):

        keys = sorted(self.sf.jumps_counter.keys())
        y = [self.sf(k) + self.sf.jumps_counter[k] for k in keys]
        xmax = [k for k in keys if k != 0]
        xmin = [k for k in keys if k != 1]
        fig, ax = plt.subplots(1, 1)
        ax.set(ylabel='signature function')
        if title is not None:
            ax.set(title=title)
        ax.hlines(y, xmin, xmax, color='blue')
        plt.savefig('sf.png')
        plt.close()


        from PIL import Image
        image = Image.open('sf.png')
        image.show()

    def step_function_data(self):
        # Transform the signature jump data to a format understandable
        # by the plot function.
        result = [(k, self.sf(k) + self.sf.jumps_counter[k])
                 for k in sorted(self.sf.jumps_counter.keys())]
        return result

    def tikz_plot(self, file_name):
        plt_sin = plot(sin(x), (x, 0, 2*pi))
        # plt_sin.show()
        plt_sin.save("MyPic.pdf")

        return
        # Draw the graph of the signature and transform it into TiKz.
        # header of the LaTeX file
        head = inspect.cleandoc(
            r"""
            \documentclass{standalone}
            \usepackage{tikz}
            \usetikzlibrary{calc}
            \begin{document}
            \begin{tikzpicture}
            """)

        body = \
            r"""
            %A piecewise linear function is drawn over the interval.
            \draw (5,0) -- (6,-4);
            %The axes are drawn.
            \draw[latex-latex] ($(0,{-4*(2/5)}) +(0pt,-12.5pt)$) --
            ($(0,{4*(2/5)}) +(0pt,12.5pt)$) node[above right]{$y$};
            \draw[latex-latex] ($({-4*(2/5)},0) +(-12.5pt,0pt)$) --
            ($({12*(2/5)},0) +(12.5pt,0pt)$) node[below right]{$x$};
            """
        tail = \
            r"""
            \end{tikzpicture}
            \end{document}
            """
        tikzpicture = re.sub(r' +', ' ', ''.join([head, body, tail]))
        tikzpicture = re.sub(r'\n ', '\n', tikzpicture)

        with open("tmp.tex", "w") as f:
            f.write(tikzpicture)

        data = self.step_function_data()
        with open(file_name, "w") as f:
            head = \
                r"""
                \documentclass[tikz]{{standalone}}
                %\usepackage{{tikz}}
                \usetikzlibrary{{datavisualization}}
                \usetikzlibrary{{datavisualization.formats.functions}}
                %\usetikzlibrary{{calc}}
                \begin{{document}}
                \begin{{tikzpicture}}
                \datavisualization[scientific axes, visualize as smooth line,
                x axis={{ticks={{none,major={{at={{, {arg0} " as \\( {val0} \\
                %]
                """.format(arg0=str(N(data[0][0] ,digits=4)), val0=str(data[0][0]))
            f.write(head)


            # f.write(", " + str(N(data[0][0],digits=4)) + " as \\(" + \
            #                          str(data[0][0]) + "\\)")
            for jump_arg, jump  in data[1:3]:
                f.write(", " + str(N(jump_arg,digits=4)) +
                        " as \\(" + str(jump_arg) + "\\)")
            f.write("}}}}\n")
            f.write("  ]\n")
            f.write("data [format=function]{\n")
            f.write("var x : interval [0:1];\n")
            f.write("func y = \\value x;\n")
            f.write("};\n")
                # close LaTeX enviroments
            tail = \
                r"""
                %};
                \end{tikzpicture}
                \end{document}
                """
            f.write(tail)


class CableSummand():
    pass





class CableSum():
    def __init__(self, knot_formula, k_vector=None, q_vector=None):

        self._knot_formula = knot_formula
        # q_i = 2 * k_i + 1
        if k_vector is not None:
            self.k_vector = k_vector
        elif q_vector is not None:
            self.q_vector = q_vector
        else:
            self.q_vector = self.get_q_vector_alg_slice(self.knot_formula)

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
        self._patt_k_list = [abs(i[-1]) for i in knot_sum]
        self._patt_q_list = [2 * i + 1 for i in self._patt_k_list]
        if any(n not in Primes() for n in self._patt_q_list):
            msg = "Incorrect q-vector. This implementation assumes that" + \
                  " all last q values are prime numbers.\n" + \
                  str(self._patt_q_list)
            raise ValueError(msg)

        self.q_order = LCM_list(self._patt_q_list)

    @property
    def patt_k_list(self):
        return self._patt_k_list

    @property
    def patt_q_list(self):
        return self._patt_q_list

    # q_order is LCM of all q values for pattern knots
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


    def __add__(self, other):
        s_formula = self.knot_formula
        o_formula = other.knot_formula
        k_vector = self.k_vector

        if self.k_vector != other.k_vector:
            msg = "k_vectors are different. k-vector preserving addition is " +\
                  "impossible."
            warnings.warn(msg)
            shift = len(self.k_vector)
            o_formula = re.sub(r'\d+', lambda x: str(int(x.group()) + shift),
                               o_formula)
            k_vector += other.k_vector

        knot_formula = s_formula[:-1] + ",\n" + o_formula[1:]
        cable = CableSum(knot_formula, k_vector=k_vector)
        s_sig = self.signature_as_function_of_theta
        o_sig = other.signature_as_function_of_theta
        shift = len(self.knot_sum)

        def signature_as_function_of_theta(*thetas, **kwargs):
            thetas = cable.parse_thetas(*thetas)
            result = s_sig(*thetas[shift:]) + o_sig(*thetas[0:shift])
            return result

        cable._signature_as_function_of_theta = signature_as_function_of_theta
        return cable

    def parse_thetas(self, *thetas):
        summands_num = len(self.knot_sum)
        if not thetas:
            return summands_num * (0,)
        if len(thetas) == 1 and summands_num > 1:
            if isinstance(thetas[0], Iterable):
                if len(thetas[0]) >= summands_num:
                    return tuple(thetas[0])
                elif not thetas[0]:
                    return summands_num * (0,)
            elif thetas[0] == 0:
                return summands_num * (0,)
            else:
                msg = "This function takes at least " + str(summands_num) + \
                      " arguments or no argument at all (" + str(len(thetas)) \
                      + " given)."
                raise TypeError(msg)
        return tuple(thetas)

    @staticmethod
    def get_q_vector_alg_slice(knot_formula):
        lowest_number = 2
        q_vector = [0] * (CableSum.extract_max(knot_formula) + 1)
        P = Primes()
        for layer in CableSum.get_layers_from_formula(knot_formula)[::-1]:
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
    def get_blanchfield_for_pattern(k_n, theta=0):

        msg = "Theorem on which this function is based, assumes " +\
              "theta < k, where q = 2*k + 1 for pattern knot T(p, q)."
        if theta == 0:
            sf = CableSum.get_untwisted_signature_function(k_n)
            return sf.square_root() + sf.minus_square_root()

        k = abs(k_n)
        assert theta <= k, msg
        results = []

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
        # return the signature function of the T_{2, 2k+1} torus knot
        k = abs(j)
        q = 2 * k + 1
        counter = Counter({(2 * a + 1)/(2 * q) : -sgn(j)
                           for a in range(k)})
        counter.update(Counter({(2 * a + 1)/(2 * q) : sgn(j)
                           for a in range(k + 1, q)}))
        return SignatureFunction(counter=counter)


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
        k_indices = re.sub(r'[k-]', '', knot_formula)
        k_indices = re.sub(r'\[\d+\]', lambda x: x.group()[1:-1], k_indices)
        k_indices = eval(k_indices)
        number_of_layers = max(len(lst) for lst in k_indices)
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
        knot_desc = self.knot_description

        def signature_as_function_of_theta(*thetas, **kwargs):
            # print("\n\nsignature_as_function_of_theta " + knot_desc)
            verbose = verbose_default
            if 'verbose' in kwargs:
                verbose = kwargs['verbose']
            thetas = self.parse_thetas(*thetas)

            untwisted_part = SignatureFunction()
            twisted_part = SignatureFunction()

            # for each cable knot (summand) in cable sum apply theta
            for i, knot in enumerate(self.knot_sum):
                ssf = self.get_summand_signature_as_theta_function(*knot)
                tp, up = ssf(thetas[i])
                twisted_part += tp
                untwisted_part += up
            sf = twisted_part + untwisted_part

            if verbose:
                print()
                print(str(thetas))
                print(sf)
            assert sf.total_sign_jump() == 0
            return sf

        signature_as_function_of_theta.__doc__ =\
                                        signature_as_function_of_theta_docstring
        return signature_as_function_of_theta

    def get_untwisted_part(self, *knot_as_k_values, theta=0):
        patt_k = knot_as_k_values[-1]
        ksi = 1/(2 * abs(patt_k) + 1)

        untwisted_part = SignatureFunction()
        # For each knot summand consider k values in reversed order,
        # ommit k value for pattern.
        for layer_num, k in enumerate(knot_as_k_values[:-1][::-1]):
            sf = CableSum.get_untwisted_signature_function(k)
            shift = theta * ksi * 2^layer_num
            right_shift = sf >> shift
            left__shift = sf << shift
            for _ in range(layer_num):
                right_shift = right_shift.double_cover()
                left__shift = left__shift.double_cover()
            untwisted_part += right_shift + left__shift
        return untwisted_part

    def get_summand_signature_as_theta_function(self, *knot_as_k_values):

        def get_summand_signture_function(theta):

            patt_k = knot_as_k_values[-1]

            # theta should not be larger than k for the pattern.
            theta %= (2 * abs(patt_k) + 1)
            theta = min(theta, 2 * abs(patt_k) + 1 - theta)

            twisted_part = self.get_blanchfield_for_pattern(patt_k, theta)
            untwisted_part = self.get_untwisted_part(*knot_as_k_values,
                                                     theta=theta)
            return twisted_part, untwisted_part
        get_summand_signture_function.__doc__ = \
            get_summand_signture_function_docsting

        return get_summand_signture_function

    def is_metabolizer(self, theta):
        # Check if square alternating difference
        # divided by last q value is integer.
        result = sum(el^2 / self.patt_q_list[idx] * (-1)^idx
                      for idx, el in enumerate(theta))
        # for idx, el in enumerate(theta):
        #     old_sum += (el^2 / self.patt_q_list[idx] * (-1)^idx)
        return result.is_integer()

    def is_signature_big_in_ranges(self, ranges_list):

        for thetas in it.product(*ranges_list):

            # Check only non-zero metabolizers.
            if not self.is_metabolizer(thetas) or not any(thetas):
                continue

            signature_is_small = True
            # Check if any element generated by thetas vector
            # has a large signature.
            for shift in range(1, self.q_order):
                shifted_thetas = [shift * th for th in thetas]
                sf = self.signature_as_function_of_theta(*shifted_thetas)
                limit = 5 + np.count_nonzero(shifted_thetas)
                extremum = abs(sf.extremum(limit=limit))
                if shift > 1:
                    print(shifted_thetas, end=" ")
                    print(extremum)
                if extremum > limit:
                    signature_is_small = False
                    break
                elif shift == 1:
                    print("*" * 10)
                    print(shifted_thetas, end=" ")
                    print(extremum)
            if signature_is_small:
                print("\n" * 10 + "!" * 1000)
                return False
        return True

    def is_signature_big_for_all_metabolizers(self):
        num_of_summands = len(self.knot_sum)
        if num_of_summands % 4:
            f_name = self.is_signature_big_for_all_metabolizers.__name__
            msg = "Function {}".format(f_name) + " is implemented only for " +\
                  "knots that are direct sums of 4n direct summands."
            raise ValueError(msg)

        for shift in range(0, num_of_summands, 4):
            ranges_list = num_of_summands * [range(0, 1)]
            ranges_list[shift : shift + 3] = \
                [range(0, i + 1) for i in self.patt_k_list[shift: shift + 3]]
            ranges_list[shift + 3] = range(0, 2)
            if not self.is_signature_big_in_ranges(ranges_list):
                return False
            else:
                print("\nOK")
        return True



def mod_one(n):
    return n - floor(n)


# CableSum.get_knot_descrption.__doc__ = \
#     """
#     Arguments:
#         arbitrary number of lists of numbers, each list encodes a single cable.
#     Examples:
#         sage: get_knot_descrption([1, 3], [2], [-1, -2], [-3])
#         'T(2, 3; 2, 7) # T(2, 5) # -T(2, 3; 2, 5) # -T(2, 7)'
#     """

CableSum.get_signature_as_function_of_theta.__doc__ = \
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
    in a dictionary self.jumps_counter.
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

CableSum.get_blanchfield_for_pattern.__doc__ = \
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

CableSum.get_summand_signature_as_theta_function.__doc__ = \
    """
    Argument:
        n integers that encode a single cable, i.e.
        values of q_i for T(2,q_0; 2,q_1; ... 2, q_n)
    Return:
        a function that returns SignatureFunction for this single cable
        and a theta given as an argument
    """
