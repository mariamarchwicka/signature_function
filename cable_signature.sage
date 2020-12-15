#!/usr/bin/env sage -python

import numpy as np
import itertools as it
import warnings
import re
from typing import Iterable
from collections import Counter
from sage.arith.functions import LCM_list
import importlib


SIGMA = 0
SIGNATURE = 1

def import_sage(module_name):

    sage_name = module_name + ".sage"
    python_name = module_name + ".sage.py"

    if os.path.isfile(sage_name):
        os.system('sage --preparse {}'.format(sage_name));
        os.system('mv {} {}.py'.format(python_name, module_name))

    if module_name in sys.modules:
        return importlib.reload(sys.modules[module_name])
    return importlib.import_module(module_name, package=None)
sig = import_sage('signature')


# #############################################################################
# 9.11 (9.8)
# 9.15 (9.9)
PLOTS_DIR = "plots"

class CableSummand:


    def __init__(self, knot_as_k_values, verbose=False):

        self.verbose = verbose
        self.knot_as_k_values = knot_as_k_values
        self.knot_description = self.get_summand_descrption(knot_as_k_values)
        self.signature_as_function_of_theta = \
                                self.get_summand_signature_as_theta_function()
        self.sigma_as_function_of_theta = self.get_sigma_as_function_of_theta()

    @staticmethod
    def get_summand_descrption(knot_as_k_values):
        description = ""
        if knot_as_k_values[0] < 0:
            description += "-"
        description += "T("
        for k in knot_as_k_values:
            description += "2, " + str(2 * abs(k) + 1) + "; "
        return description[:-2] + ")"

    @classmethod
    def get_blanchfield_for_pattern(cls, k_n, theta=0):

        msg = "Theorem on which this function is based, assumes " +\
              "theta < k, where q = 2*k + 1 for pattern knot T(p, q)."
        if theta == 0:
            sf = cls.get_untwisted_signature_function(k_n)
            return sf.square_root() + sf.minus_square_root()

        k = abs(k_n)
        assert theta <= k, msg

        results = []
        ksi = 1/(2 * k + 1)

        # print("lambda_odd, i.e. (theta + e) % 2 != 0")
        for e in range(1, k + 1):
            if (theta + e) % 2 != 0:
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

        return sig.SignatureFunction(values=results)

    @classmethod
    def get_satellite_part(cls, *knot_as_k_values, theta=0):
        patt_k = knot_as_k_values[-1]
        ksi = 1/(2 * abs(patt_k) + 1)

        satellite_part = sig.SignatureFunction()
        # For each knot summand consider k values in reversed order,
        # ommit k value for pattern.
        for layer_num, k in enumerate(knot_as_k_values[:-1][::-1]):
            sf = cls.get_untwisted_signature_function(k)
            shift = theta * ksi * 2^layer_num
            right_shift = sf >> shift
            left__shift = sf << shift
            for _ in range(layer_num):
                right_shift = right_shift.double_cover()
                left__shift = left__shift.double_cover()
            satellite_part += right_shift + left__shift
        return satellite_part

    @staticmethod
    def get_untwisted_signature_function(k=None, q=None):
        # return the signature function of the T_{2, 2k+1} torus knot

        if q is not None:
            signum = sign(q)
            q = abs(q)
            k = (q - 1)/2
        elif k is not None:
            signum = sign(k)
            k = abs(k)
            q = 2 * k + 1
        else:
            raise ValueError('k or q value must be given')

        counter = Counter({(2 * a + 1)/(2 * q) : -signum
                           for a in range(k)})
        counter.update(Counter({(2 * a + 1)/(2 * q) : signum
                           for a in range(k + 1, q)}))
        return sig.SignatureFunction(counter=counter)

    def get_summand_signature_as_theta_function(self):
        # knot_as_k_values = self.knot_as_k_values
        def get_summand_signture_function(theta):

            patt_k = self.knot_as_k_values[-1]

            # theta should not be larger than k for the pattern.
            theta %= (2 * abs(patt_k) + 1)
            theta = min(theta, 2 * abs(patt_k) + 1 - theta)

            pattern_part = self.get_blanchfield_for_pattern(patt_k, theta)
            satellite_part = self.get_satellite_part(*self.knot_as_k_values,
                                                     theta=theta)
            sf = satellite_part + pattern_part

            satellite_part.plot_title = self.knot_description + \
                                        ", theta = " + str(theta) + \
                                         ", satellite part."
            pattern_part.plot_title = self.knot_description + \
                                      ", theta = " + str(theta) + \
                                      ", pattern part."
            sf.plot_title = self.knot_description +\
                            ", theta = " + str(theta)

            return pattern_part, satellite_part, sf
        get_summand_signture_function.__doc__ = \
            get_summand_signture_function_docsting

        return get_summand_signture_function

    def get_file_name_for_summand_plot(self, theta=0):
        if self.knot_as_k_values[0] < 0:
            name = "inv_T_"
        else:
            name = "T_"
        for k in self.knot_as_k_values:
            name += str(abs(k)) + "_"
        name += "_theta_" + str(theta)
        return name

    def plot_summand_for_theta(self, theta, save_path=None):
        pp, sp, sf = self.signature_as_function_of_theta(theta)
        title = self.knot_description + ", theta = " + str(theta)
        if save_path is not None:
            file_name = self.get_file_name_for_summand_plot(theta)
            save_path = os.path.join(save_path, file_name)
        sig.SignaturePloter.plot_sum_of_two(pp, sp, title=title,
                                            save_path=save_path)

    def plot_summand_sigma(self):
        sigma = self.sigma_as_function_of_theta
        # pattern part
        th_values = list(range(abs(self.knot_as_k_values[-1]) + 1))
        y = [sigma(th)[0] for th in th_values]
        print("plot_summand_sigma")
        print(th_values)
        print(y)

        # satellite_part
        patt_k = self.knot_as_k_values[-1]
        patt_q = 2 * abs(patt_k) + 1
        ksi = 1/patt_q
        x = []
        s = self.get_untwisted_signature_function
        list_of_signatue_functions = [s(k) for k in self.knot_as_k_values[:-1]]
        for i, k in enumerate(self.knot_as_k_values[:-1][::-1]):
            layer_num = i + 1
            x.append(ksi * layer_num)
        print("\nx")
        print(x)
        print(th_values)
        print("\nx product")
        x = list(set(it.product(x, th_values)))
        x = [(a * b) for (a, b) in x]
        print(x)


    def print_sigma_as_function_of_theta(self, theta):
        if not theta:
            return

        # theta should not be larger than q for the pattern.
        patt_k = self.knot_as_k_values[-1]

        patt_q = 2 * abs(patt_k) + 1
        theta %= patt_q

        ksi = 1/patt_q

        # satellite part (Levine-Tristram signatures)
        print(3 * "\n" + 10 * "#" + " " + self.knot_description +
              " " + 10 * "#" + "\n")

        satellite_part = 0
        for layer_num, k in enumerate(self.knot_as_k_values[::-1]):

            sigma_q = self.get_untwisted_signature_function(k)
            arg = ksi * theta * layer_num
            sp = sigma_q(arg)
            satellite_part += 2 * sp

            if details and arg:
                label = "ksi * theta * layer_num = " + str(arg)
                title = self.knot_description + ", layer " + str(layer_num)
                title += ", theta = " + str(theta)
                sigma_q.plot(special_point=(mod_one(arg), sp),
                             special_label=label,
                             title=title,)

        pp = (-patt_q + 2 * theta - 2 * (theta^2/patt_q)) * sign(patt_k)
        sigma = pp + satellite_part
        print(self.knot_description + ", theta = " + str(theta))
        print("pp = " + str(pp), end=', ')
        print("satellite_part = " + str(satellite_part) + "\n")

    def get_sigma_as_function_of_theta(self):

        patt_k = self.knot_as_k_values[-1]
        patt_q = 2 * abs(patt_k) + 1
        ksi = 1/patt_q

        def sigma_as_function_of_theta(theta):
            if theta == 0:
                return 0, 0, 0

            # theta should not be larger than q for the pattern.
            patt_k = self.knot_as_k_values[-1]
            theta %= (2 * abs(patt_k) + 1)

            satellite_part = 0
            for i, k in enumerate(self.knot_as_k_values[:-1][::-1]):
                layer_num = i + 1
                sigma_q = self.get_untwisted_signature_function(k)
                sp = 2 * sigma_q(ksi * theta * layer_num)
                satellite_part += sp
            if theta:
                pp = (-patt_q + 2 * theta - 2 * (theta^2/patt_q)) * sign(patt_k)
            else:
                pp = 0
            return pp, satellite_part, pp + satellite_part
        return sigma_as_function_of_theta


class CableSum:

    def __init__(self, knot_sum, verbose=False):

        self.verbose = verbose
        self.knot_sum_as_k_valus = knot_sum
        self.knot_description = self.get_knot_descrption(knot_sum)
        self.patt_k_list = [abs(i[-1]) for i in knot_sum]
        self.patt_q_list = [2 * i + 1 for i in self.patt_k_list]

        if any(n not in Primes() for n in self.patt_q_list):
            msg = "Incorrect k- or q-vector. This implementation assumes that"\
                  + " all last q values are prime numbers.\n" + \
                  str(self.patt_q_list)
            raise ValueError(msg)
        self.q_order = LCM_list(self.patt_q_list)

        self.knot_summands = [CableSummand(k, verbose) for k in knot_sum]
        self.signature_as_function_of_theta = \
                        self.get_signature_as_function_of_theta()
        self.sigma_as_function_of_theta = \
                        self.get_sigma_as_function_of_theta()


    def __call__(self, *thetas):
        return self.signature_as_function_of_theta(*thetas)

    def get_dir_name_for_plots(self, dir=None):
        dir_name = ''
        for knot in self.knot_summands:
            if knot.knot_as_k_values[0] < 0:
                dir_name += "inv_"
            dir_name += "T_"
            for k in knot.knot_as_k_values:
                k = 2 * abs (k) + 1
                dir_name += str(k) + "_"
        dir_name = dir_name[:-1]
        print(dir_name)
        dir_path = os.getcwd()
        if dir is not None:
            dir_path = os.path.join(dir_path, dir)
        dir_path = os.path.join(dir_path, dir_name)

        if not os.path.isdir(dir_path):
            os.mkdir(dir_path)
        return dir_name

    def plot_sum_for_theta_vector(self, thetas, save_to_dir=False):
        if save_to_dir:
            if not os.path.isdir(PLOTS_DIR):
                os.mkdir(PLOTS_DIR)
            dir_name = self.get_dir_name_for_plots(dir=PLOTS_DIR)
            save_path = os.path.join(os.getcwd(), PLOTS_DIR)
            save_path = os.path.join(save_path, dir_name)
        else:
            save_path = None

        for theta, knot in zip(thetas, self.knot_summands):
            knot.plot_summand_for_theta(thetas, save_path=save_path)

        # pp, sp, sf = self.signature_as_function_of_theta(*thetas)
        # title = self.knot_description + ", thetas = " + str(thetas)
        # if save_path is not None:
        #     file_name = re.sub(r', ', '_', str(thetas))
        #     file_name = re.sub(r'[\[\]]', '', str(file_name))
        #     file_path = os.path.join(save_path, file_name)
        # sig.SignaturePloter.plot_sum_of_two(pp, sp, title=title,
        #                                      save_path=file_path)
        #
        # if save_path is not None:
        #     file_path = os.path.join(save_path, "all_" + file_name)
        # sf_list = [knot.signature_as_function_of_theta(thetas[i])[2]
        #             for i, knot in enumerate(self.knot_summands)]
        # sig.SignaturePloter.plot_many(*sf_list, cols=2)
            #     pp, sp, sf = knot.signature_as_function_of_theta(thetas[i])
            #     (pp + sp) = sp.plot
            #
            # sig.SignatureFunction.plot_sum_of_two(pp, sp, title=title,
            #                                      save_path=file_path)


        return dir_name

    def plot_sigma_for_summands(self):
        for knot in self.knot_summands:
            knot.plot_summand_sigma()

    def parse_thetas(self, *thetas):
        summands_num = len(self.knot_sum_as_k_valus)
        if not thetas:
            thetas = summands_num * (0,)
        elif len(thetas) == 1 and summands_num > 1:
            if isinstance(thetas[0], Iterable):
                if len(thetas[0]) >= summands_num:
                    thetas = thetas[0]
                elif not thetas[0]:
                    thetas = summands_num * (0,)
            elif thetas[0] == 0:
                thetas = summands_num * (0,)
            else:
                msg = "This function takes at least " + str(summands_num) + \
                      " arguments or no argument at all (" + str(len(thetas)) \
                      + " given)."
                raise TypeError(msg)
        return tuple(thetas)

    @staticmethod
    def get_knot_descrption(knot_sum):

        """
        Arguments:
            arbitrary number of lists of numbers,
            each list encodes a single cable.
        Examples:
            sage: get_knot_descrption([1, 3], [2], [-1, -2], [-3])
            'T(2, 3; 2, 7) # T(2, 5) # -T(2, 3; 2, 5) # -T(2, 7)'
        """

        description = ""
        for knot in knot_sum:
            if knot[0] < 0:
                description += "-"
            description += "T("
            for k in knot:
                description += "2, " + str(2 * abs(k) + 1) + "; "
            description = description[:-2] + ") # "
        return description[:-3]

    def get_sigma_as_function_of_theta(self, verbose=None):
        default_verbose = verbose or self.verbose

        def sigma_as_function_of_theta(*thetas, verbose=None, **kwargs):

            verbose = verbose or default_verbose
            thetas = self.parse_thetas(*thetas)
            sigma = 0
            for th, knot in zip(thetas, self.knot_summands):
                _, _, s = knot.sigma_as_function_of_theta(th)
                sigma += s
            return sigma

        return sigma_as_function_of_theta

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

            satellite_part = sig.SignatureFunction()
            pattern_part = sig.SignatureFunction()

            # for each cable knot (summand) in cable sum apply theta
            for th, knot in zip(thetas, self.knot_summands):
                pp, sp, _ = knot.signature_as_function_of_theta(th)
                pattern_part += pp
                satellite_part += sp


            sf = pattern_part + satellite_part

            if verbose:
                print()
                print(str(thetas))
                print(sf)
            assert sf.total_sign_jump() == 0
            return pattern_part, satellite_part, sf

        signature_as_function_of_theta.__doc__ =\
                            signature_as_function_of_theta_docstring
        return signature_as_function_of_theta

    def get_sign_ext_for_theta(self, thetas, limit):
        _, _, sf = self.signature_as_function_of_theta(*thetas)
        return sf.extremum(limit=limit)[1]

    def is_metabolizer(self, theta):
        # Check if square alternating difference
        # divided by last q value is integer.
        result = sum(el^2 / self.patt_q_list[idx] * (-1)^idx
                      for idx, el in enumerate(theta))
        return result.is_integer()

    def is_function_big_in_ranges(self, ranges_list, invariant=SIGMA,
                                  verbose=None):
        verbose = verbose or self.verbose
        if invariant == SIGNATURE:
            get_invariant = self.get_sign_ext_for_theta
            name = "signature (extremum)"
        else:
            get_invariant = self.sigma_as_function_of_theta
            name = "sigma value"

        for thetas in it.product(*ranges_list):

            # Check only non-zero metabolizers.
            if not self.is_metabolizer(thetas) or not any(thetas):
                continue
            #
            # cond1 = thetas[0] and thetas[3] and not thetas[1] and not thetas[2]
            # cond = thetas[0] and thetas[3] and not thetas[1] and not thetas[2]


            function_is_small = True
            # Check if any element generated by thetas vector
            # has a large signature or sigma.
            for shift in range(1, self.q_order):
                shifted_thetas = [shift * th for th in thetas]
                limit = 5 + np.count_nonzero(shifted_thetas)
                inv_value = get_invariant(shifted_thetas, limit=limit)
                abs_value = abs(inv_value)

                if verbose:
                    if shift == 1:
                        print("\n" + "*" * 10)
                        print("Knot sum:\n" + self.knot_description)
                        print("[ characters ] " + name)
                    print(shifted_thetas, end=" ")
                    print(inv_value)

                if abs_value > limit:
                    function_is_small = False
                    if invariant == SIGMA and verbose:
                        self.print_calculations_for_sigma(*shifted_thetas)
                    break
            if function_is_small:
                return False
        return True

    def print_calculations_for_sigma(self, *thetas):

        print("Calculation details for a cable sum:\n" +
              self.knot_description + "\nand theta vector: " +
              str(thetas) + "\n")

        for i, (th, knot) in enumerate(zip(thetas, self.knot_summands)):
            print("{}. {}, theta = {}".format(i + 1, knot.knot_description, th))
            if not th:
                continue
            patt_k = knot.knot_as_k_values[-1]
            q = 2 * abs(patt_k) + 1
            th %= q
            if patt_k > 0:
                print("Pattern part = pp")
            else:
                print("Pattern part = -pp")

            print("pp = -q  + 2 * theta * (q - theta)/q =")
            print("   = -{} + 2 * {}     *  ({}  - {} )/{} =".format(
                         q,      th, q,  th, q))
            print("   = -{} +      {}     *  ({}       )/{} =".format(
                         q,   2 * th, q - th, q))
            print("   = -{} + {} * {} = ".format(
                         q,   2 * th, (q - th)/ q))
            print("   = -{} +  {}  = ".format(
                         q,   2 * th * (q - th)/ q))
            print("   = {}  ".format(
                        -q + (2 * th * (q - th)/ q)))

            pp = (-q + 2 * th - 2 * (th^2/q)) * sign(patt_k)
            sigma = knot.sigma_as_function_of_theta(th)
            print("Pattern part = {} ~ {}".format(sigma[0],int(sigma[0])))
            print("Satellite part = {}".format(sigma[1]))
            print("Sigma = {} ~ {}\n".format(sigma[2], int(sigma[2])))

    def is_function_big_for_all_metabolizers(self, invariant=SIGMA):
        num_of_summands = len(self.knot_sum_as_k_valus)
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
            if not self.is_function_big_in_ranges(ranges_list, invariant):
                return False
        return True


class CableTemplate:

    def __init__(self, knot_formula, q_vector=None, k_vector=None,
                 generate_q_vector=True, slice=True, verbose=False):
        self.verbose = verbose
        self._knot_formula = knot_formula
        # q_i = 2 * k_i + 1
        if k_vector is not None:
            self.k_vector = k_vector
        elif q_vector is not None:
            self.q_vector = q_vector
        elif generate_q_vector:
            self.q_vector = self.get_q_vector(slice=slice)

    @property
    def cable(self):
        if self._cable is None:
            msg = "q_vector for cable instance has not been set explicit. " + \
                  "The variable is assigned a default value."
            warnings.warn(msg)
            self.fill_q_vector()
        return self._cable

    def fill_q_vector(self, q_vector=None, slice=True, lowest_number=2):
        self.q_vector = q_vector or self.get_q_vector(slice, lowest_number)

    @property
    def knot_formula(self):
        return self._knot_formula

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

        self.knot_sum_as_k_valus = eval(self.knot_formula)
        self._cable = CableSum(self.knot_sum_as_k_valus, verbose=self.verbose)
        self._q_vector = [2 * k_val + 1 for k_val in k]

    @property
    def q_vector(self):
        return self._q_vector
    @q_vector.setter
    def q_vector(self, new_q_vector):
        self.k_vector = [(q - 1)/2 for q in new_q_vector]

    @staticmethod
    def extract_max(string):
        numbers = re.findall(r'\d+', string)
        numbers = map(int, numbers)
        return max(numbers)

    def get_q_vector(self, slice=True, lowest_number=2):
        knot_formula = self.knot_formula
        q_vector = [0] * (self.extract_max(knot_formula) + 1)
        P = Primes()
        for layer in self.get_layers_from_formula(knot_formula)[::-1]:
            for el in layer:
                q_vector[el] = P.next(lowest_number)
                lowest_number = q_vector[el]
            lowest_number *= 4
        return q_vector

    @staticmethod
    def get_layers_from_formula(knot_formula):
        k_indices = re.sub(r'[k-]', '', knot_formula)
        k_indices = re.sub(r'\[\d+\]', lambda x: x.group()[1:-1], k_indices)
        k_indices = eval(k_indices)
        number_of_layers = max(len(lst) for lst in k_indices)
        layers = []
        for i in range(1, number_of_layers + 1):
            layer = [lst[-i] for lst in k_indices if len(lst)>= i]
            layers.append(layer)
        return layers

    def add_with_shift(self, other):
        shift = self.extract_max(self.knot_formula) + 1
        o_formula = re.sub(r'\d+', lambda x: str(int(x.group()) + shift),
                           other.knot_formula)
        return self + CableTemplate(o_formula)

    def __add__(self, other):
        knot_formula = self.knot_formula[:-1] + ",\n" + other.knot_formula[1:]
        return CableTemplate(knot_formula)


def mod_one(n):
    return n - floor(n)


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
    an object sig.SignatureFunction.

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

get_summand_signture_function_docsting = \
    """
    This function returns sig.SignatureFunction for previously defined single
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

    Returns object of sig.SignatureFunction class for a previously defined
    connected sum of len(arg) cables.
    Accept len(arg) arguments: for each cable one theta parameter.
    If call with no arguments, all theta parameters are set to be 0.
    """
#
# CableSummand.get_blanchfield_for_pattern.__doc__ = \
#     """
#     Arguments:
#         k_n:    a number s.t. q_n = 2 * k_n + 1, where
#                 T(2, q_n) is a pattern knot for a single cable from a cable sum
#         theta:  twist/character for the cable (value form v vector)
#     Return:
#         sig.SignatureFunction created for pattern signature function
#         for a given cable and theta/character
#     Based on:
#         Proposition 9.8. in Twisted Blanchfield Pairing
#         (https://arxiv.org/pdf/1809.08791.pdf)
#     """

# CableSummand.get_summand_signature_as_theta_function.__doc__ = \
#     """
#     Argument:
#         n integers that encode a single cable, i.e.
#         values of q_i for T(2,q_0; 2,q_1; ... 2, q_n)
#     Return:
#         a function that returns sig.SignatureFunction for this single cable
#         and a theta given as an argument
#     """
