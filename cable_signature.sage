#!/usr/bin/python
import collections
import numpy as np
import itertools as it
from typing import Iterable

SIGNATURE = 0
SIGMA = 1



class SignatureFunction(object):

    def __init__(self, values=None, counter=None):
        # set values of signature jumps
        if counter is None:
            counter = collections.Counter()

            if values is None:
                values = []

            msg = "Signature function is defined on the interval [0, 1)."
            assert all(k < 1 for k, v in values), msg

            for k, v in values:
                counter[k] += v
        self.cnt_signature_jumps = counter

    def sum_of_absolute_values(self):
        return sum([abs(i) for i in self.cnt_signature_jumps.values()])

    def is_zero_everywhere(self):
        return not any(self.cnt_signature_jumps.values())

    def double_cover(self):
        # to read values for t^2
        new_data = []
        for jump_arg, jump in self.cnt_signature_jumps.items():
            if jump != 0:
                new_data.append((jump_arg/2, jump))
                new_data.append((1/2 + jump_arg/2, jump))
        return SignatureFunction(values=new_data)

    def square_root(self):
        # to read values for t^(1/2)
        new_data = []
        for jump_arg, jump in self.cnt_signature_jumps.items():
            if jump_arg < 1/2 and jump != 0:
                new_data.append((2 * jump_arg, jump))
        return SignatureFunction(values=new_data)

    def minus_square_root(self):
        # to read values for t^(1/2)
        counter = collections.Counter()
        for jump_arg, jump in self.cnt_signature_jumps.items():
            if jump_arg >= 1/2 and jump != 0:
                counter[mod_one(2 * jump_arg)] = jump
        return SignatureFunction(counter=counter)

    def is_big(self):
        max = 0
        current = 0
        items = sorted(self.cnt_signature_jumps.items())
        for arg, jump in items:
            current += 2 * jump
            msg = "current = " + str(current) + ", jump = " + str(jump)
            msg += "\n" + str(self(arg))
            result = [jump for jump_arg, jump in self.cnt_signature_jumps.items() if jump_arg < mod_one(arg)]
            msg += "\nresult = " + str(sum(result))
            msg += "\narg = " + str(arg)
            assert current == self(arg) + jump, msg
            if abs(current) > abs(max):
                max = current
                # if abs(max) > 9:
                #     return max
        return max

    def __rshift__(self, shift):
        # A shift of the signature functions corresponds to the rotation.
        new_data = []
        for jump_arg, jump in self.cnt_signature_jumps.items():
            new_data.append((mod_one(jump_arg + shift), jump))
        sf = SignatureFunction(values=new_data)
        counter = collections.Counter({mod_one(k + shift) : v \
                                for k,v in self.cnt_signature_jumps.items()})
        assert SignatureFunction(counter=counter) == \
                            SignatureFunction(values=new_data)
        return SignatureFunction(counter=counter)

    def __lshift__(self, shift):
        return self.__rshift__(-shift)

    def __neg__(self):
        counter = collections.Counter()
        counter.subtract(self.cnt_signature_jumps)
        return SignatureFunction(counter=counter)

    def __add__(self, other):
        counter = copy(self.cnt_signature_jumps)
        counter.update(other.cnt_signature_jumps)
        return SignatureFunction(counter=counter)

    def __eq__(self, other):
        return self.cnt_signature_jumps == other.cnt_signature_jumps

    def __sub__(self, other):
        counter = copy(self.cnt_signature_jumps)
        counter.subtract(other.cnt_signature_jumps)
        return SignatureFunction(counter=counter)

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
        return sorted(self.cnt_signature_jumps.items(), key = lambda x: x[0])

    def step_function_data(self):
        # Transform the signature jump data to a format understandable
        # by the plot function.
        l = self.to_list()
        vals = ([(d[0], sum(2 * j[1] for j in l[:l.index(d)+1])) for d in l] +
              [(0,self.cnt_signature_jumps[0]), (1,self.total_sign_jump())])
        return vals

    def plot(self):
        # plot the signture function
        plot_step_function(self.step_function_data())

    def tikz_plot(self, file_name):
        # Draw the graph of the signature and transform it into            TiKz.
        # header of the LaTeX file

        with open(file_name, "w") as output_file:
            output_file.write("\\documentclass[tikz]{standalone}\n")
            output_file.write("\\usetikzlibrary{datavisualization,datavisualization.formats.functions}\n")
            output_file.write("\\begin{document}\n")
            output_file.write("\\begin{tikzpicture}\n")
            data = sorted(self.step_function_data())
            output_file.write("  \\datavisualization[scientific axes,visualize as smooth line,\n")
            output_file.write("  x axis={ticks={none,major={at={")
            output_file.write(", " + str(N(data[0][0],digits=4)) + " as \\(" + str(data[0][0]) + "\\)")
            for jump_arg,jump in data:
                output_file.write(", " + str(N(jump_arg,digits=4)) + " as \\(" + str(jump_arg) + "\\)")
                output_file.write("}}}}\n")
                output_file.write("  ]\n")
                output_file.write("data [format=function]{\n")
                output_file.write("var x : interval [0:1];\n")
                output_file.write("func y = \\value x;\n")
                output_file.write("};\n")
                # close LaTeX enviroments
                output_file.write("\\end{tikzpicture}\n")
                output_file.write("\\end{document}\n")


class TorusCable(object):
    def __init__(self, knot_formula, k_vector=None, q_vector=None):
        # q_i = 2 * k_i + 1
        if k_vector is None:
            if q_vector is None:
                # TBD docstring
                print("Please give a list of k (k_vector) \
                        or q values (q_vector).")
                msg = ""
                raise ValueError(msg)
            else:
                k_vector = [(q - 1)/2 for q in q_vector]
        elif q_vector is None:
                q_vector = [2 * k + 1 for k in k_vector]
        self.knot_formula = knot_formula
        self.k_vector = k_vector
        self.q_vector = q_vector
        k = k_vector
        self.knot_sum = eval(knot_formula)
        self.knot_description = self.get_knot_descrption()
        self.__sigma_function = None
        self.signature_as_function_of_theta = None

    def __get_sigma_function(self):
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

    def __get_signature_as_function_of_theta(self, **key_args):
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
                if isinstance(thetas, Iterable) and len(thetas[0]) == len_a:
                    thetas = thetas[0]
                else:
                    msg = "This function takes exactly " + str(len_a) + \
                          " arguments or no argument at all (" + str(len_t) + \
                          " given)."
                    raise TypeError(msg)

            sf = SignatureFunction()
            untwisted_part = SignatureFunction()
            # for each cable knot in cable sum apply theta
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
            a = thetas[0]
            if all(i == a or i == self.q_vector[-1] - a for i in thetas):
                print()
                print("\n" + "*" * 100)
                print(self.knot_description)
                print("one vector " + str(thetas))
                print("max sf " + str(sf.is_big()))
                print()
                assert untwisted_part.is_zero_everywhere()

            if verbose:
                print()
                print(str(thetas))
                print(sf)
            return sf
        signature_as_function_of_theta.__doc__ =\
                                        signature_as_function_of_theta_docstring
        return signature_as_function_of_theta


    @staticmethod
    def get_blanchfield_for_pattern(k_n, theta):
        if theta == 0:
            sf = TorusCable.get_untwisted_signature_function(k_n)
            return sf.square_root() + sf.minus_square_root()

        results = []
        k = abs(k_n)
        ksi = 1/(2 * k + 1)

        counter = collections.Counter()
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
        w = ([((2 * a + 1)/(2 * q), -1 * sgn(j)) for a in range(k)] +
             [((2 * a + 1)/(2 * q), 1 * sgn(j))
             for a in range(k + 1, 2 * k + 1)])
        return SignatureFunction(values=w)


    def get_knot_descrption(self):
        description = ""
        for knot in self.knot_sum:
            if knot[0] < 0:
                description += "-"
            description += "T("
            for k in knot:
                description += "2, " + str(2 * abs(k) + 1) + "; "
            description = description[:-2] + ") # "
        return description[:-3]

    def get_signature_as_function_of_theta(self, verbose=False):
        if self.signature_as_function_of_theta is None:
            self.signature_as_function_of_theta = \
                self.__get_signature_as_function_of_theta(verbose=verbose)
        return self.signature_as_function_of_theta

    def get_summand_signature_as_theta_function(self, *knot_as_k_values):
        def get_summand_signture_function(theta):
            # TBD: another formula (for t^2) description
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

    # to be removed
    def check_combinations_in_range(self, range_product):
        if self.__sigma_function is None:
            self.__sigma_function = self.__get_sigma_function()
        return self.__check_combinations_in_range(range_product, SIGMA)

    # to be removed
    def calculate_sigma(self, theta_vector):
        if self.__sigma_function is None:
            self.__sigma_function = self.__get_sigma_function()
        return self.__sigma_function(theta_vector)

    def is_sigma_for_vector_class_big(self, theta_vector):
        if self.__sigma_function is None:
            self.__sigma_function = self.__get_sigma_function()
        return self.__is_value_for_vector_class_big(theta_vector, SIGMA)


    # searching for signature == 0
    def check_for_null_theta_combinations(self, verbose=False):
        list_of_good_vectors= []
        number_of_null_comb = 0
        f = self.get_signature_as_function_of_theta(verbose=verbose)
        range_list = [range(abs(knot[-1]) + 1) for knot in self.knot_sum]
        for theta_vector in it.product(*range_list):
            if f(*theta_vector, verbose=False).is_zero_everywhere():
                list_of_good_vectors.append(theta_vector)
                m = len([theta for theta in theta_vector if theta != 0])
                number_of_null_comb += 2^m
        return number_of_null_comb, list_of_good_vectors

    # searching for sigma > 5 + #(v_i != 0)
    def __check_all_combinations_in_ranges(self, list_of_ranges,
                                            sigma_or_sign,
                                            print_results=False):
        all_combinations_pass = True
        all_bad_vectors = []
        number_of_all_good_v = 0
        for i, range_product in enumerate(list_of_ranges):
            good_v, bad_v = self.__check_combinations_in_range(range_product,
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
    def __check_combinations_in_range(self, range_product, sigma_or_sign):
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
            if self.__is_value_for_vector_class_big(vector, sigma_or_sign):
                good_vectors.append(vector)
            else:
                # print(vector)
                bad_vectors.append(vector)
        return good_vectors, bad_vectors

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
    def __is_value_for_vector_class_big(self, theta_vector, sigma_or_sign):
        [a_1, a_2, a_3, a_4] = theta_vector
        q_4 = self.q_vector[-1]
        k_4 = self.k_vector[-1]

        max_sigma = 0

        if sigma_or_sign == SIGNATURE:
            f = self.signature_as_function_of_theta
        else:
            f = self.__sigma_function
        # print(theta_vector, end="\t")

        for shift in range(1, k_4 + 1):
            shifted_theta = [(shift * a) % q_4 for a in
                             [a_1, a_2, a_3, a_4]]
            if sigma_or_sign == SIGNATURE:
                sf = f(shifted_theta)
                sigma_v = sf.is_big()
            else:
                sigma_v = f(shifted_theta)
            print(sigma_v, end=" ")
            if abs(sigma_v) > abs(max_sigma):
                max_sigma = sigma_v
            if abs(sigma_v) > 5 + np.count_nonzero(shifted_theta):
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
        if sigma_or_sign == SIGMA:
            if self.__sigma_function is None:
                self.__sigma_function = self.__get_sigma_function()
        else:
            if self.signature_as_function_of_theta is None:
                self.signature_as_function_of_theta= \
                            self.__get_signature_as_function_of_theta()

        if self.__check_all_combinations_in_ranges(list_of_ranges,
                                                sigma_or_sign,
                                                print_results=print_results):
            return True
        return False

##############################################################################
    # sigma function
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
