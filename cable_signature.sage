#!/usr/bin/python
import collections


class TorusCable(object):
    def __init__(self, knot_formula, k_vector=None, q_vector=None):
        # q_i = 2 * k_i + 1

        if k_vector is None:
            if q_vector is None:
                # TBD docstring
                print("Please give a list of k (k_vector) \
                        or q values (q_vector).")
                return None
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
        self.__signature_as_function_of_theta = None

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

    # searching for signature == 0
    def get_signature_as_function_of_theta(self, verbose=False):
        if self.__signature_as_function_of_theta is None:
            self.__signature_as_function_of_theta = \
                self.__get_signature_as_function_of_theta(verbose=verbose)
        return self.__signature_as_function_of_theta

    # searching for signature == 0
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
                msg = "This function takes exactly " + str(len_a) + \
                      " arguments or no argument at all (" + str(len_t) + \
                      " given)."
                raise TypeError(msg)

            sf = SignatureFunction()

            # for each cable knot in cable sum apply theta
            for i, knot in enumerate(self.knot_sum):
                try:
                    ssf = get_signature_summand_as_theta_function(*knot)
                    sf += ssf(thetas[i])
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
        signature_as_function_of_theta.__doc__ =\
                                        signature_as_function_of_theta_docstring
        return signature_as_function_of_theta

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
            print("Zero theta combinations: ")
            for el in list_of_good_vectors:
                print(el)
        if number_of_null_comb^2 >= number_of_all_comb:
            return number_of_null_comb, number_of_all_comb
        return None

    # check sigma for all v = s * [a_1, a_2, a_3, a_4] for s in [1, q_4 - 1]
    def __is_sigma_for_vector_class_big(self, theta_vector):
        [a_1, a_2, a_3, a_4] = theta_vector
        q_4 = self.q_vector[3]
        for shift in range(1, q_4):
            shifted_theta = [(shift * a) % q_4 for a in
                             [a_1, a_2, a_3, a_4]]
            sigma_v = self.__sigma_function(shifted_theta)
            if abs(sigma_v) > 5 + np.count_nonzero(shifted_theta):
                return True
        return False

    def __tmp_print_all_sigma_for_vector_class(self, theta_vector):
        print("\n")
        print(self.knot_description)
        print("vector = " + str(theta_vector))
        [a_1, a_2, a_3, a_4] = theta_vector
        q_4 = self.q_vector[3]
        for shift in range(1, q_4):
            shifted_theta = [(shift * a) % q_4 for a in
                             [a_1, a_2, a_3, a_4]]
            print(str(shifted_theta) + "\t\t" + \
                    str(self.__sigma_function(shifted_theta)))
        print("\n")

    def __tmp_get_max_sigma_for_vector_class(self, theta_vector):
        # print("\n")
        # print(self.knot_description)
        # print("vector = " + str(theta_vector))
        max_sigma = (theta_vector, 0)
        [a_1, a_2, a_3, a_4] = theta_vector
        q_4 = self.q_vector[3]
        for shift in range(1, q_4):
            shifted_theta = [(shift * a) % q_4 for a in
                             [a_1, a_2, a_3, a_4]]
            sigma = self.__sigma_function(shifted_theta)
            if abs(sigma) > abs(max_sigma[1]):
                max_sigma = (shifted_theta, sigma)
        assert max_sigma[1] == 0, knot_description
            # print("\n" + self.knot_description + "\t" + str(max_sigma[0]) +\
            #         "\t" + str(max_sigma[1]))
        return max_sigma[1]



    def is_sigma_for_vector_class_big(self, theta_vector):
        if self.__sigma_function is None:
            self.__sigma_function = self.__get_sigma_function()
        return self.__is_sigma_for_vector_class_big(theta_vector)

    def __get_sigma_function(self):
        k_1, k_2, k_3, k_4 = [abs(k) for k in self.k_vector]
        q_4 = 2 * k_4 + 1
        ksi = 1/q_4
        sigma_q_1 = get_untwisted_signature_function(k_1)
        sigma_q_2 = get_untwisted_signature_function(k_2)
        sigma_q_3 = get_untwisted_signature_function(k_3)

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
                    tp[i] = -q_4 + 2 * a - 2 * (a^2/q_4)
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
        q_4 = 2 * k_4 + 1
        ksi = 1/q_4
        sigma_q_1 = get_untwisted_signature_function(k_1)
        sigma_q_2 = get_untwisted_signature_function(k_2)
        sigma_q_3 = get_untwisted_signature_function(k_3)
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

    def get_number_of_combinations_of_theta(self):
        number_of_combinations = 1
        for knot in self.knot_sum:
            number_of_combinations *= (2 * abs(knot[-1]) + 1)
        return number_of_combinations

    def print_results_sigma(self, theta_vector, twisted_part):
        a_1, a_2, a_3, a_4 = theta_vector
        knot_description = self.knot_description
        q_4 = self.q_vector[-1]
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

    # searching for sigma > 5 + #(v_i != 0)
    def calculate_sigma(self, theta_vector):
        if self.__sigma_function is None:
            self.__sigma_function = self.__get_sigma_function()
        return self.__sigma_function(theta_vector)


    # searching for sigma > 5 + #(v_i != 0)
    def __check_combinations_in_range(self, range_product):
        large_sigma_for_all_combinations = True
        bad_vectors = []
        good_vectors = []
        q_4 = self.q_vector[-1]
        for vector in range_product:
            a_1, a_2, a_3, a_4 = vector
            if (a_1^2 - a_2^2 + a_3^2 - a_4^2) % q_4:
                continue
            if all(a in [1, q_4 - 1] for a in vector):
                is_all_one = True
            else:
                is_all_one = False
            if self.__is_sigma_for_vector_class_big(vector):
                good_vectors.append(vector)
                # if is_all_one:
                #     print("\nHURA" * 100)
                #     print(self.knot_description)
                #     self.__tmp_print_all_sigma_for_vector_class(vector)
                # pass
            else:
                if is_all_one:
                    self.__tmp_get_max_sigma_for_vector_class(vector)
                bad_vectors.append(vector)
                #####################################################
                if len(bad_vectors) > 8:
                    break
                ####################################################
                large_sigma_for_all_combinations = False
        return good_vectors, bad_vectors

    # searching for sigma > 5 + #(v_i != 0)
    def check_combinations_in_range(self, range_product):
        if self.__sigma_function is None:
            self.__sigma_function = self.__get_sigma_function()
        return self.__check_combinations_in_range(range_product)

    # searching for sigma > 5 + #(v_i != 0)
    def __check_all_combinations_in_ranges(self, list_of_ranges,
                                           print_results=True):
        all_combinations_pass = True
        all_bad_vectors = []
        number_of_all_good_v = 0
        for i, range_product in enumerate(list_of_ranges):
            good_v, bad_v  = self.__check_combinations_in_range(range_product)
            number_of_all_good_v += len(good_v)
            all_bad_vectors = list(it.chain(all_bad_vectors, bad_v))
            if bad_v:
                all_combinations_pass = False
            if len(all_bad_vectors) > 8:
                break
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
            if len(all_bad_vectors) < 8:
                print()
                print(all_bad_vectors)


        return all_combinations_pass

    # searching for sigma > 5 + #(v_i != 0)
    def eval_cable_for_large_sigma(self, list_of_ranges,
                                    print_results=False, verbose=False):
        if self.__sigma_function is None:
            self.__sigma_function = self.__get_sigma_function()
        if print_results:
            # print("\n\n")
            # print(100 * "*")
            # print("Searching for a large signature values for the cable sum: ")
            print(self.knot_description, end="\t\t\t")
            # print()
        if self.__check_all_combinations_in_ranges(list_of_ranges,
                                                print_results=print_results):
            return True
        return False



class SignatureFunction(object):

    def __init__(self, values=None, counter=None):
        # set values of signature jumps
        if counter is None:
            counter = collections.Counter()
            if values is None:
                values = []
            assert all(x < 1 for x, y in values),\
                "Signature function is defined on the interval [0, 1)."
            counter = collections.Counter(dict(values))
        self.cnt_signature_jumps = counter

    def sum_of_absolute_values(self):
        return sum([abs(i) for i in self.cnt_signature_jumps.values()])

    def is_zero_everywhere(self):
        return not any(self.cnt_signature_jumps.values())

    def double_cover(self):
        # to read values for t^2
        new_data = []
        for jump_arg, jump in self.cnt_signature_jumps.items():
            new_data.append((jump_arg/2, jump))
            new_data.append((1/2 + jump_arg/2, jump))
        return SignatureFunction(values=new_data)

    def square_root(self):
        # to read values for t^(1/2)
        new_data = []
        for jump_arg, jump in self.cnt_signature_jumps.items():
            if jump_arg < 1/2:
                new_data.append((2 * jump_arg, jump))
        return SignatureFunction(values=new_data)

    def minus_square_root(self):
        # to read values for t^(1/2)
        counter = collections.Counter()
        for jump_arg, jump in self.cnt_signature_jumps.items():
            if jump_arg >= 1/2:
                counter[mod_one(2 * jump_arg)] = jump
        return SignatureFunction(counter=counter)

    def __lshift__(self, shift):
        # A shift of the signature functions corresponds to the rotation.
        return self.__rshift__(-shift)

    def __rshift__(self, shift):
        new_data = []
        for jump_arg, jump in self.cnt_signature_jumps.items():
            new_data.append((mod_one(jump_arg + shift), jump))
        return SignatureFunction(values=new_data)

    def __neg__(self):
        counter = collections.Counter()
        counter.subtract(self.cnt_signature_jumps)
        return SignatureFunction(counter=counter)

    # TBD short
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
                for jump_arg, jump in sorted(self.cnt_signature_jumps.items())])
        return result

    def __repr__(self):
        result = ''.join([str(jump_arg) + ": " + str(jump) + ", "
                for jump_arg, jump in sorted(self.cnt_signature_jumps.items())])
        return result[:-2] + "."

    def __call__(self, arg):
        # Compute the value of the signature function at the point arg.
        # This requires summing all signature jumps that occur before arg.
        arg = mod_one(arg)
        cnt = self.cnt_signature_jumps
        before_arg = [jump for jump_arg, jump in cnt.items() if jump_arg < arg]
        return 2 * sum(before_arg) + cnt[arg]

def mod_one(n):
    return n - floor(n)
