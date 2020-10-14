#!/usr/bin/python
attach("cable_signature.sage")
attach("my_signature.sage")

import numpy as np


def main():
    global cable, cab_2, cab_1, joined_formula
    # self.knot_formula = "[[k[0], k[1], k[3]], " + \
    #                      "[-k[1], -k[3]], " + \
    #                      "[k[2], k[3]], " + \
    #                      "[-k[0], -k[2], -k[3]]]"

    # knot_formula = config.knot_formula
    # q_vector = (3, 5, 7, 13)
    # cab_to_update = TorusCable(knot_formula=knot_formula, q_vector=q_vector)
    # q_vector = (3, 5, 7, 11)
    # cab_to_add = TorusCable(knot_formula=knot_formula, q_vector=q_vector)
    # cab_shifted = cab_to_update.add_with_shift(cab_to_add)

    q_vector = (5, 13, 19, 41,\
                5, 17, 23, 43)

    knot_formula = "[[k[0], k[5], k[3]], " + \
                         "[-k[1], -k[3]], " + \
                         "[k[2], k[3]], " + \
                         "[-k[0], -k[2], -k[3]]]"
    cab_1 = TorusCable(knot_formula=knot_formula, q_vector=q_vector)
    knot_formula = "[[k[4], k[1], k[7]], " + \
                         "[-k[5], -k[7]], " + \
                         "[k[6], k[7]], " + \
                         "[-k[4], -k[6], -k[7]]]"
    cab_2 = TorusCable(knot_formula=knot_formula, q_vector=q_vector)
    cable = cab_1 + cab_2
    joined_formula = cable.knot_formula

def check_all_thetas(cable):
    upper_bounds = cable.last_k_list[:4]
    # upper_bounds += [0, 0, 0, 0]
    ranges_list = [range(1, i + 1) for i in upper_bounds]
    ranges_list += [range(0, 1) for _ in range(4)]


    print(ranges_list)
    for theta in it.product(*ranges_list):
        # pass
        if cable.is_metaboliser(theta):
            print("\n" * 10)
            print("!" * 100)
            for shift in range(1, cable.q_order):
                shifted_theta = [(shift * th) % cable.last_q_list[i]
                                 for i, th in enumerate(theta)]
                shifted_theta = [min(th, cable.last_q_list[i] - th)
                                 for i, th in enumerate(shifted_theta)]
                print(shifted_theta)
                sf = cable.signature_as_function_of_theta(*shifted_theta)
                extremum = abs(sf.extremum())
                print(extremum)
                if extremum > 5 + np.count_nonzero(shifted_theta):
                    print("ok")
                    break
                else:
                    print("hliphlip")
            print("!" * 100)
            print("\n" * 10)
    return


def get_q_vector(q_vector_size, lowest_number=1):
    q = [lowest_number] * q_vector_size
    P = Primes()
    next_number = P.next(lowest_number)
    for i in range(q_vector_size):
        q[i] = next_number
        next_number = P.next(4 * next_number)

    return q


    next_number = P.next(lowest_number)
    for i, q in enumerate(q_vector):
        q[i] = next_number
        next_number = P.next(lowest_number)

        q = [P.unrank(i + config.start_shift) for i in c]
        ratio = q[3] > 4 * q[2] and q[2] > 4 * q[1] and q[1] > 4 * q[0]
        if not ratio:
                # print("Ratio-condition does not hold")
            continue
        print("q = ", q)
    #     cable = TorusCable(knot_formula=knot_formula, q_vector=q)
    #     list_of_ranges = config.get_list_of_ranges(cable.q_order)
    #     if cable.eval_cable_for_large_values(list_of_ranges, SIGMA,
    #                                         verbose=verbose,
    #                                         print_results=print_results):
    #         good_knots.append(cable.knot_description)
    # return good_knots







# cab_to_update.update(cab_to_add)
# cab_to_update.update(cab_to_add)
# cab_to_update.update(cab_to_add)
# cab_to_update.update(cab_to_add)

    pass


#
# def get_blanchfield_for_pattern(k_n, theta):
#     if theta == 0:
#         sf = TorusCable.get_untwisted_signature_function(k_n)
#         return sf.square_root() + sf.minus_square_root()
#
#     results = []
#     k = abs(k_n)
#     ksi = 1/(2 * k + 1)
#
#     # print("lambda_odd, i.e. (theta + e) % 2 != 0")
#     for e in range(1, k + 1):
#         if (theta + e) % 2 != 0:
#             results.append((e * ksi, 1 * sgn(k_n)))
#             results.append((1 - e * ksi, -1 * sgn(k_n)))
#     # print("\nlambda_even")
#     # print("\nnormal")
#     results_odd = results
#     results = []
#     for e in range(1, theta):
#         if (theta + e) % 2 == 0:
#             results.append((e * ksi, 1 * sgn(k_n)))
#             results.append((1 - e * ksi, -1 * sgn(k_n)))
#     # print(results)
#     results_even_small = results
#     results = []
#     # print("reversed")
#     for e in range(theta + 1, k + 1):
#         if (theta + e) % 2 == 0:
#             results.append((e * ksi, -1 * sgn(k_n)))
#             results.append((1 - e * ksi, 1 * sgn(k_n)))
#     # print(results)
#     results_even_big = results
#
#     return results_odd, results_even_small, results_even_big
#     # return SignatureFunction(values=results)
#
# def main():
#     prim = 3
#     P = Primes()
#     for it in range(20):
#         prim = P.next(prim)
#         k_j = (prim - 1)/2
#         print(60 * "*")
#         print("k is " + str(k_j))
#         print(60 * "*")
#
#         for i in range(1, k_j + 1):
#
#             a, j, m =  get_blanchfield_for_pattern(k_j, i)
#             sf_j = SignatureFunction(j)
#             sf_a = SignatureFunction(a)
#             sf_m = SignatureFunction(m)
#             sf_jam = sf_j + sf_a + sf_m
#             assert TorusCable.get_blanchfield_for_pattern(k_j, i) == sf_jam
#             af, jf, mf =  get_blanchfield_for_pattern(-k_j, prim - i)
#             print(60 * "*")
#             print("lists")
#             print(af)
#             print(jf)
#             print(mf)
#             j = SignatureFunction(jf)
#             a = SignatureFunction(af)
#             m =  SignatureFunction(mf)
#             minus_jam = j + a + m
#             values = cmp_blanchfield_for_pattern(-k_j, prim - i)
#             print("sum of lists - 3 lists added")
#             print(sorted(jf + af + mf))
#
#             print("sum of lists - all values from Blanchfield")
#             print(sorted(values))
#
#             assert values == af + jf + mf
#             print("not equeles - sf from all values")
#
#             print(SignatureFunction(values))
#             print("not equeles - sum of sf")
#             print("sf for each list sep")
#             print("jf")
#             print(jf)
#             print("af")
#             print(af)
#             print("mf")
#             print(mf)
#             print("sum of abouve sfs")
#             print(minus_jam)
#             assert TorusCable.get_blanchfield_for_pattern(-k_j, prim - i) == \
#                 SignatureFunction(values)
#             assert TorusCable.get_blanchfield_for_pattern(-k_j, prim - i) == \
#                 minus_jam
#
#
#             # a, j, m =  get_blanchfield_for_pattern(k_j, 2 * k_j + 1 - i)
#             # j = SignatureFunction(j)
#             # a = SignatureFunction(a)
#             # m =  SignatureFunction(m)
#
#             print("*" * 100)
#             print("i is " + str(i) + ", q is " + str(2 * k_j + 1), " q - i is " + str(2 * k_j + 1 - i))
#             print("*" * 100)
#
#             ajm = sf_j + sf_a + sf_m
#             ajm_minus = a + j + m
#             print("4 times")
#             print(ajm + ajm + ajm_minus + ajm_minus)
#             print("is big")
#             print((ajm + ajm + ajm_minus + ajm_minus).is_big())
#             print()
