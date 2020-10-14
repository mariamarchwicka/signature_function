#!/usr/bin/python
attach("cable_signature.sage")
# attach("my_signature.sage")

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

def is_big_in_ranges(cable, ranges_list):
    we_have_no_problem = True
    for theta in it.product(*ranges_list):
        if all(i == 0 for i in theta):
            continue
        we_have_a_problem = True
        if cable.is_metaboliser(theta):
            # print("\n" * 10)
            for shift in range(1, cable.q_order):
                shifted_theta = [(shift * th) % cable.last_q_list[i]
                                 for i, th in enumerate(theta)]
                shifted_theta = [min(th, cable.last_q_list[i] - th)
                                 for i, th in enumerate(shifted_theta)]
                sf = cable.signature_as_function_of_theta(*shifted_theta)
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
                we_have_a_big_problem = True
                break
    if not we_have_no_problem:
        print("we have a big problem")
    return we_have_no_problem

def check_all_thetas(cable):
    upper_bounds = cable.last_k_list[:3]
    ranges_list = [range(0, i + 1) for i in upper_bounds]
    ranges_list.append(range(0, 2))
    ranges_list += [range(0, 1) for _ in range(4)]
    if not is_big_in_ranges(cable, ranges_list):
        return False
    upper_bounds = cable.last_k_list[5:8]
    ranges_list = [range(0, 1) for _ in range(4)]
    ranges_list += [range(0, i + 1) for i in upper_bounds]
    ranges_list.append(range(0, 2))
    if not is_big_in_ranges(cable, ranges_list):
        return False
    return True





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

        q = [P.unrank(i) for i in c]
        ratio = q[3] > 4 * q[2] and q[2] > 4 * q[1] and q[1] > 4 * q[0]
        if not ratio:
                # print("Ratio-condition does not hold")
            continue
        print("q = ", q)
