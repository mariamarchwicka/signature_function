#!/usr/bin/python
import os
import sys

import itertools as it
import re
import numpy as np


attach("cable_signature.sage")
attach("my_signature.sage")



# TBD: read about Factory Method, variable in docstring, sage documentation


class Config(object):
    def __init__(self):
        self.f_results = os.path.join(os.getcwd(), "results.out")

        # knot_formula is a schema for knots which signature function
        # will be calculated
        self.knot_formula = "[[k[0], k[1], k[3]], " + \
                             "[-k[1], -k[3]], " + \
                             "[k[2], k[3]], " + \
                             "[-k[0], -k[2], -k[3]]]"

        # self.knot_formula = "[[k[0], k[1], k[4]], [-k[1], -k[3]], \
        #                      [k[2], k[3]], [-k[0], -k[2], -k[4]]]"
        #
        #
        #
        # self.knot_formula = "[[k[3]], [-k[3]], \
        #                      [k[3]], [-k[3]] ]"
        #
        # self.knot_formula = "[[k[3], k[2], k[0]], [-k[2], -k[0]], \
        #                      [k[1], k[0]], [-k[3], -k[1], -k[0]]]"
        #
        # self.knot_formula = "[[k[0], k[1], k[2]], [k[3], k[4]], \
        #                      [-k[0], -k[3], -k[4]], [-k[1], -k[2]]]"
        # self.knot_formula = "[[k[0], k[1], k[2]], [k[3]],\
        #                          [-k[0], -k[1], -k[3]], [-k[2]]]"
        self.limit = 3

        # in rch for large sigma, for 1. checked knot q_1 = 3 + start_shift
        self.start_shift =  0

        self.verbose = True
        # self.verbose = False

        self.print_results = True
        # self.print_results = False

        # is the ratio restriction for values in q_vector taken into account
        self.only_slice_candidates = True
        self.only_slice_candidates = False




def main(arg=None):
    try:
        limit = int(arg[1])
    except (IndexError, TypeError):
        limit = None

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
    print(cable.is_signature_big_for_all_metabolizers())


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


if __name__ == '__main__':
    global config
    config = Config()
    if '__file__' in globals():
        # skiped in interactive mode as __file__ is not defined
        main(sys.argv)
    else:
        pass
        # main()
