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
from PIL import Image
from pathlib import Path
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



    def plot_sum_with_other(self, other,
                            save_path=None, title=''):
        tp = self
        up = other
        sf = tp + up

        fig, axes_matrix = plt.subplots(2, 2, sharey=True,
                                        figsize=(10,5))

        tp.plot(title='subplot_tp',
                subplot=True,
                ax=axes_matrix[0][1])

        up.plot(title='subplot_up',
                subplot=True,
                ax=axes_matrix[1][0],
                color='red',
                linestyle='dotted')

        sf.plot(title='subplot_up',
                subplot=True,
                ax=axes_matrix[0][0],
                color='black')

        tp.plot(title='subplot_tp',
                subplot=True,
                ax=axes_matrix[1][1],
                alpha=0.3)

        up.plot(title='subplot_up',
                subplot=True,
                ax=axes_matrix[1][1],
                color='red', alpha=0.3,
                linestyle='dotted')

        sf.plot(title='subplot_up',
                subplot=True,
                ax=axes_matrix[1][1],
                color='black',
                alpha=0.7,)

        fig.suptitle(title)

        plt.tight_layout()
        if save_path is None:
            save_path = os.path.join(os.getcwd(),"tmp.png")
        save_path = Path(save_path)
        save_path = save_path.with_suffix('.png')


        # print(save_as)

        plt.savefig(save_path)
        plt.close()
        image = Image.open(save_path)
        image.show()


        pass


    def plot(self, subplot=False, ax=None, save_as='sf',
             title="",
             alpha=1,
             color='blue',
             linestyle='solid'):

        if ax is None:
            fig, ax = plt.subplots(1, 1)

        keys = sorted(self.jumps_counter.keys())
        y = [self(k) + self.jumps_counter[k] for k in keys]
        xmax = keys[1:]
        xmin = keys[:-1]

        ax.set(ylabel='signature function')
        ax.set(title=title)
        ax.hlines(y, xmin, xmax, color=color, linestyle=linestyle, alpha=alpha)

        if subplot:
            return ax

        save_as += ".png"
        plt.savefig(save_as)
        plt.close()
        image = Image.open(save_as)
        image.show()



    def step_function_data(self):
        # Transform the signature jump data to a format understandable
        # by the plot function.
        result = [(k, self.sf(k) + self.jumps_counter[k])
                 for k in sorted(self.jumps_counter.keys())]
        return result

    def tikz_plot(self, save_as):
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
        with open(save_as, "w") as f:
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

def mod_one(n):
    return n - floor(n)

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
