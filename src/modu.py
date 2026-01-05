"""
MIT License

Copyright (c) 2025-2026 Pierre Denis

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

from __future__ import annotations
import operator
import inspect
from itertools import zip_longest, product
from math import gcd, lcm
from typing import Iterable, Iterator, Any, Callable

__version__ = "0.1.1"


class ModuError(Exception):
    pass


class Modu:
    """ A Modu instance represents a residues system for a given modulus. It provides several methods for
        arithmetic (overloading +, -, *, /, ** operators) and for set operations (overloading &, |, ~, ==, !=,
        in, not in operators).
        A Modu instance has the following attributes:
        * _modulus: integer modulus >= 1 or None to represent a finite subset of integers
        * _residues: tuple of integers giving, if _modulus is not None, the representative residues in the range
                     [0, _modulus - 1]
        * _name: string giving a name to the instance or None if no name is given
        * _format_spec: string giving output format specification
        * _table_func: function to apply to output table elements or None to use default
    """

    __slots__ = ('_modulus', '_residues', '_name', '_format_spec', '_table_func')

    # integer modulus >= 1  or None to represent a finite subset of integers
    _modulus: int | None

    # tuple of integers giving, if _m is not None, the representative residues in the range [0, _m - 1]
    _residues: tuple[int, ...]

    # string giving a name to the instance or None if no name is given
    _name: str | None

    # string giving output format specification
    _format_spec: str

    # function to apply to output table elements or None to use default
    _table_func: Callable[[int], Any] | None

    # default display format used in text environments (terminal)
    default_str_format: str = "s"

    # default display format used in LaTeX environment of Jupyter Notebook
    default_latex_format: str = "L"

    # default number of rows for table output format
    default_table_end: int = 10

    # allowed format codes
    __FORMAT_CODES: frozenset = frozenset("Llistpm+-:0123456789")

    @staticmethod
    def build(modulus: int | None = None, residues: Iterable[int] | None = None, name: str | None = None,
              format_spec: str = "", table_func: Callable[[int], Any] | None = None) -> Modu:
        if modulus is not None and (not isinstance(modulus, int) or modulus <= 0):
            raise ModuError(f"invalid modulus {modulus}: expected strictly positive integer or None for no modulus")
        if residues is None:
            return Modu(modulus, (0,), name, format_spec, table_func)
        try:
            residues = tuple(residues)
        except TypeError:
            residues = (residues,)
        if not all(isinstance(r, int) for r in residues):
            raise ModuError("all residues shall be integer")
        if name is not None and (not isinstance(name, str) or len(name) == 0):
            raise ModuError("the name shall be a non-empty string")
        _ = Modu.__decode_format_spec(format_spec)
        return Modu(modulus, residues, name, format_spec, table_func)

    def __init__(self, modulus: int | None, residues: Iterable[int] | None, name: str | None = None,
                 format_spec: str = "", table_func: Callable[[int], Any] | None = None) -> None:
        self._modulus = modulus
        if modulus is not None:
            residues = (r % modulus for r in residues)
        self._name = name
        self._residues = tuple(sorted(frozenset(residues)))
        self._format_spec = format_spec
        self._table_func = table_func

    def normalized(self) -> Modu:
        # attempt to normalize using a divisor of m as modulus
        modulus = self._modulus
        residues = self._residues
        if modulus is not None:
            if len(residues) >= 2:
                for modulus1 in range(1, modulus//2+1):
                    if modulus % modulus1 == 0:
                        if all(r1 == r2
                               for (r1, r2) in zip_longest(residues,
                                                           (b+r for b in range(0, modulus, modulus1)
                                                                for r in residues if r < modulus1))):
                            modulus = modulus1
                            residues = tuple(r for r in residues if r < modulus1)
                            break
            if len(residues) == 0:
                modulus = 1
        if modulus == self._modulus:
            return self
        return Modu(modulus, residues, self._name, self._format_spec)

    def residues(self, is_plus_form: bool = False) -> Iterable[int]:
        if is_plus_form or self._modulus is None:
            return self._residues
        m2 = self._modulus // 2
        return tuple(sorted(r if r <= m2 else r-self._modulus for r in self._residues))

    def __iter__(self) -> Iterator[int]:
        return iter(self.residues())

    def gen_expand(self, start: int = 0, end: int = 1) -> Iterator[int]:
        if self._modulus is None:
            if start != 0 or end != 1:
                raise ModuError(f"without modulus, expansion of residues requires start=0 and end=1")
            yield from self._residues
        else:
            for b in range(start*self._modulus, end * self._modulus, self._modulus):
                for r in self._residues:
                    yield b + r

    def gen_expand_table(self, start: int = 0, end: int = 1, is_plus_form: bool = True) -> Iterator[Iterator[int]]:
        if self._modulus is None:
            if start != 0 or end != 1:
                raise ModuError(f"without modulus, expansion of residues requires start=0 and end=1")
            yield iter(self._residues)
        else:
            rs = self.residues(is_plus_form)
            for b in range(start*self._modulus, end * self._modulus, self._modulus):
                yield (b+r for r in rs)

    def expand_table_str(self, start: int = 0, end: int = 1, is_plus_form: bool = True) -> tuple[tuple[str, ...], ...]:
        f = str
        if self._table_func is not None:
            def g(r):
                return str(self._table_func(r))
            if 'l' in self._format_spec or 'L' in self._format_spec:
                def f(r):
                    a = g(r)
                    if a.strip():
                        return f"\\text{{{a}}}"
                    return a
            else:
                f = g
        return tuple(tuple(f(r) for r in it)
                     for it in self.gen_expand_table(start, end, is_plus_form))

    def expand(self, start: int = 0, end: int = 1) -> Modu:
        return Modu(None, (self.gen_expand(start, end)))

    def gen_complement_residues(self) -> Iterator[int]:
        return (r for r in self.all_residues() if r not in self._residues)

    @staticmethod
    def __str_residue(r: int) -> str:
        return ('+' if r > 0 else '') + str(r)

    def get_name(self) -> str | None:
        if self._name is not None:
            return self._name
        frame = inspect.currentframe()
        while frame is not None:
            for (name, obj) in frame.f_globals.items():
                if obj is self and '_' not in name:
                    self._name = name
                    return name
            frame = frame.f_back
        return None

    def __gen_str_residues(self, is_plus_form: bool, is_minus_form: bool, latex: bool) -> Iterator[str]:
        if is_plus_form ^ is_minus_form:
            rs_iter = self if is_minus_form else self._residues
            str_rs = ('0' if r == 0 else f"{r:+d}" for r in rs_iter)
        else:
            if latex:
                pm_char = "\\pm"
            else:
                pm_char = '±'
            rs_with_signs = []
            if self._modulus is None:
                for r in self._residues:
                    if r < 0:
                        if -r in self._residues:
                            rs_with_signs.append((-r, pm_char))
                        else:
                            rs_with_signs.append((-r, '-'))
                    elif (r, pm_char) not in rs_with_signs:
                        rs_with_signs.append((r, '+' if r > 0 else ''))
            else:
                m2 = self._modulus // 2
                for r in self._residues:
                    if r <= m2:
                        if self._modulus-r in self._residues and (self._modulus % 2 == 1 or r != m2):
                            rs_with_signs.append((r, pm_char))
                        else:
                            rs_with_signs.append((r, '+' if r > 0 else ''))
                    elif self._modulus-r not in self._residues:
                        rs_with_signs.append((self._modulus - r, '-'))
            str_rs = (f"{s}{r}" for (r, s) in sorted(rs_with_signs))
        return str_rs

    @staticmethod
    def __decode_format_spec(format_spec: str) -> tuple[bool, bool, bool, bool, bool, bool, int, int]:
        if not isinstance(format_spec, str):
            raise ModuError(f"the format specification requires a string with format codes among "
                            f"{sorted(Modu.__FORMAT_CODES)}")
        unknown_codes = frozenset(format_spec) - Modu.__FORMAT_CODES
        if len(unknown_codes) > 0:
            raise ModuError(f"Unknown format specification '{next(iter(unknown_codes))}': "
                            f"requires format code(s) among {sorted(Modu.__FORMAT_CODES)}")
        is_inverted_form = 'i' in format_spec
        is_short_latex = 'l' in format_spec
        is_big_latex = 'L' in format_spec
        is_string_format = 's' in format_spec
        is_plus_form = 'p' in format_spec
        is_minus_form = 'm' in format_spec
        is_table_form = 't' in format_spec
        if is_table_form and not (is_plus_form | is_minus_form):
            is_minus_form = True
        start = 0
        end = Modu.default_table_end
        if is_table_form:
            (_, table_bounds_spec) = format_spec.split('t')
            table_bounds_spec = table_bounds_spec.strip()
            while len(table_bounds_spec) > 0 and not table_bounds_spec[-1].isdigit():
                table_bounds_spec = table_bounds_spec[:-1]
            if len(table_bounds_spec) > 0:
                try:
                    end = int(table_bounds_spec)
                except ValueError:
                    try:
                        (start_str, end_str) = table_bounds_spec.split(':')
                        start = int(start_str)
                        end = int(end_str)
                    except ValueError:
                        raise ModuError(f"invalid table format 't{table_bounds_spec}'") from None
        if (is_inverted_form and is_table_form) or (is_string_format and (is_short_latex or is_big_latex)) \
                or (is_short_latex and is_big_latex):
            raise ModuError(f"invalid format specification '{format_spec}'")
        return (is_inverted_form, is_short_latex, is_big_latex, is_plus_form, is_minus_form, is_table_form, start, end)

    def __format__(self, format_spec: str) -> str:
        (is_inverted_form, is_short_latex, is_big_latex, is_plus_form, is_minus_form, is_table_form, start, end) \
            = Modu.__decode_format_spec(format_spec)
        is_latex_format = is_short_latex or is_big_latex
        modu1 = ~self if is_inverted_form else self
        str_residues = tuple(modu1.__gen_str_residues(is_plus_form, is_minus_form, is_latex_format))
        if is_table_form:
            rows = self.expand_table_str(start, end, is_plus_form)
            if is_latex_format:
                table_header_str = " & ".join(f"{str_residue}" for str_residue in str_residues)
                rows_str = "\\\\\n".join(" & ".join(str_residue for str_residue in row)
                                         for row in rows)
                nb_cols = len(self._residues)
                out = f"\\begin{{array}}\n" \
                      f"{{{'|'.join('c' for _ in range(nb_cols))}}}\n" \
                      f"{table_header_str}\\\\\n" \
                      f"\\hline\\hline\n" \
                      f"{rows_str}\n" \
                      f"\\end{{array}}"
            else:
                if len(str_residues) == 0 or len(rows) == 0:
                    col_width = 1
                else:
                    col_width = max(max(len(r) for r in str_residues),
                                    max(len(r) for row in rows for r in row))
                table_header_str = " ".join(f"{str_residue:>{col_width}}"
                                            for str_residue in str_residues)
                rows_str = "\n".join(" ".join(f"{str_residue:>{col_width}}" for str_residue in row)
                                     for row in rows)
                out = f"mod {self._modulus}\n" \
                      f"{table_header_str}\n" \
                      f"{len(table_header_str)*'-'}\n" \
                      f"{rows_str}"
        else:
            count = len(str_residues)
            head = "\\{"
            tail = "\\}"
            delimiter = ", "
            if is_latex_format:
                if not (is_short_latex or count <= 1 or self._modulus is None) and is_big_latex:
                    head = "\\left\\{\\begin{array}{l} "
                    tail = " \\end{array}\\right."
                    delimiter = "\\\\"
            else:
                head = '{'
                tail = '}'
                delimiter = ", "
            out = delimiter.join(str_residues)
            if self._modulus is None or count != 1:
                out = f"{head}{out}{tail}"
            if is_latex_format:
                if self._modulus is None:
                    equiv_char = '='
                else:
                    equiv_char = "\\equiv"
                    out += f" \\pmod{{{self._modulus}}}"
                if is_inverted_form:
                    equiv_char = "\\not " + equiv_char
            else:
                if self._modulus is None:
                    equiv_char = '=' if not is_inverted_form else '≠'
                else:
                    equiv_char = '≡' if not is_inverted_form else '≢'
                    out += f" (mod {self._modulus})"
            name = self.get_name()
            if name is not None:
                out = f"{name} {equiv_char} {out}"
            elif is_inverted_form:
                out = f"{equiv_char} {out}"
        if is_latex_format:
            if is_table_form:
                out = f"\\begin{{array}}{{c}}\n" \
                      f"\\textbf{{mod }}{self._modulus}\\\\\n" \
                      f"{out}\n" \
                      f"\\end{{array}}"
            out = f"$ {out} $"
        return out

    def _repr_latex_(self) -> str:
        format_spec = self._format_spec
        if Modu.default_latex_format not in format_spec:
            try:
                Modu.__decode_format_spec(Modu.default_latex_format + format_spec)
                format_spec = Modu.default_latex_format + format_spec
            except ModuError:
                pass
        return format(self, format_spec)

    def __str__(self) -> str:
        format_spec = self._format_spec
        if Modu.default_str_format not in format_spec:
            try:
                Modu.__decode_format_spec(Modu.default_str_format + format_spec)
                format_spec = Modu.default_str_format + format_spec
            except ModuError:
                pass
        return format(self, format_spec)

    __repr__ = __str__

    def __rshift__(self, other: str | Callable[[int], Any]) -> Modu:
        name = self.get_name()
        if other is None or isinstance(other, str):
            return Modu(self._modulus, self._residues, name, self._format_spec + other, self._table_func)
        if callable(other):
            format_spec = self._format_spec
            if 't' not in format_spec:
                format_spec += 't'
            return Modu(self._modulus, self._residues, name, format_spec, other)
        raise ModuError(f"invalid type '{type(other).__name__}' for >> operator: requires string or callable")

    def all_residues(self) -> range:
        if self._modulus is None:
            raise ModuError(f"cannot represent result without modulus (infinitely many elements)")       
        return range(self._modulus)

    def _coerce(self, other: Any) -> Modu:
        if isinstance(other, Modu):
            return other
        return Modu.build(self._modulus, other)

    @staticmethod
    def _prepare_combinatorics(modu1, modu2) -> tuple[int, Iterator[int], Iterator[int]]:
        if modu1._modulus is None or modu2._modulus is None:
            modulus = modu1._modulus or modu2._modulus
            m1 = 1
            m2 = 1
        else:
            modulus = lcm(modu1._modulus, modu2._modulus)
            m1 = modulus // modu1._modulus
            m2 = modulus // modu2._modulus
        return (modulus, modu1.gen_expand(end=m1), modu2.gen_expand(end=m2))

    # set operations

    def __contains__(self, r: int) -> bool:
        if self._modulus is not None:
            r = r % self._modulus
        return r in self._residues

    def __invert__(self) -> Modu:
        return Modu.build(self._modulus, self.gen_complement_residues(), format_spec=self._format_spec)

    # temporary helper meta-function to create set magic methods below
    def __make_op_set(inv: int, f: Callable[[Any, Any], Any]) -> Callable[[Modu, Any], Modu]:
        def func(modu1: Modu, modu2: Any) -> Modu:
            modu2 = modu1._coerce(modu2)
            if inv:
                (modu1, modu2) = (modu2, modu1)
            (modulus, rs1_iter, rs2_iter) = Modu._prepare_combinatorics(modu1, modu2)
            return Modu(modulus, f(frozenset(rs1_iter), frozenset(rs2_iter)), format_spec=modu1._format_spec)
        func.__name__ = f.__name__
        r = 'other, self' if inv else 'self, other'
        func.__doc__ = f"returns Modu instance applying {f.__name__} function on ({r}), for operator overloading"
        return func

    __and__  = __make_op_set(0, operator.and_)
    __rand__ = __make_op_set(1, operator.and_)
    __or__   = __make_op_set(0, operator.or_)
    __ror__  = __make_op_set(1, operator.or_)

    del __make_op_set

    # arithmetic operations

    def __truediv__(self, other: Any) -> Modu:
        return self * self._coerce(other)**-1

    def __rtruediv__(self, other: Any) -> Modu:
        return self._coerce(other) * self**-1

    def __pow__(self, other: Any) -> Modu:
        if not isinstance(other, int):
            raise ModuError(f"power value shall be an integer")
        if self._modulus is None and other < 0:
            raise ModuError(f"without modulus, the power shall be a non-negative integer")
        try:
            return Modu(self._modulus,
                        (pow(r, other, self._modulus) for r in self._residues),
                        format_spec=self._format_spec)
        except ValueError:
            modulus = self._modulus
            r = None
            for r in self._residues:
                if gcd(r, modulus) > 1:
                    break
            raise ModuError(f"cannot find inverse of {r} modulo {modulus}, because {r} and {modulus} are not coprime") \
                    from None

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, Modu):
            raise ModuError(f"cannot compare {type(self).__name__} instance with {type(other).__name__} instance")
        m1 = self.normalized()
        m2 = other.normalized()
        return m1._modulus == m2._modulus and m1._residues == m2._residues

    def __mod__(self, other: int | None) -> Modu:
        if self._modulus is None:
            end = 1
        else:
            end = lcm(self._modulus, other) // self._modulus
        return Modu.build(other, (r1 % other for r1 in self.gen_expand(end=end)), format_spec=self._format_spec)

    def __pos__(self) -> Modu:
        return self.normalized()

    def __neg__(self) -> Modu:
        return Modu(self._modulus, (-r for r in self._residues), format_spec=self._format_spec)

    # temporary helper meta-function to create arithmetic magic methods below
    def __make_op_int(inv: int, f: Callable[[int, int], int]) -> Callable[[Modu, Any], Modu]:
        def func(modu1: Modu, modu2: Any) -> Modu:
            modu2 = modu1._coerce(modu2)
            if inv:
                (modu1, modu2) = (modu2, modu1)
            (modulus, rs1_iter, rs2_iter) = Modu._prepare_combinatorics(modu1, modu2)
            return Modu(modulus, (f(*r1r2) for r1r2 in product(rs1_iter, rs2_iter)), format_spec=modu1._format_spec)
        func.__name__ = f.__name__
        r = 'other, self' if inv else 'self, other'
        func.__doc__ = f"returns Modu instance applying {f.__name__} function on ({r}), for operator overloading"
        return func

    __add__  = __make_op_int(0, operator.add)
    __radd__ = __make_op_int(1, operator.add)
    __sub__  = __make_op_int(0, operator.sub)
    __rsub__ = __make_op_int(1, operator.sub)
    __mul__  = __make_op_int(0, operator.mul)
    __rmul__ = __make_op_int(1, operator.mul)

    del __make_op_int


m = mod = Modu.build
E = Modu(None, ())
O = Modu(None, (0,))
I = Modu(None, (+1,))
T = Modu(None, (-1, +1))

__all__ = ("m", "mod", "Modu", "E", "O", "I", "T")
