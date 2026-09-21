#!/usr/bin/env python3
"""Change one parameter of a generated deck (one-at-a-time sensitivity analysis).

    python set_parameter.py <deck> <name> <value>
    name: kb   microbial rate constant [M s-1]            (baseline 3.d-7)
          KO2  O2 inhibition constant [M]                 (baseline 1.d-6)
          O2   inlet dissolved O2 [M]                     (baseline 2.813d-5)
          kcal multiplier on the calcite rate constants   (baseline 1)
          kpyr multiplier on the pyrite rate constant     (baseline 1)
"""
import re, sys

deck, name, val = sys.argv[1:4]
s = open(deck).read()


def replace_once(pattern, repl):
    global s
    s, n = re.subn(pattern, repl, s, count=1)
    if n != 1: sys.exit("set_parameter(%s): expected exactly one match, found %d" % (name, n))


if name == "kb":
    replace_once(r"RATE_CONSTANT 3\.d-7", "RATE_CONSTANT " + val)
elif name == "KO2":
    replace_once(r"THRESHOLD_CONCENTRATION 1\.d-6", "THRESHOLD_CONCENTRATION " + val)
elif name == "O2":
    replace_once(r"2\.813d-5", val)
elif name in ("kcal", "kpyr"):
    mineral, factor = ("Calcite" if name == "kcal" else "Pyrite"), float(val)
    start = s.index("MINERAL_KINETICS")
    m = re.compile(r"\n    " + mineral + r"\b[^\n]*\n.*?\n    /\n", re.S).search(s, start)
    if not m: sys.exit("set_parameter: %s kinetics block not found" % mineral)
    scaled = re.sub(r"RATE_CONSTANT\s+([0-9.eEdD+-]+) mol/m\^2-sec",
                    lambda k: "RATE_CONSTANT   %.4e mol/m^2-sec" % (float(k.group(1).replace("d", "e").replace("D", "e")) * factor), m.group(0))
    s = s[:m.start()] + scaled + s[m.end():]
else:
    sys.exit("unknown parameter: " + name)
open(deck, "w").write(s)
print("set %s = %s in %s" % (name, val, deck))
