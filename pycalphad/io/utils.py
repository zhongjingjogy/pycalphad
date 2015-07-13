from sympy.printing.ccode import CCodePrinter
from sympy.printing.codeprinter import Assignment
class CustomCCodePrinter(CCodePrinter):
    
    def _print_Piecewise(self, expr):
        lines = []
        if expr.has(Assignment):
            for i, (e, c) in enumerate(expr.args):
                if i == 0:
                    lines.append("if (%s) {" % self._print(c))
                elif i == len(expr.args) - 1 and c == True:
                    lines.append("else {")
                else:
                    lines.append("else if (%s) {" % self._print(c))
                code0 = self._print(e)
                lines.append(code0)
                lines.append("}")
            return "\n".join(lines)
        
        else:
            # The piecewise was used in an expression, need to do inline
            # operators. This has the downside that inline operators will
            # not work for statements that span multiple lines (Matrix or
            # Indexed expressions).
            #print expr.args
            
            ecpairs = ["((%s) ? (%s)" % (self._print(c), self._print(e))
                    for e, c in expr.args]
            last_line = ": (0)"
            
            return ": ".join(ecpairs) + last_line + " ".join([")"*len(ecpairs)])

def ccode(expr, assign_to=None, **settings):
    return CustomCCodePrinter(settings).doprint(expr, assign_to)
    
def print_ccode(expr, **settings):
    """Prints C representation of the given expression."""
    print(ccode(expr, **settings))
