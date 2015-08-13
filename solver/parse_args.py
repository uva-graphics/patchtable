
import sys
import string

def parse_args(allowed_args=[], usage=None, allow_all=False):
    args = sys.argv[1:]
    ans = []
    kw = {}
    i = 0
    while i < len(args):
        s = args[i]
        if len(s) >= 2 and s.startswith('-') and not s[1] in string.digits:
            if (s[1:] in allowed_args and i+1 < len(args)) or allow_all:
                kw[s[1:]] = args[i+1]
            else:
                if usage is not None:
                    usage()
                raise KeyError
            i += 2
        else:
            ans.append(s)
            i += 1
    return (ans, kw)

