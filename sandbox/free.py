def use_string(wrap):
    print wrap.s
    wrap.s = None
    print wrap.s

class wrapper(object):
    pass

a = wrapper()
a.s = 'olle'
print a.s
use_string(a)
print a.s
