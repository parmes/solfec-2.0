# https://stackoverflow.com/questions/279561/what-is-the-python-equivalent-of-static-variables-inside-a-function
# Decorate function with static variables:
def static_vars(**kwargs):
    def decorate(func):
        for k in kwargs:
            setattr(func, k, kwargs[k])
        return func
    return decorate

# Example:
# @static_vars(counter = 0)
# def foo():
#   print(foo.counter)
#   foo.counter += 1
