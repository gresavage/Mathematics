__author__ = 'Tom'

def cast_as(obj, type):
    if isinstance(obj, (list, set, tuple, dict)):
        if type is list:
            return list(obj)
        elif type is set:
            return set(obj)
        elif type is tuple:
            return tuple(obj)
        elif type is dict:
            return {i: obj[i] for i in range(len(obj))}
    elif isinstance(obj, (int, float, long)):
        if type is list:
            return [obj]
        elif type is set:
            return {obj}
        elif type is tuple:
            return (obj)
        elif type is dict:
            return {0: obj}