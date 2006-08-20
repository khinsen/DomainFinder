import operator

class Variable:

    def __init__(self, name, value):
	self.name = name
	self.value = value
	self.arguments = []
	self.dependents = []
	self.function = None

    def defineFunction(self, arguments, function):
	self.function = function
	for arg in arguments:
	    if arg not in self.arguments:
		self.arguments.append(arg)
	    if self not in arg.dependents:
		arg.dependents.append(self)

    def invalidate(self):
	self.value = None
	for dependent in self.dependents:
	    dependent.invalidate()

    def setValue(self, new_value, override = 0):
	if not override and self.function is not None:
	    raise ValueError, "Variable " + self.name + " is defined by " + \
		              " a functional relation"
	if self.value != new_value:
	    self.value = new_value
	    for dependent in self.dependents:
		dependent.invalidate()

    def getValue(self):
	if self.value is not None:
	    return self.value
	if self.function is None:
	    raise ValueError, "No evaluation function defined for " + self.name
	arguments = []
	for arg in self.arguments:
	    arguments.append(arg.getValue())
	self.value = apply(self.function, tuple(arguments))
	return self.value

    def __repr__(self):
	return self.name + ': ' + repr(self.value)

    def __str__(self):
	return self.name + ': ' + str(self.value)


class DataGraph:

    def __init__(self, **kwargs):
	self._variables = {}
	for name, value in kwargs.items():
	    self._variables[name] = Variable(name, value)

    def __setattr__(self, name, value):
	if name[0] == '_':
	    self.__dict__[name] = value
	else:
	    self._variables[name].setValue(value)

    def __getattr__(self, name):
	return self.__dict__['_variables'][name].getValue()

    def defineFunction(self, variable, arguments, function):
	vars = []
	for name in arguments:
	    vars.append(self._variables[name])
	self._variables[variable].defineFunction(vars, function)

    def setValue(self, name, value):
        self._variables[name].setValue(value, 1)


if __name__ == '__main__':

    def get_x1(x2):
	return 2*x2
    def get_x2(x3):
	return x3+1

    test = DataGraph(x1 = None, x2 = None, x3 = None)
    test.defineFunction('x1', ['x2'], get_x1)
    test.defineFunction('x2', ['x3'], get_x2)
    test.x3 = 2
    print test.x1
