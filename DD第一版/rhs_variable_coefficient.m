function result = rhs_variable_coefficient(uh,A,M,b)

result = -A*uh-M*uh + b;